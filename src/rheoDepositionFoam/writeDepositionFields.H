const label iteration = label(runTime.value()+0.5);
const bool monitor = (iteration%monitorEvery == 0) || runTime.writeTime() || !runTime.running();
if (monitor)
{
    tauTotal = constEq.tauTotal();
    pPa = p*constEq.rho();
    const tmp<volSymmTensorField> polymer = constEq.tau();
    N1 = polymer().component(symmTensor::XX)-polymer().component(symmTensor::YY);
    scalar qIn=-gSum(phi.boundaryField()[inletID])/depth;
    scalar qOut=gSum(phi.boundaryField()[outletID])/depth;
    scalar massError=max(mag(qIn-qOut),mag(qIn-qHalf))/qHalf;
    const scalar residualU=equationResidual<vector>(mesh,"U");
    const scalar residualP=equationResidual<scalar>(mesh,"p");
    const scalar residual1=elastic ? equationResidual<symmTensor>(mesh,"theta_mode1") : 0;
    const scalar residual2=elastic ? equationResidual<symmTensor>(mesh,"theta_mode2") : 0;

    const vectorField outN(mesh.boundary()[outletID].nf());
    scalar outletTraction=0, outletStressScale=1e-8;
    forAll(outN, i)
    {
        const symmTensor& t = tauTotal.boundaryField()[outletID][i];
        outletTraction=max(outletTraction,mag((t & outN[i])-pPa.boundaryField()[outletID][i]*outN[i]));
        outletStressScale=max(outletStressScale,mag(t));
    }
    reduce(outletTraction,maxOp<scalar>()); reduce(outletStressScale,maxOp<scalar>());

    // Positive signed velocity points from fluid toward the respective wall.
    scalarField wallNormal(nWall,0);
    forAll(samplePoints,i) if (owner[i]==Pstream::myProcNo())
    {
        const direction axis = label(wallRows[i][0])==1 ? 0 : 1;
        wallNormal[i]=planarNormalVelocity(mesh,phi,sampleCell[i],samplePoints[i],axis);
    }
    Pstream::listCombineGather(wallNormal,plusEqOp<scalar>());
    Pstream::listCombineScatter(wallNormal);

    vectorField stress(nWall,vector::zero); // total, polymer mode 1, polymer mode 2
    label globalRow=wallStart;
    for (label w=0; w<3; ++w)
    {
        const label k=wallPatchID[w];
        const vectorField n(mesh.boundary()[k].nf());
        const vector tang = w==1 ? vector(0,-1,0) : vector(1,0,0);
        forAll(n, i)
        {
            vector& v=stress[globalRow++];
            v.x() = tang & (tauTotal.boundaryField()[k][i] & n[i]);
            if (elastic)
            {
                v.y() = tang & (mesh.lookupObject<volSymmTensorField>("tau_mode1").boundaryField()[k][i] & n[i]);
                v.z() = tang & (mesh.lookupObject<volSymmTensorField>("tau_mode2").boundaryField()[k][i] & n[i]);
            }
        }
    }
    Pstream::listCombineGather(stress,plusEqOp<vector>());
    Pstream::listCombineScatter(stress);
    FixedList<scalar,3> q(scalar(0)), wallDiff(scalar(0)), wallNorm(scalar(0));
    FixedList<scalar,3> normalDiff(scalar(0)), normalNorm(scalar(0));
    forAll(wallRows,i)
    {
        const label w=label(wallRows[i][0]);
        const scalar ds=wallRows[i][4];
        q[w]+=max(wallNormal[i],scalar(0))*ds;
        wallDiff[w]+=mag(stress[i].x()-prevWall[i])*ds;
        wallNorm[w]+=mag(stress[i].x())*ds;
        normalDiff[w]+=mag(wallNormal[i]-prevNormal[i])*ds;
        normalNorm[w]+=mag(wallNormal[i])*ds;
        prevWall[i]=stress[i].x();
    }
    scalar relWall=0, relFlux=0, relNormal=0;
    for (label w=0; w<3; ++w)
    {
        relWall=max(relWall,wallDiff[w]/max(wallNorm[w],scalar(1e-13)));
        relFlux=max(relFlux,mag(q[w]-prevQ[w])/max(mag(q[w]),qHalf*1e-5));
        relNormal=max(relNormal,normalDiff[w]/max(normalNorm[w],qHalf*1e-5));
    }
    consecutive = havePrevious && max(max(residualU,residualP),max(residual1,residual2))<resTol &&
        max(relWall,max(relNormal,relFlux))<depTol && massError<massTol &&
        outletTraction/outletStressScale<depTol ? consecutive+1 : 0;
    const bool stationary=consecutive>=stableChecks;
    prevNormal=wallNormal; prevQ=q; havePrevious=true;

    if (Pstream::master())
    {
        std::ofstream out((depOutDir/"metrics.csv").c_str(),std::ios::app);
        out.precision(16);
        out << iteration << ',' << qIn << ',' << qOut << ',' << massError << ','
            << residualU << ',' << residualP << ',' << residual1 << ',' << residual2 << ','
            << relWall << ',' << relFlux << ',' << relNormal << ',' << q[1] << ',' << q[1]/qHalf
            << ',' << outletTraction << ',' << consecutive << '\n';
        out.close();
        if (!out) FatalErrorInFunction << "Writing metrics.csv failed" << exit(FatalError);
        Info << "Steady monitor: mass=" << massError << " residual(U,p,theta1,theta2)="
             << residualU << ',' << residualP << ',' << residual1 << ',' << residual2
             << " step approach/Qhalf=" << q[1]/qHalf << " wall change=" << relWall
             << " outlet traction [Pa]=" << outletTraction << " checks=" << consecutive << nl;
    }
    if (stationary)
    {
        if (Pstream::master())
        {
            std::ofstream out((depOutDir/"wall_profile.csv").c_str());
            out.precision(16);
            out << "iteration,wall,s_m,x_wall_m,y_wall_m,x_sample_m,y_sample_m,ds_m,"
                << "u_toward_m_s,tau_total_Pa,tau_viscous_Pa,tau_mode1_Pa,tau_mode2_Pa\n";
            forAll(wallRows,i)
            {
                const WallRow& r=wallRows[i];
                out << iteration << ',' << wallNames[label(r[0])] << ',' << r[1] << ',' << r[2]
                    << ',' << r[3] << ',' << samplePoints[i].x() << ',' << samplePoints[i].y()
                    << ',' << r[4] << ',' << wallNormal[i] << ',' << stress[i].x()
                    << ',' << stress[i].x()-stress[i].y()-stress[i].z()
                    << ',' << stress[i].y() << ',' << stress[i].z() << '\n';
            }
            out.close();
            if (!out) FatalErrorInFunction << "Writing wall_profile.csv failed" << exit(FatalError);
        }
        writeCells(mesh,U,phi,depth,iteration,depOutDir);
        if (!runTime.writeAndEnd())
            FatalErrorInFunction << "Writing the converged checkpoint failed" << exit(FatalError);
        if (Pstream::master())
        {
            std::ofstream marker((depOutDir/"stationary.json").c_str());
            marker.precision(16);
            marker << "{\"status\":\"stationary\",\"iteration\":" << iteration << ",\"time\":" << iteration
                   << ",\"qHalfPerDepth_m2_s\":" << qHalf << ",\"massImbalance\":" << massError
                   << ",\"stableChecks\":" << consecutive << "}\n";
            marker.close();
            if (!marker) FatalErrorInFunction << "Writing stationary.json failed" << exit(FatalError);
        }
    }
}
