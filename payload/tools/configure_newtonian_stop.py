#!/usr/bin/env python3
"""Configure an existing stopped Gly75 case for convergence-based termination."""
import argparse
import json
import math
from pathlib import Path
import shutil
import subprocess
from datetime import datetime, timezone
from check_run import checkpoint


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--case',type=Path,required=True)
    p.add_argument('--ranks',type=int,required=True)
    a=p.parse_args()
    case=a.case.resolve()
    m=json.loads((case/'metadata.json').read_text())
    if m.get('fluid_key')!='gly75':
        raise SystemExit('Only a Gly75 case can use this helper.')
    if m.get('ranks')!=a.ranks:
        raise SystemExit('Use the original MPI rank count when resuming.')
    cp=checkpoint(case,a.ranks)
    if not math.isfinite(cp['time']) or cp['time']<=0:
        raise SystemExit('No evolved, complete checkpoint. Use a new case instead.')
    path=case/'constant/depositionProperties'
    stamp=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    backup=path.with_name(path.name+'.before_convergence_'+stamp)
    shutil.copy2(path,backup)
    values={'minimumSteadyTime':'0','steadyWindow':'0','requiredConsecutiveChecks':'3',
            'monitorInterval':'0.01','relativeTolerance':'1e-4','massTolerance':'1e-5',
            'stopWhenStationary':'true'}
    try:
        for key,value in values.items():
            subprocess.run(['foamDictionary',str(path),'-entry',key,'-set',value],check=True,
                           stdout=subprocess.DEVNULL)
    except BaseException:
        shutil.copy2(backup,path)
        raise
    print('Convergence settings applied; original dictionary saved at',backup)
    print('Resume checkpoint:',cp['time_name'],'s; three independent successful checks required.')


if __name__=='__main__': main()
