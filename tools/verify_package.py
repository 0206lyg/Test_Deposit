#!/usr/bin/env python3
"""Check physical units/model bookkeeping and shell syntax before a production run."""
from pathlib import Path
import json
import math
import subprocess
import sys
import tempfile

ROOT=Path(__file__).resolve().parents[1]

def main():
    # Inspect generated dictionaries, not an independent copy of coefficients.
    with tempfile.TemporaryDirectory(prefix='deposition_verify_') as temp:
        case=Path(temp)/'peo'
        subprocess.run([sys.executable,str(ROOT/'tools/make_case.py'),'--fluid','peo',
                        '--flow','0.5','--mesh','smoke','--case',str(case),'--ranks','1'],
                       check=True,stdout=subprocess.DEVNULL)
        metadata=json.loads((case/'metadata.json').read_text())
        assert math.isclose(metadata['flow_rate_m3_s'],.5e-6/3600,rel_tol=1e-14)
        assert math.isclose(metadata['rho_kg_m3'],998,rel_tol=1e-14)
        geom=metadata['geometry_m']
        assert math.isclose(geom['upstream_width']/geom['downstream_width'], 4)
        assert geom['upstream_length'] == geom['downstream_length'] == 900e-6
        assert math.isclose(metadata['q_half_per_depth_m2_s'],
                            metadata['flow_rate_m3_s']/(2*geom['height']))
        assert 'default steadyState;' in (case/'system/fvSchemes').read_text()
        velocity=(case/'0/U').read_text()
        assert 'type empty;' in velocity and 'type symmetryPlane;' in velocity
        assert 'type fixedGradient;' in velocity
        text=(case/'constant/constitutiveProperties').read_text()
        assert text.count('type PTTLog;')==2 or text.count('PTTLog')==2
        assert text.count('exponential')==2
        assert (case/'constant/depositionProperties').is_file()
        assert '[0 2 -2 0 0 0 0]' in (case/'0/p').read_text()
        assert '[1 -1 -2 0 0 0 0]' in (case/'0/tau_mode1').read_text()
        # Case creation must fail safely when the destination already exists.
        result=subprocess.run([sys.executable,str(ROOT/'tools/make_case.py'),'--fluid','peo',
                        '--flow','0.5','--mesh','smoke','--case',str(case),'--ranks','1'],
                       stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
        assert result.returncode!=0
    for script in [*ROOT.glob('*.sbatch'),*ROOT.glob('*.sh'),ROOT/'tools/load_environment.sh']:
        subprocess.run(['bash','-n',str(script)],check=True)
    assert (ROOT/'vendor/rheoTool/LICENSE').is_file()
    assert (ROOT/'deps/eigen-3.4.0/Eigen/Core').is_file()
    print('PASS: planar half geometry, flow per depth, steady schemes, velocity boundaries, ePTT modes, units, non-overwrite protection and shell syntax.')

if __name__=='__main__': main()
