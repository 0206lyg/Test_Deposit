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
    config=json.loads((ROOT/'config/experiment.json').read_text())
    peo=config['fluids']['peo']
    # Inspect generated dictionaries, not an independent copy of coefficients.
    with tempfile.TemporaryDirectory(prefix='deposition_verify_') as temp:
        case=Path(temp)/'peo'
        subprocess.run([sys.executable,str(ROOT/'tools/make_case.py'),'--fluid','peo',
                        '--flow','0.5','--mesh','smoke','--case',str(case),'--ranks','1'],
                       check=True,stdout=subprocess.DEVNULL)
        metadata=json.loads((case/'metadata.json').read_text())
        assert math.isclose(metadata['flow_rate_m3_s'],.5e-6/3600,rel_tol=1e-14)
        assert math.isclose(metadata['rho_kg_m3'],998,rel_tol=1e-14)
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
    print('PASS: case generation, flow units, stress/pressure dimensions, ePTT modes, non-overwrite protection and shell syntax.')

if __name__=='__main__': main()
