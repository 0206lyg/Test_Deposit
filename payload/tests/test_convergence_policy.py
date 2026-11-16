#!/usr/bin/env python3
"""Compile and exercise the actual C++ convergence policy and case defaults."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
SOURCE = r'''
#include "DepositionConvergence.H"
#include <cassert>
#include <limits>
using depositionConvergence::Tracker;
using depositionConvergence::withinTolerance;
int main()
{
    // A single passing interval and duplicate-time writes cannot stop a run.
    Tracker t;
    assert(!t.observe(false, .001, 3, 0, 0));
    assert(!t.observe(true, .011, 3, 0, 0));
    assert(!t.observe(true, .011, 3, 0, 0));
    assert(t.count() == 1);
    assert(!t.observe(true, .021, 3, 0, 0));
    assert(t.observe(true, .031, 3, 0, 0));
    // A failed physical observable resets the complete consecutive history.
    assert(!t.observe(false, .041, 3, 0, 0));
    assert(t.count() == 0);
    assert(!t.observe(true, .051, 3, 0, 0));
    // The PEO conservative minimum/window still applies when requested.
    t.reset();
    assert(!t.observe(true, .39, 3, .4, .2));
    assert(!t.observe(true, .4, 3, .4, .2));
    assert(!t.observe(true, .5, 3, .4, .2));
    assert(t.observe(true, .6, 3, .4, .2));
    // Near-zero stress/flux can converge through dimensional absolute tolerances.
    assert(withinTolerance(1e-21, 0, 1e-4, 1e-20));
    assert(!withinTolerance(2e-20, 0, 1e-4, 1e-20));
    assert(withinTolerance(.00009, 1, 1e-4, 0));
    assert(!withinTolerance(.0002, 1, 1e-4, 0));
    // Nonfinite values must never be accepted as converged.
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double inf = std::numeric_limits<double>::infinity();
    assert(!withinTolerance(nan, 0, 1e-4, 1e-20));
    assert(!withinTolerance(inf, inf, 1e-4, 1e-20));
}
'''

with tempfile.TemporaryDirectory(prefix='deposition_convergence_') as tmp:
    p = Path(tmp)
    (p/'check.C').write_text(SOURCE)
    subprocess.run(['g++', '-std=c++11', '-Wall', '-Wextra', '-Werror',
                    '-I'+str(ROOT/'src/rheoDepositionFoam'), str(p/'check.C'),
                    '-o', str(p/'check')], check=True)
    subprocess.run([str(p/'check')], check=True)
    for fluid in ('gly75','cross','peo'):
        case = p/fluid
        subprocess.run([sys.executable,str(ROOT/'tools/make_case.py'), '--fluid',fluid,
                        '--flow','0.5','--mesh','smoke','--case',str(case),'--ranks','1'],
                       check=True, stdout=subprocess.DEVNULL)
        d = json.loads((case/'metadata.json').read_text())['diagnostics']
        assert d['minimum_steady_time_s'] == (.4 if fluid == 'peo' else 0)
        assert d['steady_window_s'] == (.2 if fluid == 'peo' else 0)
        assert d['required_consecutive_checks'] == 3
        assert d['monitor_interval_s'] == .01
    bad = subprocess.run([sys.executable,str(ROOT/'tools/make_case.py'), '--fluid','gly75',
                        '--flow','0.5','--case',str(p/'bad'),'--steady-checks','1'],
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    assert bad.returncode != 0 and not (p/'bad').exists()
print('PASS: actual C++ count/window/tolerance policy and generated fluid-specific defaults')
