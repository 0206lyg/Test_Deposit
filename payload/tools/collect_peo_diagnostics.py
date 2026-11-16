#!/usr/bin/env python3
"""Collect compact diagnostic evidence; leave full field snapshots in the case."""
import argparse
import io
import json
from pathlib import Path
import tarfile


def bounded_tail(path,limit=262144):
    with path.open('rb') as f:
        f.seek(max(0,path.stat().st_size-limit))
        return f.read()


def collect(case):
    case=case.resolve()
    if not case.is_dir():
        raise ValueError('Case does not exist: '+str(case))
    archive=case/'diagnostic_bundle.tar.gz'
    info={'case':str(case),'full_field_snapshots_included':False,
          'note':'failure field snapshots are kept in the case, outside restart directories',
          'snapshot_files':[]}
    with tarfile.open(archive,'w:gz') as tar:
        def add_bytes(name,data):
            member=tarfile.TarInfo(name)
            member.size=len(data)
            tar.addfile(member,io.BytesIO(data))
        for relative in ['metadata.json','run_status.json','constant/constitutiveProperties',
                         'constant/depositionProperties','system/controlDict','system/fvSchemes',
                         'system/fvSolution','system/decomposeParDict']:
            path=case/relative
            if path.is_file(): tar.add(path,arcname=relative,recursive=False)
        for name in ['metrics.csv','stationary.json','approach_latest.csv']:
            path=case/'postProcessing/deposition'/name
            if path.is_file(): add_bytes('postProcessing/deposition/'+name,bounded_tail(path))
        for path in sorted(case.glob('log.*')):
            if not path.is_file(): continue
            add_bytes('logs/'+path.name+'.tail.txt',bounded_tail(path))
            if path.name.startswith('log.solver'):
                with path.open('rb') as f: add_bytes('logs/'+path.name+'.head.txt',f.read(8192))
                lines=[]
                with path.open(errors='replace') as f:
                    for line in f:
                        if 'EPTT_' in line: lines.append(line)
                add_bytes('logs/'+path.name+'.eptt.txt',''.join(lines[-10000:]).encode())
        roots=[case/'postProcessing/epttDiagnostics']
        roots += [p/'postProcessing/epttDiagnostics' for p in case.glob('processor[0-9]*') if p.is_dir()]
        for root in roots:
            if not root.exists(): continue
            for path in sorted(root.rglob('*')):
                if path.is_file():
                    relative=str(path.relative_to(case))
                    info['snapshot_files'].append({'path':relative,'bytes':path.stat().st_size})
                    if path.name=='failure.txt': add_bytes(relative,bounded_tail(path,1048576))
        add_bytes('CONTENTS.json',(json.dumps(info,indent=2)+'\n').encode())
    print('Diagnostic evidence:',archive)
    print('Share this compact archive; full failure field snapshots remain in the case.')
    return archive


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('case',type=Path)
    args=parser.parse_args()
    collect(args.case)
