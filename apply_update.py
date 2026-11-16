#!/usr/bin/env python3
"""Apply the audited add-on update, preserving user launch scripts and case data."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil


def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('project',type=Path,help='Existing Deposition project directory')
    args=parser.parse_args()
    base=Path(__file__).resolve().parent
    target=args.project.resolve()
    if not (target/'run_dep.sbatch').is_file() or not (target/'vendor/rheoTool').is_dir():
        raise SystemExit('Target is not the existing Deposition project.')
    manifest=json.loads((base/'UPDATE_MANIFEST.json').read_text())
    changes=[]
    for name, hashes in manifest['files'].items():
        relative=Path(name)
        if relative.is_absolute() or '..' in relative.parts:
            raise SystemExit('Invalid relative path in manifest.')
        src=base/'payload'/relative
        dst=target/relative
        if sha(src)!=hashes['updated_sha256']:
            raise SystemExit('Update payload checksum mismatch: '+name)
        if dst.is_symlink() or target not in dst.resolve().parents:
            raise SystemExit('Refusing symlink or outside-project target: '+name)
        if dst.exists():
            digest=sha(dst)
            if digest==hashes['updated_sha256']: continue
            if digest!=hashes.get('original_sha256'):
                raise SystemExit('Locally modified target; no files changed: '+name)
        changes.append((name,src,dst,dst.exists()))
    if not changes:
        print('This update is already installed.'); return
    stamp=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    backup=target/'update_backups'/stamp
    backup.mkdir(parents=True)
    for name,src,dst,existed in changes:
        if existed:
            saved=backup/name
            saved.parent.mkdir(parents=True,exist_ok=True)
            shutil.copy2(dst,saved)
    applied=[]
    try:
        for name,src,dst,existed in changes:
            dst.parent.mkdir(parents=True,exist_ok=True)
            shutil.copy2(src,dst)
            applied.append((name,dst,existed))
    except BaseException:
        for name,dst,existed in reversed(applied):
            if existed: shutil.copy2(backup/name,dst)
            else: dst.unlink(missing_ok=True)
        raise
    (backup/'applied.json').write_text(json.dumps({'update':manifest['version'],
        'files':[name for name,_,_,_ in changes]},indent=2)+'\n')
    print('Update installed:',manifest['version'])
    print('Original source backups:',backup)
    print('Existing run_dep.sbatch/build_dep.sbatch, config, results and base OpenFOAM were preserved.')
    print('After old solver jobs have stopped: sbatch build_dep_update.sbatch')


if __name__=='__main__': main()
