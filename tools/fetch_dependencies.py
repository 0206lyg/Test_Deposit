#!/usr/bin/env python3
"""Install one header-only dependency locally, checking the upstream archive hash."""
from pathlib import Path
import hashlib
import tarfile
import urllib.request

ROOT = Path(__file__).resolve().parents[1]
SHA = '8586084f71f9bde545ee7fa6d00288b264a2b7ac3607b974e54d13e7162c1c72'
URL = 'https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz'

def main():
    dest = ROOT / 'deps'
    dest.mkdir(exist_ok=True)
    archive = dest / 'eigen-3.4.0.tar.gz'
    if not archive.exists():
        temp = archive.with_suffix('.download')
        with urllib.request.urlopen(URL, timeout=120) as src, temp.open('wb') as out:
            while True:
                data = src.read(1024*1024)
                if not data:
                    break
                out.write(data)
        temp.rename(archive)
    if hashlib.sha256(archive.read_bytes()).hexdigest() != SHA:
        raise SystemExit('Eigen archive checksum mismatch; extraction refused.')
    # Validate every archive path, including links, before extraction.
    with tarfile.open(archive) as tar:
        for member in tar.getmembers():
            target = (dest/member.name).resolve()
            if dest.resolve() not in target.parents or member.issym() or member.islnk():
                raise SystemExit('Unexpected path/link in Eigen archive.')
        tar.extractall(dest)
    assert (dest/'eigen-3.4.0/Eigen/Core').is_file()
    print('Eigen 3.4.0 verified and unpacked under deps/.')

if __name__ == '__main__':
    main()
