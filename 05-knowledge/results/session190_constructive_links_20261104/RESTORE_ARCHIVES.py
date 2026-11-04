#!/usr/bin/env python3
"""Verify and restore an explicit Session190 ARCHIVES_CP*.json manifest.

Run from the repository root, for example:
  python 05-knowledge/results/<session-directory>/RESTORE_ARCHIVES.py \
    05-knowledge/results/<session-directory>/ARCHIVES_CP01.json
Supply manifests in checkpoint order, oldest first. If an original path appears
in more than one manifest, the last record is the current expected version.
Superseded records are reported. Existing identical current files are verified;
different existing files are never overwritten. Only explicitly listed paths
are read or written.
"""
import argparse, gzip, hashlib, json
from pathlib import Path, PurePosixPath

def safe(root, name):
    p=PurePosixPath(name)
    if p.is_absolute() or '..' in p.parts or str(p)!=name or name=='01-canon/MISTAKES.md':
        raise ValueError(f'Unsafe manifest path: {name}')
    target=(root / name).resolve()
    target.relative_to(root)
    return target

def verify(data, record, name):
    if len(data)!=record['bytes'] or hashlib.sha256(data).hexdigest()!=record['sha256']:
        raise ValueError(f'Byte/SHA-256 mismatch: {name}')
    if hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()!=record['git_blob_sha1']:
        raise ValueError(f'Git blob mismatch: {name}')

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('manifests',nargs='+'); p.add_argument('--repo-root',default='.')
    p.add_argument('--verify-only',action='store_true'); args=p.parse_args(); root=Path(args.repo_root).resolve()
    records = {}; superseded = []
    for manifest_path in args.manifests:
        manifest=json.loads(Path(manifest_path).read_text())
        if manifest.get('schema') != 1:
            raise ValueError(f'Unsupported manifest schema: {manifest_path}')
        local_paths = set()
        for f in manifest['files']:
            safe(root, f['path'])
            if f['path'] in local_paths:
                raise ValueError(f"Duplicate path within manifest: {f['path']}")
            local_paths.add(f['path'])
            if f['path'] in records:
                superseded.append((f['path'], records[f['path']][0], manifest_path))
            records[f['path']] = (manifest_path, f)
    count=0; restored=0
    for _, f in records.values():
        dest=safe(root,f['path'])
        if f['storage']=='direct':
            verify(dest.read_bytes(),f,f['path'])
        elif f['storage']=='gzip_chunks':
            parts=[]
            for part in f['parts']:
                data=safe(root,part['path']).read_bytes(); verify(data,part,part['path']); parts.append(data)
            packed=b''.join(parts); verify(packed,f['compressed'],f['path']+' (compressed)')
            data=gzip.decompress(packed); verify(data,f,f['path'])
            if dest.exists(): verify(dest.read_bytes(),f,f['path'])
            elif not args.verify_only:
                dest.parent.mkdir(parents=True,exist_ok=True)
                with dest.open('xb') as out: out.write(data)
                verify(dest.read_bytes(),f,f['path']); restored+=1
        else: raise ValueError(f"Unknown storage: {f['storage']}")
        count+=1
    print(f'Verified {count} current original files; restored {restored} files; superseded {len(superseded)} earlier records.')
    for name, old, new in superseded:
        print(f'Superseded {name}: {old} -> {new}')
if __name__=='__main__': main()
