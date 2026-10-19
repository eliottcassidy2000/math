#!/usr/bin/env python3
"""Reconstruct exact approved gzip bytes from authenticated ordered raw parts."""
import argparse, hashlib, json, os, tempfile
from pathlib import Path
ROOT=Path(__file__).resolve().parent

def check(ok,message):
    if not ok: raise RuntimeError(message)
def local(name):
    p=Path(name)
    check(not p.is_absolute() and '..' not in p.parts, 'Unsafe relative path')
    q=ROOT/p
    check(q.resolve().is_relative_to(ROOT), 'Path escapes package')
    check(not q.is_symlink(), 'Symlink is not an approved artifact')
    return q
def digest(path):
    h=hashlib.sha256();size=0
    with path.open('rb') as f:
        for chunk in iter(lambda:f.read(1024*1024),b''):h.update(chunk);size+=len(chunk)
    return size,h.hexdigest()
def reconstruct(item,check_only=False):
    target=local(item['original_path']);parts=item['parts'];combined=hashlib.sha256();total=0
    check([p['order'] for p in parts]==list(range(1,len(parts)+1)), 'Part order is not consecutive')
    for part in parts:
        path=local(part['path']);size,h=digest(path)
        check(size==part['bytes'] and h==part['sha256'], 'Part mismatch: '+part['path'])
        check(size<=8*1024*1024, 'Part exceeds transport limit')
        with path.open('rb') as f:
            for chunk in iter(lambda:f.read(1024*1024),b''):combined.update(chunk);total+=len(chunk)
    check(total==item['bytes'] and combined.hexdigest()==item['sha256'], 'Combined original digest mismatch')
    if target.exists():
        check(target.is_file() and digest(target)==(item['bytes'],item['sha256']), 'Refusing to overwrite mismatching original: '+item['original_path'])
        return dict(path=item['original_path'],status='ALREADY_MATCHING',bytes=total,sha256=combined.hexdigest())
    if check_only:return dict(path=item['original_path'],status='PARTS_VERIFIED',bytes=total,sha256=combined.hexdigest())
    target.parent.mkdir(parents=True,exist_ok=True)
    fd,tmp=tempfile.mkstemp(prefix='.reconstruct-',dir=target.parent)
    try:
        h=hashlib.sha256();written=0
        with os.fdopen(fd,'wb') as out:
            for part in parts:
                with local(part['path']).open('rb') as f:
                    for chunk in iter(lambda:f.read(1024*1024),b''):out.write(chunk);h.update(chunk);written+=len(chunk)
            out.flush();os.fsync(out.fileno())
        check(written==item['bytes'] and h.hexdigest()==item['sha256'], 'Parts changed during reconstruction')
        try:os.link(tmp,target)
        except FileExistsError:
            check(not target.is_symlink() and target.is_file() and digest(target)==(item['bytes'],item['sha256']), 'Refusing concurrent mismatching original')
        check(digest(target)==(item['bytes'],item['sha256']), 'Reconstructed file mismatch')
    finally:
        Path(tmp).unlink(missing_ok=True)
    return dict(path=item['original_path'],status='RECONSTRUCTED',bytes=total,sha256=combined.hexdigest())
def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--check-only',action='store_true');args=parser.parse_args()
    manifest=json.loads((ROOT/'TRANSPORT_MANIFEST.json').read_text());check(manifest['format_version']==1,'Unsupported transport manifest')
    print(json.dumps({'status':'PASS','artifacts':[reconstruct(x,args.check_only) for x in manifest['artifacts']]},indent=2))
if __name__=='__main__':main()
