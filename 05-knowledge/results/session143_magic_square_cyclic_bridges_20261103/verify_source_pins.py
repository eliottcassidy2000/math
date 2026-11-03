#!/usr/bin/env python3
"""Verify frozen Git-source snapshots offline, without opening a repository."""
import argparse,hashlib,json
from pathlib import Path

def need(ok,label):
    if not ok:raise RuntimeError(label)
def main():
    p=argparse.ArgumentParser();p.add_argument('--manifest',default='SOURCE_PINS_CHECKPOINT1.json');a=p.parse_args()
    root=Path(__file__).resolve().parent
    rows=json.loads((root/a.manifest).read_text())['sources']
    identities=set(); snapshots=set()
    for row in rows:
        rel=Path(row['snapshot']);need(not rel.is_absolute() and '..' not in rel.parts,'unsafe source path')
        f=root/rel;need(f.is_file() and not f.is_symlink(),'missing source snapshot')
        data=f.read_bytes();need(len(data)==row['bytes'],'source bytes '+str(rel))
        need(hashlib.sha256(data).hexdigest()==row['sha256'],'source SHA256 '+str(rel))
        need(hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()==row['git_blob'],'source Git blob '+str(rel))
        identities.add((row['source_commit'],row['source_path'],row['git_blob']));snapshots.add(str(rel))
    print(json.dumps({'status':'PASS','source_identities':len(identities),'source_snapshots':len(snapshots),'scope':'Offline Git blob and SHA256 identity of selected historical source bytes; no external formal build or mathematical status promotion.'},indent=2,sort_keys=True))
if __name__=='__main__':main()
