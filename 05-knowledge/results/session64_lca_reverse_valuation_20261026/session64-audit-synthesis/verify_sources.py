"""Verify the synthesis provenance snapshot and finite receipt consistency.

Run from the original session workspace, or use --base for that workspace.
The pins use original workspace locators, not a claimed publication layout.
This script does not rerun native engines or reconstruct their proofs.
"""
import argparse
import hashlib
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent

def need(ok,message):
    if not ok: raise RuntimeError(message)

def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--base',type=Path,default=HERE.parent)
    args=parser.parse_args();base=args.base
    lines=(HERE/'SOURCE_PINS.sha256').read_text().splitlines()
    files=[]
    for line in lines:
        expected,name=line.split(maxsplit=1);p=base/name
        need(p.is_file() and digest(p)==expected,'source pin: '+name)
        files.append(p)
    manifests=[]
    for p in files:
        if p.name!='MANIFEST.json':continue
        d=json.loads(p.read_text());entries=d.get('files',d)
        for name,entry in entries.items():
            q=p.parent/name
            h=entry['sha256'] if isinstance(entry,dict) else entry
            need(digest(q)==h,'manifest hash: '+str(q))
            if isinstance(entry,dict) and 'bytes' in entry:
                need(q.stat().st_size==entry['bytes'],'manifest byte count: '+str(q))
        manifests.append({'path':str(p.relative_to(base)),'files':len(entries)})
    sums=[]
    for p in files:
        if p.name not in ('SHA256SUMS','INPUTS_SHA256SUMS'):continue
        count=0
        for line in p.read_text().splitlines():
            h,name=line.split(maxsplit=1)
            q=base/name if name.startswith('session64-') else p.parent/name
            need(digest(q)==h,'checksum manifest: '+str(q));count+=1
        sums.append({'path':str(p.relative_to(base)),'files':count})
    pairs=[]
    for p in files:
        if not p.name.endswith('normal.json'):continue
        q=p.with_name(p.name.replace('normal.json','optimized.json'))
        if q not in files:continue
        need(p.read_bytes()==q.read_bytes(),'normal/optimized mismatch: '+str(p))
        d=json.loads(p.read_text())
        if isinstance(d,dict) and 'status' in d:
            expected='PASS / DIAGNOSTIC ONLY' if p.parent.name=='session64-audit-bridge' and p.name=='n10.normal.json' else 'PASS'
            need(d['status']==expected,'receipt unexpected status: '+str(p))
        pairs.append(str(p.relative_to(base)))
    a=(HERE/'lca-normal.json').read_bytes();b=(HERE/'lca-optimized.json').read_bytes()
    need(a==b==(base/'session64-root/newton-lca/normal.json').read_bytes(),'fresh LCA replay')
    need(json.loads(a)['status']=='PASS','fresh LCA status')
    root={name:digest(base/'session64-root'/name) for name in ['SYNTHESIS.md','CONNECTIONS.md']}
    print(json.dumps({'status':'PASS','source_files':len(files),'root_sha256':root,
      'checked_manifests':manifests,'checked_checksum_lists':sums,
      'normal_optimized_pairs':pairs,'fresh_lca_replay':json.loads(a),
      'scope':'Source-seal and synthesis consistency audit; no claim to rerun every long native engine.'},indent=2,sort_keys=True))
if __name__=='__main__':main()
