#!/usr/bin/env python3
"""Offline exact-byte identity of second-checkpoint producer/audit selections."""
from pathlib import Path
import hashlib,json

def need(ok,label):
 if not ok:raise RuntimeError(label)
def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def main():
 root=Path(__file__).resolve().parent;counts={}
 def verify(lane,name):
  base=root/lane;m=json.loads((base/name).read_text());rows=m['files']
  if isinstance(rows,dict):rows=[dict(path=k,**v) for k,v in rows.items()]
  for r in rows:
   f=(base/r['path']).resolve();need(f.is_relative_to(root),'manifest leaves packet')
   need(f.is_file() and not f.is_symlink(),'missing selected file')
   need(f.stat().st_size==r['bytes'] and digest(f)==r['sha256'],'frozen payload '+lane+'/'+r['path'])
  return len(rows)
 for lane,name in [('minimal_gap','manifest.json'),('constructions','manifest.json'),('independent','MINIMAL_GAP_AUDIT.json'),('independent','CONSTRUCTIONS_AUDIT.json')]:
  counts[lane+'/'+name]=verify(lane,name)
 for audit,lane in [('MINIMAL_GAP_AUDIT.json','minimal_gap'),('CONSTRUCTIONS_AUDIT.json','constructions')]:
  data=json.loads((root/'independent'/audit).read_text())
  need(digest(root/lane/'manifest.json')==data['producer_manifest_sha256'],'reviewed producer manifest '+lane)
 print(json.dumps({'status':'PASS','verified_payload_selections':counts,'scope':'Frozen second-checkpoint producer and independent-audit bytes only; all-size claims use reviewed written proofs.'},indent=2,sort_keys=True))
if __name__=='__main__':main()
