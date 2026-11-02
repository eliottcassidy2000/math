#!/usr/bin/env python3
"""Offline exact-byte identity of final source-preserving research additions."""
from pathlib import Path
import hashlib,json

def need(ok,label):
 if not ok:raise RuntimeError(label)
def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def main():
 root=Path(__file__).resolve().parent;counts={}
 def verify(lane,name,selector='files'):
  base=root/lane;m=json.loads((base/name).read_text());rows=m[selector]
  for r in rows:
   f=(base/r['path']).resolve();need(f.is_relative_to(root),'manifest leaves packet')
   need(f.is_file() and not f.is_symlink(),'missing selected file')
   need(f.stat().st_size==r['bytes'] and digest(f)==r['sha256'],'frozen payload '+lane+'/'+r['path'])
  return len(rows)
 counts['amplification_producer']=verify('path_cover_amplification','manifest.json')
 counts['amplification_independent']=verify('independent','AMPLIFICATION_AUDIT.json')
 counts['prime_independent']=verify('independent','PRIME_WITNESS_AUDIT.json','independent_files')
 prime=root/'prime80557';rows=(prime/'SHA256SUMS').read_text().splitlines();names=set()
 for line in rows:
  h,n=line.split('  ',1);f=(prime/n).resolve();need(f.is_relative_to(prime),'unsafe prime source path')
  need(digest(f)==h,'prime frozen file '+n);names.add(n)
 whitelist=set((prime/'WHITELIST.txt').read_text().splitlines());need(names|{'SHA256SUMS'}==whitelist,'prime selection closure')
 counts['prime_producer']=len(whitelist)
 a=json.loads((root/'independent/AMPLIFICATION_AUDIT.json').read_text());need(digest(root/'path_cover_amplification/manifest.json')==a['producer_manifest_sha256'],'amplification reviewed manifest')
 a=json.loads((root/'independent/PRIME_WITNESS_AUDIT.json').read_text())
 for name,h in a['producer_manifests'].items():need(digest(prime/name)==h,'prime reviewed manifest '+name)
 print(json.dumps({'status':'PASS','verified_payload_selections':counts,'scope':'Frozen final producer and independent-audit source bytes; no rerun of inherited100million prefix.'},indent=2,sort_keys=True))
if __name__=='__main__':main()
