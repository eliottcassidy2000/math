"""Lossless prime-root family package audit, isolated safe-reconstruction controls."""
from pathlib import Path
import base64,gzip,json,hashlib,tempfile,subprocess,sys,shutil
HERE=Path(__file__).resolve().parent;SRC=HERE.parent/'prime_pair_root_family'
def need(ok,msg):
 if not ok:raise RuntimeError(msg)
def main():
 names=['REASSEMBLY.md','REASSEMBLY.json','reassemble.py','full_records.json.gz.b64.part000'];pins={n:hashlib.sha256((SRC/n).read_bytes()).hexdigest() for n in names};meta=json.loads((SRC/'REASSEMBLY.json').read_text());need(len(meta['groups'])==1,'single complete record group');g=meta['groups'][0];originals={n:(SRC/n).read_bytes() for n in ('normal.json','optimized.json')};raw=b''.join(b''.join((SRC/p['file']).read_bytes().split()) for p in g['parts']);decoded=gzip.decompress(base64.b64decode(raw,validate=True));need(all(decoded==v for v in originals.values()),'full byte-for-byte original equality');need(len(decoded)==272915 and len(json.loads(decoded)['cases'])==6,'complete declared records')
 trials=[]
 for optimized in(False,True):
  for mode in ('fresh','idempotent','corrupt_chunk','late_existing_conflict','late_hash_conflict','duplicate_output','nonlocal_output'):
   with tempfile.TemporaryDirectory() as td:
    d=Path(td)
    for n in names:shutil.copyfile(SRC/n,d/n)
    if mode=='idempotent':
     for n,v in originals.items():(d/n).write_bytes(v)
    if mode=='corrupt_chunk':
     p=d/g['parts'][0]['file'];v=bytearray(p.read_bytes());v[0]=ord('A') if v[0]!=ord('A') else ord('B');p.write_bytes(v)
    if mode=='late_existing_conflict':(d/'optimized.json').write_bytes(b'preexisting conflict')
    if mode in ('late_hash_conflict','duplicate_output','nonlocal_output'):
     m=json.loads((d/'REASSEMBLY.json').read_text());second=m['groups'][0]['outputs'][1]
     if mode=='late_hash_conflict':second['sha256']='0'*64
     elif mode=='duplicate_output':second['file']='normal.json'
     else:second['file']='../escaped.json'
     (d/'REASSEMBLY.json').write_text(json.dumps(m))
    cmd=[sys.executable]+(['-O'] if optimized else[])+[str(d/'reassemble.py')];proc=subprocess.run(cmd,cwd='/',capture_output=True)
    if mode in ('fresh','idempotent'):
     need(proc.returncode==0 and json.loads(proc.stdout)['status']=='PASS','success result')
     for n,v in originals.items():need((d/n).read_bytes()==v,'complete reconstructed bytes')
    else:
     need(proc.returncode!=0,'explicit corruption/conflict failure');need(not(d/'normal.json').exists(),'no early normal output')
     if mode=='late_existing_conflict':need((d/'optimized.json').read_bytes()==b'preexisting conflict','conflicting late output untouched')
     else:need(not(d/'optimized.json').exists(),'no late output')
    trials.append({'optimized':optimized,'mode':mode,'passed':True})
 need(pins=={n:hashlib.sha256((SRC/n).read_bytes()).hexdigest() for n in names},'package source freeze')
 print(json.dumps({'status':'PASS','complete_originals':[{'file':n,'size':len(v),'sha256':hashlib.sha256(v).hexdigest()} for n,v in originals.items()],'isolated_controls':trials,'package_sha256':pins,'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2,sort_keys=True))
if __name__=='__main__':main()
