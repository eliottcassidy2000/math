"""Independent isolated lossless-package roundtrip and refusal controls."""
from pathlib import Path
import base64,gzip,hashlib,json,shutil,subprocess,sys,tempfile
R=Path(__file__).resolve().parent
FILES=['REASSEMBLY.md','REASSEMBLY.json','reassemble.py','full_records.json.gz.b64.part000']
def need(ok,msg):
 if not ok:raise RuntimeError(msg)
def digest(x):return hashlib.sha256(x).hexdigest()
def main():
 meta=json.loads((R/'REASSEMBLY.json').read_text());originals={};decoded=[]
 for g in meta['groups']:
  raw=b''.join((R/p['file']).read_bytes() for p in g['parts'])
  data=gzip.decompress(base64.b64decode(b''.join(raw.split()),validate=True))
  for item in g['outputs']:
   old=(R/item['file']).read_bytes();need(data==old,'decoded output differs from original')
   need(len(data)==item['size'] and digest(data)==item['sha256'],'decoded metadata mismatch')
   originals[item['file']]=old;decoded.append({'file':item['file'],'bytes':len(old),'sha256':digest(old)})
 checks=[]
 for opt in(False,True):
  for case in('roundtrip','corrupt_chunk','late_existing_conflict'):
   with tempfile.TemporaryDirectory() as t:
    d=Path(t)
    for name in FILES:shutil.copy2(R/name,d/name)
    if case=='corrupt_chunk':
     name=meta['groups'][0]['parts'][0]['file'];raw=bytearray((d/name).read_bytes());raw[0]^=1;(d/name).write_bytes(raw)
    if case=='late_existing_conflict':
     name=list(originals)[-1];(d/name).write_bytes(b'intentional conflicting last output\n')
    p=subprocess.run([sys.executable]+(['-O'] if opt else[])+[str(d/'reassemble.py')],cwd='/',capture_output=True)
    if case=='roundtrip':
     need(p.returncode==0,'reassembly failed')
     for name,old in originals.items():need((d/name).read_bytes()==old,'roundtrip differs')
    else:
     need(p.returncode!=0,'invalid input accepted')
     for name in originals:
      if case=='late_existing_conflict' and name==list(originals)[-1]:need((d/name).read_bytes()==b'intentional conflicting last output\n','conflict overwritten')
      else:need(not(d/name).exists(),'validation failure wrote an output')
    checks.append({'optimized':opt,'case':case,'status':'PASS'})
 return {'status':'PASS','decoded_outputs':decoded,'checks':checks,'package_pins':{name:digest((R/name).read_bytes()) for name in FILES}}
if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
