"""Independent isolated review of lossless modulo-eight packaging."""
from pathlib import Path
import base64,gzip,hashlib,json,shutil,subprocess,sys,tempfile
ROOT=Path(__file__).resolve().parent;SRC=ROOT.parent/'mod8_transfer'
def need(ok,msg):
 if not ok:raise RuntimeError(msg)
meta=json.loads((SRC/'REASSEMBLY.json').read_text());items=[x for g in meta['groups'] for x in g['parts']];outputs=[x for g in meta['groups'] for x in g['outputs']]
for group in meta['groups']:
 raw=[]
 for x in group['parts']:
  b=(SRC/x['file']).read_bytes();need(len(b)==x['size'] and hashlib.sha256(b).hexdigest()==x['sha256'],'raw chunk receipt');raw.append(b)
 data=gzip.decompress(base64.b64decode(b''.join(raw)))
 for x in group['outputs']:need(data==(SRC/x['file']).read_bytes(),'independent decoded output bytes')
def setup(root):
 for p in ['reassemble.py','REASSEMBLY.json']+[x['file'] for x in items]:shutil.copyfile(SRC/p,root/p)
def invoke(root,opt):return subprocess.run([sys.executable]+(['-O'] if opt else [])+[str(root/'reassemble.py')],capture_output=True)
records=[]
for opt in (False,True):
 with tempfile.TemporaryDirectory() as name:
  root=Path(name);setup(root);p=invoke(root,opt);need(p.returncode==0,'fresh isolated reassembly')
  for x in outputs:need((root/x['file']).read_bytes()==(SRC/x['file']).read_bytes(),'normal optimized bytes')
  need(invoke(root,opt).returncode==0,'identical output idempotence');records.append({'mode':'optimized' if opt else 'normal','roundtrip':True,'idempotent':True})
 with tempfile.TemporaryDirectory() as name:
  root=Path(name);setup(root);target=root/items[-1]['file'];b=bytearray(target.read_bytes());b[0]=ord('B') if b[0]!=ord('B') else ord('C');target.write_bytes(b)
  p=invoke(root,opt);need(p.returncode!=0 and all(not (root/x['file']).exists() for x in outputs),'late corruption must cause zero writes')
 with tempfile.TemporaryDirectory() as name:
  root=Path(name);setup(root);target=root/outputs[-1]['file'];target.write_bytes(b'CONFLICT');p=invoke(root,opt)
  need(p.returncode!=0 and target.read_bytes()==b'CONFLICT' and all(not (root/x['file']).exists() for x in outputs[:-1]),'conflicting late output must cause zero writes')
print(json.dumps({'status':'PASS','groups':len(meta['groups']),'chunks':len(items),'outputs':outputs,'restored_total_bytes':sum(x['size'] for x in outputs),'isolated_modes':records,'late_corrupt_chunk_zero_writes':True,'late_conflicting_output_zero_writes':True,'source_sha256':{p:hashlib.sha256((SRC/p).read_bytes()).hexdigest() for p in ['REASSEMBLY.json','REASSEMBLY.md','reassemble.py']+[x['file'] for x in items]}},sort_keys=True,indent=2))
