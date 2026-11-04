"""Frozen-source integrity and normal/optimized producer verifier replay in isolation."""
from pathlib import Path
import tempfile,subprocess,hashlib,json,sys
D=Path(__file__).resolve().parent;P=D.parent/'skew_polynomial'
if not P.exists():P=D.parent/'session189-skew-polynomial'
def need(c,s):
 if not c:raise RuntimeError(s)
def main():
 release=json.loads((P/'RELEASE_FILES.json').read_text());pins={}
 with tempfile.TemporaryDirectory() as tmp:
  root=Path(tmp)
  for name,rec in release.items():
   source=(P/name).read_bytes();need(len(source)==rec['bytes'] and hashlib.sha256(source).hexdigest()==rec['sha256'],'every frozen source/artifact byte matches explicit allowlist');(root/name).write_bytes(source);pins[name]=rec['sha256']
  before={n:hashlib.sha256((root/n).read_bytes()).hexdigest() for n in release}
  normal=subprocess.check_output([sys.executable,str(root/'verify.py')],cwd=root);optimized=subprocess.check_output([sys.executable,'-O',str(root/'verify.py')],cwd=root);need(normal==optimized==(P/'receipt.json').read_bytes()==(P/'receipt.optimized.json').read_bytes(),'isolated producer normal/optimized verifier agrees exactly');need(before=={n:hashlib.sha256((root/n).read_bytes()).hexdigest() for n in release},'producer replay leaves retained artifacts unchanged')
 txt=json.dumps(pins,sort_keys=True,indent=2)+'\n';p=D/'producer_source_pins.json'
 if p.exists():need(p.read_text()==txt,'frozen source pins stable')
 else:p.write_text(txt)
 return dict(status='PASS',allowlisted_files=len(release),producer_normal_optimized_receipt_sha256=hashlib.sha256(normal).hexdigest(),replay_source_sha256=hashlib.sha256((D/'replay_audit.py').read_bytes()).hexdigest(),source_pins_sha256=hashlib.sha256(p.read_bytes()).hexdigest(),scope='The documented C++ computations are independently recounted by the other audit programs; this isolated replay additionally checks the exact frozen Python verifier and every source/artifact byte.')
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
