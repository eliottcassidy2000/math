#!/usr/bin/env python3
from pathlib import Path
import tempfile,shutil,subprocess,sys,json,hashlib
P=Path(__file__).resolve().parent

def need(b,s):
 if not b:raise RuntimeError(s)
def main():
 src=Path(sys.argv[1]).resolve()if len(sys.argv)>1 else P.parent/'session189-kernel-span';base=Path(sys.argv[2]).resolve()if len(sys.argv)>2 else P.parent/'session189-odd-transitive';files=json.load(open(src/'RELEASE_FILES.json'));pins={f:hashlib.sha256((src/f).read_bytes()).hexdigest()for f in files};inherited={f:hashlib.sha256((base/f).read_bytes()).hexdigest()for f in ['verify.py','PROOF.md']}
 with tempfile.TemporaryDirectory(prefix='kernel-replay-')as d:
  D=Path(d);(D/'session189-kernel-span').mkdir();(D/'session189-odd-transitive').mkdir()
  for f in files:shutil.copyfile(src/f,D/'session189-kernel-span'/f)
  for f in ['verify.py','PROOF.md']:shutil.copyfile(base/f,D/'session189-odd-transitive'/f)
  for mode,opt in [('normal',[]),('optimized',['-O'])]:
   out=subprocess.check_output([sys.executable,*opt,str(D/'session189-kernel-span'/'verify.py')]);(P/f'producer.{mode}.stdout').write_bytes(out);r=(D/'session189-kernel-span'/'receipt.json').read_bytes();(P/f'producer.{mode}.json').write_bytes(r);need(r==(src/'receipt.json').read_bytes(),'producer receipt mismatch')
  need((P/'producer.normal.stdout').read_bytes()==(P/'producer.optimized.stdout').read_bytes(),'producer stdout modes')
 need(pins=={f:hashlib.sha256((src/f).read_bytes()).hexdigest()for f in files},'producer mutation');need(inherited=={f:hashlib.sha256((base/f).read_bytes()).hexdigest()for f in inherited},'dependency mutation')
 receipt={'status':'PASS','producer':pins,'inherited':inherited,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'artifacts':{f:hashlib.sha256((P/f).read_bytes()).hexdigest()for f in ['producer.normal.stdout','producer.optimized.stdout','producer.normal.json','producer.optimized.json']}};(P/'PRODUCER_PINS.json').write_text(json.dumps(receipt,sort_keys=True,indent=2)+'\n');print(json.dumps(receipt,sort_keys=True))
if __name__=='__main__':main()
