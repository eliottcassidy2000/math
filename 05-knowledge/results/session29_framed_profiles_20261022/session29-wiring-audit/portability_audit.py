from pathlib import Path
import hashlib,json,tempfile,shutil,subprocess
D=Path(__file__).parent;P=D.parent/'session29-profiles';f=P/'check_compression.py';new=f.read_text()
old=new.replace("from pathlib import Path; x=json.load(open(Path(__file__).parent/'transfer-results.json'))", "x=json.load(open('/workspace/shared/session29-profiles/transfer-results.json'))").replace("open(Path(__file__).parent/'compressed-matrices.json','w')", "open('/workspace/shared/session29-profiles/compressed-matrices.json','w')")
h=lambda b:hashlib.sha256(b).hexdigest()
def need(x,msg):
 if not x:raise RuntimeError(msg)
need(h(old.encode())=='4c18cba563e34cd268c45939e46c1bbfc6bd65d2b297e5ab7f187f60372c5760','exact old-byte reconstruction')
need(h(f.read_bytes())=='792d8b451e98b62b1b520d029f3dc06324ad9fc12f1016e0e34b2bcef46b51f1','new hash')
need(sum(a!=b for a,b in zip(old.splitlines(),new.splitlines()))==2,'two modified lines')
for line in (D/'SOURCE_SHA256SUMS').read_text().splitlines():
 oldhash,name=line.split();need(name=='check_compression.py' or h((P/name).read_bytes())==oldhash,'unchanged '+name)
with tempfile.TemporaryDirectory(prefix='portable-',dir=D) as td:
 t=Path(td)
 for name in ('check_compression.py','transfer-results.json'):shutil.copy2(P/name,t/name)
 runs=[]
 for flags in ([],['-O']):
  p=subprocess.run(['python',*flags,str(t/'check_compression.py')],cwd='/tmp',capture_output=True,text=True,check=True)
  need((t/'compressed-matrices.json').read_bytes()==(P/'compressed-matrices.json').read_bytes(),'output byte equality');runs.append(p.stdout)
 need(runs[0]==runs[1],'normal optimized equality')
 (D/'portability.out').write_text('PASS exact two-line path-only amendment; all other reviewed sources unchanged.\nPASS isolated scratch copy from /tmp; normal and optimized results byte-identical.\n'+runs[0])
(D/'SOURCE_SHA256SUMS').write_text('\n'.join(h((P/line.split()[1]).read_bytes())+'  '+line.split()[1] for line in (D/'SOURCE_SHA256SUMS').read_text().splitlines())+'\n')
print((D/'portability.out').read_text())
