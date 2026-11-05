from pathlib import Path
from tempfile import TemporaryDirectory
from math import prod
import hashlib,json,subprocess
P=Path(__file__).resolve().parent
expected={'ROOTED.md':'ba7e632518219eeb8d86b6db7da386a88e98317ea417c9b22e649278c132b32e','WEAK_START_ONE_ORACLE.md':'55048003e63b68ef21a51d424117145a4afaeb4bdfb43b5946d19524b6006ec0'}
for name,digest in expected.items():
 if hashlib.sha256((P/name).read_bytes()).hexdigest()!=digest:raise RuntimeError(('dependency',name))
with TemporaryDirectory(prefix='weak-ambiguity-')as tmp:
 out=[]
 for label,flags in [('optimized',['-O3']),('debug',['-O0','-D_GLIBCXX_ASSERTIONS'])]:
  binary=Path(tmp)/label
  subprocess.run(['g++','-std=c++17',*flags,str(P/'verify_weak_ambiguity.cpp'),'-o',str(binary)],check=True)
  out.append(subprocess.check_output([str(binary)]))
 if out[0]!=out[1]:raise RuntimeError('independent build modes differ')
 physical=json.loads(out[0])
controls=0;final={}
for m in range(1,51):
 total=prod((1<<j)+1 for j in range(1,m+1));deletions=set()
 for last in range(1,m+1):
  order=[j for j in range(1,m+1)if j!=last]+[last]
  h,e,R,d,n=1,1,1,1,1
  for j in order:
   L=1<<j
   h,e,R,d,n=h+L*R,h,h+L*R,h,n+j+1
  if(h,R,e,d,n+1)!=(total,total,total//((1<<last)+1),total//((1<<last)+1),2+m*(m+3)//2):raise RuntimeError(('exact recipe',m,last))
  deletions.add(d);controls+=1
 if len(deletions)!=m:raise RuntimeError('deletion distinctness')
 if m==50:final=dict(m=m,order=n+1,H=total,distinct_deletions=len(deletions))
result=dict(status='PASS',physical=physical,exact_large_recipe_controls=controls,largest_symbolic_control=final,dependencies=expected)
text=json.dumps(result,sort_keys=True,indent=2)+'\n';(P/'weak_ambiguity_certificate.json').write_text(text);print(text,end='')
