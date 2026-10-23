from pathlib import Path
import json,hashlib
from math import isqrt
D=Path(__file__).resolve().parent
S=next(p for p in (D.parent/'color-holdouts',D.parent/'session42-color-holdouts') if p.exists())
expected='6b476fd57a203f90fc1aefcd51acd043fd549dc689630d15a513711a05ab048f'
m=S/'MANIFEST.json'
if hashlib.sha256(m.read_bytes()).hexdigest()!=expected:raise RuntimeError('manifest changed')
for f in json.loads(m.read_text())['files']:
 b=(S/f['path']).read_bytes()
 if len(b)!=f['bytes'] or hashlib.sha256(b).hexdigest()!=f['sha256']:raise RuntimeError(f['path'])
F=[1,2]
while F[-1]<80:F.append(F[-1]+F[-2])
def vf(n):
 r=n;digits=[]
 for i in range(len(F)-1,-1,-1):
  if F[i]<=r:r-=F[i];digits.append(i+2)
 return (n,n.bit_count(),isqrt(n),len(digits),min(digits),sum(digits))
def ef(a,b):
 a,b=sorted((a,b));A=vf(a);B=vf(b)
 return (a,b,isqrt(a+b),A[1]+B[1],A[2]+B[2],A[3]+B[3],A[4]+B[4])
rows={n:set(map(tuple,json.loads((S/f'paths{n}.json').read_text()))) for n in range(33,37)}
W=json.loads((S/'witnesses.json').read_text())
for w in W:
 p,q=map(tuple,w['paths']);n=w['first_failure'];mod=w['q'];co=w['coefficients']
 if p==q or p not in rows[n] or q not in rows[n] or (p[0],p[-1])!=(q[0],q[-1]):raise RuntimeError('expanded carrier mismatch')
 def label(p):
  fs=map(vf,p) if w['kind']=='vertex' else (ef(a,b) for a,b in zip(p,p[1:]))
  return tuple(sum(a*b for a,b in zip(co,f))%mod for f in fs)
 if label(p)!=label(q):raise RuntimeError('expanded collision mismatch')
print(json.dumps({'status':'PASS','source_manifest_sha256':expected,'source_files':17,'expanded_witness_pairs':len(W)},indent=2))
