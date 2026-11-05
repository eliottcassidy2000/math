#!/usr/bin/env python3
"""Exact chamber recovery and recursive inverse checks; standard library only."""
from collections import deque
from functools import lru_cache
from pathlib import Path
import hashlib,json
T=(1,0,2,1)
def need(x,msg):
 if not x:raise RuntimeError(msg)
def step(s,c):
 D,E,U,d=s
 if c=='O':return(D,0,2*U,d)
 if c=='Z':return(D-E+U,d,2*U-d,d)
 if c=='J':return(D-E+U,D-E,D+U,D)
 raise ValueError(c)
def word(w):
 s=T
 for c in w:s=step(s,c)
 return s
def endpoint(D,U,d):
 if U>=2*D:return 0
 if U>D+d:return d
 return D+d-U
def lift(x):
 D,U,d=x
 return D,endpoint(*x),U,d
def project(s):return(s[0],s[2],s[3])
def reduced(x,c):
 D,U,d=x;E=endpoint(*x)
 if c=='O':return(D,2*U,d)
 if c=='Z':return(D-E+U,2*U-d,d)
 return(D-E+U,D+U,D)
def guarded(x):
 D,U,d=x
 if min(D,d)<1 or D%2!=1 or d%2!=1 or U<2:return False
 E=endpoint(*x)
 return 0<=E<=d and D-E>=d and U>=D+d-E and U>=2*d
@lru_cache(None)
def inverse(x):
 if not guarded(x):return None
 if x==(1,2,1):return ''
 D,U,d=x;E=endpoint(*x);f=U-D+E-d;p=d-E
 candidates=[]
 if E==0:
  if U%2:return None
  candidates=[((D,U//2,d),'O')]
 elif U==D:
  if(D,d)==(3,1):return 'Z'
  uw=(D-d)//2
  if d>D//3:return None
  candidates=[((d,uw,dw),'OJ')for dw in range(1,min(d,uw//2)+1,2)]
 elif f>0:
  if E!=d or f<=d or(U+d)%2:return None
  u=(U+d)//2;K=D-u
  candidates=[((K if K%2 else K+d,u,d),'Z')]
 elif f==0:
  if E<=p or p<=0:return None
  oldd=d-E if D>2*d else U-d-E
  candidates=[((d,U-d,oldd),'J')]
 else:return None
 for old,suffix in candidates:
  if old[1]>=U:continue
  prev=inverse(old)
  if prev is None:continue
  state=lift(old)
  for c in suffix:state=step(state,c)
  if state==lift(x):return prev+suffix
 return None

def main():
 cap=2000;bank={T:''};queue=deque([T])
 while queue:
  s=queue.popleft()
  for c in'OZJ':
   z=step(s,c)
   need(z[2]>s[2],'strict U growth')
   if z[2]<=cap and z not in bank:bank[z]=bank[s]+c;queue.append(z)
 triples={project(s)for s in bank}
 reduced_bank={(1,2,1)};queue=deque(reduced_bank)
 while queue:
  x=queue.popleft()
  for c in'OZJ':
   y=reduced(x,c)
   if y[1]<=cap and y not in reduced_bank:reduced_bank.add(y);queue.append(y)
 need(triples==reduced_bank,'reduced orbit identity')
 chamber={'weak':0,'strong_f0':0,'strong_fp':0};balanced=0
 for s,w in bank.items():
  D,E,U,d=s;f=U-D+E-d;p=d-E
  need(lift(project(s))==s,'endpoint recovery')
  need(f==0 or f>d,'strict f gap')
  need(E==0 or p==0 or f==0,'boundary union')
  need((E==0)==(U>=2*D),'weak chamber')
  if E>0 and f==0:need(E>p,'strong f0 separator')
  if E>0:
   need(E!=U-D,'strong deletion strict separator')
   need(d==(E if E<U-D else U-D+E),'strong deletion recovery')
  chamber['weak'if E==0 else'strong_f0'if f==0 else'strong_fp']+=1
  balanced+=int(U==D)
  iw=inverse(project(s))
  need(iw is not None and word(iw)==s,'physical inverse certificate')
  for c in'OZJ':need(lift(reduced(project(s),c))==step(s,c),'transition recovery')
 # Exhaustive triple-box comparison with a complete physical bank.
 box=0;accepted=0;guarded_count=0
 for D in range(1,82,2):
  for U in range(2,121):
   for d in range(1,min(D,U//2)+1,2):
    x=(D,U,d);box+=1;guarded_count+=guarded(x)
    iw=inverse(x);actual=x in triples
    need((iw is not None)==actual,'box membership')
    if iw is not None:accepted+=1;need(word(iw)==lift(x),'box certificate')
 controls=[(1,3,1),(13,15,7)]
 for x in controls:need(inverse(x)is None,'nonphysical control')
 witness=[word('OOJO'),word('ZOJO')]
 need(witness==[(9,0,18,1),(9,0,18,3)],'deletion witness')
 # The two strong J chamber maps have determinant +1.
 mats=[((1,1,-1),(1,1,0),(1,0,0)),((0,2,-1),(1,1,0),(1,0,0))]
 def det(m):
  a,b,c=m
  return a[0]*(b[1]*c[2]-b[2]*c[1])-a[1]*(b[0]*c[2]-b[2]*c[0])+a[2]*(b[0]*c[1]-b[1]*c[0])
 need([det(m)for m in mats]==[1,1],'unimodular J chambers')
 out={'status':'PASS','cap_U':cap,'physical_profiles':len(bank),'projected_triples':len(triples),'chambers':chamber,'balanced_profiles':balanced,'triple_box':{'D_max':81,'U_max':120,'d_max':'min(D,floor(U/2))','tested':box,'guarded':guarded_count,'physical':accepted},'nonphysical_controls':controls,'lost_deletion_witness':witness,'inverse_cache':inverse.cache_info()._asdict(),'scope':'Exact physical language and deterministic nonbalanced inverse, not static chamber realizability or all-target coverage.','proof_sha256':hashlib.sha256(Path(__file__).with_name('PROOF.md').read_bytes()).hexdigest(),'verifier_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
 data=json.dumps(out,indent=2,sort_keys=True)+'\n';Path(__file__).with_name('RESULT.json').write_text(data);print(data,end='')
if __name__=='__main__':main()
