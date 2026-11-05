from pathlib import Path
from collections import deque
from itertools import product
from math import gcd
import json,hashlib,sys,runpy,contextlib,io
P=Path(__file__).resolve().parent

def constants(L):
 word=bin(L-1)[2:]if L>1 else'';n=0;g=1
 for b in word:
  if b=='0':g+=n+(n&1)
  n=2*n+int(b)
 Q=1+int(''.join('1'if b=='0'else'0'for b in word)[::-1]or'0',2)
 return g,Q,1<<len(word)

def bank(cap):
 coeff=[None]+[constants(L)for L in range(1,cap+1)]
 S={(1,1,1,1)};queue=deque(S)
 while queue:
  D,E,U,d=queue.popleft();nexts=[(D,0,2*U,d)]
  for L in range(1,cap//U+1):
   g,q,z=coeff[L];nexts.append((g*D+L*U-E*(L&1),D-E*(L==1),q*D+z*U,D))
  for t in nexts:
   if t[2]<=cap and t not in S:S.add(t);queue.append(t)
 return S

def sig(s):
 D,E,U,d=s;return D+2*d,E+d,U-d

def N(v,w):
 D,E,U,d=v;F,G,V,f=w
 return (D-E)*(F-G)+D*V+U*F+3*U*V-d*(G+V)-f*(E+U)
def M(v,w):
 D,E,U,d=v;F,G,V,f=w
 return D*F+2*(D*V+U*F)+6*U*V-E*V-U*G-d*(G+V)-f*(E+U)

sys.path.insert(0,str(P/'source'))
from signature import recover,member
S=bank(2000);image={}
for s in S:
 a=sig(s)
 if a in image or recover(*a)!=s:raise RuntimeError('injective physical recovery')
 image[a]=s
 if gcd(*s)!=gcd(*a):raise RuntimeError('content equality')
# Independent direct O/M closure equals author's cores-then-tails inventory.
with contextlib.redirect_stdout(io.StringIO()):author=runpy.run_path(str(P/'source/verify.py'))
if S!=author['S']:raise RuntimeError('complete physical state set')
small=bank(31);expected={sig(s):s for s in small};cube=0
for A,B,C in product(range(1,94),range(1,63),range(32)):
 got=member(A,B,C)
 if got is not None and got[0][2]>31:got=None
 if(bool(got)!=((A,B,C)in expected)):raise RuntimeError('candidate cube membership')
 if got and got[0]!=expected[A,B,C]:raise RuntimeError('candidate response')
 cube+=1
literal=0
for line in(P/'literal-profiles.csv').read_text().splitlines():
 s=tuple(map(int,line.split(',')))
 if s not in S or recover(*sig(s))!=s:raise RuntimeError('literal graph source recovery')
 literal+=1
pairs=0
small=bank(41)-{(1,1,1,1)}
for v,w in product(small,repeat=2):
 D,E,U,d=v;F,G,V,f=w;a,b,c=sig(v);x,y,z=sig(w);n=N(v,w);m=M(v,w)
 if n!=(a-b)*(x-y)+a*z+c*x+3*c*z:raise RuntimeError('unimodular N quotient')
 if 4*n<15*D*F+5*D+5*F-1 or 2*m<17*D*F+3*D+3*F-1:raise RuntimeError('additive product bound')
 pairs+=1
# Independently check slack identities throughout a larger relaxed cone,
# not merely on physical samples. This controls every coefficient sign.
cone=0
states=[]
for e,r,a,b in product(range(4),repeat=4):
 d=e+r;D=2*e+r+a;U=2*e+2*r+a+b
 if d and D:states.append((D,e,U,d))
for v,w in product(states,repeat=2):
 D,E,U,d=v;F,G,V,f=w;b=U-D;c=V-F;sl=U-D+E-d;sr=V-F+G-f
 n0=6*D*F-2*(D*G+E*F)-E*G;m0=11*D*F-2*(D*G+E*F)-2*E*G
 if N(v,w)!=n0+b*(3*F-2*G)+c*(3*D-2*E)+b*c+sl*(G+V)+sr*(E+U):raise RuntimeError('N slack algebra')
 if M(v,w)!=m0+b*(7*F-3*G)+c*(7*D-3*E)+4*b*c+sl*(G+V)+sr*(E+U):raise RuntimeError('M slack algebra')
 cone+=1
sha=hashlib.sha256(json.dumps(sorted(S),separators=(',',':')).encode()).hexdigest()
print(json.dumps(dict(status='PASS',physical_quadruples=len(S),complete_state_sha256=sha,integer_cube=cube,literal_sources=literal,physical_pairs=pairs,relaxed_cone_slack_pairs=cone),sort_keys=True,indent=2))
