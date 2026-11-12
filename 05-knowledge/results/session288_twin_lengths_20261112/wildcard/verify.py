"""Exact finite skeleton proof certificates and independent literal checks.
No third-party packages required. Run from any working directory.
"""
from pathlib import Path
from itertools import product,combinations
from collections import defaultdict
import json,math
from kernel import skeleton,observable_skeleton,exp,evaluate
from common import counts,expand
HERE=Path(__file__).resolve().parent
ROWS=[6,4,8,19,39,15];O=[0,4];PORTS=[1,5]
EXPECTED={
'H':{(1,1):15,(2,1):6,(1,2):5,(2,2):2},
'W':{(1,1):14,(2,1):5,(1,2):6,(2,2):2},
'U':{(1,1):86,(2,1):72,(3,1):18,(1,2):42,(2,2):34,(3,2):8}}
def need(ok,msg):
 if not ok:raise RuntimeError(msg)
def encode(d):return [[*k,v] for k,v in sorted(d.items())]
def lift_modules(A,B):
 blocks=[[0],[1]*len(A),[2],[3],[4],[5]*len(B)];labels=sum(blocks,[])
 pos=[0]*6;within=[]
 for v in labels:within.append(pos[v]);pos[v]+=1
 out=[]
 for i,v in enumerate(labels):
  r=0
  for j,w in enumerate(labels):
   if i==j:continue
   edge= (A[within[i]]>>within[j]&1) if v==w==1 else ((B[within[i]]>>within[j]&1) if v==w==5 else (ROWS[v]>>w&1))
   if edge:r|=1<<j
  out.append(r)
 return out,[i for i,v in enumerate(labels) if v in O]
def profiles(rows):
 n=len(rows);full=(1<<n)-1;end=[[0]*n for _ in range(1<<n)];h=[0]*(1<<n)
 for s in range(1,1<<n):
  for v in range(n):
   if s>>v&1:
    t=s^(1<<v);end[s][v]=1 if not t else sum(end[t][u] for u in range(n) if rows[u]>>v&1)
  h[s]=sum(end[s])
 c2=[0]*(1<<n)
 for s in range(1,1<<n):
  t=(s-1)&s
  while t:c2[s]+=h[t]*h[s^t];t=(t-1)&s
 c3=0;t=full-1
 while t:c3+=c2[t]*h[full^t];t=(t-1)&full
 return h[full],c2[full],c3

def run():
 kernels={k:observable_skeleton(ROWS,O,PORTS,k) for k in EXPECTED}
 need(kernels==EXPECTED,'complete finite skeleton coefficients')
 d=exp(skeleton(ROWS,O,PORTS));need(d=={(1,1):1,(2,1):1,(1,2):-1},'defect exponential')
 literal=[]
 for a,b in product(range(1,5),repeat=2):
  dd,u,h=counts(*expand(ROWS,O,PORTS,[a,b]));need(dd==1+2**a-2**b,'literal defect')
  for k,z in [('H',h),('W',h-dd),('U',u)]:need(evaluate(exp(kernels[k]),a,b)==z,'literal '+k)
  literal.append([a,b,dd,u,h])
 rays=[]
 for a,b in product(range(0,7),repeat=2):
  if not a+b:continue
  grouped=defaultdict(int)
  for (x,y),c in d.items():grouped[x**a*y**b]+=c*x*y
  grouped[1]-=1;grouped={k:v for k,v in grouped.items() if v}
  need((not grouped)==(a==b),'ray classifier')
  rays.append([a,b,not grouped])
 diag=[]
 for m in range(1,41):
  U=8*6**m+10*4**m+2*3**m-4*2**m;c=U//2
  need(c%2==1 and U>0,'exact2valuation')
  diag.append([m,U,c])
 # Nontransitive physical substitution, not a merely rational Farey pair.
 A=[18,20,9,19,36,79,31];B=[14,12,8,0]
 h,c,e=profiles(A);g,d2,f=profiles(B)
 need((h,c,g,d2)==(65,846,1,14),'module profiles')
 expected_D=h*g+c*g-h*d2
 expected_H=15*h*g+6*c*g+5*h*d2+2*c*d2
 expected_U=86*h*g+72*c*g+18*e*g+42*h*d2+34*c*d2+8*e*d2
 need(expected_D==1,'unimodular condition')
 actual=counts(*lift_modules(A,B));need(actual==(expected_D,expected_U,expected_H),'literal nontransitive module seed')
 # Both substituted modules nontransitive: a distinct physical Farey witness.
 A2=[30,4,8,18,6];B2=[114,164,73,19,38,76,154,61]
 hh,cc,ee=profiles(A2);gg,dd,ff=profiles(B2)
 need((hh,cc,gg,dd)==(5,54,589,6950),'both nontransitive profiles')
 pred=(hh*gg+cc*gg-hh*dd,
       86*hh*gg+72*cc*gg+18*ee*gg+42*hh*dd+34*cc*dd+8*ee*dd,
       15*hh*gg+6*cc*gg+5*hh*dd+2*cc*dd)
 need(pred[0]==1,'both nontransitive determinant')
 literal2=counts(*lift_modules(A2,B2));need(literal2==pred,'literal17bothnontransitive')
 second={'A':A2,'B':B2,'profiles':[hh,cc,ee,gg,dd,ff],'D_U_H':list(pred),'order':17}
 # Independent ordinal collar classification, with actual profile witnesses.
 collar_checks=0
 for ha,ca,hb,cb in [(65,846,1,14),(5,54,589,6950),(1,6,1,6)]:
  need(ha*hb+ca*hb-ha*cb==1,'starting collar unit')
  x=hb*(ca+2*ha);y=ha*(cb+2*hb)
  for aa,bb in product(range(12),repeat=2):
   dd=ha*hb+x*2**aa-y*2**bb
   need((dd==1)==((aa==bb) if ha*hb==1 else (aa==bb==0)),'collar classifier')
   collar_checks+=1
 # A different direction can look locally close without being an exact ray.
 need(1+2**5-2**4!=1,'off-diagonal hostile')
 out={'scope':'Complete two-port skeleton, sixteen literal chain controls, exact ray certificate controls, literal15- and17-vertex nontransitive module seeds; no global coverage claim',
 'collar_controls':collar_checks,'skeletons':{k:encode(v) for k,v in kernels.items()},'literal_chains':literal,'ray_controls':rays,'first40_diagonal_slopes':diag,
 'both_nontransitive_witness':second,'module_witness':{'A':A,'B':B,'profiles':[h,c,e,g,d2,f],'D_U_H':list(actual),'order':len(A)+len(B)+4},
 'closed_recurrence_coefficients':{'U':[15,-80,180,-144],'H':[7,-14,8]}}
 text=json.dumps(out,indent=2)+'\n';name='REPLAY.optimized.json' if not __debug__ else 'REPLAY.json';(HERE/name).write_text(text);print(json.dumps({'status':'PASS','literal_chain_controls':len(literal),'module_witness':out['module_witness']}))
if __name__=='__main__':run()
