#!/usr/bin/env python3
"""Self-contained exact controls for scoped local-screen and fibre limitations."""
from collections import defaultdict
from functools import lru_cache
from fractions import Fraction as Q
import json

def require(ok,msg):
 if not ok:raise RuntimeError(msg)
def close(p):
 p=list(p)
 for k in range(len(p)):
  for v in range(len(p)):
   if p[v]>>k&1:p[v]|=p[k]
 return tuple(p)
def audit(p,z):
 n=len(p);require(close(p)==p,'Not transitively closed');require(all(not p[v]>>v&1 for v in range(n)),'Cycle')
 degree=[n-1-p[v].bit_count()-sum(bool(p[u]>>v&1) for u in range(n)) for v in range(n)];require(max(degree)<=7,'Degree exceeded')
 counts=[0]*5;fibres=defaultdict(lambda:[0]*4)
 outside=tuple(v for v in range(3,9) if v!=z)
 def rec(mask,word):
  if len(word)==n:
   r={v:i for i,v in enumerate(word)};c=[1,int(r[0]<r[z]),int(r[z]<r[1]),int(r[2]<r[9]),int(r[9]<r[z])]
   for k,v in enumerate(c):counts[k]+=v
   key=tuple(v for v in word if v in outside)
   for k in range(4):fibres[key][k]+=c[k]
   return
  for v in range(n):
   if not mask>>v&1 and not p[v]&~mask:rec(mask|1<<v,word+[v])
 rec(0,[])
 @lru_cache(None)
 def dp(mask):
  if mask==(1<<n)-1:return 1
  return sum(dp(mask|1<<v) for v in range(n) if not mask>>v&1 and not p[v]&~mask)
 require(dp(0)==counts[0],'Independent total mismatch')
 return degree,counts,fibres
actual=[]
for p,z in [((512,516,0,512,512,4,532,513,685,0),5),((512,516,0,4,512,528,540,513,529,0),3)]:
 d,c,f=audit(p,z);actual.append({'pred':p,'z':z,'degree':d,'counts_E_xz_zy_be_ez':c,'probabilities':[str(Q(k,c[0])) for k in c[1:]]})
require(Q(sum(actual[0]['counts_E_xz_zy_be_ez'][1:3]),actual[0]['counts_E_xz_zy_be_ez'][0])>Q(4,3),'Aggregate obstruction missing')
c=actual[1]['counts_E_xz_zy_be_ez'];require(min(Q(c[1],c[0]),Q(c[2],c[0]),Q(c[4],c[0]))>Q(2,3),'Three-orientation obstruction missing')
core=(0,4,0,0,12,28,4,1,197);z=6
p=close([core[v]|(512 if v not in [2,z] else 0) for v in range(9)]+[0]);d,c,f=audit(p,z)
w1=(3,7,4,8,5);w2=(3,7,8,4,5)
require(f[w1]==[286,217,188,94] and f[w2]==[218,149,158,78],'Fibre mismatch')
mix=[9*f[w1][k]+2*f[w2][k] for k in range(4)];require(mix==[3010,2251,2008,1002],'Mixture mismatch')
require(3*mix[1]>2*mix[0] and 3*mix[2]>2*mix[0] and 3*mix[3]<mix[0],'Mixture signs missing')
seen=set();orbits=[]
for word in f:
 if word in seen:continue
 partner=list(word);i,j=word.index(4),word.index(8)
 if abs(i-j)==1:partner[i],partner[j]=partner[j],partner[i]
 partner=tuple(partner);require(partner in f,'Invalid exchange')
 orbit=[word] if word==partner else [word,partner];seen.update(orbit)
 C=[sum(f[w][k] for w in orbit) for k in range(4)];margin=6*C[2]-15*C[3]+C[0];require(margin<=0,'Orbit deficit failed')
 orbits.append({'words':orbit,'counts':C,'margin':margin})
require(len(f)==10 and len(orbits)==8,'Orbit census mismatch')
q=(1536,1540,1024,512,1548,1564,1028,1537,1733,0,0);qd,qc,qf=audit(q,z);margin=6*qc[2]-15*qc[3]+qc[0]
require(qc[0]==5095 and margin==9379,'Eleven-point limitation mismatch')
result={'status':'EXACT_CONTROLS_PASSED','scope':'Specific local inequalities and artificial fibre mixtures only; no counterexample to Gold or balance.','actual_root_screen_obstructions':actual,'same_poset_fibre_limit':{'pred':p,'z':z,'degree':d,'actual_total_counts':c,'words':[w1,w2],'fibre_vectors':[f[w1],f[w2]],'artificial_raw_weights':[9,2],'artificial_mixture':mix},'fixed_poset_orbit_repair':{'orbits':orbits},'additional_exterior_limit':{'pred':q,'z':z,'degree':qd,'counts':qc,'positive_margin':margin}}
print(json.dumps(result,indent=2))
