"""Independent indexed source inventory, profile classification and literal graphs."""
from pathlib import Path
from functools import lru_cache
import json
HERE=Path(__file__).resolve().parent
R=json.loads((HERE/'upper-fibers-replay/VALIDATION_NORMAL.json').read_text())
@lru_cache(None)
def coeff(L):
 if L==1:return 1,1,1
 if L%2==0:
  g,Q,B=coeff(L//2);return g,Q,2*B
 x=(L+1)//2;g,Q,B=coeff(x);return g+x-x%2,Q+B,2*B
cap=5000;bank={(1,1,1,1)};queue=list(bank)
for D,E,U,d in queue:
 nxt=[(D,0,2*U,d)]
 for L in range(1,cap//U+1):
  g,Q,B=coeff(L);nxt.append((g*D-(L%2)*E+L*U,D-(L==1)*E,Q*D+B*U,D))
 for s in nxt:
  if s[2]<=cap and s not in bank:bank.add(s);queue.append(s)
bank.remove((1,1,1,1));strong={s for s in bank if s[1]}
A={s for s in strong if 4*s[2]-5*s[0]==1};B={s for s in strong if s[0]+2*s[2]-6*s[3]==1}
assert len(bank)==126910 and len(strong)==81548 and len(A)==41 and len(B)==4
Dset={s[0]for s in A};Bset={(s[2],s[3])for s in B};upper=[]
for q in bank:
 D,E,U,d=q
 if 9*D-3*E-9*U+4*d!=1:continue
 if D==U:
  assert E==d==1 and D>1 and (D-1)&(D-2)==0;label='ray'
 elif (D-2)%11==0 and (D-2)//11 in Dset:
  t=(D-2)//11;assert q==(11*t+2,5*t+1,12*t+2,6*t+1);label='A'
 else:
  assert (E-1)%5==0;dd=(E-1)//5;assert (D-4*dd-1)%4==0;uu=(D-4*dd-1)//4
  assert (uu,dd)in Bset and q==(4*uu+4*dd+1,5*dd+1,4*uu+5*dd+1,6*dd+1);label='B'
 upper.append((q,label))
assert set(upper)=={(tuple(r['profile']),r['branch'])for r in R['upper_response_one']}
for D,E,U,d in B:
 assert (E+1)%4==0;h=(E+1)//4
 assert (4*h-1,h,4*h-1,h)in bank and (D,E,U,d)==(12*h-3,4*h-1,15*h-4,7*h-2)

def build(word):
 T=[2,0];r=0
 for x in word:
  n=len(T);A=T.copy()+[0]
  for v in range(n):
   if x!='O'and v==r:A[n]|=1<<v
   else:A[v]|=1<<n
  T=A
  if x=='J':r=n
 return T,r

def reverse(T):
 n=len(T);return [((1<<n)-1)^row^(1<<i)for i,row in enumerate(T)]

def parent(word):
 V,c=build('J');W,a=build(word+'O');W=reverse(W);n=len(V);A=V.copy()+[row<<n for row in W]
 for i in range(n):
  for j in range(len(W)):
   if j==a and i in[c,0]:A[n+j]|=1<<i
   else:A[i]|=1<<(n+j)
 return A

def child(word):
 W,a=build(word+'J');W=reverse(W);other=[v for v in range(len(W))if v!=a];mp={a:0,**{v:i+5 for i,v in enumerate(other)}}
 n=len(W)+4;A=[0]*n;A[0]|=1<<1
 for u in range(len(W)):
  for v in range(len(W)):
   if W[u]>>v&1:A[mp[u]]|=1<<mp[v]
 tri=[2,3,4]
 for i,u in enumerate(tri):A[u]|=(1<<tri[(i+1)%3])|1
 for u in [1]:
  for v in tri+[mp[x]for x in other]:A[u]|=1<<v
 for u in tri:
  for v in other:A[u]|=1<<mp[v]
 return A

def counts(T):
 n=len(T);size=1<<n;P=[[0]*n for _ in range(size)];Q=[[0]*n for _ in range(size)]
 for i in range(n):P[1<<i][i]=1
 Q[1][0]=1
 for S in range(1,size):
  for u in range(n):
   if not P[S][u]:continue
   rest=T[u]&~S
   while rest:
    b=rest&-rest;rest-=b;v=b.bit_length()-1;P[S|b][v]+=P[S][u];Q[S|b][v]+=Q[S][u]
 return sum(P[-1]),sum(Q[-1][v]for v in range(n)if T[v]&1)

literal=[]
for r in R['literal_graph_checks']:
 w=r['source_word'];x=counts(parent(w));y=counts(child(w));assert list(x+y)==r['counts'];literal.append((r['branch'],x,y))
# These are outputs of the stated minimal-power pump, not minimal physical seeds.
for word,n in [('OOZZOOJJ',3823),('OZOJJJOZJ',6253)]:
 assert counts(parent(word))==(n,2)and counts(child(word))==(2*n+1,3)
 literal.append(('pump_choice',(n,2),(2*n+1,3)))
out={'status':'PASS','indexed_source_profiles':len(bank),'strong_profiles':len(strong),'upper_response_one_profiles':len(upper),'A_seed_profiles':len(A),'B_seed_profiles':len(B),'B_bridge_predecessors_confirmed':len(B),'literal_parent_child_pairs':literal,'method':'indexed M_L inventory and independent graph assembly/counting'}
(HERE/'UPPER_FIBERS_INDEPENDENT.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print(json.dumps(out,sort_keys=True))
