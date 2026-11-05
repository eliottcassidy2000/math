from pathlib import Path
from functools import lru_cache
from itertools import combinations
from math import gcd,prod
import json,hashlib
from core import inventories,constants,glue,gadget,cycles,require
from signature import *
P=Path(__file__).parent
@lru_cache(None)
def coeff(L):return constants(L)

def bank(cap):
 base=(1,1,1,1);Q=[base];seen={base}
 for D,E,U,d in Q:
  for L in range(1,cap//U+1):
   g,q,b=coeff(L);w=(g*D-(L%2)*E+L*U,D-int(L==1)*E,q*D+b*U,D)
   if w[2]<=cap and w not in seen:seen.add(w);Q.append(w)
 full=set()
 for s in seen:
  if s==(1,0,2,1):continue
  full.add(s);D,E,U,d=s;k=1
  while U<<k<=cap:full.add((D,0,U<<k,d));k+=1
 return full
cap=2000;S=bank(cap);mapped={}
for s in S:
 z=signature(s);require(recover(*z)==s,'piecewise exact recovery')
 require(z not in mapped,'physical signature collision');mapped[z]=s
 require(gcd(*s)==gcd(*z),'physical primitive content')
 D,E,U,d=s
 require(z[0]<=3*U and z[1]<=2*U and 0<=z[2]<=U,'cubic coordinate box')
small=31;target={signature(s):s for s in bank(small)};integer=0
for A in range(1,3*small+1):
 for B in range(1,2*small+1):
  for C in range(small+1):
   s=recover(A,B,C)
   if s is None or s[2]>small:got=None
   else:got=member(A,B,C)
   require((got is not None)==((A,B,C)in target),'complete bounded signature membership')
   if got:require(got[0]==target[A,B,C],'exact membership response')
   integer+=1
literal=0
for n in range(1,6):
 pairs=list(combinations(range(n),2))
 for bits in range(1<<len(pairs)):
  T=[0]*n
  for i,(a,b)in enumerate(pairs):T[a if bits>>i&1 else b]|=1<<(b if bits>>i&1 else a)
  I=inventories(T)
  for r in range(n):
   if I['s'][r]!=1:continue
   s=(I['H'],I['e'][r],I['R'][r],I['allH'][((1<<n)-1)^(1<<r)])
   require(member(*signature(s))[0]==s,'actual source graph recovery');literal+=1
# Direct rank/kernel certificates, without a symbolic-algebra dependency.
N=[[1,-1,1,0],[-1,1,0,-1],[1,0,3,-1],[0,-1,-1,0]];ker=[-2,-1,1,1]
require(all(sum(a*b for a,b in zip(row,ker))==0 for row in N),'N kernel')
J=[[1,-1,1],[-1,1,0],[1,0,3]]
det=J[0][0]*(J[1][1]*J[2][2]-J[1][2]*J[2][1])-J[0][1]*(J[1][0]*J[2][2]-J[1][2]*J[2][0])+J[0][2]*(J[1][0]*J[2][1]-J[1][1]*J[2][0])
require(det==-1,'unimodular signature form')
# Every small physical source pair, including weak sources and erased-d collisions.
pairs_checked=0
for s in bank(41)-{(1,1,1,1)}:
 for t in bank(41)-{(1,1,1,1)}:
  D,E,U,d=s;F,G,V,f=t;n=n_value(s,t);m=m_value(s,t)
  require(n==n_signature(signature(s),signature(t)),'quotient bilinear form')
  require(4*n>=15*D*F+5*D+5*F-1,'sharp additive N product bound')
  require(2*m>=17*D*F+3*D+3*F-1,'sharp additive M product bound')
  b=U-D;c=V-F;sl=U-D+E-d;sr=V-F+G-f
  n0=6*D*F-2*(D*G+E*F)-E*G
  m0=11*D*F-2*(D*G+E*F)-2*E*G
  require(n==n0+b*(3*F-2*G)+c*(3*D-2*E)+b*c+sl*(G+V)+sr*(E+U),'N positive slack identity')
  require(m==m0+b*(7*F-3*G)+c*(7*D-3*E)+4*b*c+sl*(G+V)+sr*(E+U),'M positive slack identity');pairs_checked+=1
large=0
for m in range(1,51):
 H=prod((1<<j)+1 for j in range(1,m+1));sig=set()
 for last in range(1,m+1):
  d=H//((1<<last)+1);s=(H,0,2*H,d);z=signature(s)
  require(recover(*z)==s,'unbounded ambiguity encoded');sig.add(z);large+=1
 require(len(sig)==m,'ambiguity signatures distinct')

def converse(T):return [sum(1<<j for j in range(len(T))if T[j]>>i&1)for i in range(len(T))]
def pair_graph(V,r,kind):
 U,a=glue(V,{r},[0])
 if kind=='M':U,a=glue(U,{a},[0])
 G=converse(V);order=[r]+[i for i in range(len(V))if i!=r]
 G=[sum(1<<j for j,y in enumerate(order)if G[x]>>y&1)for x in order]
 return glue(U,{r,a},G)[0]
sharp=[]
for j in range(1,4):
 V,r=glue([0],{0},gadget(1<<j));V,r=glue(V,{r},[0]);x=(1<<j)+1
 I=inventories(V);s=(I['H'],I['e'][r],I['R'][r],I['allH'][((1<<len(V))-1)^(1<<r)])
 require(s==(2*x-1,x-1,2*x,x),'physical sharp source')
 vals={}
 for kind,want in [('N',15*x*x-2*x),('M',34*x*x-6*x+1)]:
  T=pair_graph(V,r,kind);require(inventories(T)['H']==want and cycles(T)==2,'physical sharp bilinear graph');vals[kind]=dict(H=want,order=len(T))
 sharp.append(dict(j=j,x=x,source=s,values=vals))
receipt=dict(status='PASS',source_U_cap=cap,physical_response_quadruples=len(S),distinct_signatures=len(mapped),integer_membership_controls=integer,membership_U_cap=small,labelled_start_one_root_controls=literal,small_source_pair_controls=pairs_checked,unbounded_ambiguity_controls=large,rank=3,signature_matrix_determinant=det,sharp_physical_examples=sharp)
text=json.dumps(receipt,sort_keys=True,indent=2)+'\n';(P/'VALIDATION.json').write_text(text);print(text,end='')
