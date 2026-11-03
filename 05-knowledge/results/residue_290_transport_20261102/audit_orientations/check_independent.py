#!/usr/bin/env python3
"""Independent finite audit: no producer imports, exact integer/rational algebra."""
from collections import Counter
from fractions import Fraction as F
from math import gcd, comb
import json

G=0
def ck(c,s):
 global G
 G+=1
 if not c:raise RuntimeError(s)
def trim(a):
 a=list(a)
 while len(a)>1 and a[-1]==0:a.pop()
 return a
def add(a,b):
 c=[0]*max(len(a),len(b))
 for i,x in enumerate(a):c[i]+=x
 for i,x in enumerate(b):c[i]+=x
 return trim(c)
def scale(a,k):return trim([k*x for x in a])
def sub(a,b):return add(a,scale(b,-1))
def mul(a,b):
 c=[0]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  for j,y in enumerate(b):c[i+j]+=x*y
 return trim(c)
def divmodpoly(a,b):
 a=trim(a);b=trim(b);q=[0]*max(1,len(a)-len(b)+1)
 while a!=[0] and len(a)>=len(b):
  k=len(a)-len(b);c=F(a[-1],b[-1]);q[k]=c
  for j,x in enumerate(b):a[k+j]-=c*x
  a=trim(a)
 return trim(q),a
def rem(a,b):return divmodpoly(a,b)[1]
def ppow(a,k):
 r=[1]
 for _ in range(k):r=mul(r,a)
 return r
def evaluate(a,x):
 v=0
 for c in reversed(a):v=v*x+c
 return v

def cycles(mapping):
 unused=set(mapping);out=[]
 while unused:
  v=min(unused);row=[];w=v
  while w in unused:unused.remove(w);row.append(w);w=mapping[w]
  ck(w==v,'permutation cycle closure')
  out.append(row)
 return out

def char_poly_on_signed_edges(mapping,sign):
 out=[1]
 for row in cycles(mapping):
  s=1
  for x in row:s*=sign[x]
  out=mul(out,[-s]+[0]*(len(row)-1)+[1])
 return out

selected={}
for N in range(3,202,2):
 units=[r for r in range(1,N) if gcd(r,N)==1]
 mapping={r:2*r%N for r in units};orb=cycles(mapping)
 d=len(next(c for c in orb if 1 in c));hasminus=any(N-1 in c and 1 in c for c in orb)
 ck(all(len(c)==d for c in orb),('orbit length',N))
 reps=[r for r in units if r<N-r]
 qm={r:min(2*r%N,(-2*r)%N) for r in reps}
 sg={r:(1 if 2*r%N<N-(2*r%N) else -1) for r in reps}
 plus=char_poly_on_signed_edges(qm,{r:1 for r in reps})
 minus=char_poly_on_signed_edges(qm,sg)
 if hasminus:
  ck(d%2==0 and pow(2,d//2,N)==N-1,('half negative',N))
  eplus=ppow([-1]+[0]*(d//2-1)+[1],len(units)//d)
  eminus=ppow([1]+[0]*(d//2-1)+[1],len(units)//d)
 else:
  eplus=eminus=ppow([-1]+[0]*(d-1)+[1],len(units)//(2*d))
 ck(plus==eplus and minus==eminus,('signed characteristic',N))
 # Explicit all-sector action at every representative, not spectral inference.
 for r in reps:
  y=r;s=1
  for _ in range(d):s*=sg[y];y=qm[y]
  ck(y==r and s==1,('full signed return',N,r))
 if N in (7,9,145):
  selected[N]={'d':d,'oriented_lengths':[len(c) for c in orb],
    'folded_lengths':[len(c) for c in cycles(qm)],'negative_half':hasminus,
    'plus_charpoly':plus,'minus_charpoly':minus}
ck(selected[145]['oriented_lengths']==[28]*4,'145 oriented')
ck(selected[145]['folded_lengths']==[14]*4,'145 folded')
# Unique odd lifts and correct unit successor modulo290.
for r in range(145):
 lift=r if r%2 else r+145
 ck(lift%2==1 and lift%145==r,'odd lift')
 if gcd(r,145)==1:
  nxt=2*r%145;nxt=nxt if nxt%2 else nxt+145
  ck(147*lift%290==nxt,'odd lift dynamics')

phi_cache={}
def cyclo(n):
 if n not in phi_cache:
  a=[-1]+[0]*(n-1)+[1]
  for d in range(1,n):
   if n%d==0:
    a,r=divmodpoly(a,cyclo(d));ck(r==[0],('cyclotomic division',n,d))
  ck(all(F(x).denominator==1 for x in a),('integral Phi',n))
  phi_cache[n]=[int(x) for x in a]
 return phi_cache[n]
def trace_poly(n):
 p=cyclo(n);D=(len(p)-1)//2;res={j-D:c for j,c in enumerate(p)};out=[0]*(D+1)
 for k in range(D,-1,-1):
  c=res.get(k,0);out[k]=c
  for j in range(k+1):
   e=k-2*j;res[e]=res.get(e,0)-c*comb(k,j)
 ck(all(v==0 for v in res.values()),('Laurent trace reconstruction',n))
 return trim(out)

trace={}
for N,k,lam in ((7,3,8),(9,3,-8),(145,14,-16384)):
 R=trace_poly(N);y=[0,1];der=[1]
 for j in range(1,k+1):
  der=rem(mul(der,scale(y,2)),R)
  y=rem(sub(mul(y,y),[2]),R)
  ck((y==[0,1])==(j==k),('least trace return',N,j))
 ck(der==[lam],('trace derivative product',N))
 trace[N]={'degree':len(R)-1,'polynomial':R,'period':k,'multiplier':lam}
ck(trace[7]['polynomial']==[-1,-2,1,1],'R7')
ck(trace[9]['polynomial']==[1,-3,0,1],'R9')
# Exact odd-character projector in Q[zeta28].
phi28=cyclo(28)
log5={pow(2,a,5):a for a in range(4)}
log29={pow(2,b,29):b for b in range(28)}
ck(len(log5)==4 and len(log29)==28,'primitive roots for characters')
oddchars=[(a,b) for a in range(4) for b in range(28) if (a+b)%2]
ck(len(oddchars)==56,'odd characters')
for r in range(290):
 if gcd(r,290)!=1:continue
 co=[0]*28
 for a,b in oddchars:co[(7*a*log5[r%5]+b*log29[r%29])%28]+=1
 target=56*((r==1)-(r==289))
 ck(rem(co,phi28)==[target],('character projector',r))

# Every possible tail-start phase has one-signed partial coefficient sums.
def mask(n):return int(n%290==1)-int(n%290==289)
ck(sum(mask(n) for n in range(290))==0,'mask mean zero')
ck(sum(-int(gcd(n,290)==1) for n in range(290))==-112,'unit mean')
ck(sum(-int(n%290 in (1,289)) for n in range(290))==-2,'target mean')
for M in range(3,293):
 A=0;values={0}
 for j in range(M+1,M+581):
  A+=mask(j);values.add(A)
  ck(abs(A)<=1,('tail prefix bound',M,j))
 ck(values<={0,1} or values<={0,-1},('one signed prefix',M))
# Abel's bound checked exactly for independent positive decreasing weights.
for M in range(3,293):
 for power in (1,2):
  s=sum(F(mask(n),n**power) for n in range(M+1,M+581))
  # Infinite-tail bound follows analytically; finite endpoints can be bounded too
  # for these full-period cuts, whose accumulated coefficient is zero.
  ck(abs(s)<=F(1,(M+1)**power),('finite Abel control',M,power))

# Chosen three-point reversal identities in polynomial algebra.
t=[0,1];L=[F(-1,2),-1];D2=add(mul(sub(mul(t,t),[2]),sub(mul(t,t),[2])),[-2])
for N,c,expected_lam,factor in ((7,F(-7,4),1,[-1,1]),(9,F(-11,4),19,[0,1])):
 R=trace[N]['polynomial']
 difference=sub(add(mul(L,L),[c]),add(scale(D2,-1),[F(-1,2)]))
 ck(difference==mul(factor,R),('affine reversal exact',N))
 ck(8*evaluate(R,F(-1,2))==expected_lam,('image multiplier',N))
cycle=[F(-7,4),F(5,4),F(-1,4)]
for x,y in zip(cycle,cycle[1:]+cycle[:1]):ck(x*x-F(29,16)==y,'rational3 cycle')
ck(8*cycle[0]*cycle[1]*cycle[2]==F(35,8),'rational3 multiplier')
eta=F(-1,2)
ck(-eta-1==eta and -(eta*eta+7)/4==F(-29,16),'fixed eta')
ck(F(-29,16)==F(-7,4)-F(1,16),'28 plus1 parameter')

print(json.dumps({'status':'PASS','gates':G,'odd_conductors':[3,201],
 'selected_orbits':selected,'trace_controls':trace,'odd_character_count':56,
 'MRB_tail_start_phases':290,'MRB_exact_monotone_control_powers':[1,2],
 'scope':'Independent finite controls; all-order analytic statements checked by proof'},
 indent=2,sort_keys=True))
