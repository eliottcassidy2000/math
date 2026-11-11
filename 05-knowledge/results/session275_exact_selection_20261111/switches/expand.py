from carrier import *
from math import gcd
A=(6,16,2,7,13,31);B=(4,37,0,7,15,29)
N=(26,32,3,6,14,29);P=(52,41,2,5,14,28);C=(2,4,1)
I=matrix(profile((0,)))
def count(ts):
 M=I
 for T in ts:M=mul(M,matrix(profile(T)))
 return M[0][0]+M[0][3]
def literal(ts):
 T=(0,)
 for U in ts:T=fuse(T,U)
 h=hp(T)[-1];hp.cache_clear();return h
pn,pp=profile(N),profile(P)
assert pn==dict(H=37,h=9,Q=84,k=0,U=59,V=43,d=5,e=7)
assert pp==dict(H=37,h=9,Q=84,k=0,U=62,V=43,d=4,e=7)
levels=[]
for m in range(1,8):
 baseline=[B]*m+[N,C];base=count(baseline)
 def f(r):return 1244 if r==0 else (72966*64**r-396*5**r)//59
 ws=[5**i*f(m-1-i) for i in range(m)]+[6*5**m]
 assert gcd(*ws)==2 and all(ws[i]>sum(ws[i+1:]) for i in range(len(ws)))
 for bits in product((0,1),repeat=m+1):
  ts=[A if b else B for b in bits[:m]]+[P if bits[-1] else N,C]
  value=count(ts);assert value==base+sum(a*b for a,b in zip(ws,bits))
  if m==1:assert literal(ts)==value
 levels.append(dict(m=m,base=base,weights=ws,gcd=gcd(*ws)))
gates=[]
for S in [(0,),(2,0),C]:
 vals=[literal([A if x else B,S,op(A if y else B)]) for x,y in product((0,1),repeat=2)]
 p=profile(S);assert vals[3]-vals[2]-vals[1]+vals[0]==36*p['Q']+4*p['k']
 rev=[literal([op(A if x else B),S,A if y else B]) for x,y in product((0,1),repeat=2)]
 assert rev[3]-rev[2]-rev[1]+rev[0]==0
 gates.append(dict(spacer=S,counts=vals,interaction=vals[3]-vals[2]-vals[1]+vals[0],reversed_counts=rev))
# Exhaust all oriented-control words through5 slots and all Boolean assignments,
# checking exact quadratic reconstruction, not merely selected derivatives.
cubes=0
for n in range(1,6):
 for orientations in product((0,1),repeat=n):
  def h(bits):return count([(op(A if b else B) if o else A if b else B) for o,b in zip(orientations,bits)])
  z=(0,)*n;b=h(z);linear=[]
  for i in range(n):
   x=list(z);x[i]=1;linear.append(h(x)-b)
  pair={}
  for i,j in combinations(range(n),2):
   x=list(z);x[i]=x[j]=1;pair[i,j]=h(x)-b-linear[i]-linear[j]
   assert (pair[i,j]>0)==(orientations[i]==0 and orientations[j]==1)
  for x in product((0,1),repeat=n):assert h(x)==b+sum(a*t for a,t in zip(linear,x))+sum(v*x[i]*x[j] for (i,j),v in pair.items())
  cubes+=1
result=dict(parity_levels=levels,literal_gates=gates,quadratic_orientation_cubes=cubes)
Path(__file__).with_name('EXPANSION_RESULTS.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'parity_levels':len(levels),'literal_gates':len(gates),'quadratic_cubes':cubes}))
from collections import Counter
energy=[]
for m in range(1,8):
 vals=[count([A if b else B for b in bits]) for bits in product((0,1),repeat=m)]
 r=Counter(a-b for a in vals for b in vals)
 assert len(r)==3**m and sum(v*v for v in r.values())==6**m
 energy.append(dict(m=m,difference_support=len(r),energy=sum(v*v for v in r.values())))
result['signed_digit_energy']=energy
Path(__file__).with_name('EXPANSION_RESULTS.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'energy_levels':len(energy)}))
