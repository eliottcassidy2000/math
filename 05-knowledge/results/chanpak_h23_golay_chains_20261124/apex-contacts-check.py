from collections import Counter
G=sum(1<<i for i in (0,2,4,5,6,10,11))
basis=[(G<<i)|1<<23 for i in range(12)] # weight underlying generator odd
code={0}
for b in basis:code|={x^b for x in code}
assert len(code)==4096
def coeff(W):
 ob=[i for i in range(12) if W>>i&1]
 B=sum(1<<(11-i) for i in ob)
 E=sum(2 if j-i==1 else 1<<(j-i-1) for i in ob for j in ob if j>i)
 alp=sum(2**(i+1)*3**(11-i) for i in ob); gam=sum(2**(i+1) for i in ob)
 bet=8192+gam-alp
 for i in ob:
  for j in range(i+1,12):
   if not W>>j&1:
    d=j-i-1;bet+=2**(11-d)*(3**d+1)
 return [1-E,B,0],[-gam,bet,alp]
def Hcoe(w):
 D,Q=coeff(w&4095);E,R=coeff(w>>12)
 return tuple(-D[i]*E[j]+Q[i]*E[j]+3*D[i]*R[j]+4*Q[i]*R[j] for i in range(3) for j in range(3))
def evalat(C,t,u):
 A=[1,2**t,3**t];B=[1,2**u,3**u]
 return sum(C[3*i+j]*A[i]*B[j] for i in range(3) for j in range(3))
Cs={w:Hcoe(w) for w in code}
print('distinct bivariate',len(set(Cs.values())))
def same(C):return C[0],C[1]+C[3],C[2]+C[6],C[4],C[5]+C[7],C[8]
print('distinct univariate',len(set(map(same,Cs.values()))))
K=Counter(c[-1] for c in Cs.values());print('distinct leading',len(K),'maxmult',max(K.values()))
for t in [0,1,2,3,4,5,6,8,10,12,15,20,30,40]:
 V=Counter(evalat(c,t,t) for c in Cs.values())
 print(t,len(V),max(V.values()))
U=[same(c) for c in Cs.values()]
assert len(set(U))==4096
W=sum(max(c[i] for c in U)-min(c[i] for c in U) for i in range(6))
T=0
while W*3**T>=4**T:T+=1
for t in range(1,T):
    bases=[1,2**t,3**t,4**t,6**t,9**t]
    assert len({sum(x*y for x,y in zip(c,bases)) for c in U})==4096,t
print('exact all-width distinction threshold',T,'using coefficient-span',W)
from math import isqrt
P=139472229779
assert evalat(Cs[0x1007b4],1,1)==P and all(P%d for d in range(3,isqrt(P)+1,2))
assert evalat(Cs[0],1,1)==1073807359 and all(1073807359%d for d in range(3,isqrt(1073807359)+1,2))
groups={}
for w,c in Cs.items():groups.setdefault(c[-1],[]).append((w,same(c)))
for lead,A in sorted(groups.items()):
 if len(A)>1:print('equal leading',lead,[(hex(w),u[-2]) for w,u in sorted(A)])
from itertools import permutations
T=[10,4,9,18,39,15];a,b=5,1
S=Counter()
for P in permutations([0,2,3,4]):
 for w in range(1<<10):
  Q=[];x=y=0
  for g in range(5):
   if w>>(2*g)&1:Q.append(a);x+=1
   if w>>(2*g+1)&1:Q.append(b);y+=1
   if g<4:Q.append(P[g])
  if x and y and all(T[i]>>j&1 for i,j in zip(Q,Q[1:])):S[x,y]+=1
assert dict(S)=={(1,1):23,(1,2):11,(2,1):9,(2,2):4}
for t in [1,2,3]:
 for word in (0,1,31,2048,4095,1972,256):
  tags=[word>>i&1 for i in range(12)]+[0]*t
  d=sum((1<<(j-i-1))+(j==i+1) for i in range(len(tags)) for j in range(i+1,len(tags)) if tags[i] and not tags[j])+1
  q=(1<<(len(tags)+1))+sum((1<<(len(tags)-j+i))*(3**(j-i-1)+1)
       for i in range(len(tags)) for j in range(i+1,len(tags)) if tags[i] and not tags[j])
  D,Q=coeff(word)
  assert(d,q)==tuple(sum(c*b for c,b in zip(e,[1,2**t,3**t])) for e in [D,Q])
print('literal skeleton and raw slot checks verified')

assert max(K.values())==2 and len(K)==4093
assert all(len({u[-2] for _,u in group})==len(group) for group in groups.values())
assert W==17696218772230 and W*3**107<4**107
assert sum(w.bit_count()==8 and (w&4095).bit_count()==4 for w in code)==255
min_gap=None
for t in range(1,108):
    bases=[1,2**t,3**t,4**t,6**t,9**t]
    ordered=sorted((sum(x*y for x,y in zip(same(c),bases)),w) for w,c in Cs.items())
    candidate=min((y[0]-x[0],t,x,y) for x,y in zip(ordered,ordered[1:]))
    min_gap=candidate if min_gap is None else min(min_gap,candidate)
assert min_gap==(1023132,1,(3147435342621,0x18cf4b),(3147436365753,0x1e49d9))
print('exact global same-tail minimum gap:',min_gap)

# An independent subset DP for both complete modules of the prime witness.
def module_dp(controls,t):
    n=13+t
    g=[((1<<n)-(1<<(i+1))) for i in range(n)]
    g[-1]=controls
    for i in range(12):
        if (controls>>i)&1: g[i]^=1<<(n-1)
    F=(1<<n)-1
    dp=[[0]*n for _ in range(1<<n)]
    for i in range(n): dp[1<<i][i]=1
    counts=[1]+[0]*F
    for s in range(1,F+1):
        for i,z in enumerate(dp[s]):
            if not z: continue
            counts[s]+=z
            choices=g[i]&~s
            while choices:
                b=choices&-choices
                choices^=b
                dp[s|b][b.bit_length()-1]+=z
    return counts[F],sum(counts[s]*counts[F^s] for s in range(F+1))
Da,Qa=module_dp(0x1007b4&4095,1)
Db,Qb=module_dp(0x1007b4>>12,1)
assert 4*Qa*Qb+Qa*Db+3*Da*Qb-Da*Db==139472229779
