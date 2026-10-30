#!/usr/bin/env python3
"""Independent finite witnesses for the scale-aware modular observation theorem."""
import itertools,json,math

def need(b,msg):
    if not b: raise RuntimeError(msg)
def lcm(xs):
    z=1
    for x in xs:z=math.lcm(z,x)
    return z
def obs(s,d,ms):return tuple(s*d%m for m in ms)
checks=0
for s in range(1,41):
  for ms in itertools.combinations(range(1,13),3):
    L=lcm(ms);Q=L//math.gcd(L,s)
    need(Q==lcm(m//math.gcd(m,s) for m in ms),'valuation formula')
    for d in range(-4,9):
      for e in range(-4,9):
        need((obs(s,d,ms)==obs(s,e,ms))==((d-e)%Q==0),'exact kernel')
        checks+=1
for n in range(101):
  M=n+4
  L=lcm(range(1,2*M+1,2))
  need(L>2**n,'odd cutoff')
  need(math.comb(2*M,M)<=lcm(range(1,2*M+1)),'central binomial magnitude')
  need(lcm(range(1,2*M+1))%math.comb(2*M,M)==0,'central binomial divisibility')
# Finite sharp robust criterion, independently compare every pair's Hamming distance.
robust=0
for s in range(1,9):
 for ms in itertools.combinations(range(2,9),4):
  qs=[m//math.gcd(m,s) for m in ms]
  for B in range(1,13):
   vectors=[obs(s,d,ms) for d in range(B+1)]
   mind=min(sum(x!=y for x,y in zip(a,b)) for a,b in itertools.combinations(vectors,2))
   for e in [0,1,2]:
    exact=all(lcm(q for i,q in enumerate(qs) if i not in S)>B for k in range(min(2*e,4)+1) for S in itertools.combinations(range(4),k))
    need(exact==(mind>2*e),'robust deletion criterion');robust+=1
# Explicit loss and a genuinely stabilizer-sensitive Burnside example: D4 on square binary colorings.
ms=[25,50,100,200]
need(obs(100,0,ms)==obs(100,2,ms),'hostile scale')
perms=set()
for k in range(4):
 perms.add(tuple((i+k)%4 for i in range(4)))
 perms.add(tuple((k-i)%4 for i in range(4)))
X=list(itertools.product(range(2),repeat=4))
fixed=[sum(all(x[i]==x[p[i]] for i in range(4)) for x in X) for p in perms]
orbits={min(tuple(x[p[i]] for i in range(4)) for p in perms) for x in X}
need(sum(fixed)==len(perms)*len(orbits),'Burnside')
need(len(orbits)==6 and sum(fixed)==48,'D4 example')
need(len(X)//len(perms)!=len(orbits),'nonfree hostile')
for m in range(1,41):need((sum(fixed)%(8*m))//8==len(orbits)%m,'amplification')
print(json.dumps({'status':'PASS','kernel_pair_checks':checks,'robust_checks':robust,'odd_cutoff_n_through':100,'burnside_square':{'group_order':8,'points':16,'orbits':6,'numerator':48},'scaled_hostile':{'s':100,'moduli':ms,'effective_modulus':2}},indent=2))
