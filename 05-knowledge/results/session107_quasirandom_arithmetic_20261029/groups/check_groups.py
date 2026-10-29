"""Elementary checks only; no claim to construct Monster or compute its characters."""
import json, math
from fractions import Fraction
from itertools import permutations
from pathlib import Path
Mexp={2:46,3:20,5:9,7:6,11:2,13:3,17:1,19:1,23:1,29:1,31:1,41:1,47:1,59:1,71:1}
Bexp={2:41,3:13,5:6,7:2,11:1,13:1,17:1,19:1,23:1,31:1,47:1}
M=math.prod(p**a for p,a in Mexp.items());B=math.prod(p**a for p,a in Bexp.items())
new=sorted(set(Mexp)-set(Bexp));same=sorted(p for p in Bexp if Mexp[p]==Bexp[p]);inc=sorted(p for p in Bexp if Mexp[p]>Bexp[p])
N=math.prod(p**(Mexp[p]-Bexp[p]) for p in inc)
assert N==39046339332000 and M//B==194478922284018372000
assert M==B*N*math.prod(new)
assert 196883==47*59*71==1+4371+96255+96256
assert 4371==3*31*47
assert Fraction(196883,4371)==Fraction(4189,93)
f=lambda x:x*x-x+1
orb=[2]
for _ in range(4):orb.append(f(orb[-1]))
assert orb==[2,3,7,43,1807] and [x%19 for x in orb]==[2,3,7,5,2]
assert 1807==13*139
orders={q:next(k for k in range(1,q) if pow(43,k,q)==1) for q in (13,139)}
assert orders=={13:6,139:6}
def prime(n):return n>=2 and all(n%d for d in range(2,math.isqrt(n)+1))
def den(n):return math.prod(d+1 for d in range(1,n+1) if n%d==0 and prime(d+1))
assert [den(n) for n in (6,42,1806)]==[42,1806,1806]
Q={1,2,4}; arc=lambda x,y:(y-x)%7 in Q
ordinary=[];signed=[]
for perm in permutations(range(7)):
    same=all(arc(perm[x],perm[y])==arc(x,y) for x in range(7) for y in range(7) if x!=y)
    opp=all(arc(perm[x],perm[y])!=arc(x,y) for x in range(7) for y in range(7) if x!=y)
    if same:ordinary.append(perm)
    if same or opp:signed.append(perm)
affine={tuple((a*x+b)%7 for x in range(7)) for a in range(1,7) for b in range(7)}
squareaffine={tuple((a*x+b)%7 for x in range(7)) for a in Q for b in range(7)}
assert set(ordinary)==squareaffine and set(signed)==affine
# Verify slope projection is a surjective homomorphism for each group.
for slopes in (Q,set(range(1,7))):
    for a in slopes:
        for b in range(7):
            for c in slopes:
                for d in range(7):
                    composite=tuple((a*(c*x+d)+b)%7 for x in range(7))
                    assert composite==tuple(((a*c)%7*x+(a*d+b)%7)%7 for x in range(7))
result={'scope':'Elementary arithmetic and exhaustive Paley7 permutation checks; cited character degrees not recomputed','new_primes':new,'unchanged_primes':sorted(p for p in Bexp if Mexp[p]==Bexp[p]),'increased_primes':inc,'N':N,'order_ratio':M//B,'centralizer_index':M//(2*B),'degrees':{'M':196883,'B':4371},'degree_ratio':str(Fraction(196883,4371)),'productfree_density_bounds_approx':{str(d):d**(-1/3) for d in (196883,4371)},'integer_orbit':orb,'mod19_orbit':[x%19 for x in orb],'orders_of_43':orders,'bernoulli_denominator_chain':[6,42,1806,1806],'Paley7_ordinary_order':len(ordinary),'Paley7_signed_order':len(signed),'all_assertions_passed':True}
print(json.dumps(result,indent=2,sort_keys=True))
