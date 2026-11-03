from sympy.combinatorics import Permutation, PermutationGroup
from math import gcd,lcm
from collections import Counter
import json

def require(b,s):
 if not b:raise RuntimeError(s)
def comp(p,q):return tuple(p[q[i]] for i in range(len(p)))
def power(p,n):
 a=tuple(range(len(p)))
 for _ in range(n):a=comp(p,a)
 return a
def inverse(p):return tuple(p.index(i) for i in range(len(p)))
def support(p):return [i for i,j in enumerate(p) if i!=j]
a=tuple((5*i+8)%9 for i in range(9));b=list(range(9))
for x,y in [(4,2),(2,5),(5,3),(3,7),(7,4)]:b[x]=y
b=tuple(b)
c=power(comp(power(a,3),b),4)
require(c[0]==5 and c[5]==2 and c[2]==0 and len(support(c))==3,'3-cycle word')
gens=[tuple(range(9)),a,power(b,2),power(b,3),comp(a,power(b,3))]
supports=[support(comp(comp(g,c),inverse(g))) for g in gens]
require(supports==[[0,2,5],[0,6,8],[0,3,7],[0,4,7],[1,7,8]],'conjugate supports')
require(power(comp(a,power(b,3)),3)==tuple(range(9)),'odd A relation')
G=PermutationGroup(Permutation(list(a)),Permutation(list(b)))
require(G.order()==181440,'independent Schreier Sims A9')
# Product-action computation on 18 points rather than Cayley-state enumeration.
aa=Permutation([9*(1-i//9)+a[i%9] for i in range(18)])
bb=Permutation([9*(i//9)+b[i%9] for i in range(18)])
H=PermutationGroup(aa,bb)
require(H.order()==362880,'independent Schreier Sims product')
L=Permutation([(527*i+44)%72 for i in range(72)])
B=Permutation([next(j for j in range(72) if j%8==i%8 and j%9==b[i%9]) for i in range(72)])
K=PermutationGroup(L,B)
require(K.order()==362880,'actual 72 action group')
require(sorted(map(len,K.orbits()))==[9,9,18,18,18],'actual orbits')
require(Counter(L.cycle_structure)==Counter({1:2,2:11,6:8}),'cycle multiplicities')
require(all((527*x+44)%72==(34+23*(x-34))%72 for x in range(72)),'affine conjugacy')
# Controls: A alone is cyclic of order6; B alone order5; neither is A9.
require(PermutationGroup(Permutation(list(a))).order()==6,'A-only control')
require(PermutationGroup(Permutation(list(b))).order()==5,'B-only control')
# Polynomial over Z/9 preserves congruence modulo3, B does not.
require(4%3==7%3 and b[4]%3!=b[7]%3,'integer polynomial obstruction')
print(json.dumps({'status':'PASS','method':'independent Schreier-Sims on 9,18,72 points','orders':[G.order(),H.order(),K.order()],'supports':supports,'cycle_type':L.cycle_structure,'sectors':sorted(map(len,K.orbits()))},indent=2))
