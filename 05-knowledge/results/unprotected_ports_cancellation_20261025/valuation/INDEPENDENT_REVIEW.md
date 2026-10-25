# Independent review of the coefficient-characteristic cancellation jump

2026-10-25. **PASS** for odd primes p and integers r>=1, with q=p^r and
h=p^(r-1). The algebraic identity, termwise integrality, exact two depths,
and finite-logarithm specialization are valid. The concrete 50-vertex
tournament consequence also passes the independent finite audit below.
No infinite TRR existence theorem was independently reviewed here.

The cyclic extension added later was separately reviewed by the recurrence
reviewer. Its exact recursion aligns the Euler-phi coefficients and leaves
(u^(p^r)-A^(p^(r-1)))/p^r, with A=u^p+p*C_p(u). The expansion has integral
coefficients of valuation j-1-v_p(j), so exactly j=1 survives modulo p.
Iteration gives the stated all-scale formula. The unique lowest term is
at degree p^(r-1), with coefficient -1. Every other scale has strictly
higher least degree. Distinct scale supports have exponent valuation a,
so they are disjoint and give exactly1+r*(p-1) nonzero terms.

The same reviewer independently checked the common ordinary derivative:
after differentiating the cyclic formula only phi(1),phi(p) survive mod p.
The elementary abelian derivative agrees. Both polynomials are monic,
but their difference need not be; its zero constant, canceled leading
term and zero derivative put it in span{u^p,...,u^(q-p)}. No assertion
that every vector in this span is realized is made. At r=2 the exact
difference is (p-1)*C_p((1+u)^p-1). These are proof audits, supplemented
by the producer's independent direct binomial-shift controls.

## 1. Identity, valuations, and the smallest cases

For the regular elementary-abelian action, every nonidentity element has
h cycles of length p. Hence its cut-orbit polynomial is

    P(t)=((1+t)^q+(q-1)(1+t^p)^h)/q.

Put C=((1+u)^p-1-u^p)/p, an integral polynomial of degree p-1 with zero
constant term and coefficient one at u. Since q,p,h are odd,

    Q(u)=-P(-1-u)
        =[u^q+(q-1)(u^p+pC)^h]/q
        =u^q+sum_(j=1)^h a_j*u^(q-pj)*C^j,
    a_j=(q-1)*binom(h,j)*p^j/q.                     (R1)

The extraction of the j=0 term is correct: its coefficient joins the
separate u^q term to give q before division.

For 1<=j<=h,

    v_p binom(h,j)=r-1-v_p(j).

One direct proof uses binom(h,j)=(h/j)binom(h-1,j-1). Modulo p,
(1+x)^(h-1) has coefficients (-1)^i for 0<=i<h, by dividing
(1+x)^h=1+x^h by 1+x. Thus the latter binomial factor is a p-adic unit.
As q-1 is also a unit,

    v_p(a_j)=j-1-v_p(j).                            (R2)

Every a_j is an integer: its only possible denominator is a power of p,
and (R2) is nonnegative. For j>=2 and p odd it is at least one. In
particular p=3,j=3 gives valuation one, not zero. The j=1 coefficient
is exactly q-1. Reduction must be performed after the integer division;
division by q inside F_p would be invalid.

It follows term by term that

    Q mod p = u^(q-p)*(u^p-C mod p).                 (R3)

For 1<=i<p,

    C_i=binom(p,i)/p=binom(p-1,i-1)/i
        =(-1)^(i-1)/i mod p.

Thus the finite-logarithm formula is correct, and the coefficient of
u^(q-p+1) in (R3) is -1. The modular depth is exactly q-p+1.

Over Q the least degree of the j-th term in (R1) is q-(p-1)j, strictly
decreasing as j increases. Its unique minimum is h, attained at j=h.
The coefficient is

    (q-1)*p^(h-r)>0.

Here h>=r for every odd p and r>=1, so this is an integer as required.
There is no cancellation between different j terms at that lowest degree.

At r=1 there is only j=1: Q=u^p+(p-1)C, and both depths are one.
For p=3,r=2 the independently expanded polynomial is

    Q=24u^3+72u^4+96u^5+72u^6+32u^7+8u^8+u^9.

Its rational depth is three and its F_3 depth is seven. For p=5,r=2,
the rational leading coefficient is 3000 at degree five, whereas

    Q mod 5 = 4u^21+3u^22+3u^23+4u^24+u^25.

The odd-prime hypothesis is load-bearing both for the substitution signs
and for the disappearance of all j>=2 terms. The statement must not be
extended to p=2 from this proof.

## 2. Independent concrete tournament check

The sibling [positive-jet checker](../counting/check_positive_jet.py)
specifies a Cayley tournament on C_5^2. The independent implementation
below reads its connection set as input data but imports no producer code.
It refines rooted vertex blocks using only outgoing counts into each
current block. The class counts are 2,3,15,21,25, ending in singleton
classes. Every root-fixing automorphism must preserve every refinement
class, so the root stabilizer is trivial. Translations preserve the arcs
by their difference definition, giving exactly the regular group C_5^2.

After replacing each vertex by a directed two-vertex block, an independent
literal comparison of outside neighbourhoods finds exactly those 25
homogeneous pairs and no others. Their intrinsic orders prohibit internal
swaps, so the full automorphism group of the inflated tournament acts on
its 25 gaps as that same regular C_5^2. There is only one nonempty invariant
gap set, namely all gaps, and its odd cardinality makes it a valid toggle.
Thus the genuine 50-vertex example has invariant-toggle rank one,
rational cancellation depth five, and coefficient-characteristic-five
depth 21. The three statistics are different observables.

The checker also enumerates all 680 cut subsets at elementary-abelian
orders 3,5,7,9 and matches their literal orbit counts to P. Fourteen
formula cases cover p=3 with r=1..4, p=5,7 with r=1..3, and p=11,13 with
r=1..2. Each case compares direct substitution in P against (R1), checks
every prefactor valuation and integrality, and verifies the entire
polynomial modulo p rather than just its first nonzero coefficient.

## 3. Complete coefficient-characteristic spectrum

The extension to every coefficient field also passes. In characteristic
zero the depth is h. At a prime ell not dividing p(q-1), the leading
coefficient (q-1)p^(h-r) remains nonzero, so the depth remains h. At
ell=p it is q-p+1 by (R3). Finally, if ell divides q-1, then q is a unit
and the numerator in the first expression of (R1) reduces to u^q;
indeed q=1 modulo ell. Thus Q=u^q modulo ell and the depth is exactly q.
These disjoint cases exhaust all characteristics, including ell=2. The
checker below tests the entire polynomial in the last case and all depths
at eight coefficient primes across the fourteen formula cases.

## 4. Applicability to the incoming observer tournament family

Read `session62_agent_theme_bridges_20261025/clone/README.md` directly from
commit `61a9ae6734`, especially section 5. Its construction has
2*3^k+3k vertices for k>=2, with 3^k clone-pair gaps and singleton
observers. Its degree distinctions force preservation of the observer
triples, their rotations force the coordinate translations, and every
translation preserves the given arcs. The pair-separation argument makes
the intended clone pairs and singleton observers exactly the canonical
modules. Thus the full automorphism group acts regularly on the gaps as
C_3^k. It is not a regular action on the entire vertex set, which also
contains the observers.

The present formula needs regularity on gaps only, so it applies with
p=3,r=k. The resulting exact depths are 3^(k-1) over Q, 3^k-2 in
characteristic three, and 3^k in characteristic two, with invariant-toggle
rank one. At k=2 this gives a 24-vertex example with depths 3,7,9; at k=3
the 63-vertex example has depths 9,25,27. The restriction k>=2 remains in
force. This application uses the explicit observer proof, not a claim
that every elementary-abelian group is a tournament regular representation.

## 5. Reproduction

Run from this directory, optionally adding `-O` after `python`:

    python -c "from pathlib import Path; t=Path('INDEPENDENT_REVIEW.md').read_text(); exec(t.split(chr(96)*3+'python\n')[1].split(chr(96)*3)[0])"

```python
from itertools import product
from math import comb
import json

def need(ok,m):
    if not ok:raise RuntimeError(m)
def val(a,p):
    n=0
    while a%p==0:a//=p;n+=1
    return n
def multiply(a,b):
    c=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):c[i+j]+=x*y
    return c
def poly(p,r):
    q=p**r;h=q//p
    numer=[comb(q,k)+(q-1)*(comb(h,k//p) if k%p==0 else 0)
           for k in range(q+1)]
    need(all(v%q==0 for v in numer),('orbit integrality',p,r))
    P=[v//q for v in numer]
    Q=[sum((-1)**(k+1)*P[k]*comb(k,j) for k in range(j,q+1))
       for j in range(q+1)]
    C=[0]+[comb(p,j)//p for j in range(1,p)]
    power=[1];expand=[0]*(q+1);expand[q]=1
    need(all(C[j]%p==(-1)**(j-1)*pow(j,-1,p)%p for j in range(1,p)),
         'finite logarithm')
    for j in range(1,h+1):
        need(val(comb(h,j),p)==r-1-val(j,p),('binomial valuation',p,r,j))
        a=(q-1)*comb(h,j)*p**j
        need(a%q==0,'termwise integrality');a//=q
        need(val(a,p)==j-1-val(j,p),'prefactor valuation')
        power=multiply(power,C)
        for k,c in enumerate(power):expand[q-p*j+k]+=a*c
    need(Q==expand,'independent binomial change of basis')
    need(all(c>=0 for c in Q),'positive coefficients')
    need(next(i for i,c in enumerate(Q) if c)==h,'rational depth')
    need(Q[h]==(q-1)*p**(h-r),'rational leading coefficient')
    target=[0]*(q+1);target[q]=1
    for j,c in enumerate(C):target[q-p+j]-=c
    need([c%p for c in Q]==[c%p for c in target],'modular polynomial')
    need(next(i for i,c in enumerate(Q) if c%p)==q-p+1,'modular depth')
    return P,Q
records=[];characteristic_checks=0
for p,R in [(3,4),(5,3),(7,3),(11,2),(13,2)]:
    for r in range(1,R+1):
        P,Q=poly(p,r)
        records.append((p,r,p**r//p,p**r-p+1,Q[p**r//p]))
        q=p**r;h=q//p
        for ell in (2,3,5,7,11,13,17,19):
            expected=q-p+1 if ell==p else (q if (q-1)%ell==0 else h)
            need(next(i for i,c in enumerate(Q) if c%ell)==expected,
                 ('complete characteristic spectrum',p,r,ell))
            if (q-1)%ell==0:
                need([c%ell for c in Q]==[0]*q+[1], 'dividing q-minus-one polynomial')
            characteristic_checks+=1
subset_checks=0
for p,r in [(3,1),(5,1),(7,1),(3,2)]:
    V=list(product(range(p),repeat=r));index={v:i for i,v in enumerate(V)}
    perms=[[index[tuple((x+y)%p for x,y in zip(v,a))] for v in V] for a in V]
    reps={}
    for mask in range(1<<len(V)):
        key=min(sum(1<<g[i] for i in range(len(V)) if mask>>i&1) for g in perms)
        reps[key]=mask.bit_count();subset_checks+=1
    counts=[0]*(len(V)+1)
    for k in reps.values():counts[k]+=1
    need(counts==poly(p,r)[0],('literal orbit enumeration',p,r))
V=list(product(range(5),repeat=2))
S={(0,3),(0,4),(1,0),(3,0),(3,1),(3,2),
   (3,3),(3,4),(4,1),(4,2),(4,3),(4,4)}
def arc(i,j):return tuple((b-a)%5 for a,b in zip(V[i],V[j])) in S
need(all(arc(i,j)!=arc(j,i) for i in range(25) for j in range(i+1,25)),
     'base tournament')
blocks=[{0},set(range(1,25))];trace=[2]
while len(blocks)<25:
    refined=[]
    for B in blocks:
        by={}
        for v in B:
            signature=tuple(sum(arc(v,w) for w in D) for D in blocks)
            by.setdefault(signature,set()).add(v)
        refined+=list(by.values())
    need(len(refined)>len(blocks),'root stabilizer not discharged')
    blocks=refined;trace.append(len(blocks))
need(trace==[2,3,15,21,25],'independent refinement trace')
def doubled_arc(i,j):
    a,b=i//2,j//2
    return i<j if a==b else arc(a,b)
pairs={(i,j) for i in range(50) for j in range(i+1,50)
       if all(doubled_arc(i,k)==doubled_arc(j,k)
              for k in range(50) if k not in (i,j))}
need(pairs=={(2*i,2*i+1) for i in range(25)},'exact doubled homogeneous pairs')
print(json.dumps({'status':'PASS','formula_cases':records,'characteristic_checks':characteristic_checks,
    'literal_cut_subsets':subset_checks,'rooted_refinement':trace,
    'homogeneous_pair_count':len(pairs)},indent=2))
```
