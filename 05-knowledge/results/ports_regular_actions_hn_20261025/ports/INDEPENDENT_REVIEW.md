# Independent review of the nonconsecutive-port repairs

2026-10-25. Reviewed `README.md` Theorems A, B, C and `check.py`.
**PASS for the stated families.** No claim of full black-core connectivity.
The reviewer did not import the producer's checker or discovery compiler.

## Structural audit

The invariant is the literal labelled tournament, not an isomorphism class.
Three K moves flip a pair according to the complement of the XOR of the
three unordered path-edge sets. Consequently the displayed XOR identity is
exactly the required equality with K_P(T). This proves the target only
after each of the three Hamiltonian words has separately been shown legal.
The producer checks these obligations separately and does not use XOR as a
surrogate for legality.

**Theorem A.** Its four guards partition the full source family. The two
inherited cases satisfy the stated odd-prefix exclusions: h=n-5>=7 cannot
equal any excluded value for ell=1 or ell=3. In the two new cases, the
insertion predecessor h-1 remains in the unmoved block. For the triple
repair this is exactly where n>=12 is needed. All joins entering or leaving
that block use actual port information h->a and a->(h-1).

The anti-reversal witnesses have the advertised literal images. In both
new final words r(a)=6 and r(c)=2, so the known chord 2->6 in T is reversed
in B and contradicts a->c there. The triple repair's middle witness maps
2->4 to 6->5; neither pair is retained by Q1. Its first witness maps a->2
to 4->3. The single repair's two witnesses likewise compare c->a with
4->2 in T and z->b with 4->2 after Q1. None uses an unspecified middle arc.

**Theorem B.** The hypotheses min I>=5 and n-4,n-2 in I guarantee every
join and the first two blackness witnesses of the forward ledger. They
also imply h<=n-4, putting the final word within the rotation lemma. The
rotation argument needs only the complete a-row, not consecutive ports:
in each nontransitive a-orbit all vertices besides b lie below h and are
therefore incoming to a. If that orbit were an automorphism orbit, its
induced tournament would be regular, forcing orbit size three. In the
transitive case global regularity forces n=2k+3. The actual threshold
n>=4k-1 excludes all these possibilities for k>=4.

For the reverse-chord ledger, a shift by two is transitive at odd n. At
even n its even-label orbit has at most k-1 incoming neighbours of a,
because z is odd. Regularity gives n<=4k-2, excluded by the theorem.
The second move pairs outdegree k+1 with indegree k. The final reversal
product has an even cycle and cannot be a tournament automorphism: an
appropriate power of any even-order tournament automorphism would swap
the ends of an arc. No disjointness assumption on I and its reflected
image enters these arguments.

**Theorem C.** The cardinality of I is k-1 and it contains n-4,n-2 for
k>=5. Its first port h=n-k-1>=k+4 keeps the insertion predecessor outside
the removed block. The forward branch excludes n=2k+3 directly. For even
n in the reverse branch, the number of even incoming ports is k/2-1 when
k is even, and (k+1)/2 when k is odd. Regularity would require respectively
n=2k-2 or n=2k+4; odd n would require n=2k+1. Each falls strictly below
n>=2k+5. This is a uniform counting argument, not extrapolation of tests.

## Independent literal controls and sharp failure mechanisms

An independent adjacency-matrix constructor propagates only the source
anti-reflection, the directed base path, the exact incoming a-row, and
the branch guards. It fills remaining reflection orbits with seeded bits.
It then checks all actual directed joins, all positional anti-reflections,
the literal final matrix, and the unordered-edge XOR, without an affine
union-find implementation or any isomorphism/clone contraction.

The independent run passed 72 sources for A, 1,034 for B, and 64 for C.
The B universe includes every admissible port set at n=4k-1 for k=4,5,6
and both chord directions. Among its sources, 356 have overlapping port
and reflected-port sets. These tests corroborate rather than establish
the unbounded proofs.

Both intended hostile controls fail for the stated reason. At n=11 the
triple repair has a false final join 5->0; at h=4 the protected-port
repair has a false final join 3->0. Both also fail the literal target and
path-edge XOR, even though all three words are black in the corresponding
formal intermediate matrices. Thus blackness alone would have concealed
the error. These are ledger failures, not connectivity counterexamples.

The independent controls below are executable from this directory with:

    python -c "from pathlib import Path; t=Path('INDEPENDENT_REVIEW.md').read_text(); exec(t.split(chr(96)*3+'python\n')[1].split(chr(96)*3)[0])"

```python
from itertools import combinations
from collections import Counter
from random import Random
import json

N=Counter()
def require(x,m):
    if not x: raise RuntimeError(m)
def make(n,I,guards,seed):
    A=[[None]*n for _ in range(n)]
    def put(a,b,v):
        for u,w in [(a,b),(n-1-b,n-1-a)]:
            require(A[u][w] is None or A[u][w]==v,('inconsistent',n,I,guards,u,w))
            A[u][w]=v; A[w][u]=1-v
    for a in range(n-1): put(a,a+1,1)
    for v in range(1,n): put(v,0,int(v in I or v==n-1))
    for a,b in guards: put(a,b,1)
    rng=Random(seed)
    for a,b in combinations(range(n),2):
        if A[a][b] is None: put(a,b,rng.randrange(2))
    return A
def anti(A,q):
    r=dict(zip(q,q[::-1]));n=len(q)
    return all(A[u][v]==A[r[v]][r[u]] for u,v in combinations(range(n),2))
def edges(q): return {frozenset((u,v)) for u,v in zip(q,q[1:])}
def kmove(A,q):
    keep=edges(q);n=len(q)
    return [[None if u==v else (A[u][v] if frozenset((u,v)) in keep else 1-A[u][v])
             for v in range(n)] for u in range(n)]
def check(n,I,qs,guards,seed,negative=False):
    A=make(n,I,guards,seed);P=list(range(n));original=A
    require(anti(A,P),'source not blue')
    require(all(A[u][u+1] for u in range(n-1)),'source not directed')
    require(A[n-1][0]==1,'source not closed')
    missing=[];blue=[]
    for j,q in enumerate(qs):
        require(sorted(q)==P,('word failure',q))
        missing += [(j,u,v) for u,v in zip(q,q[1:]) if not A[u][v]]
        blue.append(anti(A,q));A=kmove(A,q)
    same=A==kmove(original,P)
    xor=edges(qs[0])^edges(qs[1])^edges(qs[2])==edges(P)
    if not negative:
        require(not missing and not any(blue) and same and xor,
                (n,I,guards,missing,blue,same,xor))
    return {'missing':missing,'blue':blue,'same_target':same,'xor':xor}
def ledger(n,l):
    D=list(range(3,l+3));U=list(range(l+3,n-1));z=n-1
    return [[0,2]+D+U+[z,1],D+[0,2]+U+[z,1],D+[0,1,2]+U+[z]]
def a_words(n,b):
    h=n-5;R=list(range(h,n))
    if b==0:return ledger(n,1),[(4,2)]
    if b==1:return ledger(n,3),[(2,4),(6,2)]
    if b==2:
        V=list(range(5,h))
        return [[0,1,2,3]+V+R+[4],[1,2,3]+V+[0]+R+[4],
                [1,2,3,4]+V+[0]+R],[(2,4),(2,6),(3,5)]
    U=list(range(6,h))
    return [[0,1,2]+U+R+[4,5,3],[1,2]+U+[0]+R+[4,5,3],
            [1,2,3,4,5]+U+[0]+R],[(2,4),(2,6),(5,3)]
def bc_words(n,h,f):
    if not f:return ledger(n,1),[(4,2)]
    U=list(range(4,h));R=list(range(h,n))
    return [[0,1,2]+U+R+[3],[1,2]+U+[0]+R+[3],
            [1,2,3]+U+[0]+R],[(2,4)]
for n in [12,13,14,17,20,31]:
    for b in range(4):
        qs,g=a_words(n,b)
        for seed in [1,59,233]:
            check(n,{n-5,n-3,n-2},qs,g,seed);N['A_literal_sources']+=1
for k in [4,5,6]:
    n=4*k-1;required={n-4,n-2}
    for extra in combinations([i for i in range(5,n-1) if i not in required],k-3):
        I=required|set(extra)
        for f in [False,True]:
            qs,g=bc_words(n,min(I),f)
            check(n,I,qs,g,30030+sum(extra));N['B_literal_sources']+=1
            if I & {n-1-i for i in I}:N['B_reflection_overlap_sources']+=1
for k in range(5,13):
    for n in [2*k+5,2*k+6,2*k+7,3*k+4]:
        I={n-k-1}|set(range(n-k+1,n-1))
        for f in [False,True]:
            qs,g=bc_words(n,min(I),f)
            check(n,I,qs,g,17+k);N['C_literal_sources']+=1
qs,g=a_words(11,3);h11=check(11,{6,8,9},qs,g,0,True)
qs,g=bc_words(15,4,True);h4=check(15,{4,11,13},qs,g,0,True)
require((2,5,0) in h11['missing'] and not h11['same_target'],'missing n11 hostile')
require((2,3,0) in h4['missing'] and not h4['same_target'],'missing h4 hostile')
print(json.dumps({'status':'PASS','counts':dict(N),'n11':h11,'h4':h4},indent=2))
```
