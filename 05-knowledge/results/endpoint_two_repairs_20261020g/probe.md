# Endpoint-two coordinates, a uniform open-middle repair, and finite certificates

October 20, 2026. **PROVED**: the open-middle penultimate-port repair in
section 2, for every order at least six. **FINITE-EXACT**: complete
standard-path endpoint-two coverage at orders six through nine, in section 3.
These statements concern an actual marked blue edge and its literal labelled
target. They do not assert general black-core connectivity.

## 1. Inheritance and the exact marked interface

The starting mechanisms are the two-block ledger in
[strong_blue_games_20261020f/detours.md](../strong_blue_games_20261020f/detours.md)
and the cyclic-module ledgers in
[kernels_games_detours_20261020f/detours.md](../kernels_games_detours_20261020f/detours.md).
The canonical hostile is a reversed suffix with no Hamiltonian path having
the required endpoints. Failure of this one block repair says nothing about
whether another black detour exists. The useful sidecars are the marked
exceptional endpoint ports and the actual positional reversal of each word.

Let P=(a,M,z) be a closed blue Hamiltonian word in a tournament T of order n,
with d(a)=n-3 and d(z)=2. Blue means that the reversal of P is an
anti-automorphism. There is exactly one middle vertex u with u->a, and
exactly one middle vertex v with z->v. The blue symmetry forces v=r_P(u).
Every other middle vertex satisfies a->w->z. Thus every such source is
constructed uniquely from a blue middle word by putting it into the cyclic
wrapper a->M->z->a and reversing the pairs au and vz.

With zero-based standard P=(0,...,n-1), the allowed mark is
u=2,...,n-2 and v=n-1-u. The equality u=v occurs exactly at the central
position of odd order. When they are distinct, the induced four vertices
{a,u,v,z} form TT4, in order (u,z,a,v), if u->v, and form R4 otherwise.
This is an exact interface description; degrees alone forget the mark and
the middle chords needed by a particular ledger.

The concept board is: marked ports; path-edge XOR; actual reversal witnesses;
endpoint score sums; and the distinction between a local template and an
exhaustive search. The anchor is this endpoint layer; the niche is a reusable
symbolic ledger compiler; the wildcard is recovering arbitrary block moves
from the small hostile presentations.

The compiler `compile_ledger` in [probe.py](probe.py) produces initial-arc
guards, checks their compatibility, computes the path-edge XOR, and derives
the affine endpoint score sum at each successive source. It never replaces
an anti-automorphism test by the merely necessary degree condition. The
additional [probe_symbolic.py](probe_symbolic.py) tests compatibility of a
proposed blue step with all source and ledger constraints over GF(2).

## 2. An arbitrary-block repair for an open middle

**PROVED.** Suppose n>=6, the middle is open, and its penultimate vertex is
the extra incoming port. Write

```
P=(a,b,x,D,c,z),       u=last(D),       v=x,
```

where D is a nonempty word. The open-middle condition is b->c. Then the
following three K moves are directed, black, and end exactly at K_P(T):

```
Q1=(a,b,c,z,x,D),
Q2=(b,c,a,z,x,D),
Q3=(b,x,D,c,a,z).
```

Here K_Q complements every pair except consecutive pairs of Q. No module
assumption is made about D, and no strongness assumption on K_P(T) is needed.

**Directedness and target.** Q1 uses b->c and z->x together with displayed
path arcs. Its move complements a->c and z->a, supplying c->a and a->z for
Q2; all remaining Q2 joins are retained. In the third source, b->x and
u->c have been complemented twice, while c->a and a->z were retained by
Q2. The internal D path is retained throughout. All three words are thus
directed. Their added boundary pairs bc,zx,ca,az occur twice. The surviving
XOR pairs are exactly ab,bx,x first(D), the internal D pairs,uc,cz, namely
E(P). Three K moves therefore reach the literal labelled target.

**Q1 is black.** The vertex u beats a,c,z, so d_T(u)>=3. The endpoint
score excess of Q1 over n-1 is d_T(u)-2>=1. An anti-automorphism would
force this excess to vanish.

**Q2 is black.** Put S=K_Q1(T). Directly from the marked ports and the
two Q1 edges incident to z,

```
N^-_S(z)={a,c}.
```

For n>=8 let r be the positional reversal of Q2. The vertex r(z) is x
at n=8 and a nonterminal vertex of D at larger orders. Hence a->r(z) in
T, and this pair is absent from Q1, so r(z)->a in S. Meanwhile r(a)
lies in D and is neither a nor c. Consequently r(a)->z is false in S.
The anti-automorphism equation would equate these two arcs, a contradiction.

For n=6,7, the vertices b and u each have at least three outneighbors:
b beats x,c,z, and u beats a,c,z. The endpoint excess of Q2 in S is

```
n-2-d_T(b)-d_T(u) <= n-8 < 0,
```

again excluding blue.

**Q3 is black.** Evaluate its reversal at the other endpoint B=K_P(T);
blackness is invariant across its own K move. The reversal of Q3 sends
a to x and b to z. In B the path arc a->b is retained, whereas z->x is
false, since the original exceptional pair z->x is not a P pair and was
complemented. These are the two arcs an anti-automorphism would equate.
This proves all three blackness assertions.

**Correction lineage.** A preliminary extension of the singleton-D ledger
incorrectly kept the same Q2 reversal witness after D was lengthened.
Its first failed implication was the asserted image of x under that
reversal: it changes with the block length. The proof above retains the
actual positions and instead uses the two-element set N^-_S(z). The
singleton-D witness remains valid only in its original singleton scope.

## 3. Complete finite coverage and the five residual certificates

[probe_coverage.py](probe_coverage.py) constructs every standard closed blue
endpoint-two presentation at n=6,7,8,9 from the exact middle-plus-mark
parameterization. It explicitly filters for strong mixed endpoints, without
an isomorphism or converse quotient. Every presentation in these orders
passes those filters. The generic bank consists of the present repair,
the last-port two-block repair in [last_port.md](last_port.md), the
last-port block transport, and the cyclic and second-port ledgers in
[theory.md](theory.md). At each order it tries actual directedness and
actual anti-reversal failure, including at orders below a template's
all-order proof boundary.

| Order | Complete marked universe | Template covered | Finite exceptions |
|---:|---:|---:|---:|
|6|12|9|3|
|7|64|63|1|
|8|320|319|1|
|9|3,072|3,072|0|

All five exceptions have literal three-black K certificates. They are
stored with source and target rows, every path word, explicit failed
anti equations, and full three-K search counts in
[probe_coverage.json](probe_coverage.json). In particular, the sole order-eight
last-port hostile is

```
T=(62,28,184,240,224,194,135,3),
Q1=(0,2,3,4,5,6,7,1),
Q2=(7,2,1,0,6,3,4,5),
Q3=(1,7,2,0,6,3,4,5).
```

Its three endpoint excesses are (1,2,-1), so each move is already certified
black by degrees. The complete exact search finds 252 three-K certificates
for this source; existence, rather than that count, is the consequence used.
The order-seven exception has coincident ports u=v=3, so the diagonal
interface is explicitly included rather than silently discarded.

The search in [probe_search.py](probe_search.py) enumerates every black first
word at T and every black last word at B=K_P(T). Their path-edge XOR with
E(P) determines the only possible middle undirected path. It tests both
orientations if necessary and replays each candidate. This is complete
for labelled three-K detours, not for all longer detours or quotient paths.
The finite coverage checker then independently replays every selected
certificate using sets of directed arcs rather than the producer's bitmask
K routine. All failures use explicit exceptions, so optimization cannot
remove a check.

Reproduction:

```
python 05-knowledge/results/endpoint_two_repairs_20261020g/probe_coverage.py
python -O 05-knowledge/results/endpoint_two_repairs_20261020g/probe_coverage.py
```

Both runs cover 3,468 presentations and reproduce the same compact JSON.
The initial [probe.py](probe.py) and its optional generated `probe.json`
preserve the old-bank diagnostic: nine residual presentations at orders
five through eight. That larger provenance output is regenerable and is
not an input to the final coverage check. One of the two strong mixed
order-five presentations has no labelled three-K certificate, so the
present exact order-six boundary must not be silently lowered.

## 4. Information preserved and scope

The source is a marked blue Hamiltonian word, the target is its exact K
partner, and the map is the displayed ordered triple of genuine black K
moves. The preserved predicate is literal equality of the final labelled
tournament with K_P(T); no relabelling is used to mask a failed XOR.
Forgetting the actual word destroys the port coordinates and the local
anti-automorphism witnesses. Thus this repair does not identify a degree
sequence with a metagraph vertex, and it does not infer that all black
components coincide. The cheap decisive tests are the displayed arc guards,
path-edge XOR, and the three actual reversal equations.
