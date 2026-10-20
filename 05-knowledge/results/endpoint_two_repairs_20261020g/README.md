# Three-black detours for the complete endpoint-two layer

Research continuation, October 20, 2026. Inherited base `b9cc57066630`.
This packet combines the two-block and cyclic-middle repairs, preserving
the actual exceptional endpoint arcs. No claim of literature novelty is made.

## Main result

**PROVED, with a FINITE-EXACT head at orders 6,7,8.** Let n>=6, and
let P=(a,M,z) be a closed blue Hamiltonian word in a tournament T, with

    d_T(a)=n-3,                 d_T(z)=2.

Then there are three successive directed **black** Hamiltonian moves
whose literal labelled target is K_P(T). At every n>=9 this follows
from the structural case partition below, with arbitrary interior arcs;
orders 6,7,8 are covered by the exact finite certificates. No additional
strongness or mixedness filter on the target is needed.

Consequently every such marked blue line has its two endpoints in the
same black component, and cannot be a bridge in the literal coloured
metagraph. Reading the line from its open K partner gives the equivalent
endpoint score pattern (3,n-4). These statements concern marked lines,
not arbitrary tournaments sharing a degree sequence.

**OPEN:** black-core connectivity at all sizes and arbitrary strong-to-
strong mixed blue lines. This theorem eliminates an entire endpoint
layer; it does not establish that every remaining line admits such a
presentation. No genuine bridge was found in the layer proved here.

## The retained coordinate and the exhaustive split

Here K_Q(T)=T^op xor E(Q). A word is blue iff its positional reversal
is an anti-automorphism. All ledger identities below concern labelled
arcs, before any isomorphism or converse quotient.

Because P is closed, z->a. The degree pair gives exactly one additional
middle vertex u with u->a, and exactly one middle vertex v with z->v.
Blue reflection forces v=r_P(u); every other middle vertex satisfies
a->w->z. The two ports may coincide at an odd-order middle vertex.
Write M=(b,x,...,y,c). Open middle means b->c; closed middle means
c->b. It does not mean that the induced tournament is nonstrong or strong.

The following disjoint cases exhaust all possibilities at n>=9.

| Port/middle condition | Black detour proof | Mechanism |
|---|---|---|
| u=c, x->c | [L2](last_port.md) | Two-block endpoint repair |
| u=c, c->x | [Theorem A](lift.md) | Transport an entire consecutive block |
| u!=c, closed middle, c->x | [Section 3](theory.md) | Cyclic no-chord ledger |
| u!=c,x, closed middle, x->c | [Section 4](theory.md) | Anti-product and two-neighbor obstruction |
| u=x, closed middle, x->c | [Section 5](theory.md) | Three fixed arc contradictions |
| u=y, open middle | [Theorem B](probe.md) | Move the block ending at the exceptional port |
| u not in {y,c}, open middle | [Open-middle lemma](open_middle.md) | Guarded cyclic ledger plus the actual port |

Every row supplies exactly three black K moves to K_P(T), not a
relocation inside one path fibre. Their guards are discharged by the
row hypotheses. The whole-block constructions remove the places where
the previous singleton or prescribed-endpoint constructions could fail.

## A representative repair: the penultimate port

Let P=(a,b,x,D,c,z), with D ending at u=y, and suppose b->c. Then
z beats exactly a,x. The three K words are

    Q1=(a,b,c,z,x,D),
    Q2=(b,c,a,z,x,D),
    Q3=(b,x,D,c,a,z).

The boundary pairs cancel to E(P), and the actual endpoint arcs make
each word directed. Q1 is black because its endpoint u beats a,c,z,
so its endpoint score sum is at least n. At S=K_Q1(T), the only
inneighbors of z are a,c. For n>=8, Q2 reversal sends z to the P vertex
at position n-6 and a to that at position n-5. The former beats a in S;
the latter cannot beat z. This contradicts anti-reversal. At n=6,7,
the endpoint score excess is at most n-8<0 instead.

At the final target B=K_P(T), Q3 reversal would send a->b to z->x.
The first arc is retained from P, whereas the second is false because
the original extra port z->x was complemented. Thus Q3 is black.
This is the exact repair of the tempting but false attempt to carry a
singleton reversal identity unchanged through a larger block.

## Verification and finite boundaries

The proofs were read across the three collaborating lanes and audited
against the literal positional reversals. The complete partition was
then implemented independently using directed arc sets and bit rows.

- [Finite head](probe_coverage.json): all 12,64,320 marked presentations
  at n=6,7,8 have literal three-black certificates. Five presentations
  left by the generic candidate bank have explicit certificates in the
  receipt. Every endpoint is independently checked strong and mixed;
  none of these presentations is discarded by that filter.
- [Unified arc-set replay](unified_check.json) and the
  [separate bit-row audit](theory_unified_audit.json): all 3,072 and
  28,672 presentations at n=9,10, with no strongness filter. This includes
  all 512 coincident-port presentations at n=9. Every move is directed,
  actually black, and reaches the exact labelled target.
- The arc-set replay also checks 630 deterministic presentations at
  orders 11..100. These are controls of the formulas, not the source of
  the all-order quantifier. A genuine n=8 failure of one fixed-tail word
  is retained as a hostile to extending that word's scope without proof.

Reproduction from the repository root:

```
python 05-knowledge/results/endpoint_two_repairs_20261020g/probe_coverage.py
python 05-knowledge/results/endpoint_two_repairs_20261020g/unified_check.py
python -O 05-knowledge/results/endpoint_two_repairs_20261020g/unified_check.py
python 05-knowledge/results/endpoint_two_repairs_20261020g/theory_unified_audit.py
```

The subordinate proof notes give their own checks and exact scope.
The large exploratory probe receipts are regenerable and are not needed
as proof dependencies.

## Inheritance, useful failure, and the next precise frontier

The inherited [endpoint-one theorem](../kernels_games_detours_20261020f/detours.md)
has a genuine cyclic middle module. The inherited
[two-block theorem](../strong_blue_games_20261020f/detours.md) requires a
real endpoint-preserving path after reversal, which can fail for an
antitransitive suffix. The new block transport handles that failure
without assuming that the middle is a module. The canonical hostile is
retained in [lift.md](lift.md): an obvious module expansion has a legal
macro ledger but the wrong physical ports, and even an illegal micro join.

The exact surviving substitution theorem retains **oriented port parity
at every prefix**, as well as the final target parity. Its source is a
macro K ledger; its target is a block expansion. It preserves directed
legality and the labelled endpoint iff the stated parity conditions hold.
It does not certify blackness; that needs the separate anti-involution
sidecar. This explains why a macro-edge cancellation alone was inadequate.

For a general closed blue word with scores (n-1-k,k), define its
**marked endpoint port defect** to be k-1. Its first endpoint loses to z
and a set U of k-1 middle vertices, while z beats a and r_P(U). This
coordinate retains the source, target and all exceptional incidences.
The previous cyclic-module theorem covers defect zero; the present
theorem covers defect one. The next direct target is defect two, with
scores (n-4,3) and a reflected pair of two-element port sets. These are
marked-path defects, not the black-orbit deficit or counting deficiency
used elsewhere in the repository. Coinciding numbers do not identify
those different invariants.
