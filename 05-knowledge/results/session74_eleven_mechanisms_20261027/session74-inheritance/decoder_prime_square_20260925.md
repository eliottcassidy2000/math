# Halving needs a square-class bit; the 24-to-25 transition has an endpoint-aware certificate

**Status.** PROVED (elementary proofs below; the `Q_24` argument independently
audited by the coordinating session's minor/decoder lane): the two-sheet halving identity; the dyadic layer obstruction; the
endpoint-aware nonexistence proof for `Q_24`; every Hamiltonian path of
`Q_25` uses `11--25--24`. FINITE-EXACT: independent path census for
`14 <= N <= 25`, sound endpoint-conditional certificates for all 23 possible
second endpoints of `Q_24`, and local arithmetic controls. INHERITED PROVED:
the prime doubling escape set is exactly `{2,3,11}`. No Collatz convergence
or graceful-tree implication is claimed. No novelty claim.

Session: collatz decoder, 2026-09-25, prime/square side lane. This result
uses no new literature import. The historical notes cited below contain
their own literature boundaries; this lane re-proves the claims it uses.

## Inheritance, portfolio, and board

Closest mechanism: the square-sum `49a+c,49b-c` blow-up from
[pairing transitions, section 3](procgen_brackets_20260924_pairings_transitions.md)
preserves the arithmetic edge predicate by multiplying sums by a square.
Its prime escape set is proved in section 1. The corrected near miss is the
unguarded degree-two forcing in
[square-sum Hamiltonicity, section 4](collatz_mod6_20260922_w6_square_sum_hamiltonicity.md):
graph-degree two does not certify a path-interior vertex. Canonical hostile
controls are a triangle path and the valid `Q_23` path with endpoint `22`,
which omits its available edge `22--14`. Least-used sidecar: endpoint
identity, with the nice-pair parity invariant as the inherited example of
what a valid scaling induction must remember.

Anchor: understand what halving preserves. Niche: repair the short `Q_24`
obstruction left OPEN after the correction. Wildcard: the prime `{2,3,11}`
bracket statistic. The live concept board is:

| Concept | Retained information | Cheap decisive test |
|---|---|---|
| Odd chain `q -> q+2` | ordered odd bases | does translation preserve sums? |
| Dyadic columns `2^k q` | odd base and valuation | square-class toggle under `/2` |
| Prime brackets | real scale and prime tag | doubling escape inequality |
| Square-sum paths | endpoint and coverage | `Q_23` endpoint `22` |
| Structural contraction | attachment locations | suppress new vertex `25` |

After the halving calculation the board changes: the square class is an
exact extra coordinate, while odd-base order alone is not. After the
endpoint proof it changes again: `25` repairs a specified attachment
obstruction; the decimal coincidence with a bracket endpoint is not what
proves Hamiltonicity.

## 1. An exact two-sheet halving law

Let `Q_N^(c)` have vertices `1,...,N`, no loops, and an edge `x--y` iff
`x+y=c s^2` for a positive integer `s`, where `c` is `1` or `2`.
Then

    Q_(2N)^(1)[even vertices] / 2  =  Q_N^(2),
    Q_(2N)^(2)[even vertices] / 2  =  Q_N^(1).               (H)

These are exact labeled graph isomorphisms: send `2x` to `x`. For the first
identity, `2(x+y)=s^2` forces `s=2t`, so `x+y=2t^2`. For the second,
`2(x+y)=2s^2` is precisely `x+y=s^2`. Distinct endpoints stay distinct.
Iterating gives

    Q_(4N)^(1)[multiples of 4] / 4 = Q_N^(1).

Thus doubling toggles the square-class bit; quadrupling restores it.
The smallest useful hostile is `1+3=4`, whereas `2+6=8` is not a square;
the positive control is `4+12=16`.

**Scope.** (H) preserves adjacency exactly on the induced even subgraph.
It discards odd vertices, their incident edges, and any path coverage or
endpoint obligation involving them. Consequently it is not a halving
theorem for Hamiltonian paths. A decoder must retain both the arithmetic
sheet and the attachments across the deleted layer. In particular, it is
not a map from a square-sum certificate to a Collatz certificate.

### Dyadic layer obstruction

Write `x=2^a u`, `y=2^b v`, with `u,v` odd. If `a<b`, then
`v_2(x+y)=a`. Hence a square-sum edge between distinct layers is possible
only when the smaller valuation is even. This condition is necessary,
not sufficient: `1+2=3` is the smallest hostile to sufficiency.

In particular, the seam `v_2=1` has no square-sum edges to any deeper
layer `v_2>=2`. The next allowed cross-layer minimum is `v_2=2`, which
is an exact copy of the odd-base condition after dividing by four. For
equal valuations `a=b`, the exact condition is

    a+v_2(u+v) is even, and
    (u+v)/2^v_2(u+v) is an odd square.

This explains why grouping all even numbers into one layer loses a real
constraint. The odd chain `q -> q+2` and the vertical doubling edges give
every positive integer the coordinate `(q,k)`; primality occurs only on
the odd-base axis except for `2`, while square-sum adjacency depends on
both coordinates. Translating both endpoints by two changes a sum by
four and does not generally preserve squareness (`1+3=4`, but `3+5=8`).

## 2. The prime statement is a scale theorem with an exact scope

Let `B_m=((2m-1)^2,(2m+1)^2]`, with `{1}` separate. A positive integer
`n in B_m` has a nontrivial integer multiple in the same bracket iff
`2n <= (2m+1)^2`. For `m>=3`, even the smallest possible integer fails:

    2((2m-1)^2+1) - (2m+1)^2 = 4m^2-12m+3 > 0.

The expression is `3` at `m=3` and increases thereafter. Checking the
first two brackets gives exactly

    n = 2,3,4,10,11,12,
    prime n = 2,3,11.

The phrase "a multiple below the next odd square" must mean a
**nontrivial multiple above the same bracket's lower boundary**; allowing
the multiple `p` itself makes every prime qualify. With the bracket
defined above, doubling is the decisive test and the classification is
complete, not a finite extrapolation.

This property is a comparison of multiplicative scale `2n` against an
additive square-gap boundary. It contains no edge direction, residue word,
or orbit certificate. The full admissible pairs are inherited:
`(k,p)=(2,2),(3,2),(4,2),(2,3),(3,3),(2,11)` under the right-closed
bracket convention. If "less than the next odd square" is strict, remove
`(3,3)` because `3*3=9`; the prime classification `{2,3,11}` is unchanged,
since each still has its double strictly below the upper square.

## 3. Why vertex 18 remains a leaf until 31

The square-sum neighbors of `18` are `s^2-18`, with positive results and
`s^2 != 36`. The first candidates are `25-18=7`, `36-18=18` (a forbidden
self-loop), and `49-18=31`. Thus it has exactly one neighbor through
`Q_30`, and its second neighbor appears at `Q_31`.

This is one instance of a typed diagonal hole: `x=2r^2` would have the
self-loop sum `2x=(2r)^2`. The square immediately before and after that
loop produces partners `2r^2-4r+1` and `2r^2+4r+1`, separated by `8r`.
Such a self-loop vertex has odd valuation `v_2(x)=1+2v_2(r)`, so it lives
on the first seam or its every-other-layer replicas. Larger `r` may have
earlier partners too; the identity does not imply they are leaves.

The distinction between a local arithmetic hole and an actual path
obstruction matters: `18` is an endpoint of the positive `Q_23` and
`Q_25` witnesses as well as an obstruction contributor at other sizes.

## 4. A short, endpoint-aware proof that Q_24 has no Hamiltonian path

The proof uses only displayed adjacency and certified nonendpoints.
Suppose `P` were a Hamiltonian path. Vertex `18` has only neighbor `7`,
so one endpoint is `18`; write `e` for the other endpoint.

**Endpoint gate.** If `e != 9`, degree-two vertex `9` is interior and
forces `7--9`. Together with `7--18`, this saturates vertex `7`, excluding
`2--7`. If `e` is neither `11` nor `22`, both degree-two vertices `11`
and `22` are interior, forcing `14--11` and `14--22`; this saturates `14`,
excluding `2--14`. Vertex `2` has exactly the neighbors `7,14,23`.
Therefore, if `e` were outside `{2,9,11,22}`, interior vertex `2` would
have at most the single available edge `2--23`, impossible. Thus

    e in {2,9,11,22}.                                      (E)

**Vertex 5 is saturated in every remaining case.** Interior vertex `20`
has neighbors `5,16`, so `5--20` is forced. If `e != 11`, interior
vertex `11` forces `5--11` too. If `e=11`, then `7` is saturated as
above, and interior vertex `2` must use `2--23` and `2--14`. Interior
vertex `22` forces `22--14`, saturating `14`, so edge `14--11` is
excluded. Endpoint `11` must consequently use its other edge, `11--5`.
Thus in all four cases `5` uses `20` and `11`, excluding `4--5`.

**The forced proper cycle.** Consider

    1--8--17--19--6--10--15--21--4--12--24--1.              (C)

Every vertex of this cycle is an interior vertex by (E). Its degree-two
vertices `8,17,19,10,21,24` force every displayed edge except `4--12`.
Vertex `4` has exactly the neighbors `5,12,21`; the first edge is already
excluded. Because `4` is interior, it must use `4--12` and `4--21`.
The entire eleven-vertex cycle (C) is therefore a subgraph of the path
`P`, impossible. This proves nonexistence without exhaustive path search.

**Correction lineage.** This is a repaired proof of the result whose old
short forcing proof was retracted on 2026-09-25. It does not rehabilitate
the old rule: every use of degree-two forcing above comes after a specific
endpoint exclusion. The old minimal hostile and `Q_23` endpoint hostile
remain valid and are checked by the new code. No new demonstrated mistake
was found in the corrected note; its OPEN short-proof slot is now resolved.

## 5. Vertex 25 supplies the missing attachment, with a real decoder

In `Q_25`, the new vertex `25` has exactly neighbors `11,24`, using the
sums `36,49`. If it were an endpoint of a Hamiltonian path, deleting it
would leave a Hamiltonian path of `Q_24`, ruled out above. Hence every
Hamiltonian path of `Q_25` contains `11--25--24` (in one order).

One explicit witness is

    18,7,9,16,20,5,11,25,24,12,4,21,15,10,6,19,17,8,1,3,13,23,2,14,22.

All sums are squares, and it uses every integer `1..25` exactly once.
The edge `24--25`, sum `49`, is the first possible use of that square:
the largest distinct sum in `Q_24` is `47`. In particular the new square
`49`, rather than the square `25` itself, is compulsory in every first
post-exception path.

There is an exact contraction statement:

    Hamiltonian paths of Q_25
        <--> Hamiltonian paths of Q_24 plus edge {11,24}
              that use the new edge.

Suppress `25` and mark the resulting edge `11--24`; conversely subdivide
that marked edge with `25`. The arithmetic label on the suppressed edge
is a sidecar: `11+24=35` is not a square, so an unmarked contraction
leaves the square-sum category. This is a small, fully working model for
a structural decoder: preserve exactly what tells the inverse operation
where and how to restore a valid arithmetic path.

## 6. Reproduction, controls, and stopping boundary

Run from repository root:

    python 04-computation/experiments/decoder_prime_square_20260925.py

Recorded output:
[decoder_prime_square_20260925.out](decoder_prime_square_20260925.out).
The script uses exact integer arithmetic and explicit exceptions (not
optimization-removable assertions). Explicit universes:

* Both halving sheets for every `1<=N<=128`, all vertices and adjacencies.
* All square-sum edges with endpoints at most `256`; independent
  pair-isqrt and square-reflection adjacency constructors.
* Same-bracket doubling for integers through `10000`, with completeness
  separately proved by the inequality in section 2.
* Every possible second endpoint `e!=18` for `Q_24`, with a printed
  certificate using endpoint-conditioned degree equations and forced-cycle
  detection. All 23 cases close.
* Independent plain vertex-DFS for `14<=N<=25`, using no degree-two
  forcing; counts are `0,1,1,1,0,0,0,0,0,3,0,10`. The only DFS prune
  rejects a remaining vertex isolated from all remaining vertices plus
  the current endpoint. A fixed leaf starts each search, so reversal is
  not double-counted.
* Triangle and `Q_23` with endpoints `18,22` remain undecided by local
  propagation, and their actual Hamiltonian paths exist. These are hostile
  controls for the precise earlier false implication.

The strongest survivor is an exact square-class-aware halving law and a
marked-edge decoder for the `24 -> 25` repair. The open obligation is to
combine deleted-layer attachment data with a recursive path construction;
the square-class bit alone does not preserve Hamiltonicity, and none of
these statements supplies a Collatz rank or orbit route. No further large
search is justified by this lane's local signal.
