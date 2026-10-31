# Independent audit of the positive G32 origin

**PASS — definitions and elementary claims PROVED; small-order classifications
FINITE-EXACT; Gerbicz's all-order Hamiltonicity CITED.** No mathematical correction
is required in the reviewed `origin_model/ORIGIN_MODEL.md` or `check_model.py`.
The original input hashes are in `origin_model_audit_inputs.json`.

## Typed origin and arithmetic

The original vertex label `v`, graph order `n`, and arithmetic rank `r` are
three different types. For either base `b=31` or `b=32`, the bijection
`r -> G_(b+r)` transports ordinary nonnegative-integer addition and
multiplication to `G_a +_b G_c = G_(a+c-b)` and
`G_a *_b G_c = G_(b+(a-b)(c-b))`. Therefore all semiring axioms follow at
once by conjugation; these operations are not graph union, gluing or join.
For the requested primary positive-origin convention, `G31=0` and `G32=1`.
The alternative convention `G32=0` is separately and correctly identified.
Neither introduces a new number field, negative-integer graph objects, or
changes the original vertex labels 1, 3 and 8.

## The threshold is neither a monotonicity argument nor a digit argument

A simple cycle requires at least three vertices. For orders 3..7, 8..16,
17, and 18..30, the respective degree-one witnesses 1, 8, 16, and 18 rule
out a Hamiltonian cycle. The independent auditor freshly enumerates **all**
spanning degree-two edge sets at 31, 32 and 33 without connectedness pruning.
At 31 their component sizes are exactly 11+20 and 14+17. At 32 they are
3+6+23, 3+29, 5+27, 9+23 and 32. This verifies the finite boundary and
unique G32 cycle rather than assuming Hamiltonicity is downward monotone.

The all-order implication `n>=32 -> G_n Hamiltonian` remains Gerbicz's
classical theorem, not a consequence of this finite search. Its author
statement was rechecked on [OEIS A090461](https://oeis.org/A090461), dated
2018-01-21. The pinned source-data reconstruction and its full induction
are in the separate Session124 Gerbicz packet. This lane reads that proof
and provenance; it does not claim a fresh rerun of its 17 MB external input.
The direct bibliography points to [Berlekamp–Guy's appendix](https://arxiv.org/html/2002.03705).

## Arrival and semigroup proofs

An old label adjacent to arriving `N` has the unique form `v=k²-N`.
The exact condition `1<=v<=N-1` is equivalent to `N<k²<=2N-1`, giving
`floor(sqrt(N))<k<=floor(sqrt(2N-1))`. Thus arrival degree is exactly
`floor(sqrt(2N-1))-floor(sqrt(N))`, including N=1 and all square boundaries.
The independent direct-label check covers N=1..4000.

With ordinary order addition, `S={0} union [32,infinity)` is closed.
Every 32..63 is indecomposable because two positive elements sum to at
least 64. Every larger integer reduces by 32 to this interval. Hence these
are precisely the minimal generators, conductor 32, Frobenius number 31,
and genus 31. The Apéry set relative to 32 is `{0,33,...,63}`: its elements
are the least members in the 32 residue classes. This proves uniqueness of
`n=32q+w` with `q>=0` and `w` in that set. It does **not** imply arbitrary
factorization uniqueness: `66=32+34=33+33`.

For G32->G33, arrival ports 3,16,31 have no square-sum pair. One deleted
old cycle edge would force pure subdivision by 33, so at least two are
necessary. The freshly reconstructed unique cycles realize exactly two
old deletions `{5,31}` and `{16,20}`. This is distinct from merely adding
vertex 33 to the nested underlying graph.

## Separate all-order nice-pair corollary

**PROVED.** For every original Gerbicz nice pair, closing the shorter and
longer paths yields cycles that differ by at least **two deletions of old
cycle edges**. This is not a claim about arbitrary consecutive-order cycles.

Let `delta=|E(C_n) minus E(C_(n+1))|`. The odd-order path starts at 1 and
ends at 3, so its closure contains edge 13. Both labels have odd positional
color. Common phase transfers that color to the even-order path, whose
cycle alternates positional colors; hence its closure cannot contain 13.
If delta=0 the new vertex could not have degree two without raising an old
degree. If delta=1, deleting edge ab makes only a and b deficient; adding
the new vertex requires exactly two new edges, necessarily az and bz.
Thus delta=1 is precisely pure subdivision and creates no new old-old edge.

For even shorter -> odd longer, pure subdivision cannot create old edge 13.
For odd shorter -> even longer, the sole deleted edge would have to be 13,
but no positive integer z is adjacent to both: `(z+3)-(z+1)=2` cannot be
a difference of two integer squares. Thus delta>=2 in both orientations.
The certified G43->G44 pure insertion at edge 5--20 is a hostile scope
control: it exists for general cycles and is outside this nice-pair class.

## Verification and limits

The original model checker passes in normal and -O modes with 221,649
checks each and byte-identical JSON. The independent combined audit passes
in both modes, with degree-two enumeration, direct-label arrivals,
transported-semiring, Apéry and strict-scope controls. All checks use explicit
exceptions, not `assert`. Receipt counts in the combined audit are current.
No finite range is promoted to an all-order theorem. The forest bijection,
CG digits, and scalar forest counts do not certify Hamiltonian continuation.
No source-worktree or remote write was performed; there was no access denial.
