# Independent audit: multi-root grafts and fixed-menu support

Status: **PASS for the mathematical statements below**. The independent finite
checker supports these proofs but does not establish an unrestricted Chan–Pak
coverage theorem. Producer-prose reconciliation is recorded at the end.

## 1. Scope, inheritance, and main verdict

For nonempty tournaments A and B and a subset R of V(A), put C=V(A)\R.
The graft keeps both induced tournaments and puts every arc B→R and C→B.
No strongness hypothesis is needed for the positive expansion. The strongness
criterion, growth bound, and sparse-support conclusion have their own explicit
hypotheses below.

The general maximal-block-run substitution is already proved in
`01-canon/theorems/THM-3121-path-cover-walk-content-substitution-kernel.md`.
The path-colour transform and its order-join multiplication are already proved
in `01-canon/theorems/THM-3166-falling-factorial-order-join-path-colour-transform.md`.
Thus the positive expansion and multiexponential formula should be presented
as a specialization, not as a new general substitution or diagonalization
principle. The special degree bound, exact SCC test, sharp exponential growth,
and the fixed-menu support restriction are the useful consequences audited
here. No external novelty priority is claimed.

Concept board: literal graph counts; nonempty ordered path bins; SCC endpoints;
path-colour multiplication; arithmetic support under outer ordered joins.
Canonical hostiles are an unordered profile without factorials, empty bins in
the positive expansion, nonempty proper R without the SCC condition, and equal
H(B) with different graft counts. All are retained as exact test controls.

## 2. Independent coefficient proof

Let F_B(k) count spanning covers of B by k **ordered, nonempty** directed
paths. Put F_B(k)=0 outside 1≤k≤|B|. Erase the B vertices from a Hamilton path
of the graft, retaining a permutation p of V(A). A maximal B run can occupy:

- the initial slot if the first vertex of p belongs to R;
- the final slot if the last vertex belongs to C;
- an internal slot precisely between a C vertex and an R vertex.

Every backward adjacency of p must be broken by a B run; if that adjacency
is not of type C,R, the permutation is impossible. For a compatible
permutation, let b be its number of backward adjacencies, and let g be the
number of remaining eligible slots, including the two possible endpoint
slots. Its contribution to the coefficient c_k is binom(g,k−b). Therefore

    c_k = sum over compatible p of binom(g(p), k−b(p)),
    H(graft(A,R,B)) = sum_{k≥1} c_k F_B(k).

This is a bijection: select the slots and insert the ordered paths into them.
Runs are maximal because consecutive B runs cannot occupy the same slot.
There is no division by k!, and no empty bin may be inserted. Extending the
coefficient convention to k=0 gives c_0=H(A), but that coefficient contributes
nothing when B is nonempty.

Every B run except a possible final one is followed by a distinct R vertex;
every run except a possible initial one is preceded by a distinct C vertex.
Hence c_k=0 for k>min(|R|,|C|)+1. This is sharp for every positive pair
(|R|,|C|): choose a total order beginning in R and ending in C, with as many
C,R transitions as possible. Its forward Hamilton path has exactly
min(|R|,|C|)+1 available slots.

## 3. Strongness, positivity, and a sharp growth boundary

The graft is strong if and only if the source SCC of A meets R and the sink
SCC of A meets C. Necessity follows because B can enter A only through R
and leave A only through C. For sufficiency, every A vertex is reachable
from the chosen source-root, and every A vertex reaches the chosen sink-C
vertex. Consequently any B vertex can reach every A vertex, every A vertex
can reach any B vertex, and B vertices can reach one another through A.
Internal strongness of B is unnecessary.

Under this condition c_2≥1. If A is strong, cut a spanning directed cycle
at an arc C→R to obtain a Hamilton path beginning in R and ending in C.
If A has several SCCs, concatenate a Hamilton path starting at an R vertex
in its first SCC, arbitrary Hamilton paths through the intermediate SCCs,
and a Hamilton path ending at a C vertex in its last SCC. Directed Hamilton
cycles in nontrivial strong SCCs provide the endpoint choices; singleton
SCCs cause no issue. Inserting a B run at each endpoint proves c_2≥1.

Also c_1=H(graft(A,R,K1))≥3, since this quotient is a strong tournament on
at least three vertices and the cuts of a Hamilton cycle give at least
three Hamilton paths. For every nonempty B of order n,

    F_B(1)≥1,
    F_B(2)=sum_{∅≠S⊊V(B)} H(B[S])H(B[V(B)\S])≥2^n−2.

Therefore the sharp bound is

    H(graft(A,R,B)) ≥ 2^n+1.

Equality is attained by A=T2, R its source, and B=T_n: there c_1=3,
c_2=1, F_B(1)=1 and F_B(2)=2^n−2. The weaker lower bound 2^n−2 is also
correct and suffices for the asymptotic argument.

## 4. Path colours and a literal three-mode family

Set Q_B(j)=sum_k binom(j,k)F_B(k). This allows unused colours; it is not
the same object as F_B(j). Binomial inversion gives

    F_B(k)=sum_{j=0}^k (−1)^{k−j} binom(k,j)Q_B(j).

The inherited order-join law is Q_{B1⇒...⇒Bm}(j)=prod_i Q_{Bi}(j).
Combining the two identities proves the stated finite exponential sum and
shows that the scalar graft count forgets the order of the joined factors.
It does not say that their tournaments are isomorphic.

For A=T4 with 0→1→2→3 and R={0,2}, the independent coefficient engines give

    (c_0,c_1,c_2,c_3)=(1,9,10,2),
    H(graft)=−5Q_B(1)+4Q_B(2)+2Q_B(3).

For B=C3^{⇒m}, Q_C3(j)=j^3+2j, so

    H=−5·3^m+4·12^m+2·33^m.

The direct vertex-subset Hamilton-path computation gives 99, 2709, and 78651
for m=1,2,3, on 7,10,13 vertices respectively. All three grafts are strong.
The negative coefficient in the exponential basis causes no positivity issue:
the original expansion in nonempty ordered covers is positive.

## 5. Fixed-bank support, including arbitrary outer joins

Fix finitely many compatible (A,R) shells and a fixed finite menu of k≥1
nonempty factor tournaments. Let B be any nonempty order-join word in that
menu. With H≤X, the growth bound gives |B|≤log_2(X−1), hence its word length
is O(log X). The count depends only on its k multiplicities. Consequently
there are O((log X)^k) distinct scalar counts from these grafts. The constants
may depend on the finite shell bank and menu; neither may vary with X.

Let S be the set of distinct scalar values from these strong grafts. All
s∈S satisfy s≥3, and N(t)=#{s∈S:log s≤t}≤C t^k after increasing C. Allow
arbitrary finite outer order-joins of those grafts; their scalar counts belong
to the multiplicative monoid generated by S. Counting formal multisets is
an upper bound, even when distinct multisets have the same product.
For every σ>0,

    Z(σ)=prod_{s∈S}(1−s^(−σ))^(−1),
    log Z(σ)=sum_{l≥1} (1/l) sum_{s∈S} exp(−σ l log s).

By integration of the counting function,

    sum_{s∈S} exp(−u log s)
      = u integral_0^∞ exp(−ut) N(t) dt
      ≤ C Γ(k+1) u^(−k).

Thus log Z(σ)≤C Γ(k+1)ζ(k+1)σ^(−k). The number M(X) of attainable products
up to X is at most X^σ Z(σ). Taking σ=(log X)^(−1/(k+1)) proves

    M(X) ≤ exp(O((log X)^(k/(k+1)))) = X^o(1).

This excludes complete or positive-density odd-count coverage only for the
specified fixed, one-layer shell/menu class, even after arbitrary outer
joins. It does not exclude recursive use of outputs as new graft sources,
growing menus, growing shell banks, or other constructions. It is not a
Chan–Pak impossibility result. The asymptotic statement is proved above;
the finite checker is not being cited as a proof of it.

## 6. Full-state collisions and graph-level nonisomorphism

The producer's added collision construction is also sound. Let W contain k
copies of the strong directed triangle C3 and ell copies of the strong F4
with outgoing masks (8,1,3,6). Append a fixed C3, take their ordered sum B,
and graft B onto the source-rooted edge a→c. Mark vertex 0 of the appended
triangle. Put h=3^k 5^ell and q=12^k 32^ell. The complete rooted state is

    (D,E,U,d)=(3h+12q, 3q, 6h+14q−1_(k+ell=0), h+4q).

This uses the precise conventions and proved transfer matrix in
`05-knowledge/results/session225_graft_recovery_20261108/PROOF.md`.
For B=P⇒C3, the rooted-start bin cannot contain any prefix vertex. Thus
alpha_B=3H(P)=3h; for rooted ending, both bins can carry prefix, giving
delta_B=3Q_P(2)=3q. The analogous three-bin statistic is T_B=7q because
the two ordinary bins carry the prefix. The other entries are H_B=3h,
Q_B(2)=12q, s_B=1_(P empty), e_B=h, H(B−r)=h, Q_(B−r)(2)=4q.
The source edge has state (1,0,2,1); multiplying the inherited matrix gives
(H_B+Q_B,delta_B,2alpha_B+2T_B−s_B,H(B−r)+Q_(B−r)(2)), as claimed.
No polynomial inference or arbitrary integer-state realization is used.

The nonisomorphism argument is graph-theoretic, not merely rooted. Every
vertex inside B has an internal incoming and outgoing arc because each SCC
is strong with order at least three. After adding the two anchors, its
indegree and outdegree are both at least two. Therefore a is the unique
outdegree-one vertex and c the unique indegree-one vertex. Any isomorphism
must fix both, hence induces an isomorphism on B. The ordered SCC list of B
is unique, and its sizes three and four recover W, including repetitions.
There are exactly binom(k+ell,k) resulting nonisomorphic outputs.

Equality of four-state vectors implies equality of every common future
rooted-graft-context scalar by the proved closure in Session225 and context
equivalence in `05-knowledge/results/session226_structural_continuation_20261108/graft/PROOF.md`, Section 4. This statement keeps the specified marks and
common gadget sequence. It does not cover arbitrary rerooting or tests that
inspect graph structure. The outputs are not physically identical: their
SCC-word recovery explicitly distinguishes them. A heading should call them
rooted-graft-indistinguishable rather than indistinguishable without scope.

## 7. Sharp dimension of a linear gadget summary

The additional corollary is sound for every integer r≥1. Restrict to set-root
scalar contexts with min(|R|,|C|)≤r. Their counts are rational linear
functionals of (Q_B(1),...,Q_B(r+1)), so this uniform gadget summary has
dimension r+1.

For the matching lower bound, the source-rooted edge has row (1,1), and
T3 rooted at its source has row (−1,3) in the first two Q coordinates.
Their determinant is four. For every m=2,...,r, the balanced shell T_(2m)
with its first m vertices in R adds a row supported on the first m+1
coordinates whose top coefficient is (m!)^2. These r+1 context rows form
an invertible matrix. All those shells are strong-compatible.

Evaluate them on the physical gadgets T_n for n=1,...,r+1. The Q table is
(j^n)_(j,n=1,...,r+1), a Vandermonde matrix times nonzero column/row factors,
so it is invertible. The literal context/gadget response matrix therefore
has rank r+1. Any uniform assignment of a gadget vector in Q^d, with every
context count a rational linear functional of that vector, factors this
response matrix through Q^d. Hence d≥r+1. No arbitrary integer encoding or
nonlinear evaluation is excluded. The r=0 unrestricted context class needs
only Q_B(1), but has no strong-compatible context, so the strong witness
claim is stated only for r≥1.

## 8. Independent executable scope and exact controls

`independent_multiroot_audit.py` is standard-library-only and imports no
producer code. It uses four independently specified counting methods:

1. literal vertex-subset/last-vertex Hamilton-path dynamic programming;
2. permutation/backward-adjacency cuts for ordered profiles and coefficients;
3. a repeatable-B quotient-word dynamic program for independent coefficients;
4. separate all-subset HP and fixed-root-start tables for all four rooted
   coordinates, without importing or applying the producer's transfer matrix.

All checks use explicit exceptions, not optimization-disabled assert
statements. Universes and filters are explicit:

- every labelled A of orders 1–4, every subset R, and every labelled B of
  orders 1–3: 12,078 literal grafts, with no strongness filtering;
- every labelled A of orders 1–5 and every subset R: 33,866 coefficient,
  degree and singleton-strongness cases;
- all 22,890 compatible cases in the previous universe check c_2≥1;
- all 1,099 labelled B of orders 1–5 check their Hamilton endpoint, and all
  1,098 with order at least two check the two-path lower bound;
- 25 explicit transitive colour-word constructions check sharp degree for
  1≤|R|,|C|≤5;
- 27 three-factor menu words over {C3,T1,T2} check the scalar formula and
  order blindness on literal graphs;
- the three C3-power examples above use literal graphs through 13 vertices;
- ten literal equality graphs A=T2, R its source, B=T_n, 1≤n≤10 check
  the sharp lower bound 2^n+1;
- full rooted states for every collision word with k+ell≤2 and for (k,ell)=(2,1)
  are checked by separate all-subset HP/root-start tables, including unique
  anchor degrees and reconstruction of the ordered SCC-size word;
- 25 literal context/gadget responses check the linear-summary witness
  matrix, and its four successive square minors have ranks 2,3,4,5 for
  r=1,2,3,4; three balanced-shell top coefficients are checked independently.

Named hostile consequences: forgetting factorials gives 59 instead of 99;
allowing empty bins gives 213 instead of 99; bounding degree by |R| gives a
false bound already at |R|=2,c_3=2; nonempty proper R alone need not produce
a strong graph; B=C3 and B=C3⇒K1 have equal H(B)=3 but graft counts 99 and 279.

Reproduction from this folder:

    python3 -I -B independent_multiroot_audit.py > normal.json
    python3 -I -B -O independent_multiroot_audit.py > optimized.json
    cmp normal.json optimized.json

No repository files or producer files were edited by this audit.

## 9. Producer reconciliation

The completed producer manuscript has been read. The exact expansion,
strongness criterion, sharp growth bound, finite-menu/outer-join restriction,
explicit three-mode family and full-state collision proof pass independent
audit. The additional linear-summary corollary also passes the proof audit.

Before finalization, the producer's own exact checker detected an overstrong
unequal-shore coefficient equality. The repaired statement is the valid lower
bound (m!)^2 S(M,m), with equality asserted only for equal shores. This repair
is present in the manuscript audited here. The producer also incorporated
the sharp growth bound and the independent three-mode example, and narrowed
the collision heading/claim to common rooted-graft contexts with the fixed
marking convention. No unresolved mathematical blocker remains.

Final isolated normal and optimized replays both passed, including all
added linear-rank controls, and their JSON outputs are byte-identical.
`RECEIPT.json` records the final SHA-256 hashes and exact commands.
