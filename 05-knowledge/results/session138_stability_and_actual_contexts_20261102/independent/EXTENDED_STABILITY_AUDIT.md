# Extended stability audit: first Session138 checkpoint, through order seven

**Verdict: ACCEPT, in the stated domains.** This document belongs to the
first checkpoint. Order-eight work is deliberately excluded and audited in
the separate `ORDER_EIGHT_CONTINUATION_AUDIT.md` and `n8_audit/` continuation.
The all-order optimal constant remains unknown.

## All-order multi-star rounding

The one-edge interpolation sign and minimum inequality are correct. Pair-
deleted remainders are nonempty precisely for n>=3; order2 is a genuine
counterexample to a positive variance correction. A binary uncertainty-free
core of size m>=1 supplies2^(m-1) odd-integer Pfaffian squares in each of the
two augmented-skew parity sums. Thus an uncertainty vertex cover of size c
implies `u-|v|>=2^-c`. If m=0, the inherited uniform bounds
`u+v>=2^(1-r)` and `u-v>=r2^(1-r)` give the claim for r>=1 without a fictitious
half-subset count. The empty pair-deletion case is not used.

For ordered centers x1,...,xk, preliminary rounding outside the centers
safely discards nonnegative corrections. An edge first reached at position t
has remainder cover size at most k-t (outside other endpoint), or k-t-1
(later-center other endpoint). Its unrounded probability and variance are
still their original values on every conditioning branch. The respective
assignment probabilities for a fixed edge are
`2(n-k)/[n(n-1)]` and `2(k-t)/[n(n-1)]`. The second case has zero probability
at t=k. Summing and averaging yields exactly

    Delta(A)>=E Delta(T)+gamma_(n,k)V, n>=3,1<=k<=n,
    gamma_(n,k)=4[n-k+2-(n+k+2)/2^k]/[n(n-1)].

Every ordering has the same independent product-law terminal T, so its
expectation is retained once. The k1 value2/n, consecutive-numerator
difference `-1+(n+k+1)/2^(k+1)`, and asymptotic optimum4/n all check.
The simplified k=ceil(log2 n) estimate has the correctly stated n>=4 domain.
The sharp H1 upper family is asymptotic8/n; no optimality at all orders follows.

`audit_multistar.py` supplies an independent exact Fraction replay using
Hamilton subset DP, Gaussian determinants and literal binary expectations.
It checks808 matrices (all ternary n1..4 and12 fixed-seed rational cases,
with six fractional edges, at each n5..8),808 binary-core bounds and3309
multi-star inequalities. It caches2916 independently evaluated binary gaps.
All coefficient sum/difference identities are checked for n3..200, and the
n2 failure is explicit. Normal/-O receipts match; semantic digest:
`d2d3c3c9706adc4b9730ae5be4f687888f2638a43e64ba737dd00ee06182f017`.
These finite checks corroborate, not replace, the all-order proof above.

## Local triangle variance

The degree-two Bernstein argument is sound. An active edge (index0 or1)
has variance coefficient1/4 and skew magnitude at most1/2 in the triangle
coefficient; an inactive edge has variance coefficient0 and skew magnitude1.
With two active edges, the minimum triangle coefficient is1/16; with three
it is3/16. With at most one, the variance-pair coefficient vanishes. These
minima follow from the separately affine identity4tau=1+xy+xz+yz.
The independent replay evaluates all216 orthant/index coefficients directly:
minimum0,54 zero coefficients. Thus

    tau >= W(q)=sum_(incident unordered edge pairs) q_e q_f.

Multiplying each triple inequality by q on a disjoint edge gives
`M23>=P21(q)` with no extra multiplicity. Wedges are not required to be
induced; extra graph edges do not invalidate a counted pattern. The summed
variance consequence `tau>=max(0,2V^2/n-V/4)` follows from Cauchy and q<=1/4.

## Independently reproduced six/seven graph certificates

`graph_gray_audit.cpp` uses a fresh counting method. It traverses every
labelled graph in Gray-code order and updates patterns according to the role
of the toggled edge: isolated edge of a wedge-plus-edge, either wedge arm,
or an edge of a matching. Matching increments are computed from the induced
complement after deleting its endpoints. It does not import producer code,
data, representatives or receipt values. A separate triple-edge brute count
checks25 snapshots at n6 and529 at n7. The wedge identity
`w=binom(m,2)-m2` is checked throughout. All checks are active exceptions.
Normal and undefined-behavior-sanitized executions match byte-for-byte.

* n6: all32768 labelled graphs; `11m-12m2+3m3+3p21>=0`;
  minimum0, exactly76 zero graphs.
* n7: all2097152 labelled graphs;
  `22m-24m2+6m3+8w+6p21>=0`;
  minimum0, exactly421 zero graphs.
* The graph coefficient in `P+alpha W+P21>=37V/48` is independently
  optimized to alpha=1/3. The auditor's witness3096 is P3 plus two edges,
  different from the producer's chosen mask but the same graph type.

Every pattern polynomial is affine in each individual q_e. Its endpoint
checks therefore prove the full continuous q-box inequality by tensor linear
interpolation. This is not an inference from a weighted-matrix sample.

## Sharp weak order six

The exact cycle expansion is

    delta=V+tau+C4+C5+C6-M22+M23-M24+3M33+M222.

Since a4-cycle's complementary pair has variance<=1/4, `M24<=C4/4`.
After subtracting tau+37V/48, the remaining non-cycle part is
`I=M222+M23+11V/48-M22`. Replace M23 by P21. At q-endpoints,
`192I_low=3m3+3p21+11m-12m2`, exactly the audited n6 graph residual. Hence

    delta>=tau+37V/48+3C4/4+C5+C6+3M33.

Three fair ordinal pairs give H1,u27/64,V3/4 and ratio37/48, proving
sharpness. Unreserved positive-V equality forces tau0 and H1, then the
scalar H1 minimization forces three fair pairs. The larger producer weak6
Bernstein run is corroborative, not required by this proof or independently
replayed here.

## Sharp strong order seven

The complete determinant and odd-cycle expansions in
`SHARP_SEVEN_GRAPH_CERTIFICATE.md` were checked term by term, including
M25,M34,M223 and the sign v=-det(A). Writing P=V-M22+M222 gives

    G_-=P+tau+C4+C5+C6+2C7+M23-M24+3M33,
    G_+=P+tau+C4+C5+C6+M23-M24+2M25+3M33+2M34-2M223.

Each4-cycle's complementary triple has total variance<=3/4, so
`M24<=3C4/4`. Each triangle's complementary four-set has three2-matchings,
so `M223<=3tau/16`. Both margins are therefore at least
`P+5tau/8+M23+C4/4+C5+C6+3M33`.

The audited n7 graph residual is384 times
`P+W/3+P21-37V/48`. Keeping the excess triangle reserve gives

    Delta>=37V/48+7tau/24+C4/4+C5+C6+3M33.

Three fair pairs and one singleton attain37/48. Positive-V equality in the
unreserved bound forces tau0, hence H1, and exactly those fair ordinal blocks.
The7/24 triangle reserve is not claimed optimal for actual tournaments.
The older failed graph reduction without P21 remains only a hostile control;
its negative residual is not a counterexample to this tournament theorem.

## Root independence-polynomial edge-order note

The proof in `../rounding_graph/PROOF.md` was also reviewed. For an uncertainty
supergraph J on r>=1 vertices,
`mu(J)=2^(1-r)min(i_even(J),i_odd(J))` is a valid determinant-margin lower bound:
even independent sets supply binary Pfaffians; odd ones supply observer
Pfaffians. Iterating a deterministic edge order gives exactly
`Delta>=E Delta(T)+sum_e q_e mu(J_e)`.

The sharp n4 smaller-variance-first complementary-pair argument is valid.
At n5, a three-vertex remainder with s uncertain edges has mu=1-s/4, giving
the stated support-sensitive sum of min(q_e,q_f) over disjoint edge pairs.
Neither support estimate is silently identified with the optimal bound at
larger orders.

For c disjoint uncertain pairs and t isolated vertices, the independence
polynomial is(1+2z)^c(1+z)^t. Its parity counts give mu=(3/4)^c for t>0 and
mu=(3/4)^c-4^-c for t0,c>=1. In sequential rounding with no original isolated
vertices, only the first edge incurs the imbalance term: previously rounded
pair endpoints thereafter become isolated vertices of the uncertainty graph.
Averaging orders recovers precisely the sharp H1 constants for matching
uncertainty even with arbitrary binary cross orientations.

The additional fixed-parameter parity algorithm is valid: enumerate each
independent S inside a vertex cover C; outside C is independent, so with
d=|I\\N(S)|, d>0 contributes2^(d-1) to each parity count, while d0 contributes
one only to parity|S|. This computes exact parity counts in O(2^c poly(n)).
Equality cautions: n3 has Delta=V for every weighted tournament, not only H1;
for n>=4 and one fractional pair, H1 equality permits any probability; with
at least two pairs, fairness is forced. No equality classification for the
stronger baseline-containing inequality is asserted.

## Reviewed bytes and replay scope

`REVIEWED_SOURCES.json` pins reviewed source snapshots with explicit scope.
A later mathematical edit requires renewed review; status-only changes do
not retroactively broaden this audit. `FIRST_CHECKPOINT_FILES.json` lists
only this first-checkpoint auditor's artifacts and excludes all n8 files.
No graph census beyond the stated universes, proof of all-order optimality,
unrestricted zero-rigidity closure or full-spectrum realization is claimed.
