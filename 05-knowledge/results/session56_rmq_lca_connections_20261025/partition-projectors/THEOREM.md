# Signed spectral projectors versus physical meet distributions

25 October 2026. Status: **PROVED** elementary finite-algebra statements; **FINITE-EXACT** checks at the bounds below. No novelty, group-representation, quantum, prime-identification, or LRC consequence is asserted.

## Inheritance and operation ledger

Closest proved mechanism: Session55 weighted CRT, whose labeled equality histograms multiply by meet convolution and whose upper zeta coordinates multiply pointwise. Canonical hostile: Session56 overlap signatures, components [2,5,3], where duplicating the F5 alternatives fabricates a collision-free state. Corrected near miss: a chosen partition satisfies pi meet pi = pi, but independent draws from alternative partitions need not agree. Least-used sidecar: the full signed coefficient space and its positive cone, together with the distinction between independent and coupled draws.

Concept board: named partitions; zeta coordinates; signed projectors; nonnegative physical histograms; independent product versus repeated variable; operation-preserving symmetry.

Sources are pinned local snapshots in SOURCE_PINS.json. This continuation does not modify the frozen Session55 packet or the Session56 overlap packet. Statements below have full elementary proofs rather than relying on an uncited standard-theorem invocation.

## 1. Algebra and its exact diagonalization

Let S be a nonempty finite meet-semilattice, ordered so x <= y iff x meet y = x. Full partition lattices use refinement order, with discrete bottom and indiscrete top. Work over K = Q or R. For f,g in K^S define

    (f*g)(t) = sum_{x meet y = t} f(x)g(y),
    (Zf)(r) = sum_{x >= r} f(x).

The product is bilinear, associative, and commutative because meet is. The key equivalence is

    r <= x meet y iff r <= x and r <= y.

Consequently Z(f*g)(r) = Zf(r) Zg(r). This is valid for signed coefficients; the finiteness makes rearrangement harmless.

Fix any linear extension of S. The matrix of Z is triangular with diagonal entries 1, hence invertible over Q and R. Define mu(x,x)=1 and, for x<y,

    mu(x,y) = -sum_{x <= z < y} mu(x,z),

and mu(x,y)=0 when x is not <= y. This is a triangular inverse to the zeta matrix: the displayed recurrence makes the corresponding matrix product identity; inverses of square invertible matrices are two-sided, giving also

    sum_{r <= x <= p} mu(x,p) = 1 if r=p, else 0.

For p in S set

    e_p = sum_{x <= p} mu(x,p) delta_x.

Then Z(e_p)(r) is the coordinate indicator 1_{r=p}. Therefore

    e_p*e_q = 0 (p!=q),  e_p*e_p=e_p,
    f = sum_p Zf(p)e_p.

The algebra identity is E=sum_p e_p; if S has a top, E=delta_top. A finite meet-closed subfamily of a partition lattice need not contain a top, so delta_top must not be asserted without that hypothesis.

A vector f is idempotent iff every Zf(p) is 0 or 1, since a scalar in Q or R satisfies a^2=a iff a in {0,1}. Thus all 2^|S| idempotents are precisely sums of distinct e_p. Each e_p is primitive: an orthogonal decomposition into two nonzero idempotents would split its singleton zeta support into two nonempty disjoint sets, impossible. These are algebraic projectors; positivity in zeta coordinates does not imply positivity of the physical delta coefficients.

## 2. Complete positivity obstruction, stronger than integer counts

For any coefficient vector define total mass M(f)=sum_x f(x). Then M(f*g)=M(f)M(g). Thus an idempotent has mass 0 or 1. For a nonnegative integer vector, mass zero forces zero and mass one forces a Dirac vector; conversely zero and every delta_x are idempotent.

In fact, **over the nonnegative reals the same classification holds**. Suppose p is a probability distribution and X,Y are independent with this common law. Its convolution is the law of X meet Y. Let r(x)=|{s in S:s<=x}|. This is strictly order-increasing: if x<y, the lower set of x is a proper subset of the lower set of y. Always r(X meet Y)<=r(X). If p*p=p, both variables have the same marginal law, so the nonnegative difference has expectation zero. Hence X meet Y=X almost surely. By symmetry X meet Y=Y almost surely, so X=Y almost surely. Independence now gives

    1 = P(X=Y) = sum_x p(x)^2.

But sum_x p(x)=1; if two probabilities were positive, the positive cross-term in its square would make sum_x p(x)^2<1. Thus p is Dirac. Conversely Dirac is stable. Combining mass and this theorem proves that the only nonnegative real convolution idempotents are zero and Dirac.

Independence is indispensable. If one couples Y=X instead, X meet Y=X for *every* distribution. This is the correct behavior when the same physical choice is reused; replacing that shared variable by an independent draw changes the operation.

An alternative spectral explanation is that idempotence forces every upper-set probability P(X>=r) to be 0 or 1. This does not rescue nontrivial random alternatives: the rank proof above shows that the only compatible positive law is a point mass.

## 3. Explicit Pi_3 and smallest random-mixture hostile

Order the five named states as bottom=1|2|3, a=12|3, b=13|2, c=1|23, top=123. The primitive projectors are

    e_bottom = delta_bottom,
    e_a = delta_a - delta_bottom,
    e_b = delta_b - delta_bottom,
    e_c = delta_c - delta_bottom,
    e_top = delta_top - delta_a - delta_b - delta_c + 2 delta_bottom.

For instance e_top has mass zero and negative coefficients, yet is a nonzero idempotent. The sum of all five projectors is delta_top. Only e_bottom itself is nonnegative; other nonnegative idempotents arise as sums, notably delta_p = sum_{r <= p} e_r. The 32 signed idempotents contain exactly six nonnegative ones: zero and the five Dirac vectors.

Let p=(delta_a+delta_b)/2. Since a meet b=bottom,

    p*p = delta_a/4 + delta_b/4 + delta_bottom/2 != p.

Thus a normalized average can fabricate a state absent from its own support. A coupled repetition of the same draw returns p, not p*p. This is a complete exact rational hostile to using probability normalization to justify overlapping blocks.

## 4. CRT and range-query consequences

Session55 counts arrays, not square-root choices; coordinatewise coprime CRT induces precisely the meet convolution above. Z diagonalizes multiplication but retains its ordinary multiplicative semantics. It does not turn multiplication into an idempotent operation. Repeating a block duplicates its zeta factors, changing Z(Q) to Z(Q)^2 unless Q is an actual algebraic idempotent.

A nonempty physical integer histogram is idempotent only when it contains exactly one array in exactly one state. Normalization helps only if all arrays have the same equality state: its normalized law is then Dirac, but the original integer mass still squares, so counts are not preserved. Even this exceptional support equality does not authorize reusing a physical CRT component or consuming a named port twice. Provenance is a separate constraint.

The pinned Session56 hostile [2,5,3] has true mass 112 and false overlapping mass 784, with an erroneously created zero-collision state. This packet cites that inherited finite certificate, rather than independently re-running its ring enumeration. No claim turns a repeated F5 factor into Z/25Z or F25. Sparse support and weighted multiplicity remain distinct. Independent physical component identities and fixed cell names must survive any implementation; use disjoint aggregation or a genuinely compatible-assignment join for overlap, not signed cancellation interpreted as physical alternatives.

## 5. Exact analogy boundary: S2, S3, and fixed-point marks

Label permutations act on partition states and preserve meet. For a permutation g, T_g(delta_p)=delta_{g p} extends linearly to an algebra automorphism, and sends e_p to e_{g p}. Thus the S3 action on Pi_3 fixes bottom/top and permutes its three middle projectors. This statement preserves the exact operation and is legitimate.

Forgetting labels down to block sizes is not compatible with binary meet: a and b have the same block sizes, but a meet a=a while a meet b=bottom. A simultaneous global relabeling is safe; independently choosing representatives of two orbits is not. The symmetry orbit sum e_a+e_b+e_c is a signed idempotent, not a probability law.

Do not identify this algebra with a group algebra. Under group convolution on S2=C2, the uniform distribution (delta_identity+delta_transposition)/2 *is* a non-Dirac idempotent: among four ordered products each group element occurs twice. This explicitly contradicts the meet conclusion because group multiplication is a different operation. No identification with S3 irreducibles, representation dimensions, or arithmetic primes follows.

There is a precise, limited analogy to fixed-point marks of finite G-sets: for H<=G define m_H(X)=|X^H|. Disjoint union makes these additive, while Cartesian product with diagonal action makes them multiplicative, since (X times Y)^H=X^H times Y^H. Likewise Z uses multiplicative evaluations for a different algebra, meet convolution. This paragraph proves only those evaluation identities, not an isomorphism of the Burnside ring with the partition algebra and not an equality of their primitive projectors. Any stronger transfer would require its own source, target, map, and product-preservation proof.

## 6. Exact tests and limits

Run python3 verify.py and python3 -O verify.py. The checker constructs all partitions on n=1,...,5 (1,2,5,15,52 states), computes meet by equality-mask intersection, independently builds the Mobius inverse by interval recurrence, and checks all primitive-projector products and all delta-basis diagonalization pairs. It enumerates all 32 Pi_3 idempotents and all Pi_3 probability histograms whose common denominator is 1,...,8, with exact Fraction arithmetic. Repeated representations across denominators are deliberately counted separately. It verifies independent-mixture and correlated-repetition controls and the C2 operation countercontrol.

python3 -O verify.py --self-test-failure must fail with RuntimeError. All substantive checks use explicit exceptions, never removable assertions. This finite computation corroborates the self-contained all-size proofs; it does not establish any larger arithmetic existence theorem or algorithmic speedup.
