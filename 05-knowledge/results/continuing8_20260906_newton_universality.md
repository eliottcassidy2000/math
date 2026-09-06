# Every positive Newton-circuit ratio vector has a real-rooted realization

**Status: PROVED + FINITE-EXACT; [independent audit accepted](continuing8_20260906_newton_universality_audit.md).**
The circuit map below is surjective even after fixing degree and the sum of
the positive root parameters. Every ternary sign word, including exact zero
circuits, occurs. This is a structural limit on what real-rootedness alone can
prove about the repository's Newton circuits. The fixed second moment and the
original factorial or wall-stripped coefficient row are separate constraints.

## Inheritance and the new question

The closest proved mechanism is the separated double-cluster maximal
alternation in **THM-3004**,
`01-canon/theorems/THM-3004-circuit-sign-change-cluster-law-and-classifier-refutation.md`.
The new [narrow-cluster theorem](continuing8_20260906_newton_clusters.md) fills
the interior bands with an exact factorial-curvature argument and obtains
the count `2K-3` for every fixed multiplicity profile with all entries at least
two. Its profile and width hypotheses remain essential.

The canonical hostile is the degree-five root-parameter row `(1,1,3,3,8)`,
whose circuit signs `-,+,-` defeat the two-end classifier. The corrected near
miss is the general norm/alternation dichotomy retracted in
[THM-3010's repair](continuing8_20260906_ballot_repair.md): positive reciprocal
symmetry and maximal alternation can coexist. The least-used operation here
is to invert the Newton-ratio map on **coefficients**, then certify the roots,
instead of searching root configurations for a desired word.

The live concepts are: cluster gap ratios; exact coefficient coordinates;
one-term sign domination; reversal; fixed first and second moments; and the
original coefficient consumer. Targeted searches for arbitrary circuit words,
prescribed Newton signs, circuit surjectivity, and strong coefficient
separation found no existing version of this consumer in the recovered routes.
No external priority claim is made for the elementary root-domination method.

For `N(n)=prod_i(n+r_i)`, with `d>=2` and all `r_i>0`, write

    N(n)=sum_(k=0)^d e_k n^(d-k),   h_k=e_k/binom(d,k),
    R_k=h_k^2/(h_(k-1)h_(k+1)),    1<=k<=d-1,
    C_k=R_k/R_(k-1),              2<=k<=d-1.

The circuit sign is `sgn(C_k-1)`. All roots referred to below are simple.

## 1. An exact surjectivity theorem

**Theorem.** For any `d>=2` and any prescribed positive real vector
`(c_2,...,c_(d-1))`, there is a monic polynomial `N` with `d` distinct negative
real roots, positive coefficients, and

    sum_i r_i=d,    C_k=c_k for every 2<=k<=d-1.

If the prescribed vector is rational, `N` can have rational coefficients,
with explicit rational intervals isolating its roots. Scaling all root
parameters realizes any other prescribed positive sum without changing the
circuit vector. The assertion is vacuous in its vector coordinate when `d=2`.

**Construction.** Put

    p_1=1,    p_k=product_(j=2)^k c_j,
    kappa_k=binom(d,k-1)binom(d,k+1)/binom(d,k)^2.

Choose any positive `lambda` satisfying

    lambda >= 9 max_(1<=k<=d-1) kappa_k/p_k.           (1)

Set `R_k=lambda p_k`, `h_0=h_1=1`, and recursively define

    h_(k+1)=h_k^2/(R_k h_(k-1)),
    a_k=binom(d,k)h_k,
    P(z)=sum_(k=0)^d a_k z^k,
    N(n)=n^d P(1/n).                                 (2)

All these quantities are positive. The ratio identities in (2) give the
prescribed `C` exactly, rather than asymptotically. Also `a_0=1` and `a_1=d`,
so `N` is monic and its positive root-parameter sum is `d` once its roots
are certified. For rational input, (1) permits rational `lambda`, and all
coefficients and samples below are rational.

## 2. Self-contained root certificate

We use a conservative constant nine; its optimality is neither required nor
asserted. From (1),

    a_k^2/(a_(k-1)a_(k+1))=R_k/kappa_k >=9.           (3)

Put `q_k=a_(k-1)/a_k` for `1<=k<=d`. Then
`q_(k+1)>=9q_k`. At the positive address `x_k=3q_k`, the absolute value of the
`k`th term of `P(-x_k)` strictly exceeds the sum of the other terms.

Indeed its immediately preceding term has ratio `q_k/x_k=1/3`; every
subsequent step toward the constant coefficient has ratio at most `1/3`.
The immediately following term has ratio `x_k/q_(k+1)<=1/3`, and every
later step toward the leading coefficient has ratio at most `1/3`.
Each finite tail is strictly less than

    sum_(j>=1) (1/3)^j = 1/2

of the central term. Their sum is therefore strictly less than the central
term, even if equality holds in (3). At `k=d` there is only the first tail.
Consequently

    sgn P(-x_k)=(-1)^k,   k=1,...,d,    P(0)>0.       (4)

The addresses `0<x_1<...<x_d` are strictly ordered. The intermediate value
theorem gives a root in every interval `(-x_k,-x_(k-1))`, with `x_0=0`.
Degree `d` makes these all the roots, each simple. Reversal then proves that
`N` has distinct negative real roots as claimed. The first interval for `P`
inverts to a half-line for `N`; the root-parameter sum `d` bounds every one
strictly below `d`, so it too can be replaced by a finite rational interval.
This completes the proof.

An explicit nonrecursive formula, useful as an independent verification path,
is

    h_k=lambda^(-k(k-1)/2)
        product_(j=2)^(k-1) c_j^(-(k-j)(k-j+1)/2).    (5)

## 3. Every ternary word and its equality boundary

For any word `sigma in {-1,0,1}^(d-2)`, take `c_k=2^(sigma_k)` and
`lambda=9*2^(d-2)`. Since `p_k>=2^(-(d-2))` and `kappa_k<1`, (1) holds.
Thus the construction realizes **every** circuit sign word, with zero
entries exact and with distinct roots throughout. In particular:

- Maximal alternation is possible at every degree `d>=3`.
- For `d>=3`, any number of sign changes from zero through `d-3` is possible.
- An identically zero circuit does not force coalesced root parameters.

Zero entries do not arise by rounding a small nonzero number; they are the
identities `C_k=1`. Root simplicity follows from (4) independently of those
ties. Clearing a common coefficient denominator preserves roots and
normalized ratios, but generally loses the displayed monic normalization;
the theorem claims rational coefficients, not integer root parameters.

There is a second useful construction for **strict** words: choose positive
integer gaps `g_k` with `g_k-g_(k-1)=sigma_k`, add a common positive offset,
and set `r_d=1`, `r_i=T^(sum_(j=i)^(d-1)g_j)`. For sufficiently large integer
`T`, the dominant coefficients give `C_k=T^(sigma_k)kappa_k/kappa_(k-1)`
up to an error tending to one. This gives integer root parameters. The exact
coefficient construction above is stronger concerning ties and target values;
the gap construction instead retains integer roots. They preserve different
coordinates and should not be merged into a single rational-root assertion.

## 4. The second moment is a load-bearing sidecar

Fixing the root-parameter sum is compatible with the entire positive circuit
image; fixing the sum of squares is an additional restriction. With
`S=sum r_i` and `Q=sum r_i^2`,

    e_2=(S^2-Q)/2,
    R_1 = S^2(d-1)/(d(S^2-Q)).                       (6)

Newton's inequality requires every `R_k>=1`. Since
`R_k=R_1 product_(j=2)^k C_j`, a fibre with fixed first two moments must obey

    product_(j=2)^k C_j >=1/R_1.                     (7)

Thus there is no surjectivity assertion on such a fibre. For the actual
positive root parameters `(1/2,1,3/2)`, `d=3`, `S=3`, `Q=7/2`, formula (6)
gives `R_1=12/11`. The target `C_2=1/2` would give `R_2=6/11<1`, impossible
for any positive-root-parameter polynomial with these same moments.

The repository's anchored degree-five beta polynomial has `S=13`, `Q=59`,
and hence `R_1=338/275`. The same target would give `R_2=169/275<1`.
This is an exact obstruction to transporting the unrestricted circuit theorem
into the current fixed-anchor Laurent sign model. It is a necessary-condition
failure, not a claim to classify that moment fibre or to decide its full sign.

## 5. Connection contract and remaining problem

The source is a positive vector of desired circuit ratios. The map integrates
it twice in coefficient-ratio coordinates using (2), chooses the free first
ratio using (1), and certifies real roots by (4). It preserves degree, the
entire circuit vector, positivity, simplicity, and a fixed root sum. It loses
the second moment, prescribed factorial coefficient structure, and any fixed
root-scale pattern. The essential sidecar for the anchored application is
`R_1`, equivalent to its second moment after fixing the first.

The cheap hostile (7) decides the attempted fixed-moment transfer immediately.
The constructive positive signal is stronger than another sign census: it
shows exactly which free coordinate allows all circuits. Consequently a
universal circuit-based no-return claim needs further structure. A precise
next problem is to characterize the circuit image for **fixed `R_1`** together
with any genuine interlacing or factorial constraints. Conditions (7) alone
are not asserted sufficient.

## 6. Reproduction

The standalone source `04-computation/continuing8_20260906_newton_universality.py`
uses exact fractions, constructs all `1,093` ternary words for degrees two
through eight, checks literal alternating rational root brackets, and handles
additional prescribed rational vectors and both moment hostiles. It does not
infer the universal statement from this bank. The source and independent
referee have distinct coefficient constructions.

Both normal and optimized modes pass **7,680 always-active gates**, with
byte-identical raw LF output and regenerated JSON. The complete bank is
retained by per-degree semantic hashes and explicit representative
coefficients and brackets in the certificate. No numerical root solver is
used by the producer. General LRC(14), the fixed-anchor sign problem, and
the original wall-core no-return problem remain open.
