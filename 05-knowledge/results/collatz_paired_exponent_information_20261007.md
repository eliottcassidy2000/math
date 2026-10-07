# Paired Collatz exponents: retained information and quantitative forgetting

**Status:** PROVED elementary probability and cylinder statements below;
FINITE-EXACT controls. Archimedean box recurrence, almost-sure merging, and universal
integer Collatz coverage remain OPEN. No novelty or priority claim is made.

**Reproduction:** `python 04-computation/experiments/collatz_paired_exponent_information_20261007.py`
and the same command with `python -O`; deterministic output is beside this note.

## Inheritance and scope

The anchor is the paired debt model in
[the October 6 tube/debt note](collatz_cycles_tubes_debt_walk_openai_20261006.md),
especially its corrected distinction between individual Haar exponent laws
and the joint process. The exact source-cylinder mechanism is
[THM-4512, coefficient-descent classes](../../01-canon/theorems/THM-4512-coefficient-descent-classes-one-member.md).
The hostile is the identical coupling X=Y: both marginal exponent sequences
are iid geometric, but their difference is identically zero. The corrected
near miss is inferring random-walk recurrence from marginal laws or a finite
covariance table. The least-used sidecar is **consumed binary precision**.

The live board is: legal word cylinders; dyadic source balls; retained mutual
information; finite-observation forgetting; debt returns; and digit-height
diagonals. The niche comparison is the plus/minus sheet. The wildcard is the
first-separating-coordinate device in the linked set-theory paper.

Throughout, U(n)=oddpart(3n+1) is the ordinary map on odd 2-adic integers,
with normalized odd Haar probability. Delete the countable Haar-null set
whose iterates encounter 3n+1=0. ROOT is not made absorbing in this
probabilistic model. Source-specific positive-integer certificates still
need their separate first-hit and legality guards.

## 1. Exact joint cylinder and information formulas

Let Y be odd Haar and put

\[
X=\alpha 2^vY+\beta,
\]

where alpha is an odd integer, v>=0, and beta is even if v=0 and odd if
v>=1. Thus X is uniform on one odd residue class modulo 2^(v+1).
Write a_i and b_i for the successive valuations of the Y and X orbits,
and A_s=sum_(i=1)^s a_i, B_t=sum_(i=1)^t b_i.

A prescribed valuation word w of total A has the affine expression

\[
U^{|w|}(n)=\frac{3^{|w|}n+C_w}{2^A},\qquad
n\equiv (2^A-C_w)3^{-|w|}\pmod{2^{A+1}}.
\]

This one odd residue class is exactly the legal word cylinder; its odd Haar
mass is 2^(-A). Pulling an X-word of total B back to Y gives either the empty
set or a single odd ball with E=max(0,B-v) extra binary bits. In the latter
case its mass is 2^(-E). This includes B<v, when the initial coset either
forces the entire X-prefix or makes it impossible.

Two dyadic balls are nested or disjoint. Therefore, on a compatible pair
of word cylinders,

\[
p_Y=2^{-A},\quad p_X=2^{-E},\quad
p_{Y,X}=2^{-\max(A,E)},\quad
\log_2\frac{p_{Y,X}}{p_Yp_X}=\min(A,E).                 \tag{1}
\]

For fixed word lengths s,t, this yields exactly

\[
I(W_s^Y;W_t^X)=\mathbb E\min(A_s,(B_t-v)^+),\qquad
H(W_s^Y,W_t^X)=\mathbb E\max(A_s,(B_t-v)^+).           \tag{2}
\]

The marginal entropies are 2s and E(B_t-v)^+, respectively. At equal lengths,

\[
\lim_{t\to\infty}\frac{H(W_t^Y)}t
=\lim_{t\to\infty}\frac{H(W_t^X)}t
=\lim_{t\to\infty}\frac{H(W_t^Y,W_t^X)}t
=\lim_{t\to\infty}\frac{I(W_t^Y;W_t^X)}t=2.           \tag{3}
\]

**Proof of the limits.** Under full odd Haar, the valuations are iid with
P(a=k)=2^(-k), k>=1, mean 2 and variance 2. Hence A_t/t tends to 2 in L1.
The law of X has density bounded by 2^v relative to odd Haar, so the same
L1 conclusion holds for B_t/t. Taking min and max preserves this limit.
No stationary law of the paired process is assumed.

This is lossless-history coupling at the entropy-rate level. It does not
imply a large covariance for any particular pair of exponents: two long
histories can retain essentially the same random information while
individual observations at separated times become weakly correlated.

## 2. Every finite source observation is eventually forgotten exactly

**Finite-coset lemma.** If Z is Haar on an odd residue class modulo
2^(M+1), then U^M(Z) is full odd Haar, for every M>=0.

For M>0, there are two cases. If a=v2(3r+1)<M+1, this valuation is fixed
on the coset, and U sends it to an odd coset with M-a extra bits. If
3r+1 vanishes modulo 2^(M+1), write 3Z+1=2^(M+1)T, where T is Haar on
Z_2. Its odd part is full odd Haar: conditional on each valuation, the
remaining odd unit has that law. Each noncritical step removes at least
one known bit, and full odd Haar is U-invariant. The bound M is sharp for
the all-valuation-one cylinder of length M.

Now condition Y on any one of its M-bit odd source atoms. X is then Haar
on a coset with v+M extra bits. The lemma proves:

\[
F(Y\bmod 2^{M+1})\ \text{is independent of the entire tail}\
(b_{v+M+1},b_{v+M+2},\ldots)                         \tag{4}
\]

for every bounded finite-precision observable F. Indeed, conditional on
each such source atom, U^(v+M)X has exactly the same full odd Haar law.

The retained sidecar is the deepest binary position read, not merely the
number of bits kept after compression. For example min(a_1,R) is an R-bit
observable, so it is independent of that X-tail whenever M=R. The complete
a_1 is unbounded and has no uniform finite observation depth.

## 3. A rigorous long-gap exponent-covariance bound

For s>=1, t>=v+1, put M=t-1-v. The truncated centered observable

\[
h_M=(a_s-2)\mathbf1_{\{A_s\le M\}}
\]

is measurable from Y modulo 2^(M+1): parse up to s valuation symbols
within M bits, returning zero if the bit budget is exhausted. By (4),
E[h_M(b_t-2)]=0. Since a_s and b_t both have the geometric marginal law,

\[
\operatorname{Cov}(a_s,b_t)
=\mathbb E[(a_s-2)(b_t-2)\mathbf1_{\{A_s>M\}}].
\]

That geometric law has fourth centered moment 38. Holder with exponents
4,4,2, followed by the binomial/negative-binomial identity, gives

\[
\boxed{\quad
|\operatorname{Cov}(a_s,b_t)|^2
\le\min\left\{4,\,
38\,2^{-M}\sum_{j=0}^{\min(s-1,M)}\binom Mj\right\}.
\quad}                                                       \tag{5}
\]

The cap 4 is Cauchy--Schwarz using variance 2. For fixed s, (5) decays
exponentially in t, with a polynomial factor of degree s-1. For v=2,
s=1, t=32, it gives covariance squared <=19/268435456. It is not a
uniform small bound near t=s; its dependence on both indices is essential.

### The actual Mersenne lag-one comparison

The recent debt model uses x_0=3*2^v*y_0+1, with v>=1 for oddness, and
compares its orbit with the orbit after one Y-step. Conditional on
a_0=v2(3y_0+1), let z=U(y_0). Then z is odd Haar, independent of a_0, and

\[
x_0=2^{v+a_0}z+1-2^v.
\]

Thus (1)--(5) apply conditionally with alpha=1, scale v+a_0, and
beta=1-2^v. The random initial valuation is retained; silently dropping
it would give an incorrect uniform forgetting time.

## 4. Deterministic comparison of two sources and the two sheets

Let x!=y be odd integers and h=v2(y-x). If their first k valuations under
U agree, the gap after those steps has valuation h-A_k. In particular,

\[
\text{the first k valuations agree}\quad\Longleftrightarrow\quad A_k<h.
                                                               \tag{6}
\]

To prove the cutoff, at one step compare 3x+1 and 3y+1, whose difference
has valuation h. If the first valuation is below h, the two valuations
agree and subtract from the gap precision. Otherwise the valuations
differ and their minimum is h. The first differing symbol is therefore
the first one whose cumulative A reaches or exceeds h. This proves (6)
inductively, including equality at the boundary.

For U_-(y)=oddpart(3y-1), reflection gives U_-(y)=-U_+(-y). Thus the same
cutoff compares the plus orbit of x with the minus orbit of y after
replacing v2(y-x) by v2(x+y), when x!=-y. For the same positive odd source
on the two sheets, h=1: their very first valuations differ.

These are exact source laws, not asymptotic mixing assertions for every
fixed integer. They cover the second requested orbit comparison without
identifying its different cycles or ROOT rules.

## 5. The remaining debt-return target

For a block starting after r>=v, both marginal tails are iid geometric.
Consequently the exact variance identity is

\[
\operatorname{Var}\!\left(\sum_{i=r+1}^{r+K}(a_i-b_i)\right)
=4K-2\sum_{i,j=r+1}^{r+K}\operatorname{Cov}(a_i,b_j).             \tag{7}
\]

Equation (5) controls a remote portion of this double sum. It does not
control the expanding region near comparable indices, nor correlations
of blocks after conditioning on current debt and a failed merge attempt.
The identical coupling X=Y satisfies the finite-coset and entropy theorems
and has zero cross-covariances off the diagonal, yet the left side of (7)
is zero. It is a decisive hostile to any proposed general inference from
remote decorrelation to a variance-4 walk.

For the nontrivial Mersenne affine coupling, the productive next target
is a conditional innovation or return bound retaining current source
precision, exponent imbalance L, and affine debt Delta. Even a variance
asymptotic or central limit theorem alone would not prove almost-sure
returns with a uniformly usable merge chance. Haar-almost-everywhere
merging would in turn still require the separate arithmetic transfer to
the requested integer family; it is not universal integer Collatz.

## 6. A typed connection to the linked set-theory paper

The user-supplied [paper by Zhixing You](https://arxiv.org/html/2610.07537v1)
studies GCH above a strongly compact cardinal. Its Lemma 3.1 uses the first
coordinate where two traces differ as a separator. We accept the paper's
stated results as requested and do not re-audit its proof.

The transferable device is finite: retain a separating coordinate for
each pair of unequal observations. It yields an injective code on that
finite family. Here the target predicate is distinguishing source
cylinders. The destroyed information is unretained source precision;
the needed sidecar is the largest coordinate index. One bit at position
100000 can distinguish two objects but does not provide a one-step
forgetting bound. Equation (4) pays for the position of that bit.
No large-cardinal hypothesis or cardinality theorem is asserted to prove
a Collatz probability estimate.

The decimal-window counterpart is in
[the companion note](powers_two_decimal_windows_20261007.md): many suffix
phases survive, but the full-power question concerns a narrow band where
exponent and actual decimal length agree. Both comparisons require
preserving **where** the information sits, not only how much remains.

The subsequent [overlap-kernel integration](collatz_overlap_kernel_integration_20261007.md)
improves the covariance bound using the exact saturation law, and corrects
the incoming full-past Haar-depth hypothesis. Equations (1)--(7) above
remain valid; the stronger estimate is linear in overlap probability.

The later [Terras inverse-clock bridge](collatz_terras_inverse_clock_bridge_20261007.md)
proves bounded returns of the synchronous exponent imbalance in the Haar
model, outside the merge alternative. The unresolved return obligation
is now the simultaneous control of the translation coordinate.

## 7. Finite controls and audit

The checker independently counts source residues for all word lengths
0..3 with cumulative valuation <=6 under six affine maps, including the
identity, reflection, and Mersenne-style maps. It checks each joint mass
and information ratio; every odd coset with 0..11 extra bits and the sharp
flush time; both deterministic cutoff laws for all positive odd x,y<100;
and the negative-binomial formula by independent fair-bit enumeration.
Twelve malformed-input controls reject invalid types or domains.

Normal and optimized execution perform **21,018 exact checks**. These are
controls on the algebra and finite boundaries, not finite simulations
offered as proofs of (3)--(5). The probability proofs above supply the
unbounded quantifiers. An independent agent audited the cylinder,
flushing, information, and Holder steps; no mathematical correction was
needed. Strict sign and length input guards were added during that audit.
