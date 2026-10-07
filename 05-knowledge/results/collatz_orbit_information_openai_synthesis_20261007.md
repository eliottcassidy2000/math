# Collatz orbit information, decimal powers, and accepted-paper transfers

**Current status, October 7, 2026:** elementary results below are PROVED with
FINITE-EXACT controls; the requested OpenAI paper statements are CITED and
accepted as premises for this session, without re-auditing their proofs.
Universal Collatz coverage, archimedean box recurrence/merging, the decimal
cutoff at exponent 86, and the full selected prize targets remain OPEN here.

## What changed

The main gain is an exact account of which information two Collatz orbits
share, when a finite observation is forgotten, and when an adaptive query
draws fresh randomness. The latest Terras-clock integration additionally
closes bounded returns of the exponent imbalance in the Haar model, with
the alternative that the pair has already merged. Controlling translation
debt at those returns remains the missing archimedean condition.

The assigned anchor is source-specific Collatz coverage. The niche is the
decimal digit problem. The wildcard is a group of transfers from the supplied
papers: first-separating coordinates, marked tensor representations, adaptive
witness budgets, and arithmetic positivity. The common board is **source,
precision, carry, quotient, conditioning, and height**. No chronological log
is used as the current theorem source.

## 1. Both histories retain information; distant individual exponents forget

Let Y be odd Haar, X=alpha*2^v*Y+beta an odd-output affine coupling, and
let A_s,B_t be their cumulative odd-step valuations. A compatible pair of
prefixes specifies one source ball with

\[
M=\max\{A_s,(B_t-v)^+\}
\]

extra binary bits. Its mass is exactly 2^(-M), and the information shared
by that pair of observations is min(A_s,(B_t-v)^+) bits. At equal prefix
length t, both marginal entropies, joint entropy, and mutual information
are asymptotic to 2t bits. This is an exact lossless-information connection;
it is not a model of two independent random histories.

Nevertheless a source observation using M bits is independent of the
entire X exponent tail after v+M steps. The map U removes at least one
known source bit per step; at a critical coset, its odd part becomes full
odd Haar at once. Truncating an early exponent to a finite bit budget gives

\[
|\operatorname{Cov}(a_s,b_t)|^2
\le\min\left(4,38\Pr\{A_s>t-1-v\}\right),\qquad t\ge v+1.
\]

The probability is an explicit negative-binomial tail. In particular,
for each fixed s the correlation decays exponentially at long gaps.
[Complete proof, deterministic source comparison, and exact controls](collatz_paired_exponent_information_20261007.md).

This addresses all three requested comparisons with their proper meanings:

| Comparison | Retained law | Additional information needed |
|---|---|---|
| Two Collatz source orbits | Exact dyadic joint cylinders and long-gap bounds | Their affine relation and actual source/word guards |
| Positive plus and minus maps | Reflection through signed integers; common-prefix cutoff at v2(x+y) | The sheet and ROOT convention |
| Operator, free-factor and tensor orbits in the papers | Marked multiplication, trace or tensor restrictions, as applicable | Named generators/basis and the norm or growth notion |

The same-sheet common-prefix cutoff is v2(y-x). For the same positive
source on opposite sheets the first valuations already differ. Operator
norm, free-group word growth, arithmetic complexity, and real integer
height are separate observables, even when numerical exponents coincide.

## 2. Adaptive information revelation is the next useful controller state

The finer object is the full joint prefix, not just the exponent imbalance.
Conditional on it, the current Y and X states are odd Haar cosets with
qY=M-A_s and qX=M+v-B_t known extra bits. At least one of qY,qX is zero.

For a selected coset r mod 2^(q+1), the next exponent is either already
determined by those q bits, or it is q+G with G a fresh geometric draw.
In the latter case the source precision increases by exactly G. Thus a
predictable adaptive query order exposes a stopped sequence of independent
fresh draws, interspersed with deterministic uses of previously read bits.
Conditioning on a stopping length does not preserve an iid claim.

The residue is essential. For X=6Y+1, the joint prefix X=(1), Y=(1,2)
and the prefix X=(1), Y=(2,1) have identical consumed-bit totals, but the
next X exponent is respectively 3+Geom and the fixed value 1. This is a
minimal test against a controller storing only totals.

[The adaptive precision proof and checker](collatz_lazy_precision_innovations_20261007.md)
also localize covariance to the event that the later orbit has not consumed
the earlier observation's information. This is a stronger target than a
fixed deterministic delay: the remaining correlation occupies a random
overlap region of consumed binary positions.

For fixed v and any epsilon>0, this proves exponential decorrelation as
s grows, uniformly for t>=(1+epsilon)s. The band near matching times is
the unresolved region; numerical smallness there is not used as a theorem.

### Integration with the new concurrent overlap kernel

The subsequently fetched THM-4564/4565 checkpoints provide the complementary
local law inside an overlap. Independent checking confirms the kernel
`Cov(fresh,re-read | depth D)=2-6*2^(-D)`. Combining it with our exact
negative-binomial bound improves the estimate to

\[
|\operatorname{Cov}(a_s,b_t)|\le2\Pr(\text{window overlap})
\le2\Pr(A_s>B_{t-1}-v).
\]

The [integration note](collatz_overlap_kernel_integration_20261007.md)
records the stronger bound and two needed repairs. The full joint past
already determines the comparison depth, so HYP-9218's proposed geometric
law conditioned on that full past cannot approach its target. A replacement
must keep some depth information unobserved by specifying a coarser field.
The lockstep clock bound also has an equality/cancellation exception,
realized by a small lawful integer example. Exact local identities survive.

A useful infinite family shows why richer data matter: at initial alignments
with v=4k for odd k, depth probabilities 1/3 at1 and2/3 at2 give zero
covariance, agreement1/3 and average conditional
information4/3bits, just like geometric depth. The pair is still dependent.
All three statistics use only E[2^(-D)]. The full diagonal profile
`P(F=R=j)=2^(-j)P(D>j)` recovers the whole depth law and is the sharper
diagnostic for the next experiment.

**Concurrent return advance.** The newly fetched
[THM-4569, Terras clock](../../01-canon/theorems/THM-4569-the-terras-clock-recurrence-of-two-collatz-orbits-is-unconditional.md)
gives a simple random walk on predictable disagreement times. Its core
transition/recurrence proof was independently checked here. A separate
[inverse-clock argument](collatz_terras_inverse_clock_bridge_20261007.md)
transfers recurrent count ties to infinitely many odd-event pairs whose
times differ by at most one. Thus the appropriate synchronous exponent
imbalance returns to {-1,0,1} infinitely often, unless the aligned pair
has already merged. This uses a fresh conditional success chance at
separated stopping times, not a variance fit.

**Next proof obligation.** Control the affine translation at those returns,
so the pair visits a bounded state region with a usable merge chance.
Count recurrence alone does not bound this second coordinate. The
arithmetic transfer from Haar statements to the intended integer sources
remains separately necessary. Our depth diagnostics can target that
conditional translation control; they are no longer needed merely to
justify the count's recurrence.

## 3. Decimal zero-free powers: a proved sparse residual

Write Z(N) for the number of exponents 1<=n<=N for which the decimal
expansion of 2^n contains no zero. A two-state parity weight in the exact
suffix-lifting tree proves

\[
Z(N)=O(N^\theta),\qquad
\theta=\log_5\frac{5+\sqrt{17}}2=0.9429771\ldots<1.
\]

No priority or best-bound claim is made. The proof counts four/five
children per surviving suffix, retaining the parity of its division by
2^m. A positive weight improves the naive five-child bound to a spectral
factor strictly below 5. This is an actual example of constructing a useful
weight everywhere on a complete state space.

The finite checker covers all exponents 0..100000 and finds no zero-free
example above 86. The infinite cutoff remains unproved. For every fixed
suffix width there are surviving phases; imposing any fixed zero-free
leading word still leaves a positive-density set of exponents. Fixed
windows at both ends therefore cannot settle the full-digit statement.

The exact long-gap suffix relation is

\[
v_{10}(2^{n+d}-2^n)=
\begin{cases}\min(n,1+v_5(d)),&4\mid d,\\0,&4\nmid d.\end{cases}
\]

The next target is **position**, not total surviving mass: at decimal
length m>=2 a full counterexample must have its canonical suffix phase
e=n in the short band 10^(m-1)<=2^e<10^m, with m<=e<4m. Most suffix
survivors lie far outside this band.
[Full proof, accepted entropy-dimension comparison, and census](powers_two_decimal_windows_20261007.md).

The incoming concurrent
[THM-4580, 5-adic tree and larger verification](../../01-canon/theorems/THM-4580-zeroless-powers-of-two-the-5-adic-tree-and-verification-to-1-1e11.md)
records a substantially larger finite scan, 87<=n<1.1e11, and a stronger
computer-assisted counting exponent0.93783. Those are its recorded
computations, not additional runs performed by this package; its checkpoint
was marked audit pending. Our short parity-weight proof and its bounded
independent controls are retained. Both routes still leave the infinite
cutoff and the height-sensitive diagonal open.

This connects to Collatz in a precise methodological way. A small total
exceptional mass need not exclude a designated integer; retaining the
position/height coordinate exposes the extra estimate needed.

## 4. What the requested papers contribute

The supplied list of 65 integers is Euler's **idoneal list**, associated
with quadratic orders of discriminant -4n whose class group has exponent
at most two. In this domain the genus quotient loses no class information.
The value 105 has eight form classes; this is not the class-number-one list.
All n<=1848 were checked by exact reduced-form enumeration.

The strongest arithmetic transfer found is from the accepted ordinary
two-point Liouville-correlation theorem. For each fixed legal Collatz word,
its source and target are nonproportional affine functions of the cylinder
parameter. The theorem therefore gives asymptotic frequency one half for
equal versus unequal total-prime-factor parity along that fixed cylinder.
It does not assert independent Collatz valuations or a uniform estimate
over growing words.

Serre's positivity theorem makes the desired interface sharper: one needs
nonzero modules with a finite proper intersection, plus a comparison of
that intersection multiplicity with the **designated source atom**. The
rational graph of a word intersecting the diagonal merely gives multiplicity
one even at an illegal or negative anchor. The missing work is encoding
source legality and ROOT reachability in the modules without assuming it.

Additional selections were Jacobsthal bounds, multiple mixing,
inhomogeneous Duffin--Schaeffer, entropy-rate dimension, and unitarizability.
Their usable hypotheses and failed transfers are recorded in the
[arithmetic note](collatz_openai_arithmetic_bridges_20261007.md) and
[operator note](collatz_operator_tensor_bridges_20261007.md).

The free-factor paper's most actionable mechanism is choosing future
error budgets **after** a finite witness and its sensitivity are known.
This preserves an established positive margin under refinement. It does
not create that initial margin. The Nine-Fourths paper supplies a concrete
Fourier shared-label filter and the representation split 4=3+1. Marked
tensor squaring retains both slope and carry; unrestricted tensor changes
erase them. The inherited identity lambda^2 and J^2-2 remains exact, with
the repeated legality guards retained.

The linked GCH paper supplies a first-separating-coordinate device, not a
Collatz theorem. Here it suggests preserving the maximum coordinate depth
of a compressed witness. That depth is precisely what the forgetting
theorem charges, even if only a few separator bits are stored.

## 5. Prize status and a concrete family

All six requested accessible targets were inspected. The supplied Erdős9
URL ends in `erdos-91`; the catalog's accessible Erdős9 target was used.
No full prize proof or refutation was obtained, and nothing was submitted.
The platform requires the pinned Lean target and validator/review acceptance;
partial contributions are a separate route, not an immediate prize claim.

The clearest family certificate is for Erdős23: every graph on 5n vertices
admitting a homomorphism to C5 becomes bipartite after at most n^2 deletions.
Delete the smallest adjacent class product. For a complete C5 blowup this
number is exactly min_i a_i*a_(i+1); balanced parts attain n^2. Arbitrary
triangle-free graphs need not have such a homomorphism, so the full target
is still missing. We make no novelty claim for this elementary family.

Borsuk in dimension nine does not settle the dimension-four prize.
Squarefree-value density does not by itself choose a squarefree n-2^l
for every fixed odd n. These are the same representation and quantifier
boundaries encountered in the source-positivity work.
[Exact targets, source links and scoped family proof](collatz_operator_tensor_bridges_20261007.md).

## 6. Reproduction and next priorities

Every package links its exact script and deterministic output. Normal and
optimized execution are compared; mathematical proofs provide the infinite
quantifiers, while finite controls check boundaries, carry, types and
independent enumerations. No accepted OpenAI theorem is claimed proved by
these tests. Agent audits checked the newly derived interfaces.

The seven packages pass **217,070 exact checks** in both normal and optimized
Python, with matching saved outputs. The navigation check also passes;
an inherited line-budget excess was repaired by joining wrapped lines in
the hypothesis index without changing its mathematical text.

The next priority is translation tightness at recurrent count/imbalance
returns. A coarse conditioning that leaves the overlap depth unrevealed
may support that estimate; its full depth profile is the useful diagnostic.
The orthogonal digit target is a surviving-phase bound
in the actual-height band. The positive-measure arithmetic route remains
target absence implying a controlled arithmetic collapse, so that a
nonzero small determinant could force a contradiction. These are precise
remaining obligations, not renamed assumptions of universal coverage.
