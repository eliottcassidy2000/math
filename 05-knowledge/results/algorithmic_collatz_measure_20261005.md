# Algorithmic information and arithmetic mass in Collatz certificates

2026-10-05. **CITED:** the prefix complexity and randomness framework in
Jan Reimann, [Information vs Dimension, v1](https://arxiv.org/html/2408.05121v1),
sections2.2,4.3-4.4,5.2. **PROVED, scoped:** the explicit atomic cylinder
formula, its pointwise dimensions, and the certificate distinctions below.
**FINITE-EXACT:** the accompanying checks. **OPEN:** universal Collatz
coverage, atomic survivor-mass decay, and a uniformly successful paid selector.
No literature-priority or machine-verified proof claim.

The useful connection is an exact change of measure on the SAME Collatz
binary words. Haar mass describes generic addresses. A computable atomic
measure remembers every positive integer separately. Prefix-free coding
gives a budget for searching and storing certificates. None of these
operations, by itself, supplies a missing root certificate.

## 1. Inheritance and the live questions

Closest proved mechanisms: [finite seed kernels](finite_seed_kernel_20261005.md)
and [crossroads integer carrier](crossroads_poset_20260926_integer.md), section4A.
The latter already proves that decay of positive mass assigned to every
integer is equivalent to eliminating all actual survivors. It also identifies
the loss when taking a Haar-density limit. This reduction is inherited here.

Canonical hostile: the computable signed fixed point -1 is invisible to
Haar singleton mass. Corrected near miss: mass-one or large-dimension generic
coverage is not pointwise integer coverage. Least-used sidecar: keep the
source atom and a finite verification transcript together with the code.

Anchor: certificates covering actual integers. Niche: finite Vitali-labelled
quotients and computable representatives. Wildcard: an atomic measure with
full Cantor support and nonconstant pointwise dimensions.
Board: **source atom / cylinder mass / effective complexity / affine carry /
verified stopping / quotient fibre**.

The meta-patterns used are *Classify the target over the exceptional stratum
before closing by density* and *Type every analogy and every implication* in
[META-PATTERNS](../../00-navigation/META-PATTERNS.md). The tests below distinguish
all six board coordinates instead of merging them into one entropy number.

Incoming `358bd7ee5` supplies a paid J controller and finite-pool conformal
repair in [four receipt helpers](collatz_fusion_helpers_20261005.md). Its
stock and endpoint checks are suitable certificate-verifier inputs. Its
conditional finite-pool completeness is not a global search termination
theorem. No empirical Fourier rate is used in the arguments here.

## 2. The source passage and its missing sign

The standard randomness criterion is

    K(x[:L]) >= -log2 mu([x[:L]]) - c, for every L.       (1)

The quoted paragraph in v1 section4.4 omits the minus sign; the display in
Theorem4.11 also omits the additive constant mentioned in its sentence.
Both omissions occur in the PDF as well as HTML. We use (1), not those
literal misprints. Positive linear surprisal, together with randomness,
gives positive effective dimension. The randomness hypothesis is essential.

## 3. A shared binary coordinate

Use the shortcut map T(n)=n/2 for even n and (3n+1)/2 for odd n. For m>=1 put
n=2m-1 and define x_m[j]=T^(j+1)(n) mod2. Thus the common initial odd bit is
omitted. The map continues through the usual root cycle; it is never asked
to decide whether an orbit eventually reaches1.

If a_i=v2(3U^(i-1)(n)+1), then

    x_m = 0^(a_1-1)1 0^(a_2-1)1 0^(a_3-1)1 ... .     (2)

Consequently this is exactly the renewal coordinate of the
[valuation information package](valuation_information_spectrum_20261005.md).
Each length-L binary prefix selects one residue m=r mod2^L, 0<=r<2^L.

Proof of the residue assertion: after j shortcut steps the state is
(P n+C)/2^j with odd P. Prescribing its next parity uniquely specifies
n mod2^(j+1), since P is invertible there. The first parity is fixed to1;
the remaining L choices leave exactly one odd n mod2^(L+1), equivalently
one m mod2^L. Compatible residues extend this map to a homeomorphism from
Z2 to Cantor space. It is an isometry with the usual first-disagreement
metrics. Positive m give a dense countable subset, not all Z2.

## 4. An exact computable atomic measure

Let ell(m)=2 floor(log2 m)+1 be the length of the Elias gamma code, and set

    w(m)=2^-ell(m),       mu=sum_(m>=1) w(m) delta_(x_m). (3)

The block 2^k<=m<2^(k+1) has mass2^(-k-1); hence the weights sum to1,
and the tail m>=2^K has mass2^-K. For a prefix sigma of length L, let r
be its unique m-residue in [0,2^L). Then

    mu([sigma]) = 4^-L + (w(r) if r>0, else 0).         (4)

Proof: in every block k>=L there are 2^(k-L) matching integers. Their
total mass summed over those blocks is4^-L. Below2^L there is exactly
the representative r if r>0, and none if r=0. This also covers L=0.

Equation(4) is an exact rational cylinder algorithm, not just a semicomputable
approximation. The measure is computable, purely atomic, and has full
topological support: every cylinder has positive mass. These three properties
are compatible. Its mass is concentrated on a countable dense subset.

For a fixed positive m, once 2^L>m the residue is m, so

    mu([x_m[:L]]) = w(m)+4^-L.                         (5)

Thus its local dimension is0. It is also Martin-Lof random relative to mu:
a test set of measure less than w(m) cannot contain x_m. This is a direct
atom argument, requiring no Collatz convergence assumption.

Meanwhile x_m is computable, uniformly in m. A self-delimiting description
of m and L computes its first L bits, so K(x_m[:L])<=K(m)+K(L)+O(1)=O_m(log L).
Its effective dimension is0, whether or not the actual orbit reaches1.
There is no contradiction: randomness is relative to a measure, and these
random points have atoms rather than linearly growing surprisal.

There is also an exact recursive reuse law. For q>=1 and 0<=r<2^L,

    w(q*2^L+r)=4^-L*w(q).

After separating the possible small representative r, every residue tail
is a scaled copy of the original source-index distribution. This is an
arithmetic statement about the quotient q, not a claim that shifting a
Collatz trajectory preserves this measure. It explains both the identical
tail in every cylinder and the unequal atomic hot spots.

The source-index entropy is finite and exact:

    H(M)=sum_(k>=0) 2^(-k-1)*(2k+1)=3 bits.

The entropy of the length-L prefix partition increases to3, rather than
growing linearly. Indeed, the tail indicator has probability2^-L. Off
the tail the residue determines M; on the tail, after conditioning on the
residue, the quotient q has the original entropy3. Thus
H(M | prefix)<=h(2^-L)+3*2^-L, tending to0, where h is binary entropy.
Full topological support can therefore coexist with zero entropy per bit,
positive integer atoms, and exceptional positive local dimensions.

Concurrent `63763ac5f`, [the effective-prefix construction](collatz_effective_prefix_mass_20261005.md),
uses a different complete source code, with odd-source weights
v(n)=8/(3*4^bitlength(n)). The two measures have an exact bounded comparison:
for n=2m-1,

    v(2m-1)/w(m) = 4/3 if m is a power of two, and 1/3 otherwise.

Proof: for 2^k<=m<2^(k+1), bitlength(2m-1) is k+1 at m=2^k and
k+2 otherwise, while w(m)=1/(2*4^k). Consequently every set of source
labels B satisfies mu(B)/3<=nu(B)<=4*mu(B)/3. This comparison survives
any common pushforward, including the shared trajectory coordinate. Thus
coverage limits and local dimensions are unchanged by this choice of source
code; the measures are not equal. A separate exact rational check covered
all 8192 source indices 1<=m<=8192.

## 5. Local dimensions away from the integer atoms

Equation(4) gives a small exact multifractal example without assuming an
orbit law. Let r_L be the residue of a 2-adic address m. Whenever r_L>0,
put b_L=bitlength(r_L). Then

    2*4^-b_L <= mu([x_m[:L]]) <= 3*4^-b_L.

The lower and upper local dimensions are therefore respectively
2 liminf(b_L/L) and 2 limsup(b_L/L). For a nonzero address outside the
positive integers, r_L is eventually nonzero and its bitlength is unbounded.
Infinitely many new1bits give upper local dimension2. The lower dimension
can vary with the lengths of the zero runs in the address.

Two explicit computable non-atom controls:

* m=0 is the separate boundary: source n=-1, x_m=111... and mass4^-L.
  Its local dimension is2, although its effective dimension is0. The
  cylinders [1^k] themselves form a mu-null test covering it.
* m=sum_(j>=0) 2^(2^j) has lower local dimension1 and upper local dimension2.
  Between successive new bits b_L=2^j+1; the ratio b_L/L runs from1 down
  toward1/2. Its binary word is computable, so effective dimension is0.

Local measure dimension can exceed the ambient Hausdorff dimension at
exceptional, measure-zero points. It is not the Hausdorff dimension of that
point or of its level set. We make no full multifractal-spectrum claim.

A mixture nu=(1-epsilon)*Haar+epsilon*mu, 0<epsilon<1 computable, retains
positive atoms at every arithmetic input. Those points still have local
dimension0. Mixing in a diffuse observer never makes their cylinders lose
their atomic mass. Generic and atomic measures should be separate accounting
channels, not substituted for one another in a proof.

## 6. The affine multiplier is an exact likelihood ratio

On completed valuation blocks, normalized odd Haar assigns p(a)=2^-a.
Choose instead q(a)=3/4^a, whose sum over a>=1 is1. For a word of r
valuations and total halving cost A, the product laws satisfy

    likelihood(word) = q(word)/p(word) = 3^r/2^A.      (6)

This is exactly its affine multiplier. Under the Haar valuation law these
likelihoods form a nonnegative mean-one martingale, since
sum_(a>=1) 2^-a*(3/2^a)=1. The exact proof and native guard checks are in
the [valuation information package](valuation_information_spectrum_20261005.md).
This is a change of measure on valuation words, not a probabilistic law
for an arbitrary specified integer.

Incoming `5ab42c931`, [the Syracuse information bridge](collatz_information_dimension_bridge_20261005.md),
studies the separate pushforward law on ternary output residues. Its useful
common object is exact: q(a)=p(a)^2/sum_b p(b)^2, because sum_b4^-b=1/3.
Our likelihood law is precisely the normalized quadratic word-energy tilt.
Pushing words to output residues merges their masses, so for s>=1 the
output s-moment is at least the unmerged word s-moment. This preserves an
inequality, not the individual source identities or equality of spectra.
The incoming note's coefficient/payment and numerical/asymptotic scope
are corrected separately; no claimed output-spectrum equality is imported.

If the mean valuation tends to m, the q-local dimension is
2-log2(3)/m. Thus the coefficient boundary m=log2(3) is precisely where
that local dimension equals1. This gives a structural connection between
carry-free coefficient drift and the multifractal observer. The frequency level-set
dimension h(1/m), including the inherited critical value approximately
0.949955, is a third distinct quantity; it is not the point's effective
dimension or the dimension of the nonconvergent positive-integer set.

The actual endpoint is (3^r n+B)/2^A. Strict descent requires

    (2^A-3^r)n > B.                                 (7)

Likelihood retains the coefficient but forgets B and n. Common-future
controllers also retain terminal interfaces. For native L, q/p=81/128
on its exact guard, while its affine multiplier is9/16. The factor9/8
is the output G interface. L and LG have the same source prefix1011100
yet send219 to123 and139 respectively. Equal guard mass does not mean
equal dependency or payment.

## 7. Prefix-free proof search and the missing resource

A prefix-free certificate code is useful because its parsing and total
budget are controlled. A sound machine can accept a source, finite seed
certificates, and a construction address, then output a literal checked
route. A fair dovetail finds any finite certificate in its language that
actually exists. Universal simulation costs a fixed description overhead;
it does not establish that every input has a halting certificate program.

There is an exact limitation of unbounded program-size complexity. Let R(n)
be the canonical first-hit Collatz route, defined only for rooted n. On this
domain,

    K(n,R(n)) = K(n)+O(1).                             (8)

The upper bound uses a fixed partial program that simulates n until1 and
prints the route; the lower bound projects n. Thus enormous trajectories
can have very short programs. Neither (8) nor successful compression gives
a runtime bound or proves that the partial program halts on every n.

The [atomic prefix certificate package](atomic_prefix_certificate_20261005.md)
supplies a concrete fair gamma-weighted schedule and an exact ledger. For a
finite set C of distinct sources with independently checked ROOT certificates,

    covered(C)=sum_(2m-1 in C) w(m),   residual(C)=1-covered(C).

An unverified integer cannot have its atom removed. Family overlap and
multiple descriptions are deduplicated at the source, not added as extra
coverage. Proving residual tends to0 would certify every integer; the
decay estimate remains OPEN. This is a concrete implementation of the
inherited crossroads reduction, not a new termination theorem.

## 8. Vitali and the retained fibre

The [Vitali selector package](vitali_selector_certificate_20261005.md) recovers
one precise finite use in the repo: a tournament statistic is measurable
with respect to its lambda observation exactly when it is constant on every
lambda fibre. The exact seven-vertex hostile has equal labelled lambda but
Hamiltonian-path counts109 and111. Classical Vitali non-Lebesgue-measurability
is a different statement and supplies no noncomputable Collatz oracle.

A second route is the measure/construction split in
[THM-398, lrc-reduction-to-Cprime-and-dominance-dodge](../../01-canon/theorems/THM-398-lrc-reduction-to-Cprime-and-dominance-dodge.md).
Its surviving mechanism is concrete: when no speed is divisible by the
runner count n, the single time1/n is a witness, including possible equality
cases that an open-safe-set measure misses. The multiple-of-n candidate
requires primitive speed normalization: without it n=3,S={3,6} refutes
the old strict-looseness claim. Fixed-radius danger arcs also do not form
an automatic fine Vitali cover. Both scope repairs are recorded in the
source and [HYP-2104](../hypotheses/HYP-2104-lrc-vitali-handoff.md); the
primitive candidate remains OPEN. This is an analogy about retained exact
witnesses, not a transfer of an LRC theorem to Collatz.

For any c.e. equivalence relation on positive integers, the least member
of a class is limit-computable by successively finding smaller certified
equivalents. The sequence decreases only finitely often, but its final
stabilization need not be recognisable. The least selector is computable
iff the equivalence relation is decidable. These are general scoped facts;
we do not claim that Collatz equivalence is undecidable.

In the Collatz common-future relation, reaching the minimum1 is special:
the witness is already a root certificate. An arbitrary temporarily minimal
representative is not a certified terminal. This identifies the stopping
information that a representative-only quotient would erase.

Incoming `9eb19b215`, [two-sheet receipts](collatz_two_sheet_receipts_20261005.md),
provided another exact instance: its sibling automaton BREAKs at x=7,
yet checked routes from21,3,7 give a signed boundary repair for R(3,7).
The incoming iff has been corrected. Losing a finite tracked relation is
not a proof that every possible certificate fails; the source-labelled
verifier is the interface through which a different construction can help.

| source and map | retained predicate | lost information | necessary sidecar |
|---|---|---|---|
| actual valuation word to cylinder mass | exact guard size | order, carry, endpoint | source, affine triple, native guard |
| integer to computable infinite word | all finite orbit observations | a certified time to1 | finite first-hit receipt |
| rooted family to short generation code | decoder and proof transport | literal expansion cost | runtime or bounded verification trace |
| tournament to labelled lambda | triangle-pair counts | Hamiltonian-path count | fibre-discriminating statistic |
| c.e. component to current minimum | checked equivalences | future smaller discoveries | grounded representative or paid rank |

## 9. Recommended proof target and reproducibility

Keep a generic cylinder model to propose useful guards, and an atomic
certificate ledger to judge coverage. The productive target is an adaptive
block lemma: from an unresolved source with checked guard and history, exhibit finitely many
verifiable operations that either ground its source atom or reduce a
well-founded source-labelled rank. Any expected mass gain must preserve
actual input labels, affine carries, endpoint stock, and verification cost.

The concurrent effective-prefix note's P6 gives a sharper alternate target:
positive summable weights f on odd sources other than ROOT with incoming
operator Kf<=rho*f, 0<rho<1. This forces every surviving orbit's weight to
grow geometrically beyond the finite total mass. Its converse currently
constructs f using already terminating routes, so supplies no unconditional
solution. Its P7 hostile family n_t=2^(t+1)-1, U^t(n_t)=2*3^t-1 shows why
no reweighting bounded above and below by constants times v can satisfy
this inequality. The exact comparison above transfers that obstruction to
our gamma weights w. Source mass is a coverage ledger, not already a paid
flow: this route requires an unbounded address-dependent correction or an
additional credit/time coordinate. The needed refuel bound remains OPEN.

The old source27 plateau is the first hostile: payment need not happen
at each trajectory tick. Refuel-tree families are positive controls for
reusing finite assumptions; a repeated ungrounded dependency cycle is the
negative control. The new incoming J rule is an additional candidate
operation, not an assumed universal escape.

Reproduction:

    python -X utf8 -B 04-computation/experiments/algorithmic_collatz_measure_20261005.py

[Script](../../04-computation/experiments/algorithmic_collatz_measure_20261005.py),
[output](algorithmic_collatz_measure_20261005.out).
Universe: all2047 cylinders at lengths0..10, independent residue-block
counts with five additional blocks, 257positive atom controls, six precisions
through512 for signed hostiles, 128nine-valuation renewal comparisons, and
seven sparse-address gap checks. Total10,913exact checks. The code does not
estimate K, empirical randomness, infinite coverage, or a limiting spectrum.
Normal and optimized outputs agree. All four linked proof/code packages
received independent peer audits. Each mathematical infinite statement
above has its explicit proof rather than a finite experiment as justification.
