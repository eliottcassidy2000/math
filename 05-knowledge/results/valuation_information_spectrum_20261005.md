# Valuation information, affine likelihood, and the missing payment data

2026-10-05. **PROVED** elementary coding, local-dimension formulas, and
native-guard likelihood identities below. **INHERITED/CITED** where marked;
**FINITE-EXACT** for the declared experiments; decimal entropy values are
**NUMERICAL**. No probability-one or dimension statement here proves
Collatz for a supplied positive integer.

The useful bridge is exact: a particular Bernoulli likelihood ratio equals
the carry-free affine coefficient of a completed valuation word. For the
paid controller, terminal interface bits supply explicit corrections.
The affine carry and actual source remain necessary to prove payment.

## 1. Inheritance and distinct meanings of dimension

The nearest established result is
[THM-4487, dip-spectrum-entropy-curve-and-sharpness-of-thin-divergence](../../01-canon/theorems/THM-4487-dip-spectrum-entropy-curve-and-sharpness-of-thin-divergence.md).
Its critical entropy \(h(\log_3 2)=0.949955527\ldots\) and its counting
meaning are already proved. The
[sibling dimension ladder](collatz_procgen_20260922_sibling_dimension_ladder.md)
also proves that dimension for the all-prefix multiplier barrier in the
2-adic parity model. The present mean-level sets are different sets:
recovering the same entropy function is not a new no-dip or divergence theorem.

The exact native guard and carry decoders are inherited from
[translation_phase_decoder_20261005.md](translation_phase_decoder_20261005.md).
The incoming J extension, its actual paths and its two forbidden pairs are
Propositions J1–J2 of
[collatz_fusion_helpers_20261005.md](collatz_fusion_helpers_20261005.md),
from commit 358bd7ee5. We do not reprove its funding potential.

The canonical hostile is a computable arithmetic trajectory: its prefixes
have short algorithmic descriptions even when a probabilistic code regards
them as expensive. The corrected near miss is to identify the dimension of
a set, the local dimension of a measure at a point, and the effective
dimension of that point. The least-used sidecar is the terminal guard
interface, in addition to the full affine carry.
Our board is code / Haar mass / tilted likelihood / carry / native interface /
immutable source.

**CITED context.** Reimann's
[Information vs Dimension — an Algorithmic Perspective, v1](https://arxiv.org/html/2408.05121v1)
develops prefix codes and the Cantor metric in section 2.1, effective
dimension in section 4, and local dimensions in section 5. We use these
notions with a specified measure and metric. The special Bernoulli formulas
below have direct proofs; no general multifractal identification is assumed.

## 2. The valuation code is an exact normalized 2-adic coordinate

Write \(U(x)=(3x+1)/2^{v_2(3x+1)}\) on odd 2-adic inputs wherever defined.
For a finite valuation word \(w=(a_1,\ldots,a_r)\), put
\[
 A=\sum_i a_i,\qquad
 F_w(x)=\frac{3^r x+B_w}{2^A}.
\]
Its exact odd-source cylinder is
\[
 C_w=\{x:3^r x+B_w\equiv2^A\pmod{2^{A+1}}\}.
\]
There is exactly one odd residue because \(3^r\) is a 2-adic unit. The
final oddness bit is essential. Divisibility by the full denominator gives
each earlier division, and an even intermediate nominal state would make
the next numerator odd, contradicting that division. Thus this is precisely
the actual word cylinder.

Let \(\lambda\) be normalized Haar measure on the odd 2-adics. Then
\[
 \boxed{\lambda(C_w)=2^{-A}.}                         \tag{1}
\]
Encode one valuation by \(a\mapsto0^{a-1}1\). This is a complete prefix code,
with Kraft sum \(\sum_{a\ge1}2^{-a}=1\). Its concatenation is exactly the
shortcut-map parity sequence after dropping the initial odd bit.

On odd sources use \(d_o(x,y)=2^{1-v_2(x-y)}\). The code is an isometry
onto binary sequences with infinitely many ones. Indeed if two valuation
sequences first differ after common cost \(A_0\), with distinct next
valuations a,b, their sources have
\[
 v_2(x-y)=A_0+\min(a,b).
\]
This follows from the inverse branch \((2^a u-1)/3\), with u odd, and
composition through the common prefix. Their codes agree for exactly
\(A_0+\min(a,b)-1\) bits, proving the metric identity.
Every infinite valuation sequence gives a unique point by nested cylinders.

The excluded countable set consists of \(-1/3\) and its finite odd
preimages, where a subsequent valuation is infinite. Adding their
eventually-zero codes extends the finite-cylinder coordinate to all odd
2-adics. Positive odd integers never meet this undefined state. Continuing
at 1 gives valuations \(2,2,\ldots\), not a terminated code.

## 3. Geometric valuations and an explicit multifractal spectrum

Fix \(0<z<1\) and give independent valuations the probabilities
\[
 p_z(a)=(1-z)z^{a-1}.
\]
In the unary coordinate this is simply the Bernoulli bit measure with
zero probability z and one probability \(1-z\); denote its pullback by
\(\mu_z\). For a completed word,
\[
 \mu_z(C_w)=(1-z)^r z^{A-r}.                          \tag{2}
\]
If \(A_r/r\to m\in[1,\infty)\), the binary one-frequency is \(1/m\).
This holds between valuation boundaries too: finite convergence of the
mean implies \(a_{r+1}/A_r\to0\). For \(m=\infty\), the one-frequency
is zero directly. Therefore the measure-local dimension is
\[
 d_z(m)=
 -\frac1m\log_2(1-z)-\left(1-\frac1m\right)\log_2 z,    \tag{3}
\]
where \(1/\infty=0\).

Let \(E_m\) be the set of these codes with limiting mean m.
Its Hausdorff dimension in the stated metric is
\[
 \dim_H E_m=h(1/m)
 =\frac{m\log_2m-(m-1)\log_2(m-1)}m                 \tag{4}
\]
for \(1<m<\infty\), with both endpoint dimensions zero.
For completeness, prefixes with one-frequency near \(\theta\) have at most
\((N+1)2^{N\sup h}\) possibilities; this gives the upper bound after
taking a shrinking frequency interval. Bernoulli(\(\theta\)) gives full
measure to the frequency-\(\theta\) set by the strong law, and its typical
cylinder information is \(Nh(\theta)+o(N)\). Restricting to a positive
measure subset on which this bound is uniform gives the matching
mass-distribution lower bound. The endpoint upper bounds tend to zero.
This is the standard binary frequency-set argument, also presented in
Reimann's section 4.1.

For \(z\ne1/2\), (3) is affine and one-to-one in the one-frequency, so
(3)–(4) give the spectrum restricted to points whose local-dimension
limit exists. Existence of that limit is equivalent to existence of the
bit frequency. No spectrum claim is made here for lower local dimension
at points with an irregular frequency; Reimann also considers that broader
notion.
For \(z=1/2\), every point has Haar-local dimension 1, and the sole local
dimension level has Hausdorff dimension 1. The several \(E_m\) sets must
not be called different Haar-local-dimension levels.

Under \(\mu_z\), almost every mean is \(m=1/(1-z)\), and the typical
local dimension is \(h(1-z)\). Equation (3) at an atypical point can exceed
the ambient dimension 1. This does not say the dimension of a subset
exceeds 1: it measures the decay of this particular measure at that point.

## 4. The affine coefficient is a likelihood ratio

Take \(\nu=\mu_{1/4}\). Then \(p_{1/4}(a)=3/4^a\), and (1)–(2) give
\[
 \boxed{\frac{\nu(C_w)}{\lambda(C_w)}
       =\frac{3^r}{2^A}.}                            \tag{5}
\]
This is exactly the real affine coefficient of \(F_w\).
Consequently that coefficient is a likelihood-ratio martingale on
completed valuation prefixes under odd Haar measure:
\[
 \mathbb E_\lambda\!\left[\frac3{2^a}\right]
 =\sum_{a\ge1}\frac3{4^a}=1.
\]
The equality refers to completed blocks. At an unfinished unary block,
the dropped-initial-bit convention can differ by one odd step from the
literal shortcut coefficient; keep that endpoint bit when comparing them.

In this tilted measure,
\[
 d_\nu(m)=2-\frac{\log_2 3}{m},\qquad
 \frac1r\log_2\!\left(\frac{3^r}{2^{A_r}}\right)
 \longrightarrow\log_2 3-m.
\]
The coefficient-neutral mean \(m=\log_2 3\) is exactly the measure-local
dimension level \(d_\nu=1\). Its Hausdorff dimension is the inherited
\(h(\log_3 2)=0.949955527\ldots\). This explains the shared entropy value
through an exact change of measure.

It does not identify coefficient drift with actual size drift.
The carry-free factor omits \(B_w\). For example valuations \(1,2\)
and \(2,1\) have the same mass, rank and coefficient but carries 5 and 7.
At the positive root 1, repeated valuation 2 has coefficient \((3/4)^r\)
while its carry keeps the state exactly 1. At the negative fixed point
\(-1\), repeated valuation 1 has expanding coefficient \((3/2)^r\)
while its carry again keeps the state fixed. Thus \(-1\) is a zero-size-
drift hostile, not a coefficient-neutral example.

## 5. Native controller interfaces supply the exact corrections

Use the inherited affine letters H,G,A,B,L and incoming J. Their native
prefix actions on the same binary coordinate are:

Here a guard's Haar mass means the full 2-adic congruence cylinder. It
also equals the relative density of that congruence among positive odd
integers; the countable set of its literal positive-integer members itself
has Haar measure zero.

| Letter | Action on the child's code | Relative bits / ones |
|---|---|---:|
| G | prepend 101 | 3 / 2 |
| A | prepend 1011001 | 7 / 4 |
| B | prepend 1101001 | 7 / 4 |
| H | prepend 1011110100 | 10 / 6 |
| L | replace initial 101 by 1011100 | 4 / 2 |
| J | replace initial 1 by 111101001 | 8 / 5 |

These are identities of actual common-future words, not arbitrary binary
substitutions. For H, its six-step word \(V=(1,2,1,1,1,2)\) is followed
by a valuation two greater than the child's first valuation; hence
code(V) followed by 00 is prepended. For L, the child's word begins
\((1,2,a,\ldots)\) and the source uses
\((1,2,1,1,a+2,\ldots)\), giving the stated replacement.
For J, Proposition J1 gives source word \((1,1,1,1,2,3)\) and child
word (1), giving the last row.

Read a controller word backward through these prefix interfaces.
The only incompatible adjacent letters are LB and LJ: L demands initial
101, whereas B starts 110 and J starts 111. All other native prefixes
begin with the required interface. This recovers the inherited exact
language with no probabilistic assumption.

For a nonempty legal word with composite
\[
 F(n)=\frac{3^R n+B}{2^A},
\]
let \(e_L,e_J\) indicate whether its last letter is L or J.
Its complete native source prefix has
\[
 N=A+3e_L+e_J\quad\hbox{bits},\qquad
 O=R+2e_L+e_J\quad\hbox{ones}.                        \tag{6}
\]
The additional bits are the unconsumed output interface: 101 for L
and 1 for J. Induction through the table proves (6). Hence
\[
 \lambda(\text{native guard})=2^{-N},\qquad
 \boxed{\frac{\nu(\text{native guard})}{\lambda(\text{native guard})}
  =\frac{3^R}{2^A}(9/8)^{e_L}(3/2)^{e_J}.}           \tag{7}
\]
The factors are forced by actual interface codes. They do not, by
themselves, prove the controller's credit or payment inequality.

**Loss hostile.** L and LG both have native source prefix 1011100,
and therefore exactly the same guard and both measure masses.
Their affine maps differ:
\[
 L(n)=\frac{9n-3}{16},\qquad LG(n)=\frac{81n+53}{128}.
\]
At 219 they produce 123 and 139. The terminal L factor in (7) accounts
exactly for the extra G coefficient in LG. A guard prefix alone cannot
recover a controller history or endpoint; the retained ordered translation
decoder provides the missing data.

## 6. An actionable guard-information ledger

For each proposed receipt retain its exact native cylinder, full
\((P,Q,B)\), terminal interface, actual source, and its proof/potential
labels. The information cost \(-\log_2\lambda(\text{guard})\) measures
how many binary source bits are required. It is not the same as an
algorithmic description length or a paid size decrease.
On a valid receipt the exact endpoint payment is
\[
 F(n)<n\quad\Longleftrightarrow\quad (Q-P)n>B.
\]
Neither the cylinder mass nor its likelihood ratio retains the right-hand
side or the actual source.

Overlaps must also be removed before turning a list of guards into a
coverage figure. For the **five primitive** H/G/A/B/L native guards,
the naive mass sum is \(153/1024\), whereas their union is
\[
 \{11\bmod16\}\ \cup\ \{7\bmod256\},
 \qquad \lambda=17/128.
\]
G contains H, A and L. If one instead asks for the primitive first-letter
funders H/A/B/L, their disjoint union has mass \(25/1024\).
These are scopes for those named primitive guards only: G alone is not a
paid start, and these numbers neither replace the larger paid libraries
nor add to their published coverage. The script uses exact dyadic
intersections and a disjoint union ledger rather than summing overlaps.

Formula-only membership is weaker as well: the numerical value
\((9\cdot27-3)/16=15\) decreases, but 27 is not in L's native guard.
This arithmetic equality is not an L common-future receipt.

## 7. Why the generic dimension theorem does not settle integers

Fix any positive odd integer n. Its infinite shortcut parity sequence,
continued through the root cycle if reached, is computable whether or
not its orbit reaches 1: each requested finite prefix takes finitely
many integer operations. Its unary valuation sequence is the same code.
For every prefix length N,
\[
 K(\text{code}(n)\upharpoonright N)
 \le K(n)+K(N)+O(1)=O_n(\log N).
\]
Thus both lower and upper effective dimensions of this fixed code are zero.
No convergence premise was used.

For a concrete distinction, 1 has code \((01)^\infty\) and mean 2.
It is in the frequency set \(E_2\) of Hausdorff dimension 1; its Haar-local
dimension is 1; its own effective dimension is 0.
Under \(\nu\) its measure-local dimension is
\(2-\frac12\log_2 3=1.207518749\ldots\), still with effective dimension zero.

Randomness is relative to a measure. In particular, this argument does
not assert that a computable integer code is nonrandom for every computable
measure: an atomic measure can make it random. It shows why replacing a
supplied arithmetic input by a Haar- or positive-entropy generic point
is unjustified. The separate atomic-measure analysis makes that boundary
explicit. Counts, mass bounds, and likelihood can guide a search or audit
guard costs; they do not discharge a ROOT obligation.

## 8. Reproduction and finite universe

Run:

    python -B 04-computation/experiments/valuation_information_spectrum_20261005.py
    python -B -O 04-computation/experiments/valuation_information_spectrum_20261005.py

The [script](../../04-computation/experiments/valuation_information_spectrum_20261005.py)
and [output](valuation_information_spectrum_20261005.out) check:

- all binary source-prefix bijections at depths 1–10 (2046 sources), and
  32640 independent metric pairs;
- all 2047 positive valuation words of total cost 1–11, three literal
  cylinder representatives each, and 6141 exact geometric probability tests;
- all 780 five-letter words of lengths 1–4 (695 native) by independent
  congruence intersections, source codes, and interface likelihoods;
- all 1554 six-letter words of lengths 1–4 (1316 native), including the J
  interface and both forbidden pairs, plus 32 literal J common-future joins;
- independent primitive union counting, 19 exact martingale-tail controls,
  the carry/interface/source hostiles, and numerical values of the proved
  spectrum formulas at six declared means.

The first seven items in the output before the numerical table use only
integer or rational arithmetic. Floating values illustrate formulas and
are not used to decide an inequality in a proof. Normal and optimized
execution agree; there is no claimed finite test of Kolmogorov complexity
or of an infinite frequency limit.
