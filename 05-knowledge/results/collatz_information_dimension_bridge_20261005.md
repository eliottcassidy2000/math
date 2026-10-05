> **TARGETED CORRECTION — 2026-10-05.** The incoming version conflated affine
> coefficient contraction with paid integer descent, a coefficient barrier with
> literal orbit non-descent, and numerical observations with an established
> multifractal spectrum. Those claims are repaired below. The elementary
> moment bound and entropy identity survive. The cited script uses floating
> arithmetic and a valuation cutoff; its large run was not repeated in this
> audit, and its source does not print the six local-mass values attributed to
> it. This is a targeted source/type/status audit, not a new audit of every
> inherited theorem or historical numerical claim.

# Information and dimension for the Syracuse law: exact bounds and numerical spectrum evidence

**Original session:** opus, pascal-lyapunov-20261005, fourth task, 2026-10-05.
The owner's request was to connect Jan Reimann's
[Information vs Dimension: An Algorithmic Perspective, arXiv:2408.05121v1](https://arxiv.org/abs/2408.05121v1)
to the repository's Collatz and Vitali work. The paper is **CITED** for
frequency-set dimension, effective dimension, randomness for computable
measures, and the point-to-set/multifractal framework. Its theorems retain
their stated hypotheses; they do not supply an unproved Collatz spectrum
or a convergence argument.

**Status:** **PROVED** below: the exact-cylinder/coefficient identity, the
carry-sensitive payment criterion, the moment upper bound for real
\(q\ge1\), the entropy identity, and the computable-code/non-atomic-randomness
statements. **NUMERICAL:** the inherited floating-point, cutoff-60
partition-function table. **UNREPRODUCED OBSERVATIONS:** the six reported
local-mass values. **CONJECTURED / CONDITIONAL:** equality in the moment
bounds and an identification of the formal Legendre branch with the actual
measure spectrum. P7 remains empirical support. Collatz is **OPEN**.

## 0. Recovery and the surviving connection

The relevant inherited routes are:

* [THM-4512, coefficient-descent-classes-one-member](../../01-canon/theorems/THM-4512-coefficient-descent-classes-one-member.md),
  which already requires attention to the carry and smallest member;
* [collatz_directions_20260926.md](collatz_directions_20260926.md),
  Proposition 1, for the closed **coefficient** barrier set and its
  box-counting dimension, using all-prefix word counts;
* [THM-4476, thin-divergent-orbits-reciprocal-sums-finite](../../01-canon/theorems/THM-4476-thin-divergent-orbits-reciprocal-sums-finite.md)
  and [THM-4499, thin-divergence-is-little-o-of-x-to-the-h-star](../../01-canon/theorems/THM-4499-thin-divergence-is-little-o-of-x-to-the-h-star.md),
  for separately proved arithmetic counting bounds;
* [collatz_fourier_experiments_20261004.md](collatz_fourier_experiments_20261004.md),
  especially its explicit P7 status: finite numerical support, not a
  proved uniform \(L^q\) bound;
* [valuation_information_spectrum_20261005.md](valuation_information_spectrum_20261005.md),
  for the exact odd-cylinder metric, renewal-code spectrum, native paid
  guards, and their likelihood ratios;
* [algorithmic_collatz_measure_20261005.md](algorithmic_collatz_measure_20261005.md)
  and [vitali_selector_certificate_20261005.md](vitali_selector_certificate_20261005.md),
  for the atomic-input measure and the distinction between a finite
  certificate and a global equivalence selector.

The useful comparison is exact but limited: the real affine coefficient is
a ratio of cylinder masses. The additive carry and the supplied integer
still decide payment. Similarly, merging words into residue fibres gives
an exact one-sided moment inequality; the merge does not automatically
preserve dimension or the full multifractal spectrum.

## 1. Frequency sets are not all-prefix barriers

Write the shortcut parity word of a 2-adic integer as \(b_0b_1\ldots\), and
let \(o_j=\sum_{i<j}b_i\). The inherited set is

\[
 E_{\rm coef}=\{x:o_j(x)\log_2 3-j\ge0\text{ for every }j\ge1\}.
\]

It is a coefficient barrier. For positive integers its membership prevents
strict literal descent, because every affine prefix has nonnegative carry.
The converse is not established by this definition: coefficient contraction
can occur without payment. The root \(1\), continued through its usual
cycle, already supplies a boundary example.

Eggleston's theorem concerns a **limiting frequency** set, with Hausdorff
dimension \(H(p)\). The inherited all-prefix count gives the **box**
dimension \(h_*=H(\log_3 2)\) of \(E_{\rm coef}\). These are distinct
statements sharing the same entropy constant. The all-prefix result is
not a direct identification with Eggleston's limiting-frequency set, and
this note does not establish a Hausdorff-dimension equality for the
barrier. The separate integer thin-divergence bounds also require their
own arithmetic cylinder-counting arguments; they are not an automatic
mass-distribution transfer to the integers.

## 2. Proposition 1: the coefficient ratio and exact payment

Let \(w=(a_1,\ldots,a_k)\), \(a_i\ge1\), be an exact accelerated valuation
word, let \(S=\sum_i a_i\), \(P=3^k\), \(Q=2^S\), and write

\[
 F_w(n)=\frac{Pn+B_w}{Q},\qquad B_w>0\quad(k>0).
\]

Among odd 2-adic sources, the exact word is the single residue class

\[
 Pn+B_w\equiv Q\pmod{2Q}.                                      \tag{1}
\]

Thus the ambient modulus is \(2^{S+1}\), not \(2^S\): the extra bit ensures
that the terminal output is odd, hence that the last valuation is exact.
Its mass under **normalized odd Haar measure** is \(2^{-S}\). Under Haar
measure on all of \(\mathbb Z_2\), its mass is \(2^{-S-1}\).

Every integer endpoint of this word lies in the same 3-adic cylinder

\[
 F_w(n)\equiv B_wQ^{-1}\pmod{3^k},
\]

of normalized \(\mathbb Z_3\)-Haar mass \(3^{-k}\). As an arithmetic
progression the endpoint set is dense in that cylinder. This does not
assert that the two measures live on the same space or that an arbitrary
2-adic endpoint has a canonical 3-adic value; the finite word supplies
the residue observable.

**Proposition 1 (repaired).** Coefficient contraction is exactly

\[
 \frac PQ<1
 \quad\Longleftrightarrow\quad 2^{-S}<3^{-k}.                   \tag{2}
\]

For a supplied positive integer satisfying (1), actual payment is exactly

\[
 F_w(n)<n
 \quad\Longleftrightarrow\quad (Q-P)n>B_w.                      \tag{3}
\]

If \(Q>P\), all sufficiently large positive members of the cylinder pay;
the threshold is \(n>B_w/(Q-P)\). Payment for **every** positive member
requires this inequality at its least positive representative. If
\(P\ge Q\), no positive member pays. These conclusions follow by subtracting
\(n\) from the displayed affine formula.

**Hostile.** The word \(w=(2)\) has \(P=3,Q=4,B=1\), exact sources
\(n\equiv1\pmod8\), and contracting coefficient \(3/4\). Yet \(F_w(1)=1\).
This refutes the incoming claim that coefficient contraction is equivalent
to strict descent at every member.

The condition \(2^{-S_j}\ge3^{-j}\) for every completed odd prefix defines
the corresponding coefficient barrier; it must not be renamed literal
non-descent. No real order is assigned to general 2-adic sources. A typical
Haar valuation sequence has \(S_j/j\to2\), so its coefficients contract
exponentially, but the mass comparison alone is not a paid integer receipt.

## 3. A proved moment upper bound, not an established spectrum

Let \(\mu_n\) be the word-induced law sending the first \(n\) valuations
of a normalized Haar-random odd source to \(B_w2^{-S_n(w)}\bmod3^n\).
For integer members this is the endpoint residue above. The valuations
are independent with probability \(2^{-a}\), \(a\ge1\).
Equivalently, the compatible laws are the
reductions of the 3-adic random series

\[
 Y=\sum_{j\ge1}3^{j-1}2^{-(a_1+\cdots+a_j)}.
\]

To obtain this version from the forward endpoint expression, reverse the
finite independent valuation word. Reversal preserves its law. The series
converges 3-adically, so a limiting measure \(\mu\) is well defined.

Set

\[
 Z_n(q)=\sum_c\mu_n(c)^q,\qquad
 \tau(q)=\liminf_{n\to\infty}\frac{-\log Z_n(q)}{n\log3}.
\]

For \(q>1\) write \(D_q=\tau(q)/(q-1)\). This quotient does not define
\(D_1\); at \(q=1\), \(Z_n(1)=1\) and \(\tau(1)=0\).

**Theorem 3.** For every real \(q\ge1\),

\[
 \tau(q)\le\min\{q-1,\log_3(2^q-1)\}.                          \tag{4}
\]

**Proof.** For a word \(w\) of \(n\) valuations let \(Y_n(w)\) denote its
endpoint residue. The exact word mass is \(2^{-S_n(w)}\). Since
\((\sum_i m_i)^q\ge\sum_i m_i^q\),

\[
 Z_n(q)\ge\sum_w2^{-qS_n(w)}
 =\left(\sum_{a\ge1}2^{-qa}\right)^n
 =(2^q-1)^{-n}.
\]

There are at most \(N_n=2\cdot3^{n-1}\) residues of positive mass. Convexity,
or Jensen after adding zero entries, gives
\(Z_n(q)\ge N_n^{\,1-q}\). Taking the stated liminf proves (4).
Countable word sums are justified by monotone convergence. \(\square\)

The two bounds agree at **both \(q=1\) and \(q=2\)**. For a direct check,
set \(t=2^{q-1}\); the equality is
\(2t-1=t^{\log_2 3}\), whose strict concavity after moving the power term
gives exactly the two roots \(t=1,2\) for \(t\ge1\). The dimension bound
\(q-1\) is smaller for \(1<q<2\); the word-energy bound is smaller for
\(q>2\). Thus \(q=2\) is the nontrivial crossing.

In particular \(D_q\le1\) for \(q>1\). If \(D_\infty\) is defined as the
limit of these nonincreasing generalized dimensions, (4) gives
\(D_\infty\le\log_3 2\); equality is not proved here.

On the unit residues put \(\rho_n=N_n\mu_n\). The exact relation is
\(Z_n(q)=N_n^{1-q}\mathbb E_{\rm uniform}[\rho_n^q]\).
A proved uniform bound on this expectation would imply
\(\tau(q)=q-1\). P7 supplies numerical evidence for such bounds when
\(1<q<2\), not an all-level proof. Likewise a proved subexponential second
moment would imply \(\tau(2)=1\); observed linear growth is not that proof.

**Conjectural branch.** Equality
\(\tau(q)=\log_3(2^q-1)\) for \(q\ge2\) remains conjectural. Even if this
equality were proved, identifying a Legendre transform with the actual
Hausdorff spectrum would require a multifractal-formalism argument.
Word-to-residue collisions and dimension transport cannot be discarded.

### Inherited numerical table, not finite-exact arithmetic

The original report supplied the following per-level slopes
\(-\log(Z_n/Z_{n-1})/\log3\). These are not the limiting definition of
\(\tau\), and (4) is not a finite-level slope bound. The column headed
reference branch is the conjectural expression; its \(q=0.5\) row lies
outside Theorem 3.

| \(q\) | reference branch | slope 14 | slope 15 | slope 16 | reference \(D_q\) |
|---|---|---|---|---|---|
| 0.5 | -0.5000 | -0.4978 | -0.4981 | -0.4983 | 1 |
| 1.5 | 0.5000 | 0.4855 | 0.4868 | 0.4880 | 1 |
| 1.9 | 0.9000 | 0.8546 | 0.8575 | 0.8601 | 1 |
| 2.0 | 1.0000 | 0.9430 | 0.9463 | 0.9492 | 1 |
| 2.5 | 1.4003 | 1.3586 | 1.3630 | 1.3669 | 0.9335 |
| 3.0 | 1.7712 | 1.7356 | 1.7401 | 1.7440 | 0.8856 |
| 4.0 | 2.4650 | 2.4257 | 2.4297 | 2.4332 | 0.8217 |
| 6.0 | 3.7712 | 3.7260 | 3.7291 | 3.7317 | 0.7542 |
| 8.0 | 5.0439 | 5.0038 | 5.0060 | 5.0078 | 0.7206 |

The source
[collatz_syracuse_multifractal_20261005.py](../../04-computation/experiments/collatz_syracuse_multifractal_20261005.py)
initializes NumPy floating arrays, sums valuations only through 60, and
uses floating powers and logarithms. In ideal exact arithmetic this cutoff
retains total mass \((1-2^{-60})^n\), so the omitted mass is at most
\(n2^{-60}\); this is not a floating-roundoff or spectrum-error certificate.
Rounded total mass near one does not make the run exact. The values are
retained as attributed **NUMERICAL** evidence, without a fresh large rerun
or an inference that the limiting gap closes.

### Proposition 4: the exact tilt identity and its boundary

For binary entropy \(H\),

\[
 2-H(3/4)=\tfrac34\log_2 3.                                   \tag{5}
\]

More generally
\(2-H(p)-p\log_2 3=D_{\rm KL,2}(p\Vert3/4)\ge0\), with equality only at
\(p=3/4\). Thus \((2-H(p))/(p\log_2 3)\ge1\) for \(0<p<1\), with tangency
at \(3/4\). This is an elementary proved identity.

Squaring and normalizing the source valuation weights gives

\[
 \frac{(2^{-a})^2}{\sum_{b\ge1}(2^{-b})^2}=3\,4^{-a},
\]

whose mean is \(4/3\) and whose induced binary one-frequency is \(3/4\).
For the **candidate** right-hand pressure branch, formal differentiation
gives
\(\alpha(q)=2^q\ln2/((2^q-1)\ln3)\) and the formal value
\(q\alpha(q)-\log_3(2^q-1)\). Its right limit at \(q=2\) is
\(\alpha=4/(3\log_2 3)\), with formal value
\(H(3/4)/((3/4)\log_2 3)\).
These algebraic formulas and the tilted word law do not by themselves prove
local dimensions or set dimensions after the residue quotient. In
particular the left candidate derivative at 2 is 1 and the right candidate
derivative differs; a single differentiable Legendre point is not asserted.

## 4. Periodic endpoints and finite-level concentration observations

A forward word starting at \(m\) supplies mass at its endpoint residue,
not generally at the fixed residue \(m\). The incoming blanket fixed-point
bound from the forward orbit therefore needs an endpoint hypothesis.

There is a valid bound for a periodic odd integer point \(m\). For each
\(n\), take the cyclic predecessor word of length \(n\) ending at \(m\).
If its cost is \(S_n\), then \(\mu_n(m\bmod3^n)\ge2^{-S_n}\).
If the cycle has average valuation \(\bar a\), this proves

\[
 \limsup_n\frac{-\log\mu_n(m\bmod3^n)}{n\log3}
 \le \bar a\log_3 2.                                         \tag{6}
\]

It applies to the fixed point \(-1\), the cycle through \(-5\), the cycle
through \(-17\), and the root fixed point \(1\) under accelerated iteration.
It does not transfer a transient forward orbit bound to the fixed points
5 or 7.

| point | periodic data / endpoint scope | proved upper bound (6), rounded | reported level-14 estimate |
|---|---|---|---|
| -1 | fixed, valuation 1 | 0.631 | 0.606 |
| -5 | 2-cycle, cost 3 | 0.946 | 0.873 |
| -17 | 7-cycle, cost 11 | 0.992 | 0.917 |
| 1 | fixed, valuation 2 | 1.262 | 1.049 |
| 5, 7 | transient to 1; no bound here at their fixed residues | — | 1.045, 1.108 |

The last column is **UNREPRODUCED NUMERICAL OBSERVATION** from the incoming
note. The listed script prints the largest residue mass and its location;
it does not print these six selected residue masses. These are finite-level
estimates, not established local-dimension limits. Nor does a value near
one prove a finite positive density at that point.

The terminology of an atom at finite level means a residue probability.
It must not be confused with a positive-mass atom of the limiting 3-adic
measure. In particular the conjectured estimate
\(\mu_n(-1)\asymp2^{-n}\) would tend to zero, not exhibit such an atom.
Other periodic points also have positive cylinder masses at every finite
level; \(-1\) is not the sole computable point detected by this measure.

## 5. Computable inputs: the valid boundary, not a universal no-go theorem

**Remark 5 (repaired proof).** For every fixed integer \(m\), its full
shortcut parity sequence \(x(m)\) is computable by iterating the total
integer map, without deciding eventual termination. Hence

\[
 K(x(m){\upharpoonright}L)\le K(m)+O(\log L),
\]

and both lower and upper effective dimensions are zero. This is effective
point dimension, not a claim about the local dimension of every measure
at that point or the dimension of a larger set containing it.

For a computable probability measure \(\nu\), a computable point \(x\) of
zero atomic mass is not \(\nu\)-Martin-Löf random. The proof requires more
than comparing the displayed logarithmic upper bound with a surprisal
which merely tends to infinity. For each \(j\), effectively search for
\(L_j\) with
\(\nu([x{\upharpoonright}L_j])<2^{-j}\), using computable upper
approximations to cylinder masses. Zero atomic mass guarantees termination.
The resulting uniformly effective cylinders form a Martin-Löf test
covering \(x\). This applies in particular to every computable non-atomic
measure.

It does not rule out atomic computable measures: a positive-mass atom is
random for that measure, and its local dimension is zero. The
[atomic integer measure](algorithmic_collatz_measure_20261005.md) and the
[concurrent source-code construction](collatz_effective_prefix_mass_20261005.md) make
every positive-integer parity path such an atom. Thus the statement that
randomness can never include Collatz inputs is false without its
non-atomic/zero-atom qualification.

Almost-everywhere results alone do not settle membership of a specified
integer in an exceptional set. That logical limitation is not a theorem
that no argument using dimension, moments, or randomness can be combined
with additional arithmetic information to settle it. Nor are all
multifractal statements solely statements about random points; spectra
can describe nonrandom exceptional sets. The incoming universal
methodological prohibition is withdrawn.

Finally, eventual periodicity of every positive integer orbit would
exclude divergent trajectories but would still allow a nonroot positive
cycle. It is therefore weaker than full Collatz convergence to 1.
A finite witness reaching ROOT is a different object from an oracle or
a global least-representative selector; neither is supplied merely by
invoking Vitali terminology.

## 6. Reproduction and source-status limits

The inherited commands are:

    python 04-computation/experiments/collatz_syracuse_multifractal_20261005.py 16 60
    python 04-computation/experiments/collatz_syracuse_multifractal_20261005.py 13 60

The original reported level-16 run uses roughly 43 million residues.
Its source writes a floating numerical output beside the script. This
targeted correction inspected that source without repeating the large
run or changing its numerical algorithm or historical output. Script labels
now identify the truncated floating computation and the finite-level
maximum; its formal infinite recursion is exact, while its implementation
uses truncated floating arithmetic.
Its selected-six-point local values require a separate retained
calculation to be reproducible.

## 7. Verdicts and the next proof obligation

| claim | repaired status |
|---|---|
| \(h_*=H(\log_3 2)\) occurs in both frequency entropy and inherited all-prefix counting | CITED mechanisms with different sets/metrics and proof obligations |
| coefficient contraction iff \(2^{-S}<3^{-k}\) | PROVED with normalized odd Haar |
| actual integer payment iff \((Q-P)n>B_w\), on the exact source cylinder | PROVED; coefficient contraction alone is insufficient |
| moment bound (4), with equality of its two bounds at \(q=1,2\) | PROVED for \(q\ge1\) |
| entropy identity (5) and the normalized quadratic tilt | PROVED |
| P7 and the finite partition slopes | NUMERICAL support, not uniform all-level estimates |
| actual equality of the candidate pressure branch and its Legendre spectrum | CONJECTURED / additional formalism required |
| periodic-point local upper bounds (6) | PROVED; six reported finite-level values remain unreproduced observations |
| computable parity sequences have effective dimension zero | PROVED, irrespective of convergence |
| computable zero-atom points are nonrandom for computable measures | PROVED by an effective cylinder test; atomic measures are an explicit exception |
| a universal prohibition on measure/dimension methods proving Collatz | NOT ESTABLISHED; withdrawn |
| Collatz convergence | OPEN |

For paid guards, the concrete retained interface is the affine coefficient,
the additive carry, the exact native source guard, and the supplied source.
The quadratic word tilt explains the coefficient and moment formulas.
A payment receipt still requires inequality (3) and a grounded dependency;
no spectrum estimate substitutes for those data.
