# Powers of two: decimal windows, long-gap correlation, and a sparse residual

**Status:** PROVED elementary window/count statements; FINITE-EXACT bounded
controls. The assertion that every \(2^n\), \(n>86\), contains a decimal zero
remains OPEN here. No priority claim is made for the elementary lemmas.

**Reproduction:** `python 04-computation/experiments/powers_two_decimal_windows_20261007.py`
and the same command with `python -O`. The matching `.out` is deterministic.
No floating-point calculation decides a digit, phase, count, or test result.

## 1. Inheritance and connection contract

The closest mechanism is the eventual decimal suffix cycle, together with
irrational rotation for leading digits. The hostile is exact recurrence at
large exponent gaps: separation of exponents does not imply independent
digits. The corrected near miss is “all fixed windows have survivors, hence
there are full zero-free powers”; it drops the unobserved middle and source
height. The useful sidecar is the exponent phase itself, not just the number
of surviving suffixes.

Our concept board is: full decimal words; suffix CRT phases; the parity carry
under digit extension; leading-digit rotations; and dimension versus
pointwise orbit hitting. The maps below retain the indicated local digits
and lose all other digits. Restoring an actual exponent and its decimal
length is necessary to test the original question.

[Khovanova's 2011 discussion](https://blog.tanyakhovanova.com/2011/02/86-conjecture/)
states the cutoff conjecture, records the suffix period and first seven
counts, and explains the random-digit heuristic. That heuristic is not used
as a proof. Its historical computational range is not asserted to be the
current record. We independently check the listed witness with exponent
103233492954 by computing only its last 250 digits: the rightmost 249 are
nonzero and the 250th is zero.

The related [Saye paper on ternary digits](https://arxiv.org/abs/2202.13256)
uses trailing-pattern construction for enormous finite exponent bounds.
This supports the choice of a modular search representation; its ternary
conclusions are not imported as decimal conclusions.

## 2. Exact long-gap correlation

For integers \(n\ge0\), \(d>0\), define

\[
z(n,d)=v_{10}(2^{n+d}-2^n).
\]

Then

\[
z(n,d)=\begin{cases}
\min\{n,1+v_5(d)\},&4\mid d,\\
0,&4\nmid d.
\end{cases}
\]

Indeed, the difference is \(2^n(2^d-1)\). Its 2-adic valuation is \(n\).
Modulo 5 the order of 2 is 4; for \(d=4u\), elementary lifting gives

\[
v_5(2^d-1)=v_5(16^u-1)=1+v_5(u)=1+v_5(d).
\]

One can prove the lifting step by writing \(16=1+15\), expanding for a factor 5,
and then removing a factor prime to 5. Thus arbitrary large gaps can force
arbitrarily many identical trailing digits. The formula concerns padded
decimal residues; comparing actual strings also requires that the shorter
number possesses the inspected positions.

For width \(m\ge1\), put \(P_m=4\cdot5^{m-1}\). The suffix indicator

\[
b_m(e)=\mathbf1\{2^n\bmod10^m\text{ has no zero in its }m
\text{ padded digits}\},\qquad n\ge m,\ n\equiv e\pmod{P_m},
\]

is well-defined. Its exact correlation is

\[
C_m(d)=P_m^{-1}\sum_{e\bmod P_m}b_m(e)b_m(e+d).
\]

At every gap divisible by \(P_m\), this equals the one-point density
\(\delta_m\), rather than \(\delta_m^2\). For width 2 and gap 20 these
numbers are \(9/10\) and \(81/100\). At a fixed width the sequence is periodic;
we do not claim arbitrarily long consecutive zero-free runs at that width.

## 3. Every suffix survivor has four or five children

For \(n\ge m\), the image of the suffix cycle is exactly

\[
\mathcal R_m=\{0\le r<10^m:2^m\mid r,\ 5\nmid r\}.
\]

The order of 2 modulo \(5^m\) is \(P_m\), by the same lifting calculation,
and CRT gives both the image and its cardinality \(P_m\).

Let \(r\in\mathcal R_m\) have no zero digit and write \(c=r/2^m\).
A one-digit extension \(r+d10^m\), \(0\le d\le9\), is in
\(\mathcal R_{m+1}\) precisely when \(d\equiv c\pmod2\).
The unit condition modulo 5 is unchanged. Requiring a nonzero new digit
leaves four choices if \(c\) is even, five if it is odd. The five phase lifts

\[
e+jP_m\pmod{P_{m+1}},\qquad 0\le j<5,
\]

biject to all five parity-allowed digits. This is a complete lifting rule,
not an independence approximation.

Consequently every width has surviving phases, indeed \(Q_m\ge4^m\),
where \(Q_m\) counts zero-free suffix residues. Each phase recurs at
arbitrarily large actual exponents. A fixed suffix sieve can never by itself
establish eventual disappearance of zero-free full powers.

The canonical phase may be below the start of its stable cycle. For example,
phase 1 at width 2 means residue 52, not the integer 2 or the padded string 02.
The API `stable_suffix` restores an exponent at least the width before
evaluating its modular power. This boundary is explicitly tested.

## 4. A two-state potential proves a quantitative sparsity bound

The new parity register is

\[
c'=\frac{c+d5^m}{2}.
\]

As allowed digits increase by 2, its parity alternates. An even parent has
two even and two odd children. An odd parent has either two even and three
odd children, or three even and two odd children.

Assign weight 1 to an even register and

\[
a=\frac{1+\sqrt{17}}4>1
\]

to an odd register. Put \(\lambda=(5+\sqrt{17})/2<5\). Then

\[
2+2a=\lambda,\qquad 2+3a=\lambda a,\qquad
3+2a\le2+3a.
\]

At width 1 the residues \(2,4,6,8\) have two registers of each parity,
so their total weight is \(\lambda\). Every lifting multiplies total
weight by at most \(\lambda\). Since each individual weight is at least 1,

\[
4^m\le Q_m\le\lambda^m,
\qquad
\delta_m=\frac{Q_m}{P_m}
\le\frac54\left(\frac\lambda5\right)^m\longrightarrow0.
\]

For a rational certificate, weights 10 and 13 give growth at most \(23/5\).
Their initial total 46 yields the convenient exact bound

\[
Q_m\le(23/5)^m.
\]

Now let \(Z(N)\) count positive exponents \(n\le N\) whose full powers are
zero-free. Choose the least \(m\) with \(P_m>N\). If \(n\ge4m\), then
\(2^n\ge16^m>10^m\), so its inspected suffix consists of actual digits;
also the exponent is in the stable cycle. Distinct such exponents below
\(N\) occupy distinct phases. Hence

\[
Z(N)\le4m+Q_m\le4m+\lambda^m
=O\!\left(N^{\log_5\lambda}\right),
\qquad \log_5\lambda\approx0.9429771<1.
\]

The exponent displayed in decimal is illustrative; all bounds in the code
use the rational certificate. This proves that zero-containing powers have
natural density 1 among exponents. It does not prove that the complement is
finite, much less that its final exponent is 86.

## 5. Even two fixed ends leave infinitely many survivors

Let \(\alpha=\log_{10}2\), which is irrational by prime factorization.
For any surviving suffix phase \(e\pmod{P_m}\), the sequence

\[
\{(e+P_mt)\alpha\},\qquad t=0,1,2,\ldots,
\]

is uniformly distributed modulo 1. For completeness, every nonzero Fourier
mode has mean tending to zero by a finite geometric-series formula,
because \(P_m\alpha\) is irrational; interval approximation then gives
uniform distribution.

A fixed \(k\)-digit leading word \(D\), \(10^{k-1}\le D<10^k\), occurs
on the interval

\[
I_D=[\log_{10}D-(k-1),\ \log_{10}(D+1)-(k-1)).
\]

This is the classical leading-prefix mechanism, also explained in
[Perucca's elementary note](https://antonellaperucca.net/didactics/Powers-of-2.pdf).
If \(D\) has no zero, the set of exponents with leading word \(D\) and a
zero-free suffix of width \(m\) has density exactly

\[
\delta_m\log_{10}\frac{D+1}{D}>0.
\]

Finite initial powers too short for the two windows do not affect density.
Thus inspecting any fixed number of digits at each end cannot close the
conjecture. This statement does not concern windows growing with \(n\).

## 6. Exact next target: a height-sensitive diagonal in the suffix tree

For a full power with actual decimal length \(m\ge2\), its exponent obeys

\[
m\le n<4m<P_m.
\]

The first inequality follows from \(2^n<10^n\); the second from
\(2^{4m}>10^m\); the last holds at \(m=2\) and persists by induction.
Therefore a full zero-free power must appear in the suffix tree at canonical
phase \(e=n\), with the additional exact height condition

\[
10^{m-1}\le2^e<10^m.
\]

This is an exact reformulation, apart from the finitely many one-digit
powers. The broad survivor tree retains many exponentially large phase
values; those are not full counterexamples at the inspected depth.
The next concrete target is a bound on surviving phases in this short
height band, not another count of all suffix phases or another fixed-window
correlation. Our greedy branch reaches 256 zero-free trailing digits at a
590-bit exponent, using only modular arithmetic; its full power is not
constructed and is not claimed zero-free.

## 7. Lawful use of the requested OpenAI paper catalog

The user asked that the catalog's stated results be accepted as premises,
so no proof re-audit is performed. The closest located digit-space tool is
[the entropy-rate dimension formula](https://github.com/openai/math/blob/main/preprints/The-entropy-rate-dimension-formula-for-self-similar-measures-on-the-line-September-24-2026/main.pdf).
Its stated theorem gives self-similar measure dimension from entropy rate
and contraction. Applied to the nine maps \(x\mapsto(x+d)/10\),
\(d=1,\ldots,9\), with equal weights, distinct finite compositions give
entropy rate \(\log9\), so the dimension is \(\log9/\log10<1\).
This is a statement about a measure on infinite digit sequences. A finite
zero-free power supplies a finite prefix, not an infinite zero-free decimal
expansion. No theorem in that statement forces the particular orbit
\(n\log_{10}2\) to leave all growing finite-level digit sets eventually.
The missing transfer is a pointwise, height-sensitive hitting estimate.

## 8. Finite controls, resource accounting, and scope

The script independently matches digit-DP counts against complete suffix
cycles through width 7, and runs digit DP through width 12. It checks all
survivor child transitions through width 6, the exact potential inequalities,
and the long-gap formula on \(0\le n\le40\), \(1\le d\le240\).
Counts at widths 1 through 12 are

\[
4,18,81,364,1638,7371,33170,149268,671701,3022653,13601945,61208743.
\]

The finite census covers every exponent \(0\le n\le100000\). A rolling
residue modulo \(10^{128}\) gives an exact zero witness whenever the window
contains zero. Before the power reaches 128 digits, the residue is checked
without artificial leading padding. Only survivors are expanded in full.
The census finds the known 35 positive exponents and exponent 0; the final
positive one is 86. It uses 36 full expansions, the largest exponent 86 and
largest full decimal string 26 digits. All 99,965 other inputs have explicit
modular zero positions, summarized by a deterministic SHA256 receipt in
the output. An independent full-string scan through exponent 1000 checks
the filtering boundary. No claim is made about the best published census.

The separate modular branch computes 256 suffix digits at each depth at
most, and the published witness check computes 250. The largest possible
full expansion allowed by the declared 100000 census would have 30,103
digits; Python's conversion cap is explicitly set to 40,000, not disabled.
Typed and corrupted-state hostiles include bool/float aliases, an invalid
phase, a padded-zero suffix, and a suffix violating divisibility.

Normal and optimized runs must produce identical output. There are 30,135
explicit checks; none relies on Python `assert`. The completed result is
a proved sparse residual, an exact fixed-observer obstruction, and a small
reproducible finite verification. Universal eventual zero occurrence remains
an additional obligation.
