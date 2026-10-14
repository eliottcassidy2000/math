# Continued fractions can retain every Collatz state while their exponent forgets the orbit

2026-10-14. **PROVED:** local-growth exponent theorem, faithful state/clock
codec, exact native-cell consumer, and the engineering corrections below.
**FINITE-EXACT:** the declared controls, not an estimate of an infinite
irrationality exponent. **OPEN:** convergence of every Collatz source.

## 1. Inheritance, carriers and the question

The closest proved mechanism is
[THM-4056, divisor phase compiler, §E.5](../../01-canon/theorems/THM-4056-divisor-phase-compiler-and-duffin-schaeffer-lcm-clock.md):
the irrationality exponent is a limsup involving **convergent denominators**,
not just the next partial quotient. Its hostile is a fixed finite prefix
with both exponent-two and Liouville continuations. The corrected near miss
is reading an almost-everywhere digit law as a theorem about a specified
real. The least-used sidecar here is the exact source and clock carried by
[native Collatz receipts](collatz_native_tail_recursion_20261013c.md).

The board is digit / denominator / source / clock / local growth / ROOT.
The anchor proves a universal statement about a particular Collatz code;
the niche retains source information despite a scalar quotient; the
wildcard tests how one changes the exponent by changing the coding rule.

Write an infinite regular continued fraction as

\[
\alpha=[0;a_1,a_2,\ldots],\qquad a_i\in\mathbb Z_{>0},\quad
q_{-1}=0,\ q_0=1,\ q_k=a_kq_{k-1}+q_{k-2}.
\]

**CITED classical interface:**

\[
\mu(\alpha)=2+\limsup_{k\to\infty}
\frac{\log a_{k+1}}{\log q_k}.                         \tag{1}
\]

Jonathan Sondow, *Irrationality Measures, Irrationality Bases, and a
Theorem of Jarník*, [Theorem 1, equation (2)](https://arxiv.org/pdf/math/0406300),
is the primary pin already retained in the
[CF reference sidecar](../reference/continued-fractions-khinchin-duffin-schaeffer-pins.md).
The exponent is at least two for every irrational. Formula (1) concerns
the infinite real; the checker never calls a finite ratio its exponent.

## 2. Local growth forces exponent two

**Theorem.** If a fixed real constant C≥1 satisfies
\(a_{k+1}\le C a_k\) at every step, then \(\mu(\alpha)=2\).
The same conclusion holds if the bound starts after a finite prefix.
No boundedness, convergence, average digit law, or periodicity is assumed.

Fix m≥1 and k≥m. Backward use of the local inequality gives
\(a_{k+1-j}\ge a_{k+1}/C^j\) for 1≤j≤m. Since
\(q_k\ge\prod_{i=1}^k a_i\),

\[
a_{k+1}^{m}\le C^{m(m+1)/2}q_k.                         \tag{2}
\]

The denominators grow at least as Fibonacci numbers. Divide the logarithm
of (2) by m log q_k and let k tend to infinity. The limsup in (1) is at
most 1/m for **every fixed m**; it is therefore zero. The proof also
applies to every sufficiently late m-window after a finite prefix.

There is an integer-only sufficient gate. For d≥1, k≥2d,

\[
q_k\ge C^{d(2d+1)}\quad\Longrightarrow\quad a_{k+1}^{d}\le q_k.
                                                               \tag{3}
\]

Indeed, (2) with m=2d bounds \(a_{k+1}^{2d}\) by \(q_k^2\).
The implementation uses integer C and verifies the local word before
applying the gate. A failed gate is inconclusive. The Fibonacci lower
bound gives an effective index beyond which (3) always holds. This is
an actual Diophantine control from a local growth proof, not a stopping
test for the underlying dynamics.

**The local hypothesis matters.** Set \(a_k=2^k\) when k is a power of
two and \(a_k=1\) otherwise. All digits satisfy the global envelope
\(a_k\le2^k\). Just before a spike N=2^j,

\[
q_{N-1}\le\prod_{i<N}(a_i+1)
\le2^{(N-1)+(N-1)}.
\]

Consequently (1) gives \(\mu\ge5/2\). This example has arbitrarily
large one-step jumps. An exponential height bound is not the local
ratio hypothesis of the theorem. Conversely, the unbounded synthetic
trajectory \(x_j=2^j\) does satisfy the local hypothesis and gives μ=2.

## 3. Three actual clocks, one blind exponent

Declare the map before encoding its orbit:

| Clock | Map | Uniform state-growth bound |
|---|---|---:|
| Standard C | n/2 if even; 3n+1 if odd | C(n)≤4n |
| Shortcut T | n/2 if even; (3n+1)/2 if odd | T(n)≤2n |
| Odd U | (3n+1)/2^v2(3n+1), odd n>1; U(1)=1 | U(n)≤2n |

For any positive source in the clock's domain, put
\(x_j=S^j(n)\) and \(\alpha_n=[0;x_0,x_1,\ldots]\).
Every positive orbit exists for all forward times under these definitions,
so this is an irrational real even if the orbit is hypothetical and
divergent. Theorem 2 gives

\[
\boxed{\mu(\alpha_n)=2\quad\text{for every source in all three clocks}.}
                                                               \tag{4}
\]

This is independent of the Collatz conjecture. The direct code itself is
injective: its first partial quotient is n. It is also exactly compatible
with the Gauss shift, \(1/\alpha_n-n=\alpha_{S(n)}\). The infinite
digit sequence retains the whole orbit. Passing to μ destroys all of
that distinction. In particular, μ=2 cannot distinguish a rooted orbit,
a possible different cycle, or a hypothetical divergent orbit in this code.

For the ROOT convention, standard C gives the periodic digits 1,4,2;
shortcut T gives 1,2; absorbed odd U gives 1,1,1. Absorbing U(1) is an
infinite-code convention, **not** permission to pad a strict ROOT receipt.
The separate [parity/stopping code](collatz_parity_cf_boundary_20261014a.md)
changes the observable or the coding rule; it must retain that distinction.

## 4. Faithful clock marking still does not change μ

Define

\[
b_j=jx_{j-1}+1,\qquad j\ge1.                              \tag{5}
\]

The index is the CF position. From its digits recover
\(x_{j-1}=(b_j-1)/j\) exactly; validate divisibility, positivity and
each actual transition. Thus (5) retains source, position and the whole
observed history. Its digits are at least two, avoiding the final-one
ambiguity of finite simple continued fractions. If \(x_j\le Cx_{j-1}\),
then \(b_{j+1}\le2Cb_j\). Hence the marked infinite code still has μ=2.
Adding a lossless sidecar to the real does not make the scalar μ lossless.

The finite rational \([0;2,1]=[0;3]\) is the cheap normalization hostile:
a rational approximant alone does not recover a literal finite digit
word or its clock. The production API retains the tuple, rather than
recovering it from a rational under an unstated final-digit convention.

For a supplied odd valuation word w, retain its affine carrier
\(U_w(x)=(Px+B)/Q\), P=3^length, Q=2^cost. Its exact formal odd-word
cell is

\[
x\equiv(Q-B)P^{-1}\pmod{2Q}.                              \tag{6}
\]

At every proper prefix i, impose \(P_i x+B_i>Q_i\), including the
empty prefix x>1 when w is nonempty. Final oddness and integrality in
(6) force each prescribed valuation; the cuts prevent earlier ROOT.
These inequalities give one least source on the cell, then every
nonnegative progression lift. For example, w=(2) has formal cell
x=1 mod8, but the least strict source is9: x=1 would be ROOT padding.

`OddCode` contains the exact source, word, states and marked CF digits.
Its auditor replays every valuation and recomputes every digit. `discharge`
accepts **a supplied first-hit word at that exact endpoint** and verifies
the concatenation. It never infers ROOT from μ or from a neighbouring
source. The literal positive control is

\[
7\xrightarrow{1}11\xrightarrow{1}17\xrightarrow{2}13
\xrightarrow{3}5\xrightarrow{4}1,
\]

with marked digits (8,23,52,53,26,7). Both the empty ROOT word at1
and the rejection of any post-ROOT valuation are tested.

## 5. Correct engineering of large exponents

Let a1=A≥2 and \(a_{k+1}=a_k^r\), with integer r≥2. Then
\(\log a_j=r^{j-1}\log A\). The elementary bounds

\[
\prod_{j\le k}a_j\le q_k\le\prod_{j\le k}(a_j+1)
\le2^k\prod_{j\le k}a_j
\]

show that
\(\log q_k=((r^k-1)/(r-1))\log A+O(k)\). Therefore

\[
\mu([0;A,A^r,A^{r^2},\ldots])=r+1.                       \tag{7}
\]

Thus **squaring successive digits gives μ=3, not4**; cubing gives4.
The degenerate A=1 produces all-one digits and μ=2.

For \(a_k=2^{k!}\), write \(S_k=\sum_{j=1}^k j!\). The same bounds
give \(\log_2q_k\le S_k+k\le3k!\) for k≥2. Hence
\(\log a_{k+1}/\log q_k\ge(k+1)/3\to\infty\), so μ=∞.
This example is also explicitly recorded in Sondow's Example 4; no
novelty claim is made. Finite changes to either construction do not change
the exponent. A finite computation cannot identify which infinite tail
was chosen.

These constructed reals are not logarithms of rational numbers merely
because both have a prescribed irrationality exponent. Any proposed
identification with a logarithm, a resonance ratio, or another named real
requires an independent equality theorem. Likewise, a code deliberately
switched after detecting ROOT is a valid reduction only when its exact
switch rule is retained; it is not the direct state code of §3.

## 6. Reproducible finite universe and scope

Run [the standalone checker](../../04-computation/experiments/continued_fraction_collatz_codes_20261014a.py)
with `python -B` and `python -B -O`. Its declared universe is:

- 128 sources per clock, 65 digits each: standard/shortcut sources1..128,
  odd sources1,3,..255; raw and marked local ratios, exact window bills
  at depths8,16,32,64 and widths1,2,4,8, plus gates d=1,2,3;
- every valuation word of length0..3 over1..4: 85 words and eight native
  progression lifts each; source guards, strict boundaries, continuants,
  unimodularity and CF cylinder membership;
- three power-digit engineering controls, the factorial control, a sparse
  spike hostile, the divergent synthetic local-growth control, and the
  literal supplied ROOT transport from7;
- 17 malformed-source, bool/float, wrong-clock, native-valuation,
  source-binding and ROOT-boundary hostiles.

The [saved JSON](continued_fraction_collatz_codes_20261014a.json) records
**17,222 always-active exact checks**. The infinite conclusions are the
proofs above. No finite run is presented as an asymptotic fit, a convergence
test, or a new ROOT discovery.
