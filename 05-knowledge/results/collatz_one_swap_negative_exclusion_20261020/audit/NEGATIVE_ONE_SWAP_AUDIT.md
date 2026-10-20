# Independent audit of the negative one-swap classification

20 October 2026. **PASS for the stated mechanical one-adjacent-swap family.**
The mathematical proof and the complete residual certificate have been checked
independently. The exceptional negative cycle is retained. This audit makes no
claim about all negative cycles, arbitrary parity words, several edits, or
convergence.

## 1. Scope and independently checked conclusion

Let \(0<p<q\), let \(c\) be the zero-intercept lower mechanical binary word
of length \(q\) and weight \(p\), and let \(w\) be obtained by one cyclic
adjacent swap of unequal bits. Put

\[
D=2^q-3^p,\qquad g=\gcd(|D|,B(w)).
\]

For \(D<0\), the reviewed proof establishes

\[
g=|D|\quad\Longleftrightarrow\quad(p,q)=(2,3).
\]

At the exceptional shape the source is \(011\). Its two actual outputs are
\(101\) and \(110\), with fixed points \(-7\) and \(-5\), respectively. Both
describe the cycle \(-5,-7,-10\). The primitive sparse-norm argument and the
nonprimitive repetition argument are separate and have both been checked.

The reviewed files are `../NEGATIVE_ONE_SWAP_EXCLUSION.md`,
`../certify_negative_swap_exclusion.py`, and `../negative-receipt.json`.
The independent verifier does not import any repository mathematical code,
read the production receipt, or reuse its logarithm routine.

## 2. Signed norm and critical-clock audit

The inherited cut-safe sparse-norm theorem is applicable to primitive
\(q>p\ge2\). The negative assumption gives
\[
p<q<p\log_2 3<2p,\qquad k=1,\qquad r=q-p,\qquad r/p<3/5.
\]
The last inequality follows from \(3^5<2^8\), not a decimal slope estimate.

The legality conditions remain exactly
\[
\epsilon=+1:\ a<r,\qquad
\epsilon=-1:\ a\ge p-r.
\]
They eliminate the two large-norm endpoint types. The other two endpoints
have content one. For an interior legal plus type, \(a+1\le r\), and
root-of-unity orthogonality gives
\[
\frac1p\sum_j|Q_+(z_j)|^2
=1+4^{a/p}+4^{(a+1)/p}
<1+2\cdot4^{3/5}<6.
\]
The final inequality is certified by
\(4^3<(5/2)^5\). For a legal minus type, \(p-a\le r\), so
\[
\frac1p\sum_j|Q_-(z_j)|^2
=1+4^{1/p}+4^{(p-a)/p}
<1+4^{1/3}+4^{3/5}<11/2<6.
\]
The three monomial exponents are distinct modulo \(p\) for every interior
type. Thus there are no unaccounted cross terms. AM-GM and full norm
divisibility, including prime-power depths, prove \(g<6^{p/2}\).
For \(p=2\), the only negative clock is \(q=3\), and the two legal outputs
have \(g=1\). There is no negative clock with \(p=1<q\).

For a hypothetical integer cycle, \(R=|D|/3^p\) consequently satisfies
\[
R<(2/3)^{p/2}<e^{-p/6}.
\]
The exponential inequality uses \(\log(3/2)>1/3\). At \(p\ge4\),
\(R<4/9<1/2\); hence
\[
0<\delta=p\log3-q\log2=-\log(1-R)
<\frac{R}{1-R}<2e^{-p/6}.
\]
Also \(1-R>1/2\), proving \(q=\lfloor p\log_2 3\rfloor\).
At \(p=2,3\), direct negative-length lists are \(q=3,4\), respectively,
so the same critical-clock conclusion is valid. The factor two and
the small-\(R\) argument are essential signed changes, and the proof
retains them.

## 3. External theorem, cutoff, and complete CF reduction

I inspected the primary PDF:
[Bugeaud, Mignotte and Siksek, arXiv:math/0403046](https://arxiv.org/pdf/math/0403046),
Theorem 9.4, printed page 16. Its real-field bound is for
\(\log|\Lambda|\), so it applies to the negative multiplicative form
\(\Lambda=2^q3^{-p}-1=-R\). The choices \(n=2\), field degree one,
\(B=q\), \(A_1=7/10\), \(A_2=11/10\) meet the stated height requirements.
Replacing \(2^{9/2}\) by \(23\) bounds the prefactor by \(602494200<10^9\).
No sign restriction from the positive argument has been imported.

The resulting comparison yields
\[
p<6\cdot10^9(1+\log(2p)).
\]
More precisely, \(p/(1+\log(2p))\) is increasing for \(p\ge1\), while
\(\log(2\cdot10^{12})<30\) and \(6\cdot10^9\cdot31<10^{12}\).
Thus \(p<10^{12}\).

For \(p\ge100\), the signed logarithmic bound is strong enough for
Legendre's criterion:
\[
0<\log_2 3-q/p<1/(2p^2).
\]
The independent verifier certifies
\(4\cdot100e^{-100/6}<\log2\) using an exact positive Taylor partial sum.
The left-hand expression decreases thereafter. Coprimality is present
when Legendre is used, so ordinary convergents, not semiconvergents,
exhaust the candidates.

The independent calculation finds all 18 convergents with
\(100\le p<10^{12}\), nine on each side. For **every** one, including
the sign-incompatible side, it checks the stronger absolute-gap condition
\[
|q\log2-p\log3|
>2\cdot2^{-\min(\lfloor p/6\rfloor,256)}
\ge2e^{-p/6}.
\]
The terminal convergent is
\[
(q,p)=(8573543875303,\ 5409303924479),
\]
whose denominator exceeds the cutoff. Continued-fraction denominators
increase thereafter. The production certificate uses the sign alone
to reject upper convergents and checks the lower side's gap; both
implementations therefore give the same complete exclusion.

## 4. Independent exact certificates and hostile controls

`verify_negative_independent.py` uses the ordinary positive expansion
\[
-\log(1-t)=\sum_{j\ge1}t^j/j
\]
with 300 terms and the upper tail
\(t^{301}/(301(1-t))\). It obtains \(\log2\) from \(t=1/2\) and
\(\log3=\log2+\log(3/2)\) using \(t=1/3\). The production checker instead
uses an atanh series. Every continued-fraction digit is independently
certified by matching the two exact endpoint floors. No floating-point
decision appears.

The finite residual is independently reconstructed with
\(q=\operatorname{bit\_length}(3^p)-1\), checked against the logarithm
interval, for every \(2\le p<100\). It comprises 59 primitive shapes and
3,348 actual cyclic swaps. Each carry is computed right-to-left and
cross-checked against affine composition. The only integral outputs are:

| \(p\) | \(q\) | word | \(B\) | \(D\) | verified shortcut orbit |
|---|---|---|---|---|---|
| 2 | 3 | 101 | 7 | -1 | -7, -10, -5 |
| 2 | 3 | 110 | 5 | -1 | -5, -7, -10 |

The independent hostile census is broader than the proof's residual:
it checks every primitive negative mechanical shape through \(q=160\),
including clocks not immediately below the critical slope.

- 2,881 primitive shapes.
- 113,674 actual swap occurrences, with 2,881 crossing-cut occurrences.
- The \(q\le120\) subtotal is 1,618 shapes and 47,892 swaps, reproducing
  the production subtotal.
- Only the same two integral outputs occur.
- All legality equivalences and every numerical \(g^2<6^p\) are checked.
- A fraction-free determinant computation independently evaluates 158
  distinct sparse norm moduli. There are 426 full-divisibility checks
  across actual negative-clock words with \(p\le18\).

This census corroborates, rather than replaces, the unbounded norm proof.

The genuine \(-17\) cycle is reconstructed directly from the shortcut map:
\[
-17,-25,-37,-55,-82,-41,-61,-91,-136,-68,-34.
\]
It has parity \(11110111000\), shape \((p,q)=(7,11)\), clock \(-139\),
and carry \(2363\). Its minimum Hamming distance over all rotations
from the mechanical source \(01011011011\) is four. A single adjacent
unequal swap changes exactly two positions, so this cycle is outside
the theorem's family. It is not an omitted counterexample to a claimed
general-negative-cycle exclusion; no such exclusion is made.

## 5. Unit-clock and nonprimitive boundary audits

The elementary unit-clock proof is correct. If \(3^p-2^q=1\) under
\(0<p<q\), then \(p\ge2\), \(q\ge3\), and reduction modulo eight makes
\(p=2m\) even. The two positive even factors
\[
(3^m-1)(3^m+1)=2^q
\]
are powers of two, differ by two, and have gcd two. They must be two
and four, forcing \((p,q)=(2,3)\).

For a nonprimitive source \(v^d\), \(d\ge2\), the factorization
\(D=HD_0\) has \(H>1\) and \(\gcd(H,6)=1\). A cut away from the
edited edge turns the perturbation into one internal adjacent swap.
Its carry difference is a signed \(2,3\)-unit, while the original
carry is divisible by \(H\). Consequently \(\gcd(H,B(w))=1\), and
\[
\gcd(|D|,B(w))=\gcd(|D_0|,B(w))<|D|.
\]
All relevant moduli divide \(D\), so the rotation identity legitimately
transports their gcds back to the original cut. No assumption that
\(H\) and \(D_0\) are coprime is needed: coprimality of \(B(w)\) with
\(H\) alone proves the gcd equality.

The independent verifier checks the two exact repetition identities
for every nonconstant binary base of length two through eight, every
repetition count two through five, and every unequal cyclic swap:
25,088 occurrences, including 1,016 cut crossings and 3,920 negative
clock occurrences. These finite controls support the sign-independent
proof rather than extend the primitive norm theorem illicitly.

## 6. Reproduction and disposition

Normal and optimized independent runs produce byte-identical receipts
with **1,044,780 explicit checks**. Checks are ordinary conditional
exceptions, not Python assertions that disappear under optimization.

From this directory:

```sh
python3 verify_negative_independent.py > independent-receipt.replay.json
python3 -O verify_negative_independent.py > independent-receipt.optimized.replay.json
cmp independent-receipt.json independent-receipt.replay.json
cmp independent-receipt.json independent-receipt.optimized.replay.json
python3 compare_negative_receipts.py > production-comparison.replay.json
```

`compare_negative_receipts.py` verifies the complete continued-fraction
prefix and pairs, terminal convergent, every small-shape count and
maximum content, the exception set, the \(q\le120\) count, and both
normal/optimized equalities. Its receipt records SHA-256 values for
the audited artifacts.

No mathematical correction was required after reading the full signed
proof. The initial draft's dependency links used publication-layout
paths rather than staging-layout paths; this is a packaging issue, not
a proof gap, and was reported to the author. No files were published
as part of this audit.

## 7. Addendum: the exact distance to the -17 necklace is three

The Hamming-distance control above does **not** imply that the genuine
\(-17\) cycle is reachable by two adjacent swaps. The exact minimum
number of cyclic adjacent unequal-bit swaps between the mechanical
necklace of shape \((7,11)\) and the \(-17\) parity necklace is **three**.

`verify_minus17_swap_distance.py` independently reconstructs the cycle,
then uses bit masks to enumerate the entire cyclic adjacent-swap graph
of weight-seven binary words of length eleven. It reaches all
\(\binom{11}{7}=330\) vertices and checks 1,848 directed edge occurrences.
The distance-layer sizes from \(01011011011\) are

| distance | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
|---|---|---|---|---|---|---|---|---|---|
| vertices | 1 | 8 | 29 | 60 | 82 | 90 | 43 | 14 | 3 |

For left rotations \(r=0,\ldots,10\) of \(11110111000\), the distances
are, in that order,

\[
6,\ 5,\ 5,\ 6,\ 4,\ 6,\ 5,\ 3,\ 5,\ 4,\ 4.
\]

The calculation also has an independent elementary lower-bound
certificate. For any fixed target \(v\), set
\[
b_i=\sum_{j=0}^{i}(c_j-v_j),\qquad 0\le i<11.
\]
For an arbitrary legal swap path, let \(f_i\) be the signed net number
of ones transported across the directed cyclic edge \(i\to i+1\).
Conservation gives \(f_i-f_{i-1}=c_i-v_i\), hence
\(f_i=b_i+t\) for some integer circulation \(t\). The number of swaps
is at least
\[
\sum_i |f_i|\ge\min_{t\in\mathbb Z}\sum_i|b_i+t|.
\]
A median of the integer \(b_i\) minimizes the last expression. The
independent script computes this bound for every target rotation,
obtaining exactly the eleven distances above. In particular, every
rotation requires at least three swaps.

A three-swap witness, independently reproduced by the shortest-path
search, is

\[
01011011011
\;\xrightarrow{0}\;
10011011011
\;\xrightarrow{4}\;
10010111011
\;\xrightarrow{3}\;
10001111011.
\]

Labels are zero-based swap-start positions. Each step swaps unequal
adjacent bits. The final word is the left rotation by seven of the
\(-17\) parity word; its carry is \(12649\), its clock is \(-139\), and
its integer fixed point is \(-91\). The script verifies all indicated
parities and return to \(-91\).

Rotation is a graph automorphism, so fixing \(c\) and minimizing over
target rotations is equivalent to minimizing over both necklaces.
The lower bound and explicit witness prove exact necklace distance
three, independently of interpreting the BFS as a certificate.
This supplies a genuine three-edit integer-cycle boundary; it does
not establish any general theorem about two-edit words.

Normal and optimized outputs agree with 1,894 explicit checks:

```sh
python3 verify_minus17_swap_distance.py > minus17-distance-receipt.replay.json
python3 -O verify_minus17_swap_distance.py > minus17-distance-receipt.optimized.replay.json
cmp minus17-distance-receipt.json minus17-distance-receipt.replay.json
cmp minus17-distance-receipt.json minus17-distance-receipt.optimized.replay.json
```
