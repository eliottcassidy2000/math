# Retaining excess credit across weak Collatz resets

**Status: PROVED** for the guarded families and supplied finite seed banks below.
**FINITE-EXACT** for the explicitly counted controls. Universal entry remains
**OPEN**. A policy rejection is not a nonconvergence result. This extends the
proof controller on an inherited constructive basin; it does not establish a
new basin beyond the known inverse closure of checked ROOT certificates.

The new mechanism is cumulative funding. A weak first descent may spend the
surplus from an earlier strong descent. The source-only budget uses the rate
chosen for the whole controller, rather than the maximum rate required by an
individual segment. Both the incoming source and the checked endpoint bank
remain part of the proof object.

## 1. Inheritance and the correction boundary

The closest mechanism is the per-segment guard in
[the weak-reset family](collatz_weak_reset_family_20261005.md), followed by the
[source-dependent deadline-to-measurement compiler](collatz_backward_measurement_compiler_20261005.md).
The canonical hostile is the actual first descent from 47 to 23: its length is
34, total valuation is 55, and its least individual rate is 31. The corrected
near miss is charging that rate to the entire source height even after earlier
segments have supplied usable credit. The least-used sidecar is the minimum
of the cumulative balance, which composes exactly and supports finite receipt
reuse. The live objects are **first-descent boundaries, cumulative valuation,
prefix deficit, immutable source, finite seed bank, and the measurement floor**.

The incoming proof's scalar inequality `exp(t/30) <= 2^(0.048t)` is false for
positive `t`. Instead, `ln 2 > 69/100` gives the contraction constant
`197/207`, whose reciprocal is less than `1.051`. The companion weak-reset
note and script retain that correction and make the finite seed premise
explicit. Historical floating-point census outputs are not all-input exact
deadline algorithms. No assertion here depends on the incoming unrestricted
staircase/asymptotic claims in
[expense-Diophantine](collatz_expense_diophantine_20261005.md).

The connection contract is precise: a checked actual orbit prefix maps to a
rational balance summary. It preserves the required initial credit and the
composition law. It loses the source cylinder, carries, order within a
segment, and ROOT evidence. Those remain in the actual word and seed receipts.
Words `(1,3)` and `(3,1)` have the same balance but different legal sources.

## 2. The cumulative guard and its deadline

Let \(U(n)=\operatorname{oddpart}(3n+1)\), with no outgoing ROOT edge at 1.
Fix an integer rate \(q\ge1\), an integer initial credit \(B\ge0\), and a
threshold \(Y\ge10q\). Supply a nonempty finite bank \(C\) of positive odd
seeds, each with its checked strict first-hit ROOT word. This is finite input
data, not an assumption that every number below an arbitrary threshold is
known to converge.

Start at \(n\). Stop at the first banked **running minimum**. Before that,
every segment start must be at least \(Y\). A first-descent segment runs
until its first endpoint below that start. If its length and valuation sum
are \(\ell_i,A_i\), write

\[
 L_j=\sum_{i\le j}\ell_i,\qquad A_j^*=\sum_{i\le j}A_i,\qquad
 M_j=\frac{2^{qA_j^*-L_j}}{3^{qL_j}}.
\]

At **every completed segment boundary**, require

\[
                         2^B M_j\ge1.                 \tag{1}
\]

Membership requires eventually reaching a seed in \(C\), with all these
guards. A new small unbanked running minimum is an uncovered stopping case.
The bank may contain seeds above \(Y\); a supplied seed itself is accepted
immediately, and all other hits are at segment boundaries. The experiments
use banks below \(Y\).

**Theorem.** If the bank endpoint is \(z\), and the prefix to it has length
\(L\), then

\[
 L\le \frac{207}{197}\left(q\log_2\frac nz+B\right).
                                                        \tag{2}
\]

**Proof.** Every preterminal state inside a first-descent segment is at least
that segment's start. Thus all the preterminal states of the entire prefix
are at least \(Y\), even though its final seed may be below \(Y\). The exact
product identity is

\[
 z=n\frac{3^L}{2^A}\prod_{i=0}^{L-1}
             \left(1+\frac1{3n_i}\right).
\]

Use \(\log(1+t)\le t\), and
\(\ln2>2(1/3+1/81)=56/81>69/100\), to obtain

\[
 \log_2(n/z)\ge A-L\log_2 3-\frac{L}{3Y\ln2}
 \ge\frac{L-B}{q}-\frac{10L}{207q}
 =\frac{197L}{207q}-\frac Bq.
\]

Here the second inequality uses the **final cumulative** guard. The guards
at earlier boundaries ensure the stronger funded-prefix policy; they are not
substituted for source legality. This proves (2).

Let \(c_{\min}=\min C\), and let \(b_C(n)\) be the least nonnegative integer
with \(n\le c_{\min}2^{b_C(n)}\). Every member outside the bank satisfies
\(n>z\ge c_{\min}\), so an exact source-only cap is

\[
 D(n,q,B,C)=\left\lfloor
       \frac{207(qb_C(n)+B)}{197}\right\rfloor.           \tag{3}
\]

The ratio ceiling is computed by integer bit lengths and comparison, without
a floating-point logarithm. With a bank containing 1, \(b_C(n)\) is simply
the bit length of odd \(n>1\). After a bank hit, append only that seed's
checked suffix. The resulting ROOT deadline is \(D+\tau(z)\); a precomputed
uniform version uses \(D+\max_{c\in C}\tau(c)\).

This explains the exact boundary in subtracting an endpoint height. One may
subtract \(\log_2 c_{\min}\), because the actual endpoint is in the specified
bank. One may not generally subtract \(\log_2 Y\), since the last segment
can jump below the threshold.

## 3. A finite recognizer and the reusable receipt summary

`recognize(n,q,bank,credit,threshold)` checks the finite seed words, computes
(3), and performs at most that many new odd steps. It starts a new segment
only after observing an actual first descent. At that boundary it checks
(1), tests the bank, and checks the next segment-start threshold. The only
successful leaves are the independently validated bank entries. ROOT padding
and forged valuations are rejected.

If the cap expires, the source is outside **this declared family**: any
member would already have reached its checked seed by (2). If a credit guard
fails or an unbanked small start is reached, that also rejects only this
policy. All three outcomes retain the checked actual prefix. In particular,
27 and 55 are policy misses for the small bank below despite their separately
checked convergence. No cyclic assumption is used to certify an endpoint.

For a list of segments retain

\[
 (M,m),\qquad m=\min(1,M_1,\ldots,M_j).
\]

Concatenation is the associative operation

\[
 (M,m)\star(N,u)=(MN,\min(m,Mu)).                       \tag{4}
\]

The least integral initial credit is the least \(B\ge0\) with \(2^B m\ge1\).
Thus a receipt can be cached and composed without discarding the earlier
funding. Resetting the balance to one at each boundary returns the weaker
individual-segment policy. The summary alone is not a certificate: exact
source words and the bank interface must accompany it.

At credit zero, the cumulative family contains the individual-segment family
for the same rate, threshold, and bank. Increasing the declared credit makes
the family larger and increases its explicit source deadline. Credit can be
any independently specified computable function of the source; computing the
credit from an already completed target route must not be presented as a
source-independent convergence proof. The bounded recognizer itself never
reads such a target route.

## 4. Inverse funding and a strict infinite extension

**Inverse-funding lemma.** Let \(x>1\) be an odd 3-unit. Supply a finite
checked excursion from \(x\) to a checked bank seed, ending at a first-descent
boundary. Require every preseed segment start to be at least \(Y\). Let
\((M,m)\) be its summary at the chosen rate. Choose the least exponent
\(a\ge2\) in the parity class satisfying \(2^a x\equiv1\pmod3\) such that

\[
                  \frac{2^{qa-1}}{3^q}m\ge1.            \tag{5}
\]

Then every exponent \(a+2t\), \(t\ge0\), gives

\[
                         n_t=\frac{2^{a+2t}x-1}{3}
\]

in the cumulative family with zero initial credit. The new first step has
exact valuation \(a+2t\), strictly descends to \(x\), and supplies the
balance in (5). Increasing \(t\) multiplies that balance by \(2^{2qt}\).
This proves both finite selection of the exponent and the all-height family.
An anchor divisible by 3 has no such odd inverse edge and is rejected.

For the concrete anchor 47, the checked first-descent word is

```text
(1,1,1,2,2,1,2,1,1,2,1,1,1,2,3,1,1,
 2,1,2,1,1,1,1,1,3,1,1,1,4,2,2,4,3)
```

It has length 34 and cost 55, ending at 23. The seed 23 has the independently
checked word \((1,1,5,4)\). At rate 3, (5) selects \(a=13\), hence

\[
             n_t=\frac{47\,2^{13+2t}-1}{3},\qquad t\ge0. \tag{6}
\]

The final cumulative guard is
\(2^{3a+130}\ge3^{105}\). It holds at \(a=13\) and fails at the preceding
legal phase \(a=11\). The latter needs exactly four initial binary credit
units. Every member of (6) has the strict ROOT word
\((a)\,w_{47}\,(1,1,5,4)\), of length **39**. Every one is rejected by the
individual rate-3 policy because its 47 segment requires rate 31. This proves
a strict infinite extension with the same bank and threshold.

There is also a **closed-form source-adaptive budget** for the entire
decreasing inverse ray, including its earlier phases. Test directly that
\(3n+1=47\,2^a\) for an odd \(a\ge3\); exact divisibility and a power-of-two
test recover \(a\) from the supplied source. Then

\[
                         B(n)=\max(0,37-3a)             \tag{7}
\]

is exactly the least rate-3 initial credit. Indeed, the new first boundary
already has balance \(2^{3a-1}/27>1\), and the only remaining boundary before
seed23 has balance \(2^{3a+130}/3^{105}\). The integer comparison
\(2^{166}<3^{105}<2^{167}\) proves (7), including its rounding boundary.
Every guarded source with this budget is accepted. No target iteration is
needed to compute the budget; the fixed checked excursion and seed words
are shared proof data and must still be paid for and retained. With the
singleton bank \(\{23\}\), each of the phases \(a=3,5,7,9,11\) has deadline
**42**, while its actual structural rank is 39. One less unit of credit fails
the declared policy. The phase \(a=1\) gives source31, which rises to47 and
is correctly outside this decreasing-inverse construction. Nearby integers
are not accepted merely because the budget formula would be positive.

Even fixed rate **1** can fund that excursion: the least exponent is 37,
giving another infinite phase \(a=37+2t\). This does not assert that fixed
rate 1 covers every supplied source; it constructs new sources whose strong
first edge pays the retained excursion.

At the least source 128341, the generic cumulative deadline is 57 using bank
\(\{1,23\}\), or **44** using the singleton bank \(\{23\}\). The inherited
rate-31 inequality with the actual stopping seed 23 gives the conservative
bit-length bound 1058. The structural family proof is sharper still, with
exact rank 39. These are different amounts of retained structure, not three
claims of optimality. The generic deadline57 compiles to

\[
 W(128341)\ge\frac1{142707391328010742671022200}.
\]

Its selected leaf is 2737941, and localization degree 148 gives a positive
actual measurement lower bound by the inherited universal kernel error.
The code also validates the tighter singleton-bank floor. No moment of the
target is numerically assumed or queried.

## 5. Finite controls and honest work accounting

The controls use exactly the 1009 odd requested sources 31 through 2047, rate
3, threshold 30, and the same explicitly checked bank of the 15 odds 1 through
29. Seed routes are obtained and validated within a declared audit cap128;
their maximum odd rank is 41. They are the only ROOT words supplied to the
recognizer. Each successful result is independently replayed afterward.

| policy | accepted requested sources | newly queried odd edges |
|---|---:|---:|
| individual segment, credit0 | 306 | 10029 |
| cumulative, credit0 | 367 | 11602 |
| cumulative, credit4 | 492 | 15117 |
| cumulative, credit16 | 664 | 20371 |

Seed validation and successful-route verification are separate work; the
counts are not a runtime speedup claim. The least additional source at zero
credit is 203:

\[
 203\xrightarrow{(1,2,4)}43
 \xrightarrow{(1,2,2)}37
 \xrightarrow{(4)}7.
\]

The strong first segment pays the individually weak middle segment; seed7's
checked suffix discharges the result. All policy inclusions are checked in
this common finite universe. The positive infinite theorem is (6), not an
extrapolation of these counts.

Additional controls include the all-height parameterization at 65 exponents
13 through 141, sixteen exact source-budget phases 3 through 33, nine rate-1
parameters, 64 summary associativity triples,
the exact preceding-phase credit boundary, and malformed sources, rates,
seed receipts, ROOT self-return, an unfinished excursion, a 3-divisible
anchor, and a forged summary. Explicit `check` calls remain active under
`python -O`.

## 6. Reproduction and remaining obligation

[Script](../../04-computation/experiments/collatz_amortized_excess_budget_20261005.py)
and [saved output](collatz_amortized_excess_budget_20261005.out):

```text
python 04-computation/experiments/collatz_amortized_excess_budget_20261005.py
python -O 04-computation/experiments/collatz_amortized_excess_budget_20261005.py
```

The controller now has a finite, source-aware way to retain credit instead of
discarding it when patterns change. What remains is an entry/coverage
obligation: show that a specified source reaches the checked bank while
meeting the declared cumulative budget and threshold guards, or provide
another guarded construction for it. Enlarging the budget is a mathematically
sound family extension, but no finite chosen budget is asserted to settle
all sources. The next concrete target is to combine this balance summary
with the source-only run-block parser, preserving both the actual source
guards and the original cumulative credit at each interface.
