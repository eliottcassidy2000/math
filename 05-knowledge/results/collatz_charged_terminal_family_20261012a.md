# Charged terminal families, a paid 75 exit, and a bounded band invoice

2026-10-12. **PROVED:** the algebraic compiler, first-hit guards, fixed-family
measurement floor, and source-owned finite invoice below. **FINITE-EXACT:**
the declared dispatcher comparison and band-exit census. **OPEN:** the
all-source band hypothesis and Collatz. No theorem ID is reserved here.

## 1. Inheritance and the exact change

The closest mechanism is [the terminal phase in C3](collatz_shadow_charge_20261011.md),
already implemented in the [grounded-port compiler](collatz_grounded_child_ports_20261007f.md).
The [previous terminal family](collatz_terminal_completion_20261011.md) also
supplies a positive normalized score and a measurement floor preserved under
refinement. Completing a chosen head is therefore inherited machinery, not a
new inverse basin theorem. The present addition extracts its exact score for
**any fixed head**, gives a source decoder, and applies it to the next pending
source 75 of the [charged dispatcher](collatz_charge_receipts_20261011.md).

The canonical hostile is 327: its actual odd charge is 35, while
\(\lfloor\log_2(328)\rfloor+24=32\). Thus the incoming linear charge invoice
cannot be asserted universally. The corrected near miss is to replace a fixed
source by a completed nearby source; the decoder below refuses that change.
The retained sidecars are the ordered carry, exact terminal phase, original
source, odd charge, and the supplied finite seed bank.

The live concepts are source identity, terminal phase, charge, affine carry,
refinement error, and the distinct meanings of the quadratic constants and 63.
The positive lane is the 75 family. The orthogonal lane is a total test of a
quadratic logarithmic invoice. The dynamics lane supplies a precise quotient
loss, rather than an asserted conjugacy.

## 2. A fixed-head compiler with an exact normalized score

Use \(U(n)=(3n+1)/2^{v_2(3n+1)}\) on positive odd integers, stopping at the
first 1. For a word \(w\) of positive valuation letters, write

\[
U_w(n)=\frac{Pn+B}{Q},\qquad P=3^L,\quad Q=2^A,
\qquad o=\#\{a\in w:a\text{ odd}\}.
\]

Append an even terminal valuation \(2k\). Its source must be

\[
n_k=\frac{Q4^k-Q-3B}{3P},\qquad
4^k\equiv1+3BQ^{-1}\pmod{3P}. \tag{1}
\]

The right side is a principal unit modulo \(3^{L+1}\). Since 4 generates
that group, (1) selects one exponent class modulo \(P\). The implementation
lifts one ternary digit at a time, testing the three possible next digits.
This uses no target ROOT trajectory and no heuristic discrete logarithm.

For every native \(k\ge2\), (1) supplies a positive odd source and the
strict ROOT word \(w,(2k)\). Here is the integrality and first-hit detail.
The endpoint of the head is \(h=(4^k-1)/3\ge5\). Each intermediate written
forward from integer \(n_k\) has a power-of-two denominator; written backward
from integer \(h\), it has a power-of-three denominator. Thus it is an
integer. The inverse recurrence \((2^a y-1)/3\) then makes every intermediate
odd and positive. The formula for \(n_k\) is itself odd: for nonempty heads
\(Q\) is even and \(B\) odd, while the empty-head formula is
\((4^k-1)/3\). An earlier state 1 would remain 1 under actual \(U\), which
cannot reach \(h\ge5\). There is no ROOT padding.

Let

\[
C=Q+3B-P.
\]

The empty word has \(C=0\). Appending a letter \(a\) gives
\(C_{wa}=3C_w+2^aQ_w\), so every nonempty head has \(C>0\).
Using the inherited quarter-flow score
\(f(n)=4^{-K}(3/4)^{\tau-1}\), with
\(K=\sum\lfloor(a_i-1)/2\rfloor\), gives exactly

\[
f(n_k)=\frac{4P}{Q2^o4^k},\qquad
Z(n_k)=\frac{3n_k+1}{4}f(n_k)
=2^{-o}\left(1-\frac{C}{Q4^k}\right). \tag{2}
\]

For a nonempty fixed head this increases strictly to \(2^{-o}\). For the
empty head it is constantly 1. The exported plan starts at the least native
\(k\ge2\) with \(Q4^k\ge2C\). Hence every exported source has the explicit
uniform normalized floor

\[
Z(n_k)\ge\eta_w:=2^{-o-1}>0. \tag{3}
\]

The cut in (3) pays a convenient score floor; it is not a necessary condition
for the strict ROOT certificate in (1).

`decode_source(w,n)` recovers \(4^k=(3Pn+Q+3B)/Q\), checks that it is an
even power of 2, and checks phase and cut. It returns the unique plan index
or `None`. The source is never silently replaced by a different family member.
For modular reading, reduce the numerator in (1) modulo \(3Pm\) **before**
dividing by \(3P\); this returns \(n_k\bmod m\) exactly. Literal expansion
has an explicit bit cap. The modular reader also accepts a tested index
\(10^{60}\) without expanding its source.

## 3. A new dispatcher exit and its completed subfamily

For every \(n=75+128t\), \(t\ge0\), the first two valuations are \((1,2)\).
Put

\[
a=v_2((27n+23)/8)\ge4,\qquad h=\frac{27n+23}{2^{a+3}}.
\]

Both proper prefix endpoints exceed \(n\), while

\[
4(h+1)\le n+1. \tag{4}
\]

Indeed the left side is at most \((27n+23)/32+4\), and (4) follows from
\(20n\ge476\). Thus this is an exact first descent. Its odd charge
\(\delta=1+(a\bmod2)\le2\) is paid by
\(E(n)-E(h)\), where \(E(x)=\lfloor\log_2(x+1)\rfloor\).
No claim that an arbitrary such \(h\) is rooted is made.

For \(w=(1,2)\), (1) has \((P,Q,B)=(9,8,5)\) and phase \(k=4\pmod9\).
It supplies the independently rooted subfamily

\[
n_j=\frac{2^{18j+11}-23}{27},\qquad
U_{(1,2,18j+8)}(n_j)=1,\qquad j\ge0. \tag{5}
\]

Here \(n_0=75\), every \(n_j\equiv75\pmod{128}\), and each word has odd
charge exactly 1. These sources are not all in the narrower exact
\((1,2,8)\) cylinder. Their final exponent is variable, as (5) records.
Their sharp common score floor is

\[
Z(n_j)=\frac12-\frac7{8\,2^{18j+8}}
\ge\frac{1017}{2048}. \tag{6}
\]

This proves an infinite domain from a three-letter formula; the proof does
not read a bank of individual ROOT words for the generated sources.

A separately retained next rule, derived by the parent lane, covers
\(n=123\pmod{512}\). The first four letters are \((1,2,1,2)\); all four
prefix slopes exceed 1. With
\(a=v_2((243n+319)/64)\ge3\), its first descent is
\(h=(243n+319)/2^{a+6}\). Since \(n\ge123\), the inequality
\(26n\ge1150\) proves \(2(h+1)\le n+1\). Its charge
\(2+(a\bmod2)\le3\) pays the inherited rate-three source invoice.
This second rule is not included silently in the 75-only comparison.

## 4. One measurement precision works for an entire fixed-head family

Inherit the exact readout and alias theorem [Q2](collatz_quarter_price_source_20261011.md).
For a fixed \(n=2m+1>1\), put

\[
\nu(m)=(m+1)^{-6}-(m+2)^{-6},\quad q_0=(3n+1)/2,
\quad R_d=(4^d-1)/3,\quad q_d=q_0R_d,
\]
\[
A_d=\nu(m)+4^{-d}\nu(m+q_d),\qquad
E_d=\frac9{448q_d^8}.
\]

The actual readout obeys \(A_df(n)\le M_d\le A_df(n)+E_d\).
The family formula (2) supplies an independently proved lower bound
\(L_d=A_df(n)\), so the conservative recovered floor is
\((L_d-E_d)/A_d\). It does not assume that an off-source measurement bank
already knows the target atom.

If \(Z(n)\ge\eta>0\), then
\(\nu(m)>6/q_0^7\) and \(f(n)\ge2\eta/q_0\). Therefore

\[
\frac{E_d}{A_df(n)}
<\frac9{5376\eta R_d^8}. \tag{7}
\]

Every depth with \(2688\eta R_d^8\ge9\) retains more than half the actual
score. Crucially this depth depends on the **fixed head**, not on its
terminal exponent or source height. The API finds the least such positive
depth by exact integer/rational comparisons. For (5), (7) is at most
\(8/2373\) already at depth 1, giving

\[
\frac{L_d-E_d}{A_d}>\frac{2365}{2373}f(n_j)>0
\quad(j\ge0,\ d\ge1). \tag{8}
\]

In particular this survives the nested depths \(d=2^s\). We do not say
that the moduli at every pair of consecutive depths are nested.
The inherited transfer \(W(n)\ge(9/16)f(n)^2\) also gives positive beta-mixture
floors on the same constructed domain, not on arbitrary requested sources.

## 5. Finite comparison and a new falsifiable budget target

The dispatcher comparison fixes all 15 supplied strict ROOT words at odd
sources 1 through 29 from `charge_ports_eleven_20261011.json`. It keeps the
incoming one-step, sibling, uncapped reset, two sporadic, 39 and 43 rules;
it enumerates every contracting inverse template of length at most 8, with
all positive letters allowed by the contraction inequality. There are 1,045
templates. It considers all available strict smaller dependencies, uses an
ascending pass, and chooses a deterministic minimum-charge grounded one.
Both runs have exactly the same seed bank and rule choices except (4).

On odd sources 1 through 8,191 the grounded count is **1,018 versus 1,193**,
a gain of 175; the first new source is 75 and the first remaining pending
source is 123. The separately reported third policy adds the 123 rule:
**1,193 to 1,278**, a further gain of 85, with first pending source **127**.
All charge labels are independently compared with actual
strict ROOT replay in a separately marked control-only routine, cap 1,000.
That routine is not called by the production grammar. This is a relative
gain over this bounded dispatcher, not a claim to exceed every earlier
finite ROOT census or the full unbounded inverse grammar.

The least counterexample through 327 to the proposed universal
\(o(n)\le E(n)+24\) is exactly 327, with \(o=35,E=8\).
A more flexible **OPEN hypothesis** is:

> For every odd \(n>1\), the first actual prefix reaching
> \(2(h+1)\le n+1\) exists and has odd charge
> \(\delta\le3(E(n)^2-E(h)^2)\).

This is stronger than convergence. If true, strict descent and induction
give \(o(n)\le3(E(n)^2-1)\). There is no claim that finite tests prove it.
The control universe is every odd source 3 through 1,000,001, with a
10,000-step per-source discovery cap. All reached the specified band;
the largest used length was 135. The largest exact ratio
\(\delta/(E(n)^2-E(h)^2)\) was \(26/9\), at
\(55\to23\) in 37 odd steps, with charge 26 and \(E=5\to4\).
Coefficient 1 fails 136 sources in the initial subuniverse through 100,001,
first at 27. Thus even this coarse replacement needed a genuine hostile test.

The single-source hypothesis has a **total finite tester** without assuming
convergence. Set \(B=c(E(n)^2-1)\). An even letter contracts \(n-1\) by at
least \(3/4\), and an odd letter expands it by at most 2. Let \(t\) be the
least integer with \((n-1)2^B3^t<2\,4^t\). By time \(B+t\), either ROOT
has occurred or the \((B+1)\)-st odd letter has occurred. Stop earlier at the
first band exit and check its exact invoice. If the charge budget is exceeded
before that exit, no later exit can pay it, since its allowable invoice is
at most \(B\). `band_invoice` returns this authenticated smaller dependency
or a failed invoice; failure is **not** a nonconvergence certificate. A
successful result still needs its child's ROOT proof or a general induction.

## 6. The constants name different dynamical objects

The current proved route is
[THM-4139, rational-three-cycle-order-six-lift-and-horizontal-carrier](../../01-canon/theorems/THM-4139-rational-three-cycle-order-six-lift-and-horizontal-carrier.md),
with its corrected fibre firewall
[THM-4146](../../01-canon/theorems/THM-4146-rational-three-cycle-order-six-lift-horizontal-divisor-fibre-firewall.md).
We use the elementary rational-dynamics statements, not a rejected
Jacobian-carrier identification.

* For \(x^2-29/16\), the rational cycle is
  \(-7/4\to5/4\to-1/4\to-7/4\), with derivative product \(35/8\).
  Here \(-7/4\) is a **point**.
* For \(x^2-7/4\), \(-7/4\) is the **parameter**. The critical orbit of
  zero has denominator sequence \(4,16,256,65536,\ldots\), so it is not
  that rational three-cycle. The inherited Chebyshev factor relation is
  restricted to a cubic cycle; it is not a global Collatz conjugacy.
* The elementary orders \(\operatorname{ord}_9(2)=
  \operatorname{ord}_{21}(2)=\operatorname{ord}_{63}(2)=6\) concern
  multiplication modulo odd moduli. The factorization \(63=2^6-1=3^2\cdot7\)
  is the familiar missing-new-prime example at this exponent. Neither
  statement identifies a Collatz stopping time. The integer 63 is separately
  authenticated by its literal ROOT word in the controls.

The prior [difference-diagonal packet](difference_diagonals_golden_phase_recursion_20261004.md)
already retains the ordered carry lost by a trace/norm projection.
The control words \((1,1,1,3)\) and \((1,1,2,2)\) have the same
\((P,Q)=(81,64)\), but carries 65 and 73 and different rational fixed
points \(-65/17,-73/17\). Our terminal phase and measurement floor depend
on this carry. A shared number, order, multiplier, or trace alone supplies
none of the source binding required by (1).

## 7. Reproduction and interfaces

Run from the worktree root:

```text
python -B 04-computation/experiments/collatz_charged_terminal_family_20261012a.py
python -B -O 04-computation/experiments/collatz_charged_terminal_family_20261012a.py
```

The self-contained algebra tests cover 341 heads, alphabet 1 through 4 and
length 0 through 4, three completed members each; 17 members of (5); 257
members of each coarse guard; exact arbitrary-modulus readers; source, word and phase
type hostiles; first-hit rejection; the three finite dispatchers; and the
declared band universe. Exact types reject Python bool/float aliases before
arithmetic. Normal and optimized runs reproduce **523,424 checks** and identical
LF-normalized output, SHA256
`8e6674b492764e984faccd94a2c90d342503d62707815ce86ae1aaae11b87b60`.
An independent peer proof/API audit and optimized replay passed with the same
output, including the strict first-hit phase, both coarse guards, the uniform
measurement estimate, and the total invoice/failure distinction.

No new theorem is inferred from the number of checks. The next useful target
is a proof or a counterexample for the source-owned quadratic band invoice,
or another all-height paid exit for pending 127. Both are more precise than
extending a bank of completed shadows while changing its marked source.
