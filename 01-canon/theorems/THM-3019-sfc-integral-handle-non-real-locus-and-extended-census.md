---
id: THM-3019
title: "Univariate slot-SFC: integral handle, non-real locus, an infinite 2-slot family, and a census"
status: >
  PROVED + VERIFIED-EXACT / AWAITING INDEPENDENT HOSTILE AUDIT. Four results
  on the Strong Factorial Conjecture with L(z^n) = n!. (S1) L(f^m) =
  int_0^inf f(t)^m e^{-t} dt for every f in C[z] and m >= 1 -- an analytic
  handle for a lane that has been purely algebraic. (S2) For REAL coefficient
  vectors and m EVEN, L(f^m) > 0; since any window of N >= 2 consecutive
  moments contains an even index, every common zero of a consecutive window
  is NON-REAL. This holds for every slot count N and every window k.
  (S3) For N = 2 the even-index member therefore has no real root, so
  Res(I_{k+1}, I_{k+2}) is a product over conjugate pairs of |.|^2 and is
  ALWAYS >= 0: SlotSFC_1(2) at window k is equivalent to Res > 0, never to a sign
  condition. (S4) For the 2-slot family f = a + b z the exact recurrence
  I_m = 1 + m lambda I_{m-1} holds, so I_m = I_{m+1} = 0 forces
  I_{m+1} = 1: SlotSFC_1(2) is PROVED for that family at EVERY window, unbounded.
  Census (one-way modular / rank certificates): SlotSFC_1(2) verified on 648 cells
  with p <= 8, s <= 8, k <= 8 (top exponent q <= 16, windows to 8), and
  SlotSFC_1(3) verified on 1430 cells with 0 <= p < q < r <= 12 and k <= 4 by
  Macaulay surjectivity.  The latter is now a contained control for
  THM-2836's later support-12/window-8 extension.  The unbounded
  SlotSFC_1(3) remains OPEN.
source: death-star-2026-07-31-coinC2
depends_on: []
related:
  - THM-2836
  - THM-2824
  - THM-2812
  - THM-3018
external:
  - "Edo, van den Essen: the Strong Factorial Conjecture."
script: 04-computation/sfc2_unbounded_structure_thm3019.py
output: 05-knowledge/results/sfc2_unbounded_structure_thm3019.out
script2: 04-computation/sfc3_macaulay_census_thm3019.py
output2: 05-knowledge/results/sfc3_macaulay_census_thm3019.out
---

# THM-3019 -- univariate slot-SFC: integral handle, non-real locus, and an extended census

Throughout `L : C[z] -> C` is linear with `L(z^n) = n!`, and for a support
`p_1 < ... < p_N` we write `f = sum_i a_i z^{p_i}`. `SlotSFC_1(N)` at **window
`k`** asserts that the `N` consecutive moments `L(f^{k+1}), ..., L(f^{k+N})`
have no common nonzero zero.  This is the local notation fixed by
MISTAKE-350 for an `N`-monomial restriction inside ambient `SFC(1)`; it is
not the original ambient conjecture `SFC(N)`.

## 1. (S1) The functional is an integral (PROVED)

Since `int_0^inf t^n e^{-t} dt = n!`, linearity gives

```text
L(F) = int_0^inf F(t) e^{-t} dt      for every F in C[z],
so   L(f^m) = int_0^inf f(t)^m e^{-t} dt.                            (1)
```

Explicitly `L(f^m) = sum_{|alpha| = m} (m choose alpha) a^alpha (alpha . p)!`,
the multinomial form used in the lane. (1) is the analytic handle: the SFC
lane (THM-2812/2824/2836/2849/2854) works with Macaulay matrices and graded
ranks; (1) adds positivity and asymptotics.

## 2. (S2) Common zeros of a consecutive window are non-real (PROVED)

If all `a_i` are **real** then `f` is a real polynomial and, for `m` even,
the integrand in (1) is `>= 0` and not identically zero, so

```text
L(f^m) > 0     for every real a != 0 and every even m.               (2)
```

Any window of `N >= 2` consecutive moments contains an even index. Hence:

**Corollary.** For every slot count `N >= 2`, every support, and every window
`k`, a common zero of `L(f^{k+1}), ..., L(f^{k+N})` must have **non-real**
coefficient vector (up to scaling). Equivalently the real locus of the SFC
system is empty, always.

## 3. (S3) The `N = 2` resultant is never negative (PROVED)

For `N = 2`, scale `a_1 = 1` (the case `a_1 = 0` is a monomial, where
`L(f^m) = a_2^m (p_2 m)! != 0`), put `s = p_2 - p_1`, `lambda = a_2/a_1`:

```text
L(f^m) = a_1^m I_m(lambda),   I_m(lambda) = sum_j C(m,j) (p_1 m + s j)! lambda^j,
```

a real polynomial of degree `m` with **positive** coefficients. By (2),
`I_m` has no real root when `m` is even. One of `k+1, k+2` is even; taking
the resultant from that side writes it as `lc^{...}` times a product over
conjugate root pairs of `|I(.)|^2`. Hence

```text
Res(I_{k+1}, I_{k+2}) >= 0 always, and SlotSFC_1(2) at window k  <=>  Res > 0.  (3)
```

This explains the census below, in which every computed resultant sign is
`+1`: a *negative* resultant is impossible, so the only failure mode is exact
vanishing.

## 4. (S4) An infinite family, proved at every window (PROVED)

For the 2-slot family `f = a + b z` (`p_1 = 0`, `s = 1`), `g(t) = 1 + lambda t`
has `g' = lambda` constant, so integrating (1) by parts,

```text
I_m = [-g^m e^{-t}]_0^inf + m int g^{m-1} g' e^{-t} dt
    = 1 + m lambda I_{m-1},        I_0 = 1.                          (4)
```

(Verified symbolically for `m = 1..8`.) If `I_m(lambda) = I_{m+1}(lambda) = 0`
then (4) at `m+1` gives `I_{m+1} = 1 + (m+1) lambda I_m = 1`, contradiction.
So **SlotSFC_1(2) holds for `f = a + b z` at every window `k`, with no bound on
`k`** -- the first unbounded-window family in this lane. Equivalently
`I_m(lambda) = m! lambda^m e_m(1/lambda)` with `e_m` the truncated
exponential, and `e_{m+1} = e_m + x^{m+1}/(m+1)!` forces a common root to be
`x = 0`, where `e_m(0) = 1`.

## 5. Extended censuses (FINITE-EXACT, one-way certificates)

**SlotSFC_1(2).** `Res(I_{k+1}, I_{k+2}) != 0` certified modulo large primes
(nonzero mod one prime is a sound one-way certificate) for

```text
0 <= p <= 8,  1 <= s <= 8,  0 <= k <= 8      -- 648 cells, 0 failures,
```

i.e. top exponent `q = p + s <= 16` and windows to `8`.

**SlotSFC_1(3).** For `f = a z^p + b z^q + c z^r`, `L(f^m)` is the form
`M_m(a,b,c) = sum_{i+j+l=m} m!/(i!j!l!) (pi+qj+rl)! a^i b^j c^l`. If the
Macaulay map `S_{D-d_1} (+) S_{D-d_2} (+) S_{D-d_3} -> S_D` with
`d_i = k+i`, `D = sum(d_i - 1) + 1`, is **surjective** over a field, the ideal
contains `S_D` and the projective variety is empty; surjectivity mod a prime
implies it over `Q` (rank only drops under reduction). Certified for

```text
0 <= p < q < r <= 12,  0 <= k <= 4          -- 1430 cells, 0 failures.
```

The `SlotSFC_1(2)` box extends the earlier two-slot tests in both support and
window.  The three-slot box originally extended THM-2836's 9/6 core in the
support direction; it is now contained in THM-2836's later 12/8 extension and
serves as an independently generated sub-box.

## 6. Scope

(S1)-(S4) are proofs; section 5 is a finite census with one-way certificates
and proves nothing outside its box. **The unbounded SlotSFC_1(3) remains OPEN**, as
does unbounded SlotSFC_1(2) outside the family of section 4. What (S2) contributes
to the open problem is that the search may be restricted to the non-real
locus at every slot count and window, and what (S3) contributes is that for
`N = 2` only exact vanishing can occur, never a sign change -- so a proof of
SlotSFC_1(2) needs a strict-positivity mechanism, not a discriminant sign analysis.
Section 4 shows such a mechanism exists when the leading factor `g'` is
constant; the obstruction to generalising is that for `s >= 2` or `p >= 1`,
`u g' = c g` has no polynomial solution `u`, so integration by parts leaves
the family `K(alpha, beta) = int t^alpha (1 + lambda t^s)^beta e^{-t} dt`
rather than closing on the `I_m`.

Note the type distinction with THM-3018: that file concerns `n` variables
with `L(x^alpha) = alpha!` and homogeneous `f`; this one is the repo's
one-variable, term-count-indexed SFC. They agree only at `n = 1`.
