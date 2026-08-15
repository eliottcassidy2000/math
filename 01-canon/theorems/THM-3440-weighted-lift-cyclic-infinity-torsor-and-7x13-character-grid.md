---
id: THM-3440
title: "Weighted-lift cyclic infinity torsor and the exact 7x13 character grid"
status: >
  RESERVED / PROVISIONAL PROOF CANDIDATE / INDEPENDENT AUDIT REQUIRED.
  The degree-n weighted lift of THM-3438 has a totally ramified geometric
  place of index n over Q=infinity.  Its local branch monodromy is an n-cycle.
  At n=91, CRT turns this branch torsor into an exact C7 x C13 Fourier carrier
  and supplies a typed H^1(C13,F13) boundary bridge after a generator gauge.
  It does not identify LRC amplitudes with Keller amplitudes, prove the LRC
  7-tensor-13 bispectrum nonzero, construct a physical current, or describe a
  finite Jelonek component.
source: codex2 weighted-lift reciprocal-Eisenstein derivation, 2026-08-15
depends_on:
  - THM-3438-weighted-lift-keller-degree-spectrum
related:
  - THM-3431-valuative-persistence-multiset-and-lrc-jc-boundary
  - THM-3437-derived-boundary-jet-euler-conservation-and-prufer-recovery
script: 04-computation/jc_weighted_lift_infinity_inertia_thm3440.py
output: 05-knowledge/results/jc_weighted_lift_infinity_inertia_thm3440.out
---

# THM-3440 -- weighted-lift cyclic infinity torsor and the exact 7x13 character grid

**RESERVED / PROVISIONAL PROOF CANDIDATE / INDEPENDENT AUDIT REQUIRED.**

## 1. Reciprocal Eisenstein law

Use the weighted lift of THM-3438.  Write

```text
R(w)=integral_0^w p(s) ds=sum_(j=1)^n r_j w^j,
n=deg(p)+1,                 r_n!=0,
T(w)=R(w)-Pw+cQ.                                      (1)
```

The source function field is the target field with one root of `T` adjoined.
Fix `P` generically and take the target valuation `Q=infinity`, with
uniformizer

```text
t=Q^(-1).
```

Put `s=w^(-1)`.  Reciprocal reversal and clearing by `t s^n` gives

```text
H(s)=t s^n T(s^(-1))
    =c s^n+t(sum_(j=1)^n r_j s^(n-j)-P s^(n-1)).       (2)
```

Over the DVR `k(P)[[t]]`, the leading coefficient `c` is a unit, every lower
coefficient is divisible by `t`, and the constant coefficient `t r_n` is not
divisible by `t^2`.  Thus `(2)` is Eisenstein.  There is one geometric place
above `Q=infinity`, its ramification index is exactly `n`, and `s` is a
uniformizer there.  Characteristic zero makes the ramification tame.

After strict henselization (or analytically after fixing a generic complex
value of `P`), tame inertia is cyclic and acts on the `n` local branches as a
single `n`-cycle.  Equivalently, the branch set over an oriented punctured
`t`-disk is a `C_n` torsor.  This is local geometric monodromy; it does not say
that the global degree-`n` field extension is cyclic or Galois.

For the explicit maps `F_n` of THM-3438 one may take the target line

```text
(A,B,C)=(Q,P,1),
```

because `P=BC` and `Q=AC^2`.  Reconstruction from `w` is valid near this
place: `w` tends to infinity and the leading term of `p(w)` keeps `gamma`
nonzero.  Hence the Eisenstein place belongs to the Keller cover itself, not
only to an auxiliary resolvent.

## 2. The degree-91 branch grid

Choose `n=91=7*13`.  THM-3438 supplies a determinant-one polynomial Keller
map `F_91:A^3->A^3`, and Section 1 gives its marked local branch torsor

```text
B_91 ~= C_91.                                           (3)
```

Once an oriented inertia generator `tau` is chosen, CRT gives

```text
C_91 -> C_7 x C_13,
j     |-> (j mod 7,j mod 13).                           (4)
```

The normalized factor generators inside `C_91` are

```text
(1,0) <-> tau^78,             (0,1) <-> tau^14,         (5)
```

since `78=1 mod 7, 0 mod 13` and `14=0 mod 7, 1 mod 13`.
Consequently, over any splitting field containing `zeta_91`, the branch
Fourier space factors exactly:

```text
k[C_91] ~= k[C_7] tensor_k k[C_13],
C_91^hat ~= C_7^hat x C_13^hat.                         (6)
```

Thus every `(alpha,beta)` character in a `7 tensor 13` spectral table has a
literal branch character on this Keller boundary torsor.  Relabelling by
`(4)` preserves character multiplication, convolution, Fourier support, and
the assertion that a specified coefficient is zero or nonzero.

## 3. The typed H1 bridge

The characteristic mismatch in THM-3431/3437 rules out a nonzero additive
map from the LRC class `H^1(C_13;F_13)` into a characteristic-zero Tor module.
The torsion object `(3)` avoids that type error.  Given

- a marked LRC deck generator `sigma_13`, and
- the oriented Keller inertia generator `tau`,

define

```text
iota_13:C_13(LRC)->C_91(Keller),
iota_13(sigma_13)=tau^14.                               (7)
```

Its image is the CRT `C_13` factor.  Pullback along `(7)` identifies the
one-dimensional classes

```text
H^1(<tau^14>;F_13) ~= H^1(C_13;F_13).                   (8)
```

This is an explicit local monodromy-class bridge.  Changing either generator
rescales `(8)` by a unit of `F_13`; without the orientation/marking sidecar
there is only a projective class.  Formula `(8)` is not THM-3437's homological
`Tor_1`, not de Rham flux, and not a physical LRC current.

## 4. Connection and loss ledger

| field | exact content |
|---|---|
| source | the degree-91 weighted inverse equation and its `Q=infinity` place |
| target | a marked `C_91` branch torsor, CRT grid, and local `H^1(C_13;F_13)` class |
| map | reciprocal Eisenstein, oriented inertia, then CRT `(4)--(8)` |
| preserved | cyclic order, branch ancestry around the loop, character products, Fourier support, and zero/nonzero under pure relabelling |
| destroyed | source coordinates away from the punctured disk, finite Jelonek geometry, LRC word weights, response amplitudes, positivity, and physical time |
| needed sidecar | an amplitude/intertwiner identifying an LRC response array with a Keller branch function |
| cheapest hostile | reverse the inertia generator: the class rescales, showing why a generator-free scalar equality is illegal |
| decisive next test | construct or obstruct a natural amplitude on `B_91` whose `7 tensor 13` coefficient equals the LRC contraction |

## 5. Boundaries

- The place is over target infinity.  It is not, by itself, a component of the
  finite Jelonek set and proves no componentwise discriminant-parity law.
- The inertia permutation is an `n`-cycle.  No claim is made that the global
  monodromy group is `C_n`; generically it may be much larger.
- CRT supplies a faithful representation carrier, not an equality of
  amplitudes.  Therefore the LRC `7 tensor 13` bispectrum nonvanishing remains
  exactly as open as before this theorem.
- Degree `91` is numerically factorable.  The local `91`-cycle does not alone
  decide whether the global map is composition-irreducible.
- Nothing here proves LRC(14), `JC(2)`, a D5 physical-current map, or a global
  comparison of boundary cohomology theories.

## 6. Exact companion

The companion must verify the reciprocal identity, Eisenstein valuations,
ramification index, CRT exponents `(5)`, all `91` character addresses, and
hostile cases where `r_n=0`, `c=0`, or a generator is forgotten.  Until that
audit and an independent proof review are complete, this file remains a
provisional candidate.
