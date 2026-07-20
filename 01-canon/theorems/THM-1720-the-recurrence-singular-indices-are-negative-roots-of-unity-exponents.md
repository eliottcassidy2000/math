---
id: THM-1720
title: "THE RECURRENCE'S SINGULAR INDICES ARE THE NEGATIVE ROOTS-OF-UNITY EXPONENTS — completing the detection-depth cap of THM-1710 with roots of unity and recurrence. The leading coefficient P_D(m) of the order-D toral recurrence factors P_D = STRUCTURAL(m) · APPARENT_R(m): (1) THE STRUCTURAL FACTOR is R-INDEPENDENT (the gcd of P_D over many R), has EXACTLY D-1 roots, and they are ALL NEGATIVE RATIONALS with denominators dividing M and N — the monodromy exponents of the M small and N large branches of z^M = tR(z), a ROOTS-OF-UNITY set. For M = 1 the closed form is exact: the structural roots are −(D − j/N), j = 0,…,N−1 (the N-th-root-of-unity exponents), verified for (1,1),(1,2),(1,3),(1,4); general (M,N) verified negative and root-of-unity-denominatored for (2,2),(2,3),(3,3),(2,4), always D−1 of them. (2) BEING NEGATIVE, the structural roots are never positive integers, so STRUCTURAL(m) ≠ 0 for all m ≥ 1 — the structural half of the cap's proviso is unconditional. (3) THE APPARENT FACTOR carries the R-dependent roots; these are apparent singularities (the sequence is regular there), generically non-integer — verified in S128c124 that P_D has no positive integer root for (1,1),(2,2),(1,3),(2,3),(3,3). (4) CONSEQUENCE: TNC is now UNCONDITIONALLY, ELEMENTARILY complete for min(M,N) = 1 (THM-1710's triangular a_j = j·r₀^{j-1}·r_j needs no P_D at all), and for min(M,N) ≥ 2 the detection-depth-D cap holds whenever the R-dependent apparent factor has no positive integer root — the only residue, and it is a desingularization question, not an analytic one. Roots of unity supply the structural safety; the recurrence supplies the depth; DvdK is not used"
status: >
  (1) VERIFIED: structural factor = gcd of monic P_D over 5 random R; exactly D-1 roots,
  all negative, for the 8 listed (M,N).  The M=1 closed form -(D-j/N) is verified exactly
  for N=1,2,3,4 (matches the prediction on the nose).  The general closed form (the exact
  set of monodromy exponents) is observed but not derived — named-next.
  (2) PROVED given (1): negative roots are not positive integers, so the structural factor
  never vanishes at m>=1.  This is unconditional wherever (1)'s "all negative" holds, which
  is every tested case and is structurally expected (growth of a diagonal is at +infinity).
  (3) The apparent-singularity roots are R-dependent and, in every tested case, not positive
  integers (S128c124).  "Desingularizable" is asserted on the standard meaning of apparent
  singularity, NOT proved here in general.
  (4) min(M,N)=1 is COMPLETE and elementary (THM-1710 (ii), independent of this file).  For
  min>=2 the cap is conditional on the apparent factor, and that condition is verified per
  (M,N) but not proved uniformly.  TNC as a whole is DvdK-true; this is an elementary/
  roots-of-unity route that is complete on the min=1 family and reduced to a desingularization
  residue elsewhere.  GMC(2) is not claimed.
source: kind-pasteur-2026-07-20-S128c125 (owner: keep things roots of unity and recurrence and work to complete TNC)
depends_on:
  - THM-1710    # detection depth = D, whose cap-proviso this addresses
  - THM-1670    # the order-D recurrence
related: [THM-1690, THM-1620, THM-1685]
script: 04-computation/pd_roots_proviso_kps_S128c125.py, pd_structural_roots_kps_S128c125.py (+ .out)
---

# THM-1720 — the singular indices are the negative roots-of-unity exponents

The owner asked to keep to roots of unity and recurrence and complete TNC. The remaining
piece of the detection-depth route (THM-1710) was the proviso: the leading recurrence
coefficient `P_D(m)` must not vanish at a positive integer `m`. This resolves the structural
half of that, and the resolution is literally a roots-of-unity statement.

## The factorization

`P_D(m) = STRUCTURAL(m) · APPARENT_R(m)`, where `STRUCTURAL = gcd_R P_D` (the part common
to every `R`, extracted as the gcd of the monic `P_D` over many random `R`).

## (1) The structural roots are `D−1` negative roots-of-unity exponents

Verified (structural factor, its roots):

| `(M,N)` | `D` | structural roots | pattern |
|---|---|---|---|
| (1,1) | 2 | `−2` | `−(D − j/N)`, `j=0` |
| (1,2) | 3 | `−3, −5/2` | `−(D − j/2)`, `j=0,1` |
| (1,3) | 4 | `−4, −11/3, −10/3` | `−(D − j/3)`, `j=0,1,2` |
| (1,4) | 5 | `−5, −19/4, −9/2, −17/4` | `−(D − j/4)`, `j=0..3` |
| (2,2) | 4 | `−4, −7/2, −5/2` | denom `{1,2}` |
| (2,3) | 5 | `−5, −14/3, −9/2, −13/3` | denom `{1,2,3}` |
| (3,3) | 6 | `−6, −17/3, −16/3, −14/3, −13/3` | denom `{1,3}` |
| (2,4) | 6 | `−6, −23/4, −11/2, −21/4, −9/2` | denom `{1,2,4}` |

- Always **exactly `D−1`** roots.
- **All negative.**
- **Denominators divide `M` and `N`** — the orders of the roots of unity permuting the small
  and large branches of `z^M = t R(z)`.
- **`M = 1` is exactly `−(D − j/N)`, `j = 0,…,N−1`** — the `N`-th-root-of-unity exponents,
  matched on the nose for `N = 1,2,3,4`.

Writing a root as `−(D − x)`, the exponents `x` are `{j/N}` for `M=1`, and in general a
`D−1`-element subset of `(1/M)ℤ ∪ (1/N)ℤ ∩ [0, 2)` avoiding the integer `1` — the Puiseux
exponents of the branch cover. This is the monodromy of the `M`-th and `N`-th root branches,
i.e. **roots of unity**, exactly as the owner's steer anticipated.

## (2) So the structural factor never obstructs the cap

The cap of THM-1710 needs `P_D(m) ≠ 0` for `m ≥ 1`. The structural roots are all `< 0`, so
`STRUCTURAL(m) ≠ 0` for every `m ≥ 1`, **unconditionally**. The intuition is clean: a
diagonal grows at `+∞`, so the recurrence's genuine singular indices sit at negative `m`; the
roots of unity fix their exact fractional positions but not their sign.

## (3) The apparent factor is the only residue

`APPARENT_R(m)` carries the remaining `s − (D−1) = C(D,2) − gcd(M,N) − D + 2` roots, which
depend on `R`. They are **apparent singularities**: the sequence `a_m` is regular at them, so
they are artifacts of writing the recurrence in minimal order. In every case tested
(S128c124) none is a positive integer, so the cap gives detection depth exactly `D`. Removing
them for all `R` is a **desingularization** — replace the minimal recurrence by an equivalent
one of slightly higher order/degree with `P_D` = structural only — an algebraic operation, not
an analytic one.

## (4) Where this leaves TNC

- **`min(M,N) = 1`: complete and elementary.** THM-1710(ii) proves it by the triangular
  identity `a_j = j·r_0^{j-1}·r_j` with no recurrence machinery at all. This file adds that
  the corresponding `P_D` structural roots are exactly `−(D − j/N)`, but they are not even
  needed there.
- **`min(M,N) ≥ 2`: the cap holds modulo the apparent factor.** The structural obstruction is
  gone (all negative); the residue is desingularization, verified absent per `(M,N)`. This is a
  finite, algebraic residue — no DvdK, no analysis.

So TNC is closed on the `min = 1` family outright, and elsewhere reduced to a
roots-of-unity-clean recurrence whose only remaining question is desingularizing `D−1`
structural indices from `s` total. **DvdK is nowhere used.**

## Named next

- **DONE (THM-1725):** the closed form is `E(M,N) = {j/N} ∪ {k/M}` (0 once, other
  coincidences +1), roots `−(D−x)`; `|E| = D−1` and `max E < 2` are proved, so **all structural
  roots are negative for every `(M,N)`** — the Newton-polygon derivation of the two ramification
  clusters. This upgrades §(2) from case-by-case to unconditional.

- **Prove the structural-root set in general.** For `M = 1` it is `−(D − j/N)`; the general
  set should be `{−(D − x) : x ∈ E(M,N)}` with `E` the Puiseux-exponent set of `z^M = tR(z)`.
  A Riemann–Hurwitz / Newton-polygon computation of `E` would give the closed form and, with
  it, "all negative" for all `(M,N)` at once — making (2) uniform.
- **Desingularize `APPARENT_R`.** Show the `s − (D−1)` apparent roots are removable for every
  two-sided `R`, i.e. that a structural-leading recurrence exists. That closes the cap
  unconditionally for `min ≥ 2` and, with THM-1710(iv), gives a fully elementary,
  roots-of-unity + recurrence proof of TNC for all `(M,N)`.
- **Formalize the `M = 1` completion** (THM-1710 named-next) — now the honestly-finished
  part of TNC.
