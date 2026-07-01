---
id: HYP-3767
title: SIGN THEORY ON THE OPEN CRUX -- the antipodal involution iota (a->-a = t->1-t = complement, THM-584) reduces the multi-metric-sheaf coverage crux to an iota-ORBIT set-cover, halving it, and identifies every witness as an iota-symmetric hole. Since ||x|| is EVEN, the danger zone of each speed s, D_s(r)={a: sa in {-r..r} mod D}, is iota-SYMMETRIC, so a rotation a is covered IFF -a is: coverage is a question over the orbit space (Z/D)/iota. (P) PARITY LEMMA (PROVED): for ODD D, iota is fixed-point-free on nonzero residues, so #lonely rotations at radius r is EVEN and lonely rotations come in +-PAIRS; a witness = an uncovered iota-ORBIT. The q-witness (radius 0) = the iota-fixed hole {0}; the k- and (n+q)-witnesses (radius 1) = the iota-pair {+1,-1} (verified: the (n+q)-witness hole is exactly {+1,-1} mod n+q). (REDUCTION) each speed covers <=r nonzero orbits at radius r; the crux 'danger covers the site' = an orbit set-cover over the (D-1)/2 nonzero orbits by n-1 speeds -- HALF the rotation space, the natural ODD/EVEN (sign) grading. (BRIDGE) the mod-2 / iota-equivariant structure is the LRC instance of the Borsuk-Ulam 'odd index' (OPEN-Q-108, THM-584): the crux = an iota-equivariant degree over the orbit space vanishes. HONEST: this SHARPENS (orbit halving + parity 'within-one-orbit') but does NOT close -- radius>=2 witnesses genuinely exist (overlaps beat the naive orbit-count: e.g. drop-11 lonely at D=41 r=2, D=53 r=3), so the orbit set-cover is not trivially full
status: MIXED (proved lemma + reduction + honest bridge). (P) PARITY PROVED (iota fixed-point-free on (Z/D)*\{0} for odd D => #lonely even; even-D caveat: a=D/2 is iota-fixed) + VERIFIED odd D=3..59 all sets. iota-symmetry of D_s(r) EXACT (||-x||=||x||). Witness=iota-orbit VERIFIED ((n+q)-hole={+1,-1}). REDUCTION to orbit set-cover is exact (iota-symmetry). BRIDGE to Borsuk-Ulam is aspirational (the per-modulus parity is trivially even, so the useful invariant is the GLOBAL iota-equivariant degree over all moduli = OPEN-Q-108, open). HONEST NEGATIVE: no radius<=1 localization (radius>=2 witnesses exist; the naive (n-1)r>=(D-1)/2 orbit-count fails on overlaps).
source: klein-2026-06-30-S55
depends_on:
  - THM-584    # complement R = the ANTIPODAL map (iota); the sign involution
  - HYP-3766   # the multi-metric witness sheaf (this is its sign-theory refinement)
related:
  - THM-523    # q-witness (radius 0 = the iota-fixed hole {0})
  - OPEN-Q-108 # the odd index as a Borsuk-Ulam degree (the global iota-equivariant aspiration)
  - HYP-3538   # R-eigenspace organizing principle (iota-even/odd split)
  - HYP-3741   # witness hierarchy (radius-r = the metrics)
results:
  - 04-computation/sign_theory_lonely_parity_klein.py
  - 05-knowledge/results/sign_theory_lonely_parity_klein.out
  - 05-knowledge/results/sign_theory_parity_general_klein.out
---

# HYP-3767 — sign theory on the open crux

## The sign involution
The crux (HYP-3766): prove the multi-metric danger sheaf COVERS the modulus site — no lonely point
survives all metrics. Sign theory enters through the **antipodal involution** `iota: a -> -a` on
rotations mod `D` (equivalently `t -> 1-t`, the complement map = THM-584's `R`). Because the distance is
even, `||(-a)s/D|| = ||as/D||`, each speed's **danger zone**
```
    D_s(r) = { a mod D : sa mod D in {-r,...,r} } = s^{-1}{-r,...,r}
```
is `iota`-symmetric. Hence **`a` is covered iff `-a` is covered**: coverage is not a question about
rotations but about their `iota`-orbits.

## (P) The parity lemma (proved)
> **Lemma.** For ODD `D`, `iota` is fixed-point-free on the nonzero residues (`a=-a => 2a≡0 => a≡0`,
> excluded), so the lonely set `L(D,r) = { a≠0 : all runners avoid {±1,...,±r} }` is `iota`-invariant
> with no fixed point: **`#L(D,r)` is EVEN**, and lonely rotations occur in `±` pairs `{a,-a}`.

Verified odd `D=3..59` (prime and composite), all sets. **Even-`D` caveat:** `a=D/2` is an `iota`-fixed
point, so the parity there is `#(free part)` `+` `[D/2 lonely]` and can be odd. So the clean parity is an
**odd-modulus** phenomenon — the same odd/even split that runs THM-584/HYP-3538.

## Every witness is an `iota`-symmetric hole
A witness is an uncovered `iota`-orbit:
- the **`q`-witness** (radius 0) is the `iota`-FIXED hole `{0}` (resonance: `0` uncovered);
- the **`k`-witness** and the **`(n+q)`-witness** (radius 1) are the `iota`-PAIR `{+1,-1}`. Verified: the
  `(n+q)`-witness hole is exactly `{1, n+q-1}` mod `n+q` (the residues `±1`), the pair the vacated
  speeds `q` (dropped) and `n` (out of range) fail to cover (HYP-3766). So "the pair `{q,n}` vacates
  `±1`" IS "the `iota`-orbit `{±1}` is uncovered."

## The reduction: coverage = an `iota`-orbit set-cover
By `iota`-symmetry, the crux over odd `D` is an **orbit set-cover on the quotient `(Z/D)/iota`**: the
`(D-1)/2` nonzero orbits must be covered by the `n-1` speeds, each of which covers `<= r` orbits at
radius `r` (`s^{-1}{±1},...,s^{-1}{±r}`). This **halves** the space — the natural sign (odd/even)
grading. The parity lemma is the statement that you cannot be "off by one rotation": you are off by
whole orbits, so a count that reaches "all but `<1` orbit" is complete.

## (BRIDGE) the Borsuk-Ulam odd index
The `iota`-equivariant danger cover is the LRC instance of OPEN-Q-108 / THM-584's program: the
"odd index" as a **Borsuk-Ulam degree** of the antipodal action. The crux — coverage — is the
statement that a certain `iota`-equivariant degree over the orbit space **vanishes**. HONEST: the
*per-modulus* parity is trivially `0 mod 2` (always even), so it is a constraint, not a detector; the
useful invariant is the *global* `iota`-equivariant degree across all moduli (the sheaf's `iota`-Cech
class), which is exactly OPEN-Q-108 and remains open. Sign theory says the crux LIVES in
`iota`-equivariant cohomology; it does not yet compute the class.

## Honest scope (what sign theory does NOT do)
No radius-`<=1` localization: radius-`>=2` witnesses genuinely exist (drop-`11` at `n=14` is lonely at
`D=41, r=2` and `D=53, r=3`), because danger zones OVERLAP and the naive orbit-count `(n-1)r >= (D-1)/2`
is an upper bound on coverage, not a lower one. So the orbit set-cover is not trivially full; the parity
+ orbit halving SHARPEN the crux (work over `(Z/D)/iota`, need only rule out uncovered orbits, parity
kills the off-by-one) but do not close it.

## Net
The antipodal `iota` reduces the sheaf-coverage crux to an `iota`-orbit set-cover over `(Z/D)/iota`
(sign halving); the parity lemma (odd `D`: `#lonely` even, `±`-paired) makes every witness an
`iota`-symmetric hole (the `(n+q)`-witness = the orbit `{±1}`), and forbids off-by-one coverage; the
whole structure is the LRC face of the Borsuk-Ulam odd index (OPEN-Q-108/THM-584), with the crux =
vanishing of a global `iota`-equivariant degree. This sharpens the crux to a signed/equivariant
covering problem on half the space — an honest step, not a closure.
