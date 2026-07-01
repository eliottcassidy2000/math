---
id: HYP-3768
title: B2/E2/DEDEKIND on the hard residual -- the covering-min construction's Dedekind sum has a PROVED closed form s(n,Phi6(n)) = -(Phi6-1)/(12*Phi6) -> -1/12 = zeta(-1) = -B2/2 (the E2 anomaly constant), so the covering-min IS the iota-EVEN Eisenstein/E2 'bulk', and the genuinely hard residual is the iota-ODD genus cusp form (HYP-3587). Three bridges: (1) B2=1/6 unifies the witness-sum heartbeat zeta(2)=pi^2*B2 (HYP-3746), E2=1-(4/B2)sum sigma_1 q^n, and Dedekind reciprocity's B2/2=1/12. (2) The SAWTOOTH ((x)) is iota-ODD (((-x))=-((x))) and ||x||=1/2-|((x))|, so loneliness (all ||st|| large) <=> every runner st near the HALF-INTEGER (the iota-fixed point) -- Dedekind's atom = the S55 sign atom. (3) DEDEKIND RECIPROCITY s(h,k)+s(k,h)=-1/4+(B2/2)(h/k+k/h+1/hk) is the TWO-MODULUS gluing of the multi-metric sheaf (HYP-3766); the covering-min descends in ONE step because Phi6==1 mod n (the Farey-neighbor condition HYP-3732), giving the closed form. The E2 quasimodular ANOMALY = the residual obstruction; genus X_0(2p)=0,0,1,2,2 (p=3,5,7,11,13) grows as the Eisenstein bulk stays constant (=#cusps-1=3): n=14 => genus-1 cusp form f_14 (elliptic curve 14a) = the iota-odd degree of OPEN-Q-108
status: MIXED (one PROVED closed form + bridges + honest identification). PROVED (sympy-verified + reciprocity derivation): s(n,Phi6(n)) = -(Phi6-1)/(12 Phi6) via ONE Dedekind reciprocity step using Phi6==1 mod n and s(1,n)=(n-1)(n-2)/(12n); limit -1/12=zeta(-1). VERIFIED: B2 unification (exact); ||x||=1/2-|((x))| + ((x)) iota-odd (exact); s(a,D) at bindings share the modulus denominator. HONEST NEGATIVE: the Dedekind sum at the binding does NOT separate tight/non-tight (AP & GW same binding => same s; not a tightness detector). The residual=genus-cusp identification is a synthesis of HYP-3587/3041 + S55, NOT a proof; reciprocity-as-gluing is a mechanism, not yet a termination proof of the crux.
source: klein-2026-06-30-S56
depends_on:
  - HYP-3767   # sign theory: iota-odd degree = the residual (the sawtooth is the iota-odd atom)
  - HYP-3587   # genus X_0(2p) = local-global gap = cusp-form obstruction (the residual's dimension)
related:
  - HYP-3746   # witness-sum zeta(2) heartbeat = pi^2 B2
  - HYP-3766   # multi-metric sheaf (reciprocity = its two-modulus gluing)
  - HYP-3732   # covering-min Farey neighbor of 1/(n-1), D==1 mod (n-1) ~ Phi6==1 mod n
  - HYP-3041   # cusps of X_0(14) = Klein four-group; genus jump at N=14
  - THM-584    # complement = antipodal iota (the sign involution)
results:
  - 04-computation/dedekind_e2_lrc_residual_klein.py
  - 05-knowledge/results/dedekind_e2_lrc_residual_klein.out
  - 05-knowledge/results/dedekind_e2_separation_klein.out
---

# HYP-3768 — B2/E2/Dedekind on the hard residual

## The one constant: B2 = 1/6
`B2` (the second Bernoulli number) runs the whole modular package AND the LRC witness sum:
- **witness-sum heartbeat** `zeta(2) = pi^2/6 = pi^2 * B2` (the Farey-grid reach density `1/zeta(2)`, HYP-3746);
- **weight-2 Eisenstein** `E2 = 1 - (4/B2) sum sigma_1(m) q^m = 1 - 24 q - ...`;
- **Dedekind reciprocity** constant `1/12 = B2/2`.
One Bernoulli ties the runner side to the modular side.

## The sawtooth is the iota-odd atom (bridge to S55)
The Dedekind sawtooth `((x)) = x - floor(x) - 1/2` (`0` at integers) satisfies:
- `((-x)) = -((x))` -- **iota-ODD** (the sign atom of HYP-3767);
- `||x|| = 1/2 - |((x))|` -- so a **lonely** observer (all `||st||` large) has every runner `st` near the
  **half-integer `1/2`** = the iota-FIXED point. Loneliness clusters the runners at the sign-fixed point;
  Dedekind sums correlate the same sawtooths that measure loneliness.

## PROVED: the covering-min's Dedekind sum closed form
> **Lemma.** At the covering-min construction binding `(D=Phi6(n)=n^2-n+1, a=n)`,
> `s(n, Phi6(n)) = -(Phi6-1)/(12 Phi6) = -n(n-1)/(12(n^2-n+1))`, and `-> -1/12 = zeta(-1) = -B2/2`.

**Proof (one reciprocity step).** `Phi6 = n^2-n+1 ≡ 1 (mod n)` (since `n^2-n=n(n-1)`), so `gcd(n,Phi6)=1`
and `s(Phi6, n) = s(1, n) = (n-1)(n-2)/(12n)`. Dedekind reciprocity
`s(n,Phi6)+s(Phi6,n) = -1/4 + (1/12)(n/Phi6 + Phi6/n + 1/(n Phi6))` then gives, after simplification
(sympy-verified), `s(n,Phi6) = -(Phi6-1)/(12 Phi6)`. ∎

So the covering-min construction's Dedekind sum is a **single Farey-neighbor reciprocity step** away from
the base `s(1,n)`, and its limit is exactly the **E2 anomaly constant** `-1/12`. The covering-min is the
`iota`-EVEN **Eisenstein/E2 "bulk"**; its arithmetic (the Dedekind sum) is the E2 anomaly.

## Reciprocity = the sheaf gluing (the descent)
Dedekind reciprocity `s(h,k) <-> s(k,h)` is the Euclidean/continued-fraction step -- the **two-modulus
GLUING** of the multi-metric danger sheaf (HYP-3766) and the three-gap renormalization (HYP-3762). The
covering-min glues in ONE step (`Phi6 ≡ 1 mod n`, the Farey-neighbor condition HYP-3732); a general
binding descends its CF `[n-1, a(n)]`. This is the arithmetic realization of the sheaf's CRT gluing.

## The genuinely hard residual = the genus cusp form
The E2/Eisenstein object is the COVERED, `iota`-EVEN bulk; the residual is the `iota`-ODD cusp-form
obstruction (HYP-3587): for `n=2p`, `genus X_0(2p) = 0,0,1,2,2` (`p=3,5,7,11,13`) GROWS while the
Eisenstein bulk (`#cusps-1=3`) stays constant. **`n=14` => genus 1 => the weight-2 cusp form `f_14` =
elliptic curve `14a`** = the `iota`-odd degree of OPEN-Q-108 / the danger sheaf's antipodal Cech class
(HYP-3767). So:
> covering-min construction  =  E2/Eisenstein bulk  (Dedekind sum `-> -1/12`, iota-EVEN, COVERED);
> the hard residual          =  the genus cusp form `f_14` (iota-ODD, the uncovered degree).

## Honest scope
PROVED: the `s(n,Phi6)` closed form (reciprocity). VERIFIED: the B2 unification, the sawtooth bridge.
HONEST NEGATIVE: the Dedekind sum at the binding does NOT separate tight from non-tight (AP and GW share
the binding `(14,1)` hence `s=13/14`), so it is not a tightness detector. The residual=genus identification
is a synthesis (HYP-3587/3041 + the S55 sign theory), and reciprocity-as-gluing is a mechanism, NOT a
termination proof of the coverage crux. What is new and solid: the covering-min IS the E2/Eisenstein
Dedekind object (closed form, `-> -1/12`), which LOCATES the hard residual precisely as the complementary
genus cusp form.

## Net
Through B2 (`= 1/6`), the LRC witness sum (`zeta(2)=pi^2 B2`), the weight-2 Eisenstein `E2`, and Dedekind
reciprocity are one object; the iota-odd sawtooth is their shared atom and the S55 sign atom. The
covering-min construction's Dedekind sum is `-(Phi6-1)/(12 Phi6) -> -1/12 = zeta(-1)` (PROVED, one
reciprocity step via `Phi6≡1 mod n`) -- it is the `iota`-even Eisenstein/E2 bulk. The genuinely hard
residual is therefore the complementary `iota`-odd **genus cusp form** (`n=14`: curve `14a`), the
Borsuk-Ulam odd index of OPEN-Q-108 -- now identified as a specific weight-2 newform, not just an abstract
obstruction.
