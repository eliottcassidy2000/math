---
id: THM-1655
title: "TNC PROVED FOR ALL BINOMIAL R, AND FOR EVERY R WITH A UNIQUE MINIMAL CHARGE-REPRESENTATION OF 0. (0) SIGN FIX to THM-1635: G(t) = +t (log Pi)', not minus (verified on R=1+u: t(log Pi)' with Pi=t/(1-t) gives exactly 1+t+t^2+...). The criterion Pi=ct is unchanged. (1) CT(Lambda^m) = [u^0] Lambda^m = sum over m-term charge-representations of 0, where Lambda = u^{-N}R has charges {k-N : r_k != 0}. CT(1) = r_N, so r_N = 0 (Lambda has NO constant term) is the first necessary condition -- the charge-0 condition of THM-1535, recovered. (2) BINOMIAL THEOREM, PROVED: for R = r_0 + r_d u^d (d=M+N, two terms), 0 has a UNIQUE minimal representation a*(-N)+b*M=0 with a=M/g, b=N/g, g=gcd(M,N), at m_0=(M+N)/g, and CT(m_0) = C(m_0; a,b) r_0^a r_d^b != 0 -- a single product, no cancellation. So NO binomial is a nullcone element: TNC holds for all two-term R at every bidegree. Verified 16/16 (N,M in 1..4) with exact predicted CT. (3) GENERAL SUFFICIENT CONDITION: if 0 has a UNIQUE minimal charge-representation, CT(m_0) is one nonzero product, so TNC holds -- no all-order argument needed. (4) THE RESIDUAL IS EXACTLY NON-UNIQUE MINIMAL REPRESENTATIONS. First trinomial where it occurs: charges {-2,1,4} (R=r_0+r_3 u^3+r_6 u^6, N=2), m_0=3 with reps {-2,-2,4} and {1,1,-2}, giving CT(3) = 3 r_0 (r_0 r_6 + r_3^2). TUNED r_0=1,r_3=1,r_6=-1 makes CT(3)=0 -- but CT(6) = -30 != 0, so R=1+u^3-u^6 is STILL not a nullcone element. Vertex cancellation starts, the branch-product/nondegenerate-saddle argument (THM-1635) finishes"
status: PROVED (1),(2),(3); (4) VERIFIED with explicit witness. Binomial and unique-minimal cases are closed at all bidegrees.
author: opus-2026-07-20-S419
depends_on: [THM-1635 (branch product; sign fixed here), THM-1625, THM-1615, THM-1535 (charge-0), THM-1530]
corrects: THM-1635 (sign of the G = t(log Pi)' identity)
---

# THM-1655 — TNC for binomials, and the unique-minimal-representation criterion

## 0. Sign fix to THM-1635

THM-1635 wrote `G(t) = −t·(log Π)′`. **The sign is wrong; it is `+`.** From
`R(u_i)/S(u_i) = −t·d(log u_i)/dt` and `G = −Σ_i R(u_i)/S(u_i)`,

```
G(t) = +t · d/dt log Π(t).
```

Verified on `R = 1+u`, `N=1`, where `Π = t/(1−t)`: `t(log Π)′ = 1/(1−t) = 1+t+t²+⋯ = Σ CT tᵐ`.
The **criterion `Π = ct` is unaffected** (both signs give `Π = ct ⟺ G` constant).

## 1. The charge picture, exact

`Λ = u^{−N}R = Σ_{k=0}^{M+N} r_k u^{k−N}` has **charge set** `𝒞 = {k−N : r_k ≠ 0} ⊆
[−N, M]`, containing `−N` (from `r_0 ≠ 0`) and `M` (from `r_{M+N} ≠ 0`). Then

```
CT(Λ^m) = [u^0] Λ^m = Σ_{ (c_1,…,c_m) ∈ 𝒞^m,  Σ c_i = 0 }  ∏ r_{·} / (symmetry) .
```

**First condition.** `CT(Λ^1) = [u^0]Λ = r_N`, so a nullcone element needs `r_N = 0` — `Λ`
has no constant term. This is exactly the THM-1535 charge-0 condition, recovered from the
`m=1` moment.

## 2. The binomial theorem (PROVED)

> **Theorem.** For `R = r_0 + r_d u^d` (`d = M+N`, `r_0, r_d ≠ 0`) and any `N ≥ 1`, `Λ =
> u^{−N}R` is **not** a nullcone element. Hence **TNC holds for every binomial `R`, at every
> bidegree.**

*Proof.* The charge set is `{−N, M}`. A representation of `0` is `a·(−N) + b·M = 0` with
`a+b = m`, i.e. `Na = Mb`. With `g = gcd(N,M)` the minimal solution is `a = M/g`, `b = N/g`,
`m_0 = (M+N)/g`. It is **unique** (two charges, one linear relation). Therefore

```
CT(Λ^{m_0}) = \binom{m_0}{a} r_0^{a} r_d^{b} ≠ 0
```

— a single product of nonzero factors. So some moment is nonzero. ∎

**Verified 16/16** for `N, M ∈ {1,…,4}`, predicted `CT(m_0)` matching exactly (e.g. `(N,M) =
(2,3)`: `m_0 = 5`, `CT = 720`; `(3,4)`: `m_0 = 7`, `CT = 15120`), with all earlier moments
zero.

## 3. The general sufficient condition (PROVED)

> **Corollary.** If `0` has a **unique** minimal representation as a sum of charges from `𝒞`,
> then `CT(Λ^{m_0})` is a single nonzero product and TNC holds — **with no all-order
> argument.**

This closes, besides all binomials, every `R` whose charge geometry gives a unique minimal
`0`-representation. Examples verified: charges `{−2,1,2}`, `{−2,−1,2}`, `{−2,1,3}`,
`{−3,1,2}`, `{−3,−1,2}` — all forced-TNC at `m_0 ∈ {2,3}`.

The mechanism is "Newton-polytope vertices cannot cancel" (THM-1625 §5) made exact: a unique
minimal representation *is* a vertex of the `m_0`-dilate with no competitor.

## 4. The residual, with an explicit witness

The only way to evade §3 is a **non-unique minimal representation**, where the competing
representations can cancel under tuned coefficients. The first trinomial where non-uniqueness
occurs is charges `{−2, 1, 4}`:

```
R = r_0 + r_3 u^3 + r_6 u^6 ,  N = 2 ,   m_0 = 3 ,
reps of 0:  {−2,−2,4} and {1,1,−2}
CT(Λ^3) = 3 r_0 (r_0 r_6 + r_3^2)         (two terms — CAN cancel)
```

**Tuning** `r_0 = 1, r_3 = 1, r_6 = −1` gives `r_3^2 = −r_0 r_6`, so `CT(Λ^3) = 0`: the
leading obstruction cancels. **But it is not a nullcone element:**

```
R = 1 + u^3 − u^6 ,  N = 2 :   CT[1..11] = 0,0,0,0,0,−30,0,0,126,0,0
```

`CT(Λ^6) = −30 ≠ 0`. **Vertex cancellation starts at `m_0 = 3` and the next order (`m = 6`)
finishes the job.** This is exactly the THM-1635 phenomenon: the dominant saddles of
`1 + u^3 − u^6` are nondegenerate (`|g″| = 8.873 ≠ 0`), so the generating function is
genuinely singular and some `CT ≠ 0`. TNC holds for this `R` too — but the proof needs the
branch-product/singularity argument, not the unique-minimal shortcut.

## 5. Status of the TNC (updated)

| case | status | by |
|---|---|---|
| `M = 0` (extreme weight) | PROVED | THM-1530 B |
| `min(M,N)=1`, `(2,2)`, `(2,3)` | PROVED | THM-1595 |
| **binomial `R`, any bidegree** | **PROVED** | **§2** |
| **unique minimal `0`-rep, any bidegree** | **PROVED** | **§3** |
| non-unique minimal rep | ⟺ branch product `Π ≠ ct` (THM-1635); nondegenerate-saddle cases hold (§4) | open in general |

## 6. Next

1. **Prove: non-unique minimal representations still cannot cancel to all orders.** §4's
   witness shows leading cancellation is real but the tail survives; the branch-product
   `Π = ct ⟹ R` monomial (HYP-8470) is the all-order form. The sufficient condition §3 has
   now peeled off everything *except* the genuine cancellation locus, which is a measure-zero
   coefficient variety — the residual is coefficients tuned onto `{CT(m) = 0 ∀m}`, an
   infinite intersection that §4 suggests is empty for non-monomial `R`.
2. **Second-order at a cancelled vertex.** When the unique-minimal shortcut fails, the *next*
   representation level gives a new CT; is it ever simultaneously cancellable? §4 says no in
   the first case (`m=6` survives) — a general "two consecutive obstruction levels cannot
   both vanish" would finish TNC.

## Verification

`04-computation/tnc_binomial_theorem_opus_S419.py` (16/16 bidegrees, exact `CT(m_0)`),
`tnc_unique_min_rep_opus_S419.py` (the sufficient condition; the first non-unique trinomials),
`tnc_vertex_cancellation_opus_S419.py` (the `1+u^3−u^6` witness: `CT(3)=0`, `CT(6)=−30`,
nondegenerate saddles), `tnc_large_branch_product_opus_S419.py` (the sign fix). Outputs in
`05-knowledge/results/`.
