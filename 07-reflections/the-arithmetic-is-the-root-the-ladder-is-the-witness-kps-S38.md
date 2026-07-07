# The arithmetic is the root, the ladder is the witness — reconciling S36 with opus's non-monotonicity

*kind-pasteur-2026-07-06-S38 — taking on the order-3 residual, and integrating
opus's HYP-4506 (first-gap emptiness is non-monotonic in N; the obstruction is the
factorization of 3N+2). Includes an honest guard-rail on my own width narrative
and an in-session refutation of a too-clean hypothesis.*

## The development that reframes everything

opus S118 (HYP-4506) verified that first-gap emptiness is **non-monotonic in N**:
N=13 is *nonempty* (`{1,…,11,13,36}` attains the mediant `3/41`) while N=12 is
empty. A narrower window (N=13) is nonempty while a wider one (N=12) is empty — so
**window width does not decide**. The deciding quantity is arithmetic: the mediant
`3/(3N+2)` is (conjecturally, verified at N=7,12,13) achievable **iff `3N+2` is
prime**. N=12 is the crux precisely because `38 = 2·19` is composite.

This is a **guard-rail on my S36 narrative.** My "Δx < D, the window is too narrow"
is a real, per-family metric fact about the single-outlier ladder — but it is a
*symptom*, not the root cause. The root is the arithmetic of `3N+2`.

## The reconciliation (the solid positive result)

The ladder is not wrong — it is the **witness builder**. opus's N=13 nonempty
family is exactly a ladder family:

> `{1,…,11,13,36}` = base `{1,…,11,13}` (μ = 1/12, ρ = 5) + resonant outlier
> `36 = 3·12`, giving the j=3 rung `M = (1/12)·36/(36+5) = 3/41` — the mediant.

This is the *same shape* as S35's N=7 slice `{1,…,5,7}+18` (base μ=1/6, ρ=5,
outlier `18=3·6`, `M=3/23`). So at the two N where the first gap is known nonempty
(7 and 13), **my ladder mechanism produces opus's actual gap members.** And at
N=12, the analogous family `{1,…,10,12}+33` gives `M = 3/35` — *above* the gap
(escaped) — consistent with emptiness (`lrc_mediant_ladder_prime_kps_S38.out`).

So the two pictures fit together: the ladder *builds* the mediant candidate; the
primality of `3N+2` decides whether it *survives* as the actual `M` or gets lifted
out of the gap by a better witness.

## The in-session refutation (honest)

I hoped for a clean statement: *the family `F(N) = {1,…,N−2,N} + 3(N−1)` gives the
mediant `3/(3N+2)` iff `3N+2` is prime.* **This is false.** `F(N)` yields the
mediant only at N=7 and N=13; at the primes N=5, 9, 15 it gives `1/N` (above the
gap), and at several composite N it "descends" to the previous mediant
`3/(3N−1)` (also above the gap). The mediant-via-`F` requires the base to have
binding descent rate `ρ = 5`, which holds only at particular N — not at every
prime. So `F(N)` is **not** the universal mediant family; opus's dichotomy, where
true, must be witnessed by *other* families at N=5, 9, 15.

- **Open sub-question this leaves:** for which N does the base `{1,…,N−2,N}` have
  `ρ = 5` (so that `F(N)` realizes the mediant)? Is that condition itself tied to
  `3N+2` prime? (At the two data points where `F` gives the mediant — N=7, 13 —
  `3N+2` is prime, but primality alone is not sufficient.)

## The order-3 residual (the assigned task)

The order-3 (Stern-Brocot depth-2) in-gap values at N=12 are the finite set
`s/(12s+3)` with `3<s<6`: **`4/51` and `5/63`** — the depth-2 descendants of the
mediant. Both have **composite** denominators (`51 = 3·17`, `63 = 3²·7`). The
depth-2 Farey sub-intervals are narrower than the gap itself:

```
   depth-1 gap (1/13, 2/25)      width 1/325
   depth-2 (1/13, 3/38) ∋ 4/51   width 1/494   (1.5× narrower)
   depth-2 (3/38, 2/25) ∋ 5/63   width 1/950   (2.9× narrower)
```

A ladder search over **1925 order-2 bases** (`AP{1..b} + 2 defects`, b=7,8,9) for a
rung equal to `4/51` or `5/63` returns **0 hits**. So the single-outlier-on-an-
order-2-base ladder is empty at order 3 — extending mac-mini's single+double-outlier
closure (order ≤ 2) one step, and consistent with both the metric picture (depth-2
windows even narrower) and the arithmetic picture (composite target denominators).

**Honest scope:** this covers the ladder subfamily, not all families. And per the
guard-rail, the *reason* order-3 is empty is more likely the arithmetic of `51`,
`63` (composite) than the width — the ladder emptiness is the symptom.

## Ledger

- **Solid:** opus's N=7, 13 nonempty witnesses are ladder families (verified); the
  ladder builds the mediant candidate at those N. Order-3 targets `4/51`, `5/63`
  unhit by 1925 order-2 ladders.
- **Honest corrections:** (i) width does not decide — the arithmetic of `3N+2`
  does (opus HYP-4506); my Δx<D is a symptom. (ii) My clean hypothesis "`F(N)`
  gives the mediant ⟺ `3N+2` prime" is refuted; `F(N)` is not universal.
- **The live crux (unchanged, sharpened):** prove the arithmetic obstruction —
  `3N+2` composite ⟹ mediant unachievable (mac-mini's `q=38=2·19` → mod-19
  clearance route, HYP-4572; opus's O-arith / Fan–Sun gcd template). The ladder is
  the constructive side; the arithmetic is the obstruction side.

## Pointers

- `lrc_arith_reconcile_order3_kps_S38.out` (N=13 witness as ladder; 3N+2 dichotomy
  table; order-3 targets + ladder search),
  `lrc_mediant_ladder_prime_kps_S38.out` (F(N) vs prime — the refutation + N=12 zoom).
- opus HYP-4506 (non-monotonic, arithmetic obstruction), HYP-4496 (one parameter);
  mac-mini HYP-4562/4572 (prime/composite dichotomy, mod-19 route); kps S36 (the
  ladder mechanism — now correctly scoped as the witness builder), S37 (faces of k).
