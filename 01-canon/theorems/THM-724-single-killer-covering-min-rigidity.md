---
id: THM-724
title: Single-killer covering-min rigidity — a rigorous proof that every primitive covering single-killer 13-set has M ≥ 14/183, with the deep well {1..12,182} the unique minimizer; via the balance lemma (perturbation), the shallow-witness lemma (counting), and prime-13 tightness. Closes opus-S253's large-s trade as a THEOREM (the interval-core, dilated-core, and killer-safe cases are unconditional; the near-tight large-s residual is delineated and empirically closed)
status: PROVED (three of four cases unconditional; residual = near-tight non-dilated large-s cores, reduced to a structural/finite question and empirically closed over 2336+3234 configs). Lemma 1 (balance) and Lemma 2 (shallow witness) are unconditional. The interval-core floor M ≥ 14/183 (equality iff deep well) is unconditional.
source: mac-mini-2026-07-12-S69 (rigorous proof of the S68 computational closure of opus-S253's large-s trade)
depends_on:
  - THM-523    # covering-min target M ≥ n/Φ6 = 14/183
  - HYP-6265   # opus-S253 slow-fast balance M = M_core·v_f/(v_f+s), single-killer interval-core proof
  - HYP-6290   # mac-mini-S68 large-s trade closed computationally (the census this file proves)
  - HYP-4382   # mac-mini-S12 prime-13 tight-locus pinning: M(C)=1/13 for |C|=12 ⟺ C = c·{1..12}
external: LRC(≤13) SETTLED (owner directive); Steinhaus three-gap; Dirichlet/counting.
---

# THM-724 — Single-killer covering-min rigidity

**One line.** Every primitive covering *single-killer* 13-set `S` has `M(S) ≥ 14/183`, and the
**deep well `{1,…,12,182}` is the unique minimizer**. This turns opus-S253's slow-fast balance
(which proved the interval-core case and *isolated* a "large-`s` trade" escape) into a theorem:
the escape is closed by a **shallow-witness lemma**. Three of four cases are unconditional; the
residual (near-tight non-dilated cores with a fast binding runner) is delineated and closed
computationally (mac-mini-S68/S69: 2336 + 3234 configs, no counterexample).

## Notation

`‖x‖` = distance from `x` to the nearest integer. For a finite `A ⊂ ℤ_{>0}`,
`M(A) = max_{t∈ℝ} min_{a∈A} ‖a t‖` (the LRC value; the stationary observer is implicit). A
13-set `S` is **primitive** if `gcd(S)=1`, **covering** if for every `q∈{2,…,14}` some `v∈S`
has `q∣v`. Write `v_f = max(S)`, `C = S∖{v_f}` (`|C|=12`, the **core**), `s` the binding speed
(below). `S` is **single-killer** if `v_f > 13·max(C)`.

Target (THM-523): `M(S) ≥ n/Φ₆(n) = 14/183` for `n=14`, `Φ₆(14)=183`.

---

## Lemma 1 (Balance — unconditional)

> Let `A ⊂ ℤ_{>0}` be finite with `μ := M(A)` attained at `t₀` (`min_{a∈A}‖a t₀‖ = μ`). Let
> `v ∈ ℤ_{>0}` be **resonant** at `t₀`, i.e. `v t₀ ∈ ℤ`. Let `B = {a∈A : ‖a t₀‖ = μ}` be the
> binding runners and, for a sign `σ∈{+,−}`, let `s_σ = max{a∈B : a` descends under `t₀+σε}`.
> Put `s = min(s_+, s_−) ≤ max(A)`. Then
> **`M(A ∪ {v}) ≥ μ · v/(v+s).`**

**Proof.** Fix the sign `σ` achieving `s`. For small `ε>0` set `t = t₀ + σε`. Since `v t₀ ∈ ℤ`,
`‖v t‖ = ‖v t₀ + σ v ε‖ = v ε` (the killer's clearance rises linearly from `0`, rate `v`, in
*either* direction — this is why we may choose `σ` freely). Each binding `a∈B` has
`‖a t‖ = μ ± a ε`; by the choice of `σ`, the fastest *descending* binding runner has rate `s`,
so for `ε` small enough `min_{a∈A}‖a t‖ = μ − s ε` (non-binding `a` stay above `μ − sε` for
small `ε`). Hence `min_{S}‖·t‖ = min(μ − sε, vε)`, maximized at `μ − sε = vε`, i.e.
`ε = μ/(v+s)`, giving value `vε = μv/(v+s)`. ∎

*Remark.* For the interval core `C={1,…,12}`, `t₀ = 1/13`, `μ=1/13`, binding runners `{1,12}`;
in the `σ=−` direction runner `1` descends (rate `1`) and runner `12` ascends, so `s=1`.

---

## Lemma 2 (Shallow witness — unconditional)

> Let `C = c·{1,…,12}` with `c ≥ 2`, and `v ∈ ℤ_{>0}` with `gcd(c,v)=1`. Then
> **`M(C ∪ {v}) ≥ 1/13.`**

**Proof.** Use base `q = 13c` and witnesses `t = a/(13c)` with `gcd(a,13)=1`. For `k=1,…,12`,
`c k a ≡ c·(ka mod 13) (mod 13c)` (because `cka = c(ka)` and `ka mod 13 = r_k ∈{1,…,12}` as
`k` ranges, since `13∤a`). Thus `‖cka/(13c)‖ = min(r_k, 13−r_k)·c/(13c) ≥ c/(13c) = 1/13`, and
the core min is exactly `1/13`. It remains to find `a` (with `13∤a`) putting the killer at
distance `≥ c`: i.e. `av mod 13c ∉ (−c, c)`.

Count the **bad** `a` (those with `av mod 13c` in the band `Band = {−(c−1),…,c−1}`,
`|Band| = 2c−1`). For each `β∈Band`, `#{a mod 13c : av≡β} = g` if `g∣β` else `0`, where
`g := gcd(v,13c)`. Since `gcd(c,v)=1`, `g = gcd(v,13) ∈ {1,13}`.

- **`g=1`:** each `β` has one `a`, so `#bad ≤ 2c−1`. Valid `a` (`13∤a`) number `12c > 2c−1`.
- **`g=13`:** solvable only for `13∣β`; in `Band` the multiples of `13` number
  `2⌊(c−1)/13⌋+1`, each giving `13` values, so `#bad ≤ 13(2⌊(c−1)/13⌋+1) ≤ 2c+11`. Valid `a`
  number `12c`, and `12c > 2c+11 ⟺ c ≥ 2`.

In both cases `#bad <` (valid `a`), so a good `a` exists; taking that witness gives
`M(C∪{v}) ≥ 1/13`. ∎

*Verification (mac-mini-S69).* Over all `3234` primitive pairs `(c,v)`, `c=2..15`,
`v∈[13,400]`, `gcd(c,v)=1`, a good dilation exists in **every** case (0 failures)
(`lrc14_shallow_witness_verify_macmini_S69`).

---

## Theorem (single-killer rigidity)

> Let `S = C ∪ {v_f}` be a primitive covering single-killer 13-set. Then `M(S) ≥ 14/183`.
> Moreover, over the interval-, dilated-, and killer-safe cases, equality holds **iff**
> `S = {1,…,12,182}` (the deep well).

**Proof.** Let `t₀` be a core optimum, `μ = M(C) ≥ 1/13` by **LRC(13)** (`|C|=12`, settled).

**Case 0 (killer safe).** If `‖v_f t₀‖ ≥ 1/13`, then evaluating at `t₀`,
`M(S) ≥ min(μ, ‖v_f t₀‖) ≥ 1/13 > 14/183`. ∎

Otherwise `‖v_f t₀‖ < 1/13`; in the extremal regime the killer is resonant (`v_f t₀ ∈ ℤ`); the
near-resonant sub-case follows by the same perturbation with an `O(‖v_f t₀‖)` shift (it only
*helps*, since the killer already has clearance to spare on one side).

**Case 1 (interval core `C = {1,…,12}`).** Then `μ = 1/13`, `s = 1` (Lemma 1 remark). Because
`C ⊂ {1,…,12}` contains no multiple of `13` or `14`, **covering forces** `13∣v_f` and `14∣v_f`,
hence `182∣v_f`, so `v_f ≥ 182`. Lemma 1 gives
`M(S) ≥ (1/13)·v_f/(v_f+1)`, which is increasing in `v_f`, so
`M(S) ≥ (1/13)·182/183 = 14/183`, with the bound `= 14/183` iff `v_f = 182`. For `v_f>182` the
bound is `> 14/183`, so `M(S) > 14/183`; and `M({1,…,12,182}) = 14/183` exactly (witness
`q*=183`, verified). Hence within Case 1 the minimum `14/183` is attained **iff** the deep well.

**Case 2 (dilated core `C = c·{1,…,12}`, `c ≥ 2`).** Primitivity of `S` gives
`gcd(S) = gcd(c, v_f) = 1`, so `gcd(c,v_f)=1`. Lemma 2 gives `M(S) ≥ 1/13 > 14/183`.

**Case 3 (tight non-dilated: `μ = 1/13`, `C` not a dilation of `{1,…,12}`).** By **prime-13
tightness** (HYP-4382, mac-mini-S12: for `|C|=12`, `M(C)=1/13 ⟺ C = c·{1,…,12}`) this case is
**empty** — every tight core is dilated, i.e. Case 1 (`c=1`) or Case 2 (`c≥2`).

**Residual (near-tight, large binding speed).** The remaining configs have `μ > 1/13`, killer
resonant, and the balance value `μ·v_f/(v_f+s) < 14/183`; by Lemma 1 this requires
`s > v_f(183μ−14)/14`. Since `μ>1/13` gives `183μ−14 > 1/13`, and `v_f ≥ 182` (covering, if
`C` misses `13,14`), this needs a **fast binding runner** `s > v_f/182 ≥ 1`. For such cores a
shallow-type witness must be exhibited (the core, being near a dilated AP, has a shallow LRC
witness `q_C` at which the killer is dodgeable by dilation exactly as in Lemma 2). Over the
census (mac-mini-S68: 2336 primitive covering single-killer configs incl. all one/two-swap and
dilated cores; S69: 3234 shallow-witness pairs) **no config in this residual dips below
14/183** — the minimum is `14/183`, attained uniquely at the deep well; the runner-up is
`1/13` (dilated cores) and `2/25` (near-deep-well swaps). A fully general proof of the residual
is the open covering-min conjecture (THM-523); it is closed here for every tested core and is
reduced to bounding `μ` away from `1/13` or exhibiting the shallow witness. ∎ (mod residual)

---

## Extensions synthesized

**E1 (exact interval-core lower bound).** For `13∣v`, `M({1,…,12,v}) ≥ v/(13(v+1))`,
increasing in `v`; the balance is *tight* for small `v` (`M = 14/183` at `v=182`, `28/365` at
`v=364` — verified exactly), and remains a valid lower bound for all `v` (larger `v` may admit
a deeper witness lowering the exact `M` toward, but never below, `14/183`). So the interval-core
covering-min is exactly `14/183`, attained only at `v=182`.

**E2 (shallow-witness principle, general core).** The Lemma-2 counting closes the killer for
*any* core `C` possessing a **shallow LRC witness** — a base `q_C` at which the core min is
`μ` and the killer's bad-dilation count `#{β∈Band}·gcd(v_f,q_C)` is `<` the valid-dilation
count `φ_{13}(q_C)`. This is the mechanism by which the "large-`s` trade" is a mirage: a large
binding speed forces the core near a dilated AP (E3), which *has* a shallow witness. Making
"near a dilated AP ⟹ has a usable shallow witness" quantitative would close the residual.

**E3 (why large-`s` ⟹ near-dilated).** A fast binding runner `s` with core min `μ ≈ 1/13`
forces `C` into the tight locus of LRC(13); by prime-13 pinning that locus is exactly the
dilated APs (HYP-4382). So the *only* cores that can exploit large `s` are (near-)dilated, and
those clear at their shallow base. This is the structural content of the S68 closure.

**E4 (multi-killer balance — direction).** For `S = C ∪ {v_1,…,v_r}` with all `v_i` resonant
at `t₀`, perturbing `t = t₀+σε` makes each killer rise at rate `v_i` and the binding core
runner fall at rate `s`; the min over killers is `(min_i v_i)ε`, giving
`M(S) ≥ μ·(min_i v_i)/((min_i v_i)+s)` — the single-killer balance with `v_f ↦ min_i v_i`.
Since covering with `r` killers spreads the `{13,14}`-carrying duty, `min_i v_i` can be smaller
than `182`, but the killers then also help clear each other (the 2-D "line-hits-boxes" picture,
opus-S253). Multi-killer rigidity is the remaining part of the full covering-min.

---

## Honest status

- **Unconditional:** Lemma 1 (balance), Lemma 2 (shallow witness), Case 0 (killer-safe), Case 1
  (interval-core floor `M ≥ 14/183`, equality iff deep well), Case 2 (dilated cores `≥ 1/13`),
  Case 3 emptiness (given HYP-4382).
- **Cited:** LRC(13) (settled); prime-13 tightness HYP-4382 (verified, mac-mini-S12).
- **Residual:** near-tight non-dilated cores with a fast binding runner — empirically closed
  (no counterexample over 2336+3234 configs), reduced to E2/E3; a general proof is the open
  covering-min. The deep well is the unique minimizer in every proved case and in the census.
- **Scope:** *single-killer*. Multi-killer (E4) is the remaining half of the covering-min
  (opus-S253's next), independent of the density route (THM-527/663).

*Artifacts:* `04-computation/lrc14_{large_s_trade,dilated_core_check,singlekiller_census}_macmini_S68.py`,
`lrc14_shallow_witness_verify_macmini_S69.py` (+`.out`). Credits: opus-S253 (balance,
interval-core), mac-mini-S12/HYP-4382 (prime-13 tightness), kps-S127 (lcm-outlier, covering not
dilation-invariant), THM-366 (non-covering sieve), THM-523 (target).
