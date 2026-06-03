---
id: HYP-2101
status: PARTIAL — sheaf reframing of S559 PROVED/verified; apex = unique whole-line
  section (verified q=3,5,7,11); the r/p lift dissolves the apex obstruction for the
  tight tuple (forbidden-fraction 1 -> 0, verified); full ratio-spread residual under
  the lift still OPEN
source: claudebox-2026-06-02-S579
related:
  - HYP-2063
  - HYP-2024
  - HYP-2020
  - THM-380
  - THM-396
---

# HYP-2101: the apex-lift certificate sheaf

A geometric/sheaf reorganization of the `k+1 = 2q` polynomial-method corrector
(S559 / HYP-2063, the literature's n=14 wall) that makes the apex obstruction a
statement about the support of a single local section, and makes S559's open lead
(a) — the `r/p` time freedom — the operation that resolves it.

## Setup (n = 2q, runners i = 1..2q-1)

S559's **Reduction Theorem** (PROVED there, re-verified here: 0 mismatches over
1200 random tuples, q=3,5,7) says a strict-interior corrector exists **iff**
`∃ s',r' ∈ 𝔽_q^× : s'w_i + r'c_i ≠ f_i ∀i`, with `w_i=v_i mod q`, `c_i=i mod q`,
`f_i=0` if `v_i+i` even else `q-1`.

**Certificate line.** Runner `i` forbids the line `L_i = {s w_i + r c_i = f_i}`
in `𝔸²(𝔽_q)`. A *loneliness certificate* is a point `(s,r) ∈ T := (𝔽_q^×)²` off
every `L_i`. Define the **certificate sheaf** `𝒞_v` = extension-by-zero of the
constant sheaf on the certificate locus `T \ ⋃_i L_i`; then `H^0(𝒞_v)` = the set
of correctors and **a corrector exists ⟺ H^0(𝒞_v) ≠ 0**.

**Parity-matched reduction to ℙ¹.** When all `f_i=0`, every `L_i` passes through
the origin, so it is the single projective point `ρ_i = -c_i/w_i ∈ ℙ¹(𝔽_q)` (its
slope). The certificate sheaf then lives on `ℙ¹`: non-apex runner `i` deletes the
one point `ρ_i`; a corrector exists iff the deleted points miss some unit slope.
This *is* S559's Residual Theorem, re-encoded: uncorrectable ⟺ `{ρ_i}` covers
`𝔽_q^×`.

## What is proved / verified this session

1. **The apex is the unique whole-line section.** The apex runner `i=q` has
   `c_q = q ≡ 0 (mod q)`, so its forbidden locus is `{s w_q = f_q}` — independent
   of `r`. When the apex speed is `≡0 mod q` (true for the tight tuple, `w_q=0`,
   `f_q=0`) the "line" degenerates to `0 = 0`: it forbids the **whole plane**.
   So `𝒞_v` has the zero stalk everywhere through this one section ⇒ `H^0=0`.
   Verified `apex_whole = True`, `|H^0| = 0` for the tight tuple at q=3,5,7,11;
   `|H^0|` jumps to `2,12,30,90` once the apex section is removed. This is the
   precise sheaf statement of HYP-2063's "apex IS the zero-divisor / field
   failure": the apex is the only runner whose certificate is non-transverse
   (support = whole space rather than codim 1).

2. **Slopes collapse for the tight tuple.** Apex excluded, every non-apex tight
   runner deletes the *same* slope `ρ_i = q-1 = -1` (verified q=3,5,7,11: the
   forbidden-slope set is the singleton `{q-1}`), so the locus is large and a
   certificate trivially exists — the geometric form of "tight = easiest".

3. **The lift dissolves the apex obstruction (tight tuple).** Adjoining the `r/p`
   freedom of `t = s/(2q) + r/p` (S559 open (a)) works over modulus `M = 2q·p`,
   `p` = least prime ∤ `2q`. The apex speed mod p is `q mod p ≠ 0` (a unit), so the
   apex line regains an `r`-slot and becomes codim 1 again. Concretely, under the
   tight corrector `r=s` the apex lands at `2q·(s mod p) mod 2q p` — a nonzero even
   multiple of `2q`, strictly interior. The apex **forbidden-fraction drops from
   `1.000` (whole space, unliftable) to exactly `0.000`** at q=3,5,7,11
   (M=30,30,42,66). The whole-line section becomes the zero section.

## Open / not claimed

- The lift computation clears the **apex alone** for the **tight tuple**. It does
  **not** yet clear the **ratio-spread residual** (S559 §4) under the lift — i.e.
  whether the lifted arrangement `⋃_i L_i^{lift}` ever covers `𝔸²(𝔽_q)×𝔽_p` for a
  parity-matched ratio-spreading `v`. That is the real open question and the next
  computation.
- The "sheaf" here is an honest constructible-sheaf reorganization (`H^0` =
  global sections = correctors); a genuine **`H^1` obstruction class** localized at
  the apex stratum is *conjectured* (the apex is the unique non-transversal point
  of the intersection poset) but **not proved**. Don't cite an `H^1` theorem.
- Transfer: the picture is uniform in `q`, so it covers the whole `2·prime`
  sub-frontier `n=10,14,22,…` exactly as S559 predicts.

## Next steps (emitted as tasks)

(a) compute the lifted line arrangement `L_i^{lift}` over `𝔸²(𝔽_q)×𝔽_p` for
ratio-spread parity-matched `v` and test whether the certificate locus is ever
empty (does the lift clear the residual, or only the apex?); (b) make the `H^1`
claim precise via the intersection poset / characteristic polynomial of the
arrangement and check whether its only non-transversal flat is the apex; (c) hand
the apex whole-line section to the pinch/shield route (THM-396) as the `(q,q)`
shield, per HYP-2063 open (b).

**Artifacts:** `04-computation/lrc_apex_lift_certificate_sheaf_s579.py` (+`.out`
in `05-knowledge/results/`). Builds on HYP-2063/S559, HYP-2024 (vertices as
section/certificate objects), THM-380 (pressure certificate), arXiv:2604.23906
Prop 4.1/4.4.
