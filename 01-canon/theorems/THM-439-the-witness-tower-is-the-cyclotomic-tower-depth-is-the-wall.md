# THM-439 — The LRC witness tower is the cyclotomic (abelian) tower automatically; its substance is the cyclotomic DEPTH q*(S), and the Abel–Ruffini analog is the Bonferroni tower, not the field tower

**Status:** PROVED (rationality ⟹ witness tower ⊆ ℚ^ab; the depth definition) + VERIFIED (depth
stratification + the per-config-finite / uniformly-unbounded dichotomy, n=5,6,7,8 in-window) +
HONEST DEFLATION of part of HYP-2303. Sharpens HYP-2303/THM-436 into a precise statement.
**Source:** opus-2026-06-07-S704 (developing the radical tower ↔ witness tower, user request). Builds
on THM-420/430 (witness hierarchy), THM-403 (cyclotomic worry-set), THM-406 (Vitali wall / no finite
Bonferroni), THM-436/HYP-2303 (the Galois solvability lens), S700 (residual R(n)), Kronecker–Weber.

## (A) The witness tower is cyclotomic — automatically (PROVED)

> For integer speeds, `M(S)=max_t min_i ‖v_i t‖` is the maximum of a continuous **piecewise-linear**
> function whose breakpoints (`j/v_i`, `j/(2v_i)`) and tent-crossings (`j/(v_i±v_j)`) are all
> **rational**. Hence the optimal witness `t*` is **rational**, every runner position `v_i t*` is a
> **root of unity** `ζ_q^{v_i m}` (`t*=m/q`), and the entire witness apparatus lies in
> `ℚ(ζ_q) ⊆ ℚ^{ab}`, with `Gal(ℚ(ζ_q)/ℚ)=(ℤ/q)^×` **abelian**.

> **Consequence (honest deflation of HYP-2303).** There is **no non-abelian witness obstruction**: the
> field/Galois "solvable tower" is *trivially* solvable because `t*∈ℚ` always. By Kronecker–Weber the
> witness tower is *exactly* the abelian (= cyclotomic) part — but that reach is automatic, not a
> theorem about LRC difficulty. The "unsolvable-quintic" analogy is therefore **not** at the
> field-of-the-witness level (always abelian); it lives one level up (D).

## (B) The substance: the cyclotomic DEPTH q*(S) (VERIFIED)

> Define the **cyclotomic depth** `q*(S) = min{ q≥2 : max_m min_i ‖v_i m/q‖ ≥ 1/n }` — the smallest
> cyclotomic modulus whose ticks already clear the floor `1/n`. The witness hierarchy (THM-420/430) is
> the **magnitude stratification** of `q*`:
> `clock (q*≤n) ⊂ sub-shell (n<q*<2n−1) ⊂ shell (q*=2n−1) ⊂ super-shell (q*>2n−1)`.
> Verified depth distribution (gcd-1 configs):
> ```
>   n=5 (B≤16): clock 1185 | sub-shell 338 | shell 159 | super-shell 63    max q*=17
>   n=6 (B≤13): clock  888 | sub-shell 280 | shell  78 | super-shell 35    max q*=17
>   n=7 (B≤11): clock  401 | sub-shell  43 | shell  18 |  —                max q*=13
>   n=8 (B≤10): clock   90 | sub-shell  30 |  —        |  —                max q*=14
> ```
> **The S700 residual `R(n)` (divisibility-complete ∧ shell-free) NEVER lands at the clock level** —
> it is entirely sub-shell/shell/super-shell (n=5: 138/46/47; n=6: 54/8/14). So `R(n)` is exactly the
> *positive-depth* core: the configs the cyclotomic floor (`q*≤n`) cannot certify.

## (C) Per-config finite, uniformly unbounded — the dichotomy (VERIFIED)

> `q*(S) < ∞` for **every** config (n≤8 in-window): each config is individually certified loose by a
> finite cyclotomic tick — **LRC holds here, constructively, inside ℚ^{ab}**. But `sup_S q*(S)` is
> **unbounded**, growing with the speed bound:
> ```
>   n=6: max q* = 10,11,11,11,12,14,17,17   for B = 6..13
>   n=7: max q* = 11,13,13,13,13,14,19,19,21 for B = 7..15
> ```
> So there is **no uniform finite cyclotomic certificate** — no single depth certifies all configs.

## (D) The true Abel–Ruffini analog is the Bonferroni tower (the reframe)

> LRC(n) `⟺` `q*(S)<∞` for all `S` (every config cyclotomically/abelian-solvable). Since `t*` is
> always rational (A), the *per-config* answer is always "solvable" exactly when `M(S)≥1/n`; a
> counterexample (`M<1/n`) would be a config where **no** rational tick clears `1/n` — a
> **combinatorial/measure** failure, not a non-abelian *field*. Therefore:

> **The Abel–Ruffini "no finite tower" mirror (THM-436(4)) is NOT the field tower (automatically
> abelian) but the BONFERRONI tower (THM-406).** The wall is the **non-uniformity** of (C): the
> cyclotomic depth `q*(S)` is finite per config but unbounded over all configs, exactly as the
> covering-depth inclusion–exclusion terminates for any finite truncation yet cancels to all orders
> with no uniform Bonferroni certificate. *Each quintic has a finite splitting field; there is no
> uniform radical formula* `↔` *each config has a finite `q*`; there is no uniform cyclotomic depth.*
> **The Vitali wall is the unbounded depth, not any single config.**

## Proofs / verification

**(A)** Standard: max of integer-breakpoint PL functions is at a rational; `e^{2πi m/q}∈μ_q⊂ℚ^{ab}`;
Kronecker–Weber. ∎
**(B)(C)** Exact rational computation, `04-computation/lrc_cyclotomic_witness_tower_s704.py` (+`.out`):
`M(S)` exact, `q*(S)` exact, magnitude buckets, residual cross-classified, `max q*` vs `B`.
**(D)** Reframe combining (A)(C) with THM-406 (imported); an identified correspondence.

## Scope / honesty

- (A) is rigorous and **deflates** the HYP-2303 reading "tower-certifiable ⟺ solvable monodromy / a
  counterexample is non-abelian": the witness field is *unconditionally* abelian, so that reading was
  vacuous. The corrected substance is the **depth** `q*` (B,C) and the **non-uniformity = wall** (D).
- (B)(C) are verified in-window only; `max q*` growth is observed, not proven unbounded in general
  (though THM-406 says no uniform certificate exists). Resolves no open LRC case.
- Last session's loose "depth of tower = derived length of local monodromy" is **corrected**: the
  tower is uniformly **abelian** (derived length 1, cyclotomic); the meaningful depth is the modulus
  `q*` (finite per config, unbounded uniformly), and the Abel–Ruffini analog is the Bonferroni order.

**Artifacts:** `04-computation/lrc_cyclotomic_witness_tower_s704.py` (+`.out`). Reflection
`07-reflections/the-cyclotomic-witness-tower-and-the-depth-wall-s704.md`. New: **HYP-2309**. Builds on
THM-420/430, THM-403, THM-406, THM-436/HYP-2303, S700, Kronecker–Weber.
