---
source: opus-2026-07-11-S230
status: A CORRECTION + a redirect of the last LRC(14) rigorous gap. The "pair-sum shallow lemma" (kps's
  remaining gap, flagged "LRC-hard") is FALSE unconditionally — it degrades with diameter. The correct
  clean-ruler instrument is a BOUNDED SMALL MODULUS (diameter-free), which always exists (verified to
  Vmax 5000). The dissociation/anti-concentration content, correctly applied, lands on bounded moduli, and
  the cleanest case (primes q∈{17,19,23}, band {0,±1}) has a PROVED exact criterion.
tags:
  - lrc14
  - clean-ruler
  - pair-sum
  - bounded-modulus
  - anti-concentration
  - dissociation
  - THM-707
  - correction
---

# The clean ruler is a bounded modulus, not a pair-sum

**opus-2026-07-11-S230.** Owner asked me to take on the pair-sum shallow lemma via the
dissociation/anti-concentration angle. Doing so honestly **refutes the lemma as stated** and **redirects the
instrument**. This is the last rigorous gap of LRC(14), so the redirect matters.

## The setup (kps THM-707)

`hB5` — the single non-foundational LRC(14) obligation — reduces to: every **residual** covering family (13
speeds) has a **clean ruler** `q` (`liveCount(q) ≥ 1` and `maxBand(v,q) ≤ 5`, where
`bandCount(v,q,p) = #{i : v_i·p mod q ∈ danger arc [0,q/14)∪(13q/14,q)}`). kps's chosen instrument was the
**pair-sum** `q = v_a + v_b`, and kps cont.32 flagged "some pair-sum has `maxBand ≤ 5`" as the remaining
rigorous gap — "LRC-hard, overwhelmingly evidenced."

Crucially, kps's own architecture only needs the clean ruler **after the THM-701 peel**: `hB5 ⟸ [peel far
elements to a bounded-spread core] + [clean ruler on the core]`.

## The finding 1 — pair-sum shallow is false unconditionally

`q = v_a + v_b ≈ 2·Vmax` scales with the diameter. The expected number of 6-clusters in the danger arc (an
arc of length `q/7`) grows `∝ Vmax` (a second-moment / Poisson count: `Σ_{|T|=6} #{p} ≈ C(13,6)·q/7^6 ∝ q`).
So `min`-over-pairsums `maxBand` **climbs with diameter and crosses 6**. Verified on prime-rich dissociated
residuals: pair-sum fails (`min`-over-pairsums `maxBand ≥ 6`) for **0, 1, 29 of 120** families as
`Vmax ≈ 60, 250, 900`.

These are **genuine dissociated residuals** (primitive, longest-AP ≤ 7, prime-rich), *not* the
coarse-reducible (dilated/detuned) families kps correctly excluded. The subtlety kps's framing missed:
**THM-701's peel bounds the *ratio* (no far element ⟺ `w < 90191·Σe'`), not `Vmax`.** A primitive
bounded-*ratio* core can still have large `Vmax` (e.g. `{N, …, N+12}`), and pair-sums fail there. So the
pair-sum instrument does **not** cover THM-701's actual base case.

## The finding 2 — the bounded modulus is the right instrument

Cleanness at `q` depends **only on `{v_i mod q}`**, so a **bounded** `q` is **diameter-free**. And a bounded
clean ruler always exists: every dissociated prime-rich residual has a clean ruler with **`q ≤ 37`**, and the
bound **does not grow** with diameter (max smallest-clean-`q` = 37 over `Vmax` up to **5000**, 0 failures). A
single fixed 8-element set `{8,17,19,23,27,29,35,37}` is a clean ruler for **100%** of a 2000-family diverse
pool (best single modulus `q=37` covers 69%).

So the clean ruler on THM-701's bounded-spread (possibly large-`Vmax`) core **exists as a bounded modulus** —
the pair-sum was simply the wrong family. This closes the instrument gap: use `q ∈ {8,…,40}`, not `v_a+v_b`.

## The finding 3 — the anti-concentration content, made exact (PROVED)

Where does dissociation enter? At the cleanest moduli. For **prime `q ∈ {17,19,23}`** the danger set is
**exactly `{0,1,q−1} = {0,±1}`** (true for every `q ∈ [15,28]`). Writing `z = #{i : v_i ≡ 0 (q)}` and the
**±-folded multiplicities** `c_r = #{i : v_i ≡ ±r (q)}` for `r = 1,…,(q−1)/2`:

> **Proposition (small-prime clean-ruler criterion — PROVED, elementary).**
> `maxBand(v,q) = z + max_r c_r` and `liveCount(v,q) = 2·#{r : c_r = 0}` if `z=0` (else `0`).
> Hence `q ∈ {17,19,23}` is a **clean ruler** for `v` ⟺ **`z = 0` ∧ `max_r c_r ≤ 5` ∧ some `c_r = 0`.**

*Proof.* `q` prime ⟹ `x ↦ px` bijects `ℤ/q` fixing `0`, so `#{v_i p ≡ 0} = z` for all `p`, and
`v_i p ≡ ±1 ⟺ v_i ≡ ±p^{-1}`. Thus `bandCount(v,q,p) = z + c_{⟨p^{-1}⟩}` (fold to `{1,…,(q−1)/2}`). As `p`
runs over `{1,…,q−1}`, `p^{-1}` runs over all nonzero residues, so `⟨p^{-1}⟩` hits each fold-class exactly
twice. Therefore `{bandCount(p)} = {z + c_r}` (each value twice); read off `maxBand` and `liveCount`. ∎
(Verified: 0 mismatches / 90000.)

**`max_r c_r ≤ 5` = no antipodal ±-pair mod `q` holds ≥ 6 of the 13 speeds** — an **anti-concentration
(Sidon-type)** condition. A dissociated family satisfies it for *most* small primes; the dilated AP `7·{1..13}`
(all speeds ≡ 0 mod 7, hence hugely concentrated mod many `q`) is exactly what violates it — and is peeled by
the dilation/detuned branch, not a residual. So dissociation is **essential and non-modular** (it is a
ℤ-condition; the residues mod `q` of an arbitrary set can be anything — the dilated AP proves the lemma is
false without dissociation).

## The redirected last gap (cleanly posed, diameter-free)

The rigorous finish of hB5's hard core is no longer "pair-sum shallow" but:

> **Bounded-modulus clean-ruler lemma.** Every dissociated prime-rich primitive 13-family has a clean ruler
> `q ∈ {8,…,40}`. Equivalently (sufficient): some prime `q ∈ {17,19,23,29,31,37}` has no antipodal ±-class
> holding ≥ 6 speeds and an empty ±-class.

This is an **inverse-additive / anti-concentration** statement — the honest home of the "dissociation" the
whole clean-ruler route was implicitly leaning on. It is diameter-free (unlike pair-sums), exact-criterion
backed for the `{0,±1}`-band primes, and verified with a hard bound (`q ≤ 37`, `Vmax` to 5000). Proving it in
full is the remaining work; the `{17,19,23}` criterion is the formalizable rung.

## Net

- **Retract** the pair-sum shallow lemma as an unconditional statement — it is diameter-dependent and false
  for large `Vmax` even at bounded ratio (so it does not cover THM-701's base case).
- **Adopt** the bounded small modulus as the clean-ruler instrument — diameter-free, always present (verified),
  with a **proved** exact criterion for the `{0,±1}`-band primes `{17,19,23}`.
- The last LRC(14) gap is an **anti-concentration lemma on residues mod a bounded prime**, not an
  equidistribution bound on an unbounded pair-sum modulus.

→ THM-707 (kps, the clean-ruler reduction), THM-701 (the peel — bounds ratio not Vmax), THM-709/712 (the
prime ruler `q ≤ 13`, band = all nonzero — this extends it to the `{0,±1}` band at `q ∈ {17,19,23}`),
opus-S225 (the earlier "bypass the non-truncatable series" lesson — same spirit: the rigor is finite/exact,
here a bounded-modulus occupancy condition).

---

## Addendum (opus-S231) — following the redirect: adversarial robustness and the honest proof limit

Pushing the math further (owner: "follow the redirect as far as you can"), four refinements:

**1. The decomposition — shallow is the robust core, live/composites the rare hard part.** Over the prime
window `{17,19,23,29,31,37}` on 1200 random dissociated prime-rich residuals (Vmax 40–5000): SHALLOW (some
prime `maxBand ≤ 5`) **1200/1200**; LIVE (some prime has a live multiplier) **1198/1200**; CLEAN (both, same
prime) **1196/1200**; CLEAN in the full window `{8..40}` (composites) **1200/1200**. So on *random* families
the shallow half is essentially free at a prime; the last few need composites for the *live* half.

**2. The prime window is NOT adversarially sufficient — composites are essential.** Hill-climbing to maximize
heavy-count finds families **heavy at all 6 window primes** with **longest-AP = 2** (fully dissociated). So
longest-AP dissociation does *not* bound the prime-window shallow disjunction. The `{17,19,23}` criterion of
S230 is therefore a *partial* tool (handles most random families), **not** a standalone closer — composite
moduli (8, 16, 25, 27, 35, …) are genuinely required.

**3. The full bounded-modulus lemma IS adversarially robust and diameter-free.** Hill-climbing *within the
prime-rich core* to maximise the smallest clean modulus: the adversarial max smallest-clean-`q` is **stable at
47–59** across Vmax ceilings 300 → 64000 — it does **not** grow with Vmax (the earlier 37→49→59 creep was
search variance, not diameter scaling). **0** families (adversarial, Vmax to 56000) lack a clean ruler ≤ 600.
So: *every prime-rich primitive family has a clean ruler at a bounded modulus (≤ ~60 in extensive adversarial
search), independent of diameter.* The bound is a bit larger than the random-family `{8..40}`, but absolute.

**4. The diameter is already bounded (LEM-010, PROVED).** A covering family with `Vmax > 3^12` has an explicit
good period (lonely time) by Dirichlet pigeonhole, so the genuine residual has `Vmax ≤ 3^12`. Hence the
bounded-modulus lemma is a **finite check** in principle (though `3^12` is far beyond exhaustion).

**The honest limit.** The lemma is robustly true, **diameter-free**, finitely-checkable, and exact-criterion
backed. Its full unconditional proof is exactly `[SHALLOW: some bounded modulus has no scaled q/7-arc holding
≥ 6 speeds — anti-concentration]` + `[LIVE: the same modulus has a live multiplier — bounded-denominator
loneliness = LRC content]`. Both are genuinely LRC-adjacent (kps: the live half is "LRC-equivalent"), and
the prime-window shallow half can be *adversarially defeated* (finding 2), so neither half is elementary. **The
value of the redirect is not a proof but the correct posing**: a bounded, diameter-independent modulus (with
composites) in place of an unbounded pair-sum — turning "does an unbounded equidistribution hold" into "does
a bounded-modulus occupancy hold," which LEM-010 makes finite. The clean-ruler route's hardness is real and
concentrated in these two bounded-modulus statements, no longer entangled with the diameter.

Files: `lrc14_bounded_modulus_adversarial_opus_S231.py` (+`.out`). → LEM-010 (the Dirichlet diameter bound),
HYP-6015 (updated).
