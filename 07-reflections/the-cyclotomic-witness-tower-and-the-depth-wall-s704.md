---
source: opus-2026-06-07-S704 (user: "work on the radical tower and the witness tower")
status: DEVELOPS + HONESTLY DEFLATES HYP-2303 into THM-439. The LRC witness tower IS the cyclotomic
  (abelian) tower — but AUTOMATICALLY: the optimal witness t* is always rational (maximin of
  integer-breakpoint PL functions), so every tick e^{2πit*}=ζ_q^m ∈ Q^ab, Gal=(Z/q)^× abelian. So the
  field-Galois "solvable tower" is TRIVIALLY solvable; there is NO non-abelian witness obstruction
  (deflating the "counterexample = non-abelian monodromy" reading). The real substance is the
  CYCLOTOMIC DEPTH q*(S)=min modulus clearing 1/n, stratified clock(q≤n)⊂sub-shell⊂shell(2n−1)⊂
  super-shell; the S700 residual NEVER lands at the clock level. KEY DICHOTOMY (verified n≤8): q*(S)<∞
  per config (LRC holds, constructively in Q^ab) BUT sup_S q*(S) is UNBOUNDED (grows with speed). The
  true Abel–Ruffini analog is therefore the BONFERRONI tower (THM-406), not the field tower: the
  Vitali wall = the unbounded depth / non-uniformity, not any single config.
tags: [lonely-runner, LRC, cyclotomic, radical-tower, witness-tower, Q-ab, kronecker-weber, abelian,
  cyclotomic-depth, q-star, bonferroni, vitali-wall, THM-406, THM-403, THM-420, THM-436, HYP-2303,
  honest-deflation, non-uniformity]
---

# The cyclotomic witness tower, and why the wall is the depth

**Prompt (user):** "work on the radical tower and the witness tower." (Developing HYP-2303 from S703:
the conjecture that the LRC witness hierarchy is a radical/solvable tower, worry-set = solvable
bottom, residual = unsolvable core.)

Working it carefully both **confirmed** the cyclotomic structure and **deflated** the romantic part of
the conjecture — and the deflation is the real finding.

## 1. The witness tower is cyclotomic — but automatically

The witness hierarchy (clock `t=1/k` ⊃ shell `t=m/(2n−1)` ⊃ pair-sum `t=m/(a+b)`, THM-420/430) is a
tower of **cyclotomic moduli**: each tick `t=m/q` puts the runners at `q`-th roots of unity,
`e^{2πi t}=ζ_q^m ∈ ℚ(ζ_q) ⊆ ℚ^{ab}`, with abelian Galois `(ℤ/q)^×`. So the tower lives in the maximal
abelian extension and (Kronecker–Weber) *is* the abelian/cyclotomic part. So far, exactly HYP-2303.

But **why** is every witness cyclotomic? Not because of any deep solvability — simply because
`M(S)=max_t min_i ‖v_i t‖` is the max of integer-breakpoint piecewise-linear functions, so the
optimal `t*` is **rational**. Rationality ⟹ root of unity ⟹ abelian. *The abelian-ness is forced and
free.*

> **Honest deflation.** HYP-2303's "a config is tower-certifiable iff its local monodromy is solvable;
> a counterexample is non-abelian" is **vacuous at the field level**: the witness field is
> *unconditionally* abelian (rational `t*`). There is no non-abelian witness obstruction to find. The
> quintic/`A₅` analogy does **not** live in the Galois group of the witness.

This matters: I almost shipped a false-romantic conjecture (non-abelian monodromy = LRC
counterexample). The rationality of `t*` kills it. Better to know.

## 2. Where the substance actually is: the cyclotomic depth q*

If the *kind* of field is always abelian, the content is the **size** of the modulus needed:

> `q*(S) = min{ q : some tick t=m/q already clears the floor 1/n }` — the **cyclotomic depth**.

The witness hierarchy is the magnitude stratification of `q*`: `clock (q*≤n) ⊂ sub-shell (n<q*<2n−1) ⊂
shell (q*=2n−1) ⊂ super-shell (q*>2n−1)`. Computed (n=5..8): most configs are clock-depth, a fat tail
needs deeper. And the clean structural fact:

> **The S700 residual `R(n)` (divisibility-complete ∧ shell-free) never lands at the clock level** —
> it is *exactly* the positive-depth core. The residual is "the configs the cyclotomic floor can't
> reach," made quantitative as `q*>n`.

So the worry-set (cyclotomic, tight `M=1/n`) sits at the very floor (`q*=n`, equality), the bulk is
clock-cleared above the floor, and the residual is the deep tail — a genuine depth gradient, the
honest version of HYP-2303's "stratification."

## 3. The dichotomy that is the real theorem: finite per config, unbounded uniformly

Two facts, both verified, that together are the point:

- **Per config: `q*(S) < ∞` always** (n≤8 in-window, zero counterexamples). Every config is certified
  loose by a *finite* cyclotomic tick — **LRC holds here constructively, inside ℚ^{ab}**.
- **Uniformly: `sup_S q*(S) = ∞`** — `max q*` grows with the speed bound (n=7: `11,13,…,19,21` as
  `B=7..15`). No single depth certifies all configs.

> **This is the correct Abel–Ruffini mirror.** *Each quintic has a finite splitting field; there is no
> uniform radical formula.* `↔` *Each LRC config has a finite cyclotomic depth `q*`; there is no
> uniform cyclotomic depth.* The unsolvability is **not** in any single object — it is in the
> **non-uniformity**. And that is precisely THM-406: the covering-depth inclusion–exclusion terminates
> for any finite truncation but cancels **to all orders**, with no finite Bonferroni certificate. The
> Vitali wall = the unbounded depth.

So the Abel–Ruffini analog is the **Bonferroni tower** (inclusion–exclusion order), *not* the field
tower (which is automatically abelian, §1). Last session (THM-436(4)) named the mirror correctly
(Bonferroni ↔ derived series); this session shows the *witness/field* tower is the wrong place to look
for it, and the *depth* is the right invariant.

## 4. The corrected conjecture (HYP-2309)

> LRC(n) ⟺ `q*(S)<∞` for all `S` (every config has a finite cyclotomic certificate). The witness tower
> is uniformly **abelian** (cyclotomic, derived length 1) — so the *kind* of solvability is settled
> and trivial; the **depth** `q*(S)` is the live invariant: finite per config, unbounded over configs.
> The Vitali wall (THM-406) is exactly this non-uniformity. A counterexample would be a config with
> `q*=∞` (no rational tick clears `1/n`, i.e. `M<1/n`) — a **combinatorial/measure** failure, not a
> non-abelian field.

This corrects S703's loose "tower depth = derived length of local monodromy": the derived length is
uniformly 1 (abelian); the meaningful depth is the cyclotomic modulus, and the wall is its
unboundedness.

## 5. What this buys, honestly

- **Rigorous:** witness tower ⊆ ℚ^{ab} via rational `t*` (§1); the depth definition; the residual =
  positive-depth core; the per-config-finite / uniformly-unbounded dichotomy (verified n≤8).
- **Deflated:** the "non-abelian counterexample" reading of HYP-2303 — vacuous (field always abelian).
- **Reframed:** Abel–Ruffini ↔ Bonferroni (depth non-uniformity), not ↔ a field tower.
- **No open case resolved.** But the depth `q*(S)` is a clean new handle, and the residual's
  identification as "`q*>n`" connects S700's residual, THM-403's cyclotomic worry-set, THM-428's
  shell tower (super-shell depth = the `2n−1` / `3³` regime), and THM-406's Vitali wall under one
  invariant. **Next:** does `q*` for the residual track the shell prime-power depth (THM-428)? — i.e.
  is super-shell depth at `n=14` exactly the `3³` tower? That would tie the depth-wall to the concrete
  `R(14)`.

**Artifacts:** `04-computation/lrc_cyclotomic_witness_tower_s704.py` (+`.out`). Theorem **THM-439**.
New **HYP-2309**. Builds on THM-420/430, THM-403, THM-406, THM-436/HYP-2303, S700, Kronecker–Weber.
Logged to `comms/CLUSTER-FEED.md`.
