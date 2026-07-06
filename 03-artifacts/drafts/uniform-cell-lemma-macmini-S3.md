# The Uniform Cell Lemma — the k-stratification apparatus at every Farey cell

**mac-mini-2026-07-06-S3 (HYP-4252).** Generalizes HYP-4232's k-stratification
(and HYP-4242's binder parity + witness determinism) from the mediant cell
(3, 38) to EVERY cell of THM-622's reduction. Verification (0 violations,
436 families, 268 distinct attained cells):
`04-computation/lrc_uniform_cell_lemma_macmini_S3.py`.

## Setting

W a finite set of distinct positive speeds with M(W) = c/q in lowest terms,
attained (grid attainment, THM-592) at t* = m/q* in lowest terms. Write
d_Q(x) for the distance of x to 0 in Z/Q. Margins: d_{q*}(v·m) ≥ c·q*/q for
all v ∈ W, with equality exactly for the binders.

**The lemma never uses in-gapness** — it is a general attainment-structure
lemma; its application to (G) is the case c/q ∈ (1/13, 2/25), q ∈ (12.5c, 13c)
(THM-622), where W ranges over hypothetical 12-element attainers.

## The lemma

**(0) Grid divisibility: q | q*.** The binder's grid distance q*·(c/q) is an
integer, so q | c·q*; gcd(c, q) = 1 forces q | q*. Write q* = qk.

**(i) Binder divisibility: k | v for every binder.** v·m ≡ ±ck (mod qk) gives
v·m ≡ 0 (mod k); gcd(m, qk) = 1 gives gcd(m, k) = 1, so k | v.

**(ii) Scaling identity: d_{qk}(kx) = k·d_q(x).** Hence quotient binders
v′ = v/k sit at exactly c mod q: the full mod-q level-c structure descends
to the quotients — every cell recurses onto its own k = 1 base.

**(iii) Kink pair: both binder signs occur, and q | v′ + w′.** At a maximizer
both an up-binder (v′m ≡ +c) and a down-binder (w′m ≡ −c mod q) exist (else a
small shift of t* raises the min); their quotient sum satisfies
(v′ + w′)m ≡ 0 (mod q), so q | v′ + w′. (v = w allowed: the single-runner
peak case 2ck = qk.)

**(iv) BINDER UNITS (new): gcd(v′, q) = 1 for every quotient binder.**
If p | gcd(v′, q) then p | v′m ∓ c and p | v′m force p | c, contradicting
gcd(c, q) = 1. — At q = 38 this is HYP-4242's both-odd parity (p = 2 case)
PLUS 19 ∤ v′; in general it collapses the k = 1 pair shapes to the
**φ(q)/2 unit pairs {a, q − a}, gcd(a, q) = 1**. Shape counts:
(3,38): 9 ✓ (matches HYP-4242); (4,51): 16; (5,63): 18; (5,64): 16; (6,77): 30.

**(v) WITNESS DETERMINISM (new, general): m ≡ ±c·(v′)⁻¹ (mod q).** By (iv)
v′ is a unit mod q, so the binder pins the witness dilation class up to the
global sign (t* ↦ 1 − t*). Every other runner's mod-q residue constraint is
then EXPLICIT per shape: for k-multiples, u′ ∉ ±v′·c⁻¹·{0, 1, …, c−1} (the
2c − 1 forbidden residues). — Cleaner proof of HYP-4242's determinism: from
q | v′(m − m′), coprimality (iv) cancels v′ directly; no prime-by-prime split.

**(vi) The k-reduction is CELL-INDEPENDENT.** Any in-gap attainer has
M < 2/25, i.e. is a no-2/25-point family, so kps's cluster-gcd ladder
(HYP-4217/4227) applies verbatim at every cell: with S = non-multiples of k,
- |S| = 0 ⟹ k | gcd(W) = 1 (primitivity) ⟹ k = 1;
- 1 ≤ |S| ≤ 6 ⟹ k ≤ 50·Σ_S|w| / (25 − 4|S|) — height-bounded;
- 7 ≤ |S| ≤ 10 — the residual, now with (iv): the ≤ 5-element quotient side
  carries a UNIT pair {a, q − a} mod q.

**Consequence.** THM-622's entire cell list — not just the mediant — now has
the same finite anchor structure: per cell (c, q), the attack is
(φ(q)/2 unit-pair shapes) × (forced witness class ±c·a⁻¹) × (explicit
forbidden residues) × (the ladder-bounded k strata + the |S| ≥ 7 residual).
The only dimension the apparatus does not bound is the cell index itself
(the c-tail), which remains the open frontier of (G).

## Non-theorem logged (sign-error trap)

A tempting mod-4 kill of cells with 4 | q, c odd ("binders odd ⟹ pair sum
≡ 2 mod 4 ≠ 0 mod 4") is FALSE: the down-binder satisfies w′ ≡ −v′ (mod q),
so v′ + w′ ≡ 0 (mod 4) automatically — odd+odd ≡ 2 (mod 4) only when
v′ ≡ w′ (mod 4). Caught before claiming; do not re-derive (the analogous
trap exists at every 2-power cell: 5/64, 7/88, 9/116, …).

## Probe (five cells, anchored templates)

8,079 families across cells (3,38), (4,51), (5,63), (5,64), (6,77): each
containing a unit pair {a, q−a} as actual elements, a 2..12 covering core
clearing level c at the forced m, random clearing completions (heights ≤ 90):
**zero in-gap realizations**; nearest approach to any cell ≈ 0.05 (all
templates land at M ≥ 0.127 — the THM-622 quantization seen at every cell);
zero families below 1/13 (the LRC(13) citation silently re-verified, 8,079×).

## Status

(0)–(vi): PROVED (elementary; verified 0 violations on 436 families spanning
268 distinct attained cells (c,q), k-spectrum reaching k = 19).
Lean: LRCUniformCell.lean (this session) — the parametric core
(dInt_scale_cell, binder_dvd_cell, grid_div_cell, pair_sum_dvd_cell,
binder_unit_cell, witness_determined_cell).

-> HYP-4232 (the (3,38) instance), HYP-4242 (parity/determinism at 38),
HYP-4217/4227 (the ladder), THM-592 (grid attainment), THM-622 (the cell
frame), OPEN: the c-tail; the |S| ≥ 7 residual.
