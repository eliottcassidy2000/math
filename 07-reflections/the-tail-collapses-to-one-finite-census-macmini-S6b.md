# The tail collapses to one finite census: why (A) is a corollary, not a lemma

**mac-mini-2026-07-06-S6b (HYP-4302).**  A work-skipping reframe of the (G)
tail, prompted by the owner's "Newman-shaped covering impossibilities —
reframe, skip work."  Builds on S4's J-K reduction, opus-S98's residue bridge,
and the sibling/kps/my 7-spread work (which it subsumes).  Verification:
`05-knowledge/results/lrc_accumulation_frombelow_macmini_S6b.out`.

## The set-up (S4's reduction)

(G) [the second-value gap is empty: no covering 12-family has M ∈ (1/13, 2/25)]
was reduced to **(A)** no coupled proper 2-subtorus U ⊆ (ℝ/ℤ)¹² has
M(U) ∈ (1/13, 2/25], **plus (C)** a finite 1-dim census.  Three lanes grew up
to attack (A): kps's `torus_split_rung` (≤6-lifted), the sibling's support-6
kill, my 7-spread census (≥3-class).  All aimed to prove a **2-D covering
impossibility**.

## The reframe: (A) is not a covering lemma — it's a limit

The Newman-shaped impossibility here is **trivial once you stop thinking in
2-D**.  A 2-torus U = ⟨r, ℓ⟩ arises as the lift-limit of the 1-D families
v⁽ᴺ⁾ = r + N·ℓ (J-K Section 3: the 1-tori ⟨v⁽ᴺ⁾⟩ equidistribute into U,
so M(v⁽ᴺ⁾) → M(U)).  Two facts:

1. **Sub-torus containment (trivial).**  v⁽ᴺ⁾ = r + N·ℓ ∈ span_ℝ(r, ℓ), so the
   1-torus ⟨v⁽ᴺ⁾⟩ ⊆ U.  Hence
   `M(v⁽ᴺ⁾) = max_{⟨v⁽ᴺ⁾⟩} min-dist ≤ max_U min-dist = M(U).`
   *No G-K, no accumulation theorem — just "a subset's max is ≤ the whole's."*
2. **Convergence (J-K lift-limit, already assumed by S4).**  M(v⁽ᴺ⁾) → M(U).

Together: **M(U) = sup_N M(v⁽ᴺ⁾), approached from below.**  So if
M(U) ∈ (1/13, 2/25], then for all large N, M(v⁽ᴺ⁾) ∈ (M(U)−ε, M(U)] ⊆ (1/13,
2/25] — **in-window 1-D families at unbounded height N.**

## Why that is impossible (opus's bridge closes it)

opus-S98's `margin_of_residue_witness` (GREEN): margin(v, a/q) = margin(v mod q,
a/q) at each modulus q.  So **"M(v) ≥ 2/25" is a residue-class property** — it
holds iff some q ≤ 50 witness clears (margin ≥ 2/25), a fact about v's residues
mod q, not v's height.  The Q50/template census establishes this for every
residue class (finitely many mod each q ≤ 50).  Hence M(v) ≥ 2/25 for **every**
covering family at **every** height — in particular every v⁽ᴺ⁾ — contradicting
the unbounded-height in-window (M(v⁽ᴺ⁾) < 2/25) families forced above.
(The endpoint M(U) = 2/25 falls too: then M(v⁽ᴺ⁾) ≤ 2/25 and ≥ 2/25 force
M(v⁽ᴺ⁾) = 2/25 for all large N — infinitely many exact-2/25 attainers on one
lift line — but the attainer species is finite along any line (S4: M_s = 2/25
only at s=1, then retreats to 2/17).)

> **Therefore M(U) ∉ (1/13, 2/25].  (A) holds — including the endpoint 2/25 —
> as a COROLLARY of (the finite census is gap-empty) + (the J-K lift-limit).**

No 2-D covering impossibility is needed.  The 7-spread census (mine), the
support-6 kill (sibling), and `torus_split_rung` (kps) are **good independent
confirmation** — and remain the *direct* proof for anyone who distrusts the
J-K preprint — but they are **not required** for the tail.

## What this skips, and what remains

- **DROPPED:** the entire ≥3-class covering lemma; and S4's G-K dependency
  (F1 used G-K for the 1/13 endpoint; the sub-torus argument replaces it and
  handles the whole half-open window, endpoint included).
- **REMAINS (the single obligation):** the finite Q50/template census is
  gap-empty — one finite, decidable object (kps's `TemplateDichotomy`
  predicate; my S55 Q50 census 511,947 survivors + spectroscopy are its
  verified core).  opus-S98 already said "only the Q50 census remains open";
  this reframe makes it *literally* the only thing, by showing (A)/(C)/rays
  are one quotient-by-residues object (opus's phrase, now precise).

## The three lanes are one

| lane | was | now |
|------|-----|-----|
| CRT rays | two_band_transport (opus, formal) | a residue-family class in the census |
| coupled 2-tori (A) | 7-spread / support-6 / rung | corollary of the 1-D census (this) |
| 1-D cells (C) | k-stratification chain | the finite census itself |

All three are the statement **"the finite residue-family census (mod q ≤ 50)
has no member with M ∈ (1/13, 2/25)."**  The tail is that single check.

## Honest ledger

- **Trivial/formal:** sub-torus containment; opus's residue bridge.
- **Preprint (already assumed by S4, not new):** J-K lift-limit convergence
  M(v⁽ᴺ⁾) → M(U).  For canon, either cite with the standard caveat or re-prove
  the equidistribution (Weyl on the 2-torus — elementary).
- **The open crux:** the finite census / template dichotomy gap-emptiness.

The reframe doesn't prove the tail; it **removes (A) from the obligation list**
and points all remaining effort at the one finite census.  That is the
work-skip: three lanes → one check, and a 2-D covering lemma → a one-line
sub-torus inequality.
