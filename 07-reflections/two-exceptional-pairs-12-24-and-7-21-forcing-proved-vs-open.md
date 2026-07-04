# Two exceptional pairs: {12,24} and {7,21} — the same shape, one proved, one open

*kind-pasteur-2026-07-04. The owner asked for analogies between {12,24} (the tight LRC coverers
= AP, GW) and {7,21} (the forbidden Hamiltonian-path counts). The literal bridge is false — but
the structural one is exact, and it points at a proof technique, not a coincidence.*

## What the two pairs are

- **`{12, 24}`** — in the single-swap coverer family `{1,…,11,13,X}` (`12 ∣ X`), these are the
  only `X` with `M = 1/14` (tight): `X = 12` is the AP, `X = 24` is GW. Characterizing them —
  the tight-locus rigidity — is **open** (Perarnau–Serra). My HYP-4076 formula `M = k/(12k+5)`
  proves it *for this family*: `X ≥ 36 ⟹ M > 1/14`, so the tight `X` are exactly `{12, 24}`.
- **`{7, 21}`** — the only two `H`-values (Hamiltonian-path / OCF counts) that **no** tournament
  achieves. Characterizing them — the H-gap — is **proved**: `H = 7` by THM-029 (three
  pairwise-conflicting 3-cycles force a fourth), `H = 21` by THM-079 (a 464-line six-way block
  of every `α₁ + 2α₂ = 10` decomposition).

## The literal bridge is false

The obvious guess — that the tight LRC families *map* to the forbidden-H tournaments — is wrong.
The AP at its tight time sits on the 13 nonzero 14th-roots; its rotational tournament (`i→j`
iff `(j−i) mod 14 ∈ {1,…,6}`) has `H = 2 641 713`, a Paley-scale count, `≡ 4 mod 7`. GW's tight
config has a residue collision (`24 ≡ 10 ≡ 10`, missing `12`), so it isn't even a simple
tournament there. The tight families are maximally-connected, high-`H`; they are the *opposite*
of the `H`-hole. So `{12,24}` and `{7,21}` are not the same object in two languages.

## But the shape is identical

Both are two-element exceptional sets of the form **`base · {1, m}`**:

| pair | base | multiplier set | status | forcing mechanism |
|---|---|---|---|---|
| `{12, 24}` | `12 = n−2` | `{1, 2}` (GW = double) | **open** rigidity | residue: single-swap must cover `q=12`; only `12, 24` tight |
| `{7, 21}` | `7 = n/2` | `{1, 3}` | **proved** H-gap | cycle: `α₁=3` forces a 4th; `α₁+2α₂=10` six-way-blocked |

Each is "the exceptional configuration is small, and beyond it a combinatorial forcing kicks
in." For `H`: past `α₁ = 3` you cannot avoid a fifth cycle; every larger blocked count is a
forced overflow. For LRC: past `X = 24` you cannot stay on the 14th-roots; `M = k/(12k+5)` rises
monotonically to `1/12`, the coverer forced loose. The pair is small **because** the forcing is
tight — the same reason in both problems.

## The resonances (poetic, flagged as such)

- **`84 = 12·7 = 4·21`.** The first *covering* member of the residue-liar family,
  `{1,…,11,13,84}`, is `12 · 7` (tight-base × heptagon) and `4 · 21` (four × the top forbidden
  `H`). It blocks the 14-grid (`84 = 6·14`) and is lonely at `37/89`, `89 = F₁₁`. The two
  obstructions' bases multiply to the one number where the census's hardest family turns
  covering. I do not have a theorem here; it is a genuine numeric confluence.
- **Multipliers `2 · 3 = 6`.** GW doubles (`24 = 2·12`), the H-gap triples (`21 = 3·7`), and
  `2·3 = 6` is the order of `14` mod `Φ₆(14) = 183` — opus's Eisenstein lever (HYP-4047). The
  two forcing multipliers compose to the cyclotomic order that governs the deep-well
  covering-min. Again a confluence, not a proof.

## The real lead: the proof technique transfers

The useful content is not a bridge object but a **method**. The H-gap is the settled instance
of exactly the problem the LRC rigidity is stuck on — "prove a two-element exceptional set is
*all* of the exceptional set" — and it was settled by **exhaustive decomposition-blocking**
(THM-079: enumerate every way to reach the count `10`, show each is forced to overflow). My
HYP-4076 is the *single-swap* LRC analogue of the `H=7` case (THM-029): one obstruction, closed
by an explicit forcing (the residue table). The open LRC rigidity — all lift and coverer
combinations — is the analogue of the *full* H-gap, and THM-079's six-way block is the template:
enumerate the finitely many residue/lift combinations that could stay tight, and show each is
forced loose by a `k/(12k+5)`-type formula. The tournament side already walked this exact path
to the end. That is the transferable lesson: **rigidity proofs in this family are
decomposition-blocking proofs**, and one of them (the H-gap) is finished.

## Honest placement

Not a theorem, a structural reading. The literal tournament-`H` bridge between `{12,24}` and
`{7,21}` is false. What is true: they are the *same shape* — small exceptional sets forced by a
combinatorial obstruction, `base·{1,m}` with base `∈ {n−2, n/2}` — one proved (the H-gap, by
exhaustive blocking) and one open (the rigidity). The proved one is the method template for the
open one, and my residue-liar formula is the first block of that template on the LRC side. The
`84 = 12·7 = 4·21` and `2·3 = 6` confluences are real numbers and flagged as not-yet-theorems.

---
*Linked: [[the-residue-liar-family-closes-by-formula-fibonacci-in-the-denominator]] (the single-swap
block), [[the-tight-locus-is-the-arithmetic-progression]] (S37). THM-029/THM-079 (H-gap, proved),
HYP-4070/mac-mini (GAP-A), HYP-4078/kps (single-swap rigidity), HYP-4047/opus (Eisenstein order
6). Script: `lrc14_tight_vs_Hgap_analogy_kps.py` (AP tournament H = 2 641 713, bridge falsified).*
