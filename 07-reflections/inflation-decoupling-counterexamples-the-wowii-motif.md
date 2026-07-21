# The inflation/decoupling counterexample motif (the WOWII-103 lesson)

*opus-2026-07-20-S437. Prompted by the owner sharing google-deepmind/formal-conjectures PR #4482
(disproof of Written-on-the-Wall II Conjecture 103).*

## What the WOWII-103 counterexample is

Conjecture 103 (an *auto-generated* graph inequality) coupled the **independence number** `α` to
the **average eccentricity** and the **largest induced bipartite subgraph** `b`. It was killed by
an explicit **11-vertex graph**: a triangle with four pendant leaves on each of two of its
vertices, giving `α = 9`, `b = 10`, avg eccentricity `30/11` — so the conjectured bound evaluates
to `8` and the inequality reads the false `9 ≤ 8`. Found by **exhaustive search over all 2¹¹
subsets** and **certified in Lean** (`decide +native`).

Two techniques, one lesson:

- **Technique (construction): inflation via pendants.** The leaves *pump one invariant* (`α`:
  leaves are mutually independent) while leaving a *coupled* invariant fixed (eccentricity: leaves
  sit next to the core). The conjecture assumed the two moved together; the leaves **decouple**
  them.
- **Technique (certification): exhaustive small-case + Lean.** The invariants are finite and
  machine-checked.

> **The lesson: a conjectured inequality between two invariants dies to a construction that
> *decouples* them.** "Inflation" (adding sacrificial structure that moves one invariant and not
> the other) is the generic way to decouple.

## The repo already runs this motif — three times

**1. The LRC extremal family is a leaf-inflation (verified).** The second tight family
`{1,…,11,13,24}` is `{1,…,13}` with the core speed `12` **replaced by the sacrificial "leaf"
`24 = 2·12`**. Both have `M = 1/14` at `t* = 1/14`, because `24·(1/14) = 12/7` reproduces the same
danger comb as `12`. The leaf pumps the *speed set* without changing the *binding invariant* `M` —
**exactly the WOWII pendant.** The `(1/14, 3/41)`-gap witness `{1,…,11,13,36}` is another leaf.
So the repo's LRC extremal-family archaeology (THM-1230/1235) is inflation-counterexample hunting
under another name.

**2. THM-1820 is the repo's own decoupling disproof.** HYP-8600 conjectured "H-extremal =
3-cycle-extremal (the character/Paley tournament)." THM-1820 refuted it by showing `H` (Schur-
concave) and the 3-cycle count (Schur-convex, `= \binom n3 − Σ\binom{s_i}2`) are **maximised at
opposite strata** — the concentrated-spectrum vs the balanced one. That is precisely WOWII-103's
move: a conjectured coupling of two invariants, broken by exhibiting they *decouple*.

**3. GMC(4) and the Sym³ counterexample are small-explicit disproofs.** The owner's `P =
(1+Z₂)(W₁(1−Z₁)+W₂)` (6 terms) disproving GMC(4), and the Sym³ `π|X` (THM-1770) disproving a naive
Keller-injectivity hope, are the same genre — a *small, checkable* object breaking a clean
conjecture — differing only in that they are algebraic rather than found by subset enumeration.

## The transferable method, made explicit

| WOWII step | repo instrument (already present) |
|---|---|
| **inflation** (pendant leaves) | LRC: `2×` speed / near-multiple; JC: append an identity coordinate (THM-1300 raises dim 3→n); tournament: a dominating/dominated vertex |
| **exhaustive small-case** | the repo's sweeps (59048 nullcone, 32768 switching classes, 2242028→8156 witness-branch) |
| **Lean certification** | `TournamentH7` — the repo *already* Lean-certifies tournament invariants; it could certify the extremal LRC invariants (`M = 1/14` at `t*`, the `(D,s)` data) by `decide`/`native_decide` exactly as the WOWII PR certified `α=9, b=10, ecc=30/11` |

## Concrete new targets this unlocks

1. **Lean-certify the LRC extremals (cheap, high-value).** The WOWII PR is a template: encode
   `{1,…,13}` and `{1,…,11,13,24}` as finite speed sets, and `native_decide` that
   `M = 1/14` at `t* = 1/14` (a finite rational check over the pair-sum denominators — THM-401
   modulus `2n−1 = 27`). This gives a **machine-checked** anchor for the whole LRC(14) program,
   matching the repo's ≤13-runner citation standard with an in-house certificate.

2. **Inflation-hunt the H-extremiser (open since THM-1820).** `H` is Schur-concave and the
   maximiser drifts off Paley for large `n`. Apply pendant-inflation to a concentrated-spectrum
   core (add a near-dominating vertex) and exhaustively search small `n` for the `H`-maximiser's
   structure — the WOWII 2¹¹-enumeration is directly imitable on tournaments (`2^{C(n,2)}` is
   feasible to `n ≤ 6`, and switching-class reduction (THM-474) cuts it 64× at `n=7`).

3. **Decoupling audit of the repo's invariant inequalities.** Every conjectured inequality
   between two tournament/LRC invariants (H vs c₃, score-spread vs diameter, `D` vs `s` in the
   located-maximizer) is a WOWII-103 candidate: try to *decouple* the two by inflation before
   trying to prove it. THM-1820 shows this catches real errors (it caught mine, HYP-8600).

## The meta-point

WOWII conjectures are *machine-generated*; this repo generates *its own* conjectures (the HYP
index) at scale, and several have been refuted the same way (HYP-8230 numerology, HYP-8450
cyclotomic, HYP-8600 H-character, MISTAKE-156 the AP dilate). **The healthy reflex is to attempt
an inflation/decoupling counterexample *first*, exhaustively and (where finite) in Lean, before
investing in a proof.** The WOWII-103 PR is an external, Lean-certified instance of exactly the
discipline the repo's MISTAKES ledger keeps re-teaching: *stress a proposed coupling against a
decoupling construction before believing it.*

## Verification

`04-computation/inflation_wowii_motif_opus_S437.py` — the `{1,…,11,13,24}` leaf-inflation of
`{1,…,13}` (both `M=1/14`); the decoupling framing. The WOWII invariants (`α=9, b=10, ecc=30/11`)
are from PR #4482's Lean-certified enumeration.
