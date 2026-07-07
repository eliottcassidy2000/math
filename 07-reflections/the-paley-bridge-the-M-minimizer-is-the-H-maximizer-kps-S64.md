---
source: kind-pasteur-2026-07-07-S64
status: PROVED (THM-640) + honest negatives on the geometric cutoffs. Owner directive:
  "pair statistics among integer families can be considered tournaments by deciding
  creatively a meaningful binary cutoff."
tags:
  - lonely-runner
  - tournaments
  - paley
  - binary-cutoff
  - THM-126
  - self-complementary
  - cross-domain
---

# The Paley bridge: the M-minimizer IS the H-maximizer

**kind-pasteur-2026-07-07-S64 (HYP-4927, THM-640).** The owner asked to turn LRC pair
statistics into tournaments by *creatively deciding a binary cutoff — as custom or convoluted
as I can imagine*. I designed four cutoffs; the arithmetic one closed a loop the project has
circled since THM-126, and the geometric ones failed informatively.

## The design goal that paid off

I did not want just *a* tournament. I wanted a cutoff whose **reversal symmetry equals
tournament complementation**, so that the LRC palindromic-extremizer phenomenon (S62) and the
tournament self-complementary theory would become the *same* statement. Two cutoffs achieve
it — the **QR-mod-p rule** and its **CRT-mod-14** cousin — and the QR one turned out to be an
identity I should have predicted:

> **`T_p^{QR}(\{0,1,…,p−1\}) = ` the Paley tournament `T_p`** (arrow `i→j` iff
> `(v_j−v_i) mod p ∈ QR_p`, `p ≡ 3 mod 4`). Exact, not up to iso.

Because the observer-inclusive LRC(p) family is `{0,1,…,p−1} = ℤ/p` and the cutoff is the
Paley circulant's own rule. Trivial once seen — and that triviality is the point.

## Why it is not trivial

The **same family** `{0,…,p−1}`:
- **minimizes `M`** (loneliness): at `t = 1/p` its orbit is the `p`-th roots of unity, the
  maximally equidistributed configuration — THM-126's spectral-flatness minimizer; and
- **maximizes `H`** (Hamiltonian paths): `T_p` is the canonical Rédei maximizer,
  `H(T_7)=189`, `H(T_11)=95095`.

Minimizing loneliness and maximizing Hamiltonicity are *opposite-feeling* extremal problems
on *different objects* (a real orbit vs a tournament). The bridge says they are the **same
object** viewed through two functors: `orbit ← family → tournament`. THM-126 already knew the
two extremal principles coincide *spectrally* (flat `|λ|` ⟺ spread orbit); the binary cutoff
makes it a **construction you can carry out with your hands** — draw the QR arrows on the
runners and you have drawn the H-maximizer.

The founding intuitions of the two halves of this project turn out to be one fact:
- *tournament side:* "a tournament is complete pairwise binary relations; the Hamiltonian
  path is the gauge that reveals its symmetry."
- *LRC side:* "the pairwise data is forced featureless" (exact uniformity, S59).
The reconciliation: **the arithmetic that forces the pairwise data featureless is the QR
field's**, and reading that field as an orientation *is* the Paley construction. The shared
fixed point is `{0,…,p−1}` — `p` runners for `p` players, the owner's "`n` events, `n−1`
edges, one Hamiltonian path."

## The reversal ↔ complement corollary (palindromes are SC)

The cutoff intertwines `V ↦ (max+min)−V` with `T ↦ T^{op}` (because `−QR_p = ` non-residues
when `p ≡ 3 mod 4`). So **palindromic families map to self-complementary tournaments**, and
the palindromic AP maps to `T_p`, which is self-complementary (`x ↦ gx`). The S62
palindromic-extremizer conjecture and THM-024's SC theory are now two faces of THM-640.
(Independently the same day, mac-mini-S44 pinned the palindrome conjecture exactly — valley
arrangement minimizes, parity obstruction — and monad-S5 built centering/attention runner
tournaments: three agents converging on the owner's cutoff directive.)

## The composite-14 shadow — the hardness of LRC(14), retold

`14 = 2·7` has no residue field, so the CRT-mod-14 cutoff cannot orient every pair by a QR
rule: the AP `{1,…,13}` maps to a **near-regular but not regular** tournament (`c3=88`), and
**42 of 78 pairs carry non-unit differences** (divisible by 2 or 7) — the arcs the cutoff
cannot make Paley. That is klein-S151's "`k+1` must be prime" seen from inside the bridge:
**LRC(14) is hard exactly because the QR-cutoff tournament of its AP cannot be a regular
Paley tournament.** A fresh, concrete way to *see* the composite obstruction.

## The honest failures (they teach the shape)

The two *geometric* cutoffs — the ones I hoped would read the density floor — both collapse
to the **transitive** tournament:
- `T_good` (leader-on-`Good_{1/7}`) is a consensus total order (leader-on-average always is).
- `T_witness` (semicircle order at the max-gap time) is transitive because the witness time
  is a *clustering* time (clustering is what makes the gap wide).

So a single orientation read from one time or one average cannot separate tight from loose —
converging with monad-S5's negative and extending HYP-4767. The floor's information is
genuinely a *measure* over all times (the density), not a snapshot; no static tournament
sees it. That is a real constraint on what any tournament-language attack can hope to do:
the bridge relates *extremal principles*, not *difficulty separators*.

## Where this points

- **THM-126 is now constructive** — worth re-examining whether the Rédei/OCF `H`-machinery
  (deletion-contraction, the independence polynomial at 2) can be pulled back along the cutoff
  to say something about `M`, or the covering/sieve.
- **The composite-14 tournament defect** (42 unorientable pairs) is a concrete finite object;
  quantifying its irregularity vs the sieve's hard residues (klein-S151) may make "why 14 is
  hard" *measurable* rather than merely structural.
- **CRT catalog item (S63 #1) advanced**: the CRT-mod-14 cutoff *is* the bi-`(mod 2, mod 7)`
  factorization applied to the pair statistics — the tournament is the visible form of the
  question "does the S–T counting object factor over `ℤ/2 × ℤ/7`?"

## Ledger

- Files: `lrc_runner_tournaments_kps_S64.py` (the four cutoffs, T_good/T_witness negatives),
  `lrc_paley_bridge_kps_S64.py` (+outs; the identity, observer-inclusive theorem, SC, CRT-14).
- Canon: THM-640. Builds on THM-126 (Paley flatness), THM-024 (SC anti-aut), HYP-4877 (step
  gauge), S59 (pair featurelessness); converges mac-mini THM-639/S44, monad-S5 HYP-4937.
- Does NOT prove LRC(14). It welds the two projects' extremal principles by an explicit
  construction and re-derives the composite-14 hardness in tournament language.
