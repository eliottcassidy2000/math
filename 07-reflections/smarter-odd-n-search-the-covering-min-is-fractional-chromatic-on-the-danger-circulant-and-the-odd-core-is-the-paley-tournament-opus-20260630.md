# Smarter odd-n covering-min search + its repo connections: M(S)≥m/q ⟺ the speeds {sa mod q} are an INDEPENDENT SET in the circulant "danger graph" C(q; ±1..±(m−1)), so the covering-min is a FRACTIONAL-CHROMATIC / independent-set problem (search the danger-graph independent sets that are also covering — a constraint-satisfaction, not brute force); and the off-cusp odd core is a CIRCULANT on Z_p whose apex gap is its Cayley spectral gap (THM-590) — the apex-7 IS the QR/Paley TOURNAMENT (connection {1,2,4}, 14 directed 3-cycles, the OCF's odd cycles), so the LRC apex and the tournament side meet at the circulant

*opus-2026-06-30. Owner: smarter search for the odd-n realizability + connections with tournaments / repo
structures; poke around. Poked: the LRC covering-min is a chromatic/independent-set problem on a circulant,
and its apex (the off-cusp odd core) is the same Z_p circulant as the Paley tournament — the bridge.*

## The smarter search: covering-min = fractional-chromatic on the danger circulant
The brute-force enumeration of covering sets is hopeless (and greedy witness-first picks bad multiples). The
right reformulation:
> `M(S) ≥ m/q` ⟺ every speed avoids the `m`-ball of `0` at witness `a/q` ⟺ **`{s·a mod q : s∈S}` is an
> INDEPENDENT SET in the circulant "danger graph" `C(q; ±1,…,±(m−1))`** (vertices `Z/q`, edges joining
> residues within distance `< m`). So the covering-min is: *the densest packing of a COVERING set into a
> danger-graph independent set* — a **constraint-satisfaction / fractional-chromatic** problem, not a search
> over all sets. This is the classical lonely-runner ↔ circular-chromatic link, specialized to covering sets.
**Search ideas (poked):**
- **A. Witness-first / band (done right).** Fix `(q,a,m)`; the allowed speed-residues are `a^{-1}·{m,…,q−m}`.
  ENUMERATE covering sets with residues in that band (not greedy) — the densest is the covering-min at `q`.
  The band restricts residues mod `q`, collapsing the search.
- **B. Descent / apex.** The off-cusp covering-min has a proper **odd core `O ⊊ Z_p`** (a circulant); minimize
  over `O` using the apex gap `g(O)` (the Cayley spectral gap, THM-590) as the per-level density bound, then
  lift through the 2-adic descent.
- **C. Circulant-tournament.** The odd core is a Z_p connection set; the **OCF odd-cycle count** is an
  invariant — search connection sets by their spectral / odd-cycle profile.
- **D. Fractional-chromatic / IP.** Solve the covering-min as an integer program: maximize the witness gap
  subject to (independent in the danger circulant) ∧ (covering `{2,…,n}`). LP-relax for bounds. Far smarter
  than enumeration.

## The connection: the apex IS a circulant — and IS the Paley tournament
Computed apex gaps `g(O) = min_{k≠0}|Σ_{x∈O}ω^{kx}|²` (Cayley spectral gaps) on `Z_7`:
| `O` | `g(O)` | what it is |
|---|---|---|
| doublet `{0,1}` | `0.198` | the descent attractor (`C_7` cycle slack) |
| `{1,3,5}` | `0.308` | the next descent core |
| **`{1,2,4}` (QR)** | **`2.000`** | **the Paley TOURNAMENT connection set** |
| `{3,5,6}` (NQR) | `2.000` | the complement (reversed Paley) |
| **`{0,1,5}` (covmin-7 odd core)** | **`2.000`** | the off-cusp covering-min's circulant — same gap as Paley |
| `Z_7` (full) | `0.000` | the cusp (apex measure vanishes) |
> **The apex-7 IS the Paley/QR tournament.** The QR set `{1,2,4}` is the Paley tournament's connection set; its
> **14 directed 3-cycles** are the OCF's odd cycles (`H(Paley_7)` odd, Rédei). And the **odd core of the
> odd-n covering-min** (`{0,1,5}` at n=7) sits at the *same* apex-gap level `2.0` as the Paley tournament — an
> off-cusp circulant on `Z_7`. So the LRC covering-min's apex and the tournament side **meet at the `Z_7`
> circulant**: the LRC odd core, the Paley connection set, the OCF odd cycles, and the apex spectral gap are
> one structure. (This is the project's `apex-7 = the C_7 = the Paley/circulant odd-cycle core`, now seen from
> the covering-min side.)

## What the connections buy (the bridge, both ways)
- **Tournament → LRC:** the OCF / odd-cycle machinery (the apex gap = the Cayley spectral gap, the Paley
  circulant, Rédei parity) is the **measure column** of the LRC; the odd-n covering-min's odd core is a
  circulant whose spectral gap THM-590 controls. The smarter search B/C uses this — minimize over circulant
  odd cores by their spectral/odd-cycle profile.
- **LRC → tournament:** the LRC's "danger graph independent set" / circular-chromatic view (search idea D) is
  a graph-coloring frame that the tournament metagraph also lives in (the metazeta, the chromatic number of
  `G_n`). The covering-min IP and the metagraph chromatic problem are the same species.
- **The observer** unifies: the marked origin (LRC) and the marked vertex (tournament) both read the circulant
  `Z_p` apex — the odd core / the connection set — through its spectrum (the apex gap) and its odd cycles.

## Honest status
- **Poked/computed (opus):** the covering-min = independent-set-on-danger-circulant (chromatic) reformulation;
  apex gaps of `Z_7` cores/connection sets (QR/Paley `{1,2,4}` and covmin core `{0,1,5}` both `=2.0`); Paley
  tournament 3-cycle count `=14` (the OCF odd cycles); 4 smarter-search ideas (band, descent/apex,
  circulant-tournament, fractional-chromatic/IP).
- **NOT yet done (honest):** none of A–D has been run to convergence to PIN the exact odd-n covering-min
  (e.g. n=9's `4/33`); the greedy witness-first found suboptimal sets. The reformulation is the deliverable;
  the IP/band-enumeration is the next implementation.
- **The bridge:** the LRC apex (off-cusp odd core, a `Z_p` circulant) = the Paley/QR tournament structure (the
  OCF odd cycles, the Cayley spectral gap THM-590). Tournaments and runners meet at the circulant `Z_p`.

Related: both-…cycle-spectral-ramanujan (the apex = Cayley spectral, Ramanujan), reconciliation-…IS-the-cusp
(even/odd, off-cusp), the-observer-on-the-tournament-side (the marked vertex), the-metazeta (chromatic);
THM-590 (apex gap), mac-mini HYP-3594 (C_p), the Paley/QR tournament work; OPEN-Q-108.
