# LRC(14) residue-orbit tournament maps forbid the regular tournament

**kps-S2-wf (2026-06-17). Theme: residue-orbit tournament generation.**

## The finding in one line

The **band-depth orbit-majority map (M4)** — a genuinely non-transitive tournament map
built from LRC(14) mod-14 residue data — **forbids the unique H-maximizer**, the regular
tournament `(2,2,2,2,2)` with `H=15` at `n=5`, **residue-exhaustively (PROVED)**. The
companion lonely-tau half-plane map (M3) only forbids it at small speed ranges
(`1..24`) but **DOES realize it** at larger speeds (e.g. `S=(3,8,16,21,26)`), so M3's
forbiddenness was a low-speed artifact. **M4 is the real find**: its forbidden classes
are provably unreachable because M4 depends ONLY on residues mod 14.

## Why this is non-obvious

The dead "overtaking/position-snapshot" map (`i->j` iff `frac(v_i tau) > frac(v_j tau)`)
is always transitive: `frac` is a total order at a single time, so `H=1` always —
no information. The way out is to **aggregate over the unit-orbit `(Z/14)*` or over the
set of lonely optima**, turning a single total order into a **majority/Condorcet
tournament over many "voters"** (the 6 units, or the several optimal taus). With >=3
voters, McGarvey-type cycles appear, so `H>1` is achievable. Both M3 and M4 clear the
non-triviality gate (verified: `H>1` realized, 8-11 of 12 iso classes at `n=5`).

## The structural mechanism (M4, the clean one)

M4 mod 14 vertex = runner; arc `i->j` iff runner `i` has greater **orbit-depth
majority** than `j`, where depth `d_i(a) = min(v_i a mod 14, 14 - v_i a mod 14)` over
units `a in {1,3,5,9,11,13}`. The depth multiset of a residue depends only on `gcd(v,14)`,
giving exactly **four depth-types**:

| type (depth multiset) | residues | role |
|---|---|---|
| `(0,0,0,0,0,0)` | `{0}` | the **parked** runner (section 0 forever) |
| `(1,1,3,3,5,5)` | coprime `{1,3,5,9,11,13}` | shallow band |
| `(2,2,4,4,6,6)` | even `{2,4,6,8,10,12}` | mid band |
| `(7,7,7,7,7,7)` | `{7}` | the **deepest/loneliest** band |

The **type tournament is strictly transitive**: `7 > even > coprime > 0`. All cross-type
arcs follow this fixed linear order; the only freedom is the **within-block** 6-voter
majority among coprime residues (or among even residues), plus a speed tie-break inside
a single residue. This transitive 4-type skeleton biases everything toward low-`H`
(transitive-leaning) tournaments. Residue-exhaustive enumeration (all residue-5-multisets
mod 14, all tie-break orders) reaches only **8 of 12** classes; the regular tournament —
which needs maximal balance impossible to manufacture from one 6-residue block over the
transitive skeleton — is **provably unreachable**.

## Why it matters for LRC

The regular/near-regular tournaments are the **maximally balanced** classes (the
`H`-maximizers, the Paley `T_5`). A counterexample to LRC(14) would be a maximally
"balanced" section occupancy — every band equally contested, no runner dominating.
The residue-orbit map says exactly that maximal-balance class **cannot arise** from the
mod-14 depth geometry: the parked-runner / midpoint-7 / even / coprime hierarchy is a
**forced transitive spine** with only block-local cyclic freedom. This is a tournament-
language echo of the covering obstruction (THM-523): the mod-14 structure is too rigid
(too transitive) to realize the perfectly-balanced configuration a counterexample needs.

## Honesty / extent

- M4 forbidden set `{H9(1,1,2,3,3), H13(1,2,2,2,3), H15(1,2,2,2,3), H15-regular}`:
  **residue-exhaustive = PROVED** (and cross-checked over 139,246 actual LRC 5-sets).
- M3 forbidden set shrinks with speed range (`1..16`: 4 forbidden; `1..20`: 2;
  `1..24`: just the regular class) but the regular class **IS realized** at larger
  speeds (`S=(3,8,16,21,26)` gives the regular tournament under M3). So **M3 forbids
  nothing exhaustively** — it is a slow-to-saturate map, NOT a genuine constraint.
- The genuine, provable forbidden class lives ONLY in M4 (residue-determined). M4's
  forbidden regular is airtight because M4 is a function of residues mod 14 alone.
- Dead maps under this theme: **M1 orbit-majority** (always transitive, `H=1`),
  **M5 crossing-parity** (realizes everything, no constraint), **M2 QR-character**
  (forbidden set is a modulus/residue-pool artifact: shrinks to 0 as the prime grows).

Scripts: `04-computation/lrc14_tourmap_residue-orbit_kps-S2-wf.py`,
`lrc14_tourmap_M4_structure_kps-S2-wf.py`, `lrc14_tourmap_confirm_kps-S2-wf.py`,
`lrc14_tourmap_M3only_kps-S2-wf.py`. Outputs in `05-knowledge/results/`.
