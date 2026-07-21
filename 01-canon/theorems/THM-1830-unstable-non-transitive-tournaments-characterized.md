---
id: THM-1830
title: "THE UNSTABLE NON-TRANSITIVE TOURNAMENTS ARE A TRANSITIVE SKELETON WITH ONE 3-CYCLE ATOM, AND EXACTLY ONE OF THEM IS BLUE (self-complementary) — at ODD n only. Characterizing the GIT-unstable non-transitive tournaments of THM-1825 (which exist from n=7). (1) STRUCTURE: an unstable tournament has a char_A-root of multiplicity > n/2, which by Lemma B is a Jordan block; a REDUCIBLE tournament has char_A = ∏ char(SCC) (block-triangular), a singleton SCC contributes x (eigenvalue 0) and a 3-cycle SCC contributes x³−1 = (x−1)(x²+x+1) (no 0-eigenvalue; verified strongly-connected SCCs of size 3,4,5 have no 0-eigenvalue). So the 0-multiplicity = #singleton SCCs, and unstable ⟺ #singletons > n/2. Since a non-trivial SCC needs ≥ 3 vertices, for 7 ≤ n ≤ 12 the ONLY way is exactly (n−3) transitive singletons + one 3-cycle atom, giving char_A = x^{n−3}(x³−1) and 0-multiplicity n−3 > n/2 ⟺ n > 6 — which is exactly WHY they first appear at n = 7. (2) STRONGLY-CONNECTED TOURNAMENTS ARE NEVER UNSTABLE: exact check of 95661 strongly-connected n=7 tournaments gives max root multiplicity 3 = (n−1)/2 < n/2, attained by Paley-7 with char_A = (x−3)(x²+x+2)³ (SEMISTABLE, its non-Perron eigenvalues share Re = −½ — a numpy real-part cluster faked multiplicity 6, caught by exact factoring). So unstable non-transitive ⟹ reducible. (3) COUNT + BLUE: the iso classes = the n−2 ranks of the 3-cycle atom in the SCC order; the complement involution reverses that order (rank ↦ (n−3) − rank), so a class is SELF-COMPLEMENTARY (BLUE / grid-symmetric / SC) iff the atom is CENTERED, which happens for exactly ONE rank iff n is ODD, and for NONE iff n is EVEN. Verified: n = 7 → 5 classes, exactly 1 blue (score [0,1,3,3,3,5,6]); n = 8 → 6 classes, 0 blue. (4) POINTS OF SYMMETRY: every unstable non-transitive tournament has automorphism group exactly Z₃ — the 3-cycle atom's rotation, a LOCAL point of symmetry, the rest being rigid — and the unique blue one additionally carries the global complement-symmetry. So these tournaments sit NEAR the transitive nullcone vertex (a large nilpotent block) with only a local Z₃, the OPPOSITE pole from opus-S434's symmetric-intransitive Paley/j=0 (global symmetry, semistable); the single centered/blue class at odd n is the one point where a local atom's symmetry becomes global"
status: >
  (1),(3),(4) VERIFIED by construction at n = 7 (5 classes, char_A = x⁴(x³−1), |Aut| = 3, one
  SC) and n = 8 (6 classes, 0 SC), and the SCC-product / centering arguments are exact.  The
  "no 0-eigenvalue in size-3,4,5 strongly-connected SCCs" is a 3000-sample check per size
  (size-3 is exact: the 3-cycle char is x³−1).
  (2) EXACT over 95661 strongly-connected n=7 tournaments (integer char poly + factor): max
  multiplicity 3.  The Paley-7 artifact is recorded as a method note (cluster only after the
  exact check).
  SCOPE: the clean form "(n−3) singletons + one 3-cycle" is proved for 7 ≤ n ≤ 12; at n ≥ 13
  two 3-cycle atoms fit ((n−6) singletons + two 3-cycles, 0-mult n−6 > n/2 ⟺ n > 12), so the
  characterization broadens there — NAMED as open.  Also: only the λ = 0 (nilpotent-skeleton)
  case is characterized; a reducible tournament with an SCC contributing a repeated NONZERO
  integer eigenvalue is not covered, and none occurs at n = 7 (all 5 have λ = 0).
source: kind-pasteur-2026-07-20-S128c133 (owner: characterize unstable non-transitive tournaments; blue iso classes; points of symmetry)
depends_on:
  - THM-1825    # transitivity is the deepest nullcone point; unstable non-transitive exist from n=7
related: [THM-1810, THM-1820]   # opus: the symmetric-intransitive (Paley/j=0) opposite pole
script: 04-computation/unstable_nontransitive_construct_kps_S128c133.py (+ .out)
---

# THM-1830 — the unstable non-transitive tournaments, characterized

THM-1825 found that from `n = 7` the GIT-unstable cone contains non-transitive tournaments.
Here is exactly what they are, and where the "blue" and "points of symmetry" live.

## (1) Structure: a transitive skeleton with one 3-cycle atom

An unstable tournament has a `char_A`-root of multiplicity `> n/2`, which by Lemma B (THM-1825)
is a **Jordan block** (geometric multiplicity is `≤ ⌊n/2⌋`). A **reducible** tournament splits
into strongly-connected components in a transitive order, and

> `char_A = ∏_SCC char(SCC)`  (block upper-triangular).

A **singleton** SCC contributes `x` (eigenvalue `0`); a **3-cycle** SCC contributes
`x³ − 1 = (x−1)(x²+x+1)` — no `0`-eigenvalue (and no strongly-connected SCC of size 3, 4, 5 has
a `0`-eigenvalue). So the multiplicity of `0` equals the **number of singleton SCCs**, and

> **unstable `⟺` `#singletons > n/2`.**

A non-trivial SCC needs `≥ 3` vertices, so for `7 ≤ n ≤ 12` the only option is `(n−3)`
singletons plus **one 3-cycle atom**:

> `char_A = x^{n−3}(x³−1)`,  `0`-multiplicity `= n−3 > n/2 ⟺ n > 6`.

That inequality is exactly why they first appear at **`n = 7`**.

## (2) Strongly-connected tournaments are never unstable

Exact factoring over **95 661 strongly-connected `n = 7` tournaments** gives max root
multiplicity **3** `= (n−1)/2 < n/2`, attained by **Paley-7**, `char_A = (x−3)(x²+x+2)³`
(semistable). Its six non-Perron eigenvalues all have `Re = −½`, so a real-part cluster faked
multiplicity 6 — caught by the exact check. So **unstable non-transitive `⟹` reducible.**

## (3) Count, and the blue class

The iso classes are the `n−2` **ranks** of the 3-cycle atom in the SCC order. The complement
`T ↦ T^op` **reverses** that order, `rank ↦ (n−3) − rank`, so a class is **self-complementary
(BLUE / grid-symmetric / SC)** iff the atom is **centered** — one rank iff `n` is **odd**, none
iff `n` is **even**:

| `n` | unstable non-transitive classes | blue (SC) |
|---|---|---|
| 7 | 5 | **1** (score `[0,1,3,3,3,5,6]`) |
| 8 | 6 | **0** |

The blue one is the fixed point of the complement involution — the atom sitting at the exact
middle of the transitive skeleton.

## (4) Points of symmetry, and the two poles

Every unstable non-transitive tournament has automorphism group **exactly `Z₃`** — the 3-cycle
atom's rotation, a **local** point of symmetry embedded in an otherwise rigid transitive frame.
The unique blue class (odd `n`) additionally carries the **global** complement-symmetry.

So these tournaments live **near the transitive nullcone vertex** (`x^{n}` with the vertex
replaced by `x^{n−3}(x³−1)` — a large nilpotent block lightly perturbed by one cycle), carrying
only a local `Z₃`. That is the **opposite pole** from opus-S434's most-intransitive Paley/`j = 0`
tournaments, which have **global** symmetry and are semistable. The single centered/blue class
at odd `n` is the one point where a local atom's symmetry becomes global — the bridge between the
"near-transitive, locally-symmetric" and the "self-complementary" worlds.

## Named next

- **`n ≥ 13`:** two 3-cycle atoms fit (`(n−6)` singletons + two 3-cycles, `0`-mult `n−6 > n/2`),
  so the characterization broadens; enumerate the atom-multiset strata and their blue count.
- **Nonzero-integer unstable eigenvalues:** a reducible tournament whose SCC contributes a
  repeated nonzero integer eigenvalue would be unstable off `λ = 0`; none occurs at `n = 7`, but
  the general possibility is uncharacterized.
- **The blue parity** (`1` blue at odd `n`, `0` at even `n`) is a clean instance of the project's
  SC/complement parity lore — worth linking to the blue-line count formulas.
