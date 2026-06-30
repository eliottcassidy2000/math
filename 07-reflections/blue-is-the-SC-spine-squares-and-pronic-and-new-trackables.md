# Blue is the SC spine: squares & pronic, the blue=SC identity, and new things to track

*klein-2026-06-29-S13. Hunting patterns in the blue (grid-symmetric tilings = the odd-boundary, S675b) and black (Eulerian) structure, and defining new invariants. Two clean patterns fell out, one of them sharp; and the genus/blue synthesis needs a correction.*

## Pattern 1 — the blue count is `2^{square}` (odd) / `2^{pronic}` (even)

`BLUE(n) = #grid-symmetric tilings = 2^{(#orbits of the grid involution on the m tiles)}`. The exponent
`e(n)` has a clean closed form (verified n=3..9):

> **`e(n) = k^2` for odd `n = 2k+1` (SQUARES); `e(n) = k(k-1)` for even `n = 2k` (PRONIC).**

So `BLUE(n) = 2,4,16,64,512,4096,65536` for `n=3..9`, exponents `1,2,4,6,9,12,16`. The blue tilings form a
sub-cube `Fix(grid) ≅ Q_{e(n)}` — a "self-anti-diagonal" core of dimension a square (odd n) or pronic (even
n). The square/pronic alternation is the staircase's anti-diagonal reflection counting its own orbits;
it is the cleanest closed form in the black/blue story (sharper than CLAUDE.md's fraction-exponent
`e − m = 0,-1,-2,-4,-6,-9`).

## Pattern 2 (sharp) — BLUE = the SC spine, exactly

Per-iso-class (verified n=3..6): the number of classes containing **≥1 blue (grid-symmetric) tiling**
equals the **self-complementary count `SC(n)`**:

| n | classes | SC | blue-classes | pure-blue | mixed | pure-black (= NS) |
|---|---|---|---|---|---|---|
| 3 | 2 | 2 | **2** | 2 | 0 | 0 |
| 4 | 4 | 2 | **2** | 1 | 1 | 2 |
| 5 | 12 | 8 | **8** | 3 | 5 | 4 |
| 6 | 56 | 12 | **12** | 2 | 10 | 44 |

> **A class contains a grid-symmetric tiling ⟺ it is self-complementary (SC).** Equivalently: pure-black
> classes = the NS (non-self-complementary) classes exactly; blue lives precisely on the SC spine.

The structural reason: for tournaments transpose = complement (adjacency-transpose = arc-reversal), and the
grid involution (anti-diagonal staircase reflection) is the tiling-level transpose; a grid-symmetric tiling
is a transpose-fixed = self-complementary tournament. So blue support = SC = the spine = the `R`-fixed
classes (THM-584's diagonal `{T,S}`). This sharpens CLAUDE.md's "transpose-self never pure-black" to the
exact count identity `blue-classes = SC`.

## The correction to the genus/blue synthesis (be honest)

S12 said "blue = the obstruction = genus = #blue floor atoms." With blue-support = SC now known, that needs
care: `SC(n) = 2,2,8,12,88` is NOT the genus `0,0,1,2,2`. So blue is NOT the obstruction *count*. The honest
picture:

- **Blue support = the SC spine** — the *arena* where the odd-boundary lives (it is `R`-fixed / σ-even).
- **S675b's "odd boundary" = the ODD multiplicity** of grid-sym tilings within each SC class (over an SC
  class the grid-fixed count is odd; over NS it is 0). The boundary-ness is a *parity of multiplicity*, not
  a class count.
- **The genus / cusp-form obstruction = the binding SUBSET** within this structure (the doublet at
  `N=14`), not the whole blue support. The genus counts the *independent binding atoms*; the blue/SC spine
  is where they sit.

So: blue tells you WHERE the obstruction lives (the SC spine, with odd grid-sym multiplicity); the genus
tells you HOW MANY independent binding atoms there are (1 for LRC14 = the doublet). Different invariants on
the same arena. (This is itself a polysemy lesson — "blue" the support vs "the obstruction" the binding
subset are not the same size.)

## New things to track (defined)

1. **`BLUE(n)`** — the grid-symmetric tiling count `= 2^{e}`, `e = k^2`(odd)/`k(k-1)`(even). The
   square/pronic exponent is the trackable signature; `Fix(grid) ≅ Q_{e(n)}` is the blue sub-cube.
2. **The blue=SC identity**: `#{classes with a grid-sym tiling} = SC(n)`; `pure-black = NS`. Track whether
   it persists (verified n≤6; structurally = transpose=complement).
3. **Blue multiplicity `μ♦(C) = #grid-sym tilings in class C`** — `0` if `C` is NS, ODD if `C` is SC
   (S675b's odd-boundary). Track its distribution over the SC spine (the "odd-boundary spectrum"). The
   pure-blue split (`μ♦(C) = #tilings(C)`) vs mixed is `2,1,3,2` pure-blue (n=3..6).
4. **The black Eulerian cycle rank** (S675b's open question, defined here, uncomputed): `b_1` of the black
   (non-grid-sym) line-lift carrier — "stable cycle basis / component rank." A trackable for the bulk side.
5. **The square/pronic alternation** as a parity diagnostic: odd `n` (apex prime exists, Paley, the
   LRC-relevant case) → square exponent; even `n` → pronic. The LRC apices (`n` odd at the apex level) sit
   on the SQUARE side. Worth watching against the genus/nu_2 parity facts (HYP-3586).

## The shape of it

The black/blue split has a clean skeleton: blue is `2^{square/pronic}` tilings sitting **exactly on the SC
spine** (R-fixed / transpose-self), carrying an **odd multiplicity** that is S675b's boundary; black is the
even/Eulerian bulk on the NS sea. The genus-obstruction is a small binding subset of the blue spine (the
doublet at `N=14`), not the spine itself. So the new clean trackables are `BLUE(n) = 2^{sq/pronic}`, the
`blue=SC` identity, and the odd-boundary multiplicity `μ♦` — and the correction that blue marks the *arena*,
the genus marks the *count*, of the obstruction.

See HYP (this session), [[merged-line-parity-even-odd-s675b]] (codex, black/blue),
[[three-evens-the-royle-sandwich-and-the-genus-is-the-odd-boundary]] (HYP-3591),
[[the-genus-is-the-local-global-gap-and-the-one-master-dichotomy]] (HYP-3587), THM-584 (SC = R-fixed spine),
THM-578 (the doublet binding atom).
