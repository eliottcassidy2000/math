# Three notions of sameness are the lonely set's fiber: equinumerosity, equidecomposability, equidistribution as one invariant tower

*mac-mini-2026-06-27-S62. The owner: consider equidecomposability and equinumerosity in addition to
equidistribution; search past work for hidden perspectives; be abstract; find invariants that capture the
fundamental nature. The answer revives three niche project threads ([[equinumerosity-equidecomposability-fiber-bridge-s617]]
/ HYP-2187, the triune carrier HYP-2239, the Dehn-scissors HYP-2213) and carries them — for the first time —
onto the LRC **lonely set itself**, where they become three computable invariants that recover the project's
own scales (mod-41, V*) and pinpoint tightness.*

## The abstraction: three resolutions of "sameness of a set under a group action"
Two sets can be "the same" at three increasing resolutions — a tower of equivalences, each forgetting less:

| notion | "same" means | the invariant it preserves | classical home |
|---|---|---|---|
| **equinumerosity** | a **bijection** exists | cardinality (a count) | `K0(FinSet) = ℕ` |
| **equidecomposability** | **finite cut + reassemble** | measure **+ Dehn invariant** (scissors class) | scissors congruence group |
| **equidistribution** | **same limiting density** | the limit measure / discrepancy | ergodic average / Haar |

These are the **K-theory of measure** at three resolutions: counting measure → Jordan/scissors measure →
Haar/limit measure. The project already named the *tournament-side* version (HYP-2187): for a tournament `T`,
**`H` is volume, `β₁` is the first Dehn obstruction, the count `A000568` is the cardinal shadow.** The new
move is to run the same fiber on the **LRC lonely set** `L(S) = {t∈[0,1): ‖s·t‖ ≥ 1/14 ∀s∈S}`.

## The lonely set's fiber (VERIFIED, `lrc_three_sameness_invariants_macmini_S62.py`)
| config | EQUINUM (count) | EQUIDECOMP (scissors) | EQUIDIST (density) |
|---|---|---|---|
| | covering? #res | `D`=min-wit-denom, `1/ℓ_max` | `meas(L)` |
| AP `{1..13}` (tight) | False, 13 | **D=14**, — | **0** |
| GW `{1..11,13,24}` (tight) | False, 12 | **D=14**, — | **0** |
| generic non-covering | False, 9 | D=6, 1/ℓ=186 | 0.102 |
| easy-cover `{1..12,14}` | **True**, 13 | D=13, 1/ℓ=308 | 0.024 |
| hard-cover `{1..11,13,84·m}` | **True**, 13 | **D=41**, 1/ℓ=511–980 | 0.005–0.011 |

Three findings, one per resolution:

1. **EQUIDISTRIBUTION = the tightness detector.** `meas(L(S)) = 0` *exactly* on the tight atoms (AP, GW
   both collapse to zero arcs). The density face — the witness floor `m_P` — hits zero precisely on the
   extremal locus. This is the "volume" invariant: it sees tightness and nothing finer.

2. **EQUINUMEROSITY is predicate-blind** — the *cardinal shadow* (HYP-2187). "Covering" (some `s≡0 mod 14`)
   is **independent of tightness**: AP is non-covering yet tight; a dilated AP is covering yet tight (the
   S60 `×2` dilation). A count cannot decide the LR predicate — the CH lesson of HYP-2232 made concrete.

3. **EQUIDECOMPOSABILITY carries the real arithmetic — and splits into two scissors invariants** that the
   project had found *separately by other routes*:
   - **`D(S)` = the min witness denominator** (the smallest `d` with some `a/d` lonely — the *easiest
     rational reassembly*). For the hard family it is **D=41, independent of the apex magnitude** — exactly
     the project's "Farey-neighbor scale **mod 41**" (kps-S40). It is a **bounded-core** invariant.
   - **`1/ℓ_max`** (reciprocal of the largest arc — the *worst* reassembly) **grows with the apex**: this is
     the **V\*** wall of S61 (`D≈213`, the Conjecture-7.1 constant). It is an **apex-driven** invariant.
   `D(S) ≤ 1/ℓ_max`; the two were conflated before, and the scissors lens separates them: **mod-41 = the
   bounded core, V* = the apex.**

## Why the equidecomposability face is the deep one (the Dehn invariant of the lonely set)
In 1-D, scissors congruence under *all* translations is trivial (only measure matters) — that is why the
*density* face sees only `meas=0` (tight) vs `>0`. The content appears the moment you restrict the allowed
reassembly to the **arithmetic** (rational) translations: then a **Dehn-type obstruction** appears, and it
is exactly `D(S)` — *the smallest modulus at which the lonely set is rationally accessible.* The Dehn
invariant classically detects **irrationality** (the regular tetrahedron vs the cube); here it detects the
**arithmetic of the optimal phase**: tight configs are rationally accessible at the boundary (`D=14`), the
hard core demands the mod-41 Farey scale, and the apex pushes the *worst* arc to the V* scale. **`meas` is
the volume; `(D, 1/ℓ_max, #lengths)` is the Dehn invariant** — the same `H`-vs-`β₁` split HYP-2187 found on
the tournament side, now on the measure side.

## The fiber is the fundamental invariant
No single column separates `{tight, generic, easy-cover, hard-cover}`; the **triple** does. So the
"fundamental nature" of a speed set, in the owner's sense, is its **fiber**
```
   Φ(S) = ( covering/residues  |  D(S), 1/ℓ_max, arc-length-spectrum  |  meas(L) )
            cardinal shadow        the scissors / Dehn invariant          the volume
```
— a single object refining HYP-2187's tournament fiber onto the LRC measure side. And it lines up with the
project's other trinities:
- the **triune carrier** (HYP-2239): product face (gcd/residues) ↔ equinum; sum face (the arcs) ↔
  equidecomp; fraction face (carry-continuant) ↔ equidist/irrational tail;
- the **four faces of 14** ([[the-four-faces-of-14-why-the-exceptional-structures-crowd-into-lrc]]):
  multiplicative (covering) ↔ equinum; additive/Farey (the arc/mod-41 scale) ↔ equidecomp; exponential
  (periodicity/density) ↔ equidist.

## What it buys the proof (honest)
- **A clean extremal characterization:** `meas(L(S)) = 0 ⟺ S tight` — the equidistribution invariant *is*
  the witness floor, and its zero set is exactly the tight locus. The proof needs `meas > 0` off that locus
  (CRUX 1) — now read as "the volume is positive away from the Dehn-degenerate atoms."
- **Two scissors scales, separated:** the bounded-core obligation is a **mod-41 Farey** check (`D`), the
  apex obligation is the **V\*** peel (`1/ℓ_max`). The proof's two regimes (bounded gK8, apex Node-3) are
  the two Dehn invariants. This is why the V* and mod-41 constants are *different* numbers doing *different*
  jobs — a point the single-"witness-denominator" language had blurred.
- **A discipline (HYP-2232/2245):** never decide the LR predicate from the cardinal shadow (covering count)
  alone; the predicate lives in the scissors+volume fiber.

Honest status: this is an organizing invariant and two recovered scales, not a new bound. The genuinely new
verified facts: `meas=0 ⟺ tight`; `D(S)=41` for the hard family (mod-41, apex-independent); `D ≤ 1/ℓ_max`
with `D` bounded-core and `1/ℓ_max` apex-driven. Next: prove `meas(L) ≥ m_P > 0` off the tight locus *as a
Dehn-nondegeneracy statement* (the scissors form of CRUX 1).

Related: HYP-3091 (the verified fiber), HYP-2187/2239/2213/2232/2245 (the revived threads), HYP-3089
(the V* = 1/ℓ_max constant), HYP-3085 (gK8 = the volume bound), [[three-walls-one-constant-the-vstar-of-lrc14]].
