# The Apex-7 unification: one obstruction across four worlds, and what it does (and doesn't) prove

*kind-mendel-2026-06-22-S6. Owner inspiration: balanced cuts (w=3) vs a single unbalanced weight (w=1)
supplying H=49,75 = the apex tile; tournament↔even-graph equinumerosity; {7,21} forbidden H; Kuratowski/
Wagner; Beurling–Selberg/trig. Pulled main repeatedly; the team converged hard on my S5 (mac-mini S34/S36
HYP-2876/2879 = "user Ideas 1/2/3"; kps S31e E_7 odd holes). This consolidates the unification with
independent verification + one proved fragment + an honest proof-status, and pushes the proof thread.*

## The hint, decoded and verified

mac-mini-S36 (HYP-2879) is the owner's hint precisely: **"weight = atom count."** H is multiplicative
over strong components (Moon; `H=∏H(C_i)`), so the "atoms" are strongly-connected tournaments and the
H-spectrum is their multiplicative closure. "Single unbalanced weight w=1" = a **single strong core**
(one atom, the most unbalanced condensation); "balanced cut w=3" = a 3-component condensation.

**Independently verified this session** (`04-computation/lrc14_apex_atom_hspectrum_kindmendel.py`,
`..._prime_coverage...`):
- **49 and 75 are strong-core atoms at n=7 and n=8** (49,75 ∈ spectrum, strongly-connected-realized).
- **The apex tile is the strong-connectivity switch**: the transitive tournament has `H=1` (not strong,
  fully reducible); flipping *only* the apex tile `(0,n-1)` (the source-sink arc) gives
  **`H = 1 + 2^{n-2}`** (33 at n=7, 65 at n=8), strongly connected — mac-mini's "single-tile = 1+2^d
  hypotenuse." The forced atom **49 sits at apex-flip distance 2** (e.g. apex + tiles {(2,0),(6,1)} at n=7).

## PROVED fragment — why 49 is an irreducible atom (and it's the unique such)

> **Proposition.** Any tournament with `H=49` has a single non-trivial strong component, with `H=49`
> (i.e. 49 is an *atom* value, never a product).

*Proof.* `H` is multiplicative over strong components (Moon). The factorizations of `49=7^2` into
integers `>1` are `7·7` and `49`. By THM-029, **`H=7` is forbidden** for every tournament, so no strong
component has `H=7`; hence `7·7` is unrealizable. So a tournament with `H=49` must have its non-trivial
H concentrated in one strong component of value 49. ∎

This is the cleanest face of the apex-7 obstruction: **49 = 7² exists *only* as a single irreducible
unit, precisely because 7 is forbidden.** And `7` is the *unique forbidden prime* H-value (`21=3·7` is
composite), so `49=7²` is the *unique* "forced-from-a-forbidden-prime" atom. The apex tile is what
creates the single strong core that carries it.

## The four worlds, unified by the prime 7 = apex of 14

| world | the apex-7 obstruction |
|---|---|
| **H-spectrum** | `H=7` forbidden (THM-029, Kuratowski K5 = 5 overlapping odd cycles in Ω); `49=7²` exists *only* as an irreducible atom (proved above) |
| **Tiling model** | the **apex tile** `(0,n-1)` = source-sink arc: unique boundary tile, max range, H-coefficient `2^{n-2}`, the cut/cycle hinge that creates strong connectivity |
| **Even-graph metagraph** | `E_7` loses chordality with **odd holes = 1496 C₅ (H=7 pentagons) + 196 C₇ (apex heptagons = 14²)** (kps S31e) — the same odd-cycle obstruction |
| **LRC(14)** | the apex arc = **observer's loneliness gap** (lonely ⟺ apex clears 1/n both sides); covering forces a multiple of 14 ⟹ **D≥15 floor** (apex-7, proved); `14=2·7` composite is why the polynomial method (needs a large prime factor) misses it |

**The organizing principle — the 2·7 difficulty factorization.** LRC(14)'s hardness factors as its
modulus does, `14 = 2 · 7`:
- the **2** = the *even-graph* structure: kps-S36 (HYP-2867) proved GOOD depends only on the
  *complement-even* part of the cluster, reducing the uniform floor to complement-symmetric clusters;
- the **7** = the *apex-prime odd-cycle* obstruction: forbidden `H=7`, `E_7` odd holes, `49=7²` atom,
  the `D≥15` apex floor.

A proof should plausibly *dispatch the 2 first* (even-graph / complement-even reduction) and then face
the irreducible **7** (the single atom / apex). This is exactly why `n=14=2·7` is the first open case:
the polynomial method wants a large prime factor, but 14's prime factor `7` is small *and* it is the apex
prime where the odd-cycle obstruction becomes irreducible.

## Proof push — the over-determination route, and its honest ceiling

The team's live route (mac-mini HYP-2878): LRC failure ⟺ `N(S,D)=0` ⟺ the unsafe sets cover `Z/D` (a
covering system); a covering set should cover only *few* primes, so a witness prime exists. **Verified**:
random covering sets cover only 3–11 of the 24 primes in `[15,120]`; even the **loosest** `{1..11,13,84}`
covers 17/24 but is still witnessed at 7 primes (incl. `D=41`); the (non-covering) AP `{1..13}` covers
**all** 24 (witnessed only at `D=14`). So the AP/covering split is real and covering sets are
prime-certified *in the tested range.*

**Honest ceiling (verifier note).** For primes `p > max speed`, the residues *are* the speeds, so "S
covers `Z/p`" becomes the *continuous* non-loneliness condition — i.e. **the large-prime over-determination
is literally Node-3 (decorrelation) in disguise.** Therefore the finite-certificate / covering-system /
character-sum / over-determination formulations are all *equivalent*, and all reduce to:
> **(compact reduction: worst case has bounded speeds ≤ V*) + (a finite check over bounded cores).**
They sharpen and unify the *small-modulus* content but do **not** bypass the compact-reduction crux
(= Node 2, "consec maximizes / spreading decorrelates"). The finite basis `{83,89,21}` and the apex-7
`D≥15` floor are genuine, but "the finite certificate closes LRC(14)" would be an overstatement until the
compact reduction is proved.

## New leads (many, as requested)
1. **Dispatch the 2 then the 7**: do the even-graph/complement-even reduction (HYP-2867) to land on a
   purely apex-7 problem, then attack the 7 via the atom / `D≥15` / `E_7`-odd-hole structure.
2. **Forced-atoms `p²` for forbidden prime `p`**: 7 is the *only* forbidden prime, so `49` is the *only*
   forced atom of this kind — a genuinely distinguished value. Track where `49` sits in the winding
   tournament `T(x)` of an LRC config: does a *lonely* `x` force `H(T(x))` to avoid the `7`-divisible
   atom shells, or to hit `49`?
3. **Lonely ⟺ apex-tile state of `T(x)`**: loneliness = apex arc opens = (winding tournament) observer
   is a source (THM-381) = apex tile "down"/master-cycle structure. Make this an exact tile-state
   criterion and read the LRC threshold off the apex tile directly.
4. **Beurling–Selberg on the apex factor**: the character-sum main term `(6/7)^13` has exponent 13 =
   #runners; factor it by atoms — if the speeds split into resonance atoms, `N(S,D)` should factor, and
   the apex atom is the binding one.
5. **`E_7` heptagon = 14²**: verify/exploit the `196 = 14²` coincidence (kps) — is the count of apex
   heptagons in `E_n` always `(2n)²`-like, and does it count something in the LRC(2n) floor?

## Status
No node closed; but the apex-7 unification is now *verified* across all four worlds, with one clean proved
fragment (`49` forced atom), the `2·7` organizing principle, and an honest statement of the proof ceiling
(the reformulations ≡ compact reduction + finite check). Best next step toward a proof: the even-graph
"dispatch the 2" reduction, landing on the irreducible apex-7 core.

→ HYP-2880 (new), HYP-2879/2878/2876/2872/2867/2864/2856/2847, THM-029/079/519/527, kps-S31e (E_7),
OPEN-Q-106/108.
