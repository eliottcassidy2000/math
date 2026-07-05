# LRC(14) formalization — state, validity, and what must be pursued
*mac-mini-2026-07-05-S58. A synthesis + validity audit, written to get the formalization into its best state.
Owner directive: consider the fleet's recent work and its validity; see the big picture; challenge assumptions.*

## 1. The canonical route (clean) — and the dead one (unmarked)

**CANONICAL, sorry-free:** `lrc14_of_dichotomy_and_corner` (LRCHcompSurface.lean, **0 sorry**):
```
LRC(14)  ⟸  LRCUpTo13 (citation)  +  TightLooseDichotomy  +  CornerLonely
```
- `LRCUpTo13` is a **named citation hypothesis** (not a sorry, not an axiom) — clean.
- The two open predicates are the **entire remaining mathematics**. Everything above them is a theorem.

**DEAD, and correctly marked:** `lrc14_of_template_and_corner` (LRCTemplateSurface.lean) rests on
`TemplateDichotomy` = "tight OR a **2/25-witness at denominator s ≤ 50**". **kps-S11 (HYP-4137, MISTAKE-110)
REFUTED this**: an explicit ratio-12 Fin-13 family at height ~10²² passes every hypothesis yet has no witness
below s=53. So `lrc14_of_template_and_corner` is a **dead reduction** (a true implication whose hypothesis is
now known false; it can never be discharged). **kps already added a ⚠️ CORRECTION block to the file header
(lines 1–14)** pointing to the refutation and to the canonical surface — so the hygiene is handled. *(I initially
flagged this as an unmarked bug; that was my error — I read only the legacy header line, not the correction block
above it. Lesson to myself: read the full header before claiming a hygiene defect.)*

**Bottom line:** the formal reduction of LRC(14) is in good shape — clean, to two named predicates, with the
dead Q50 surface already flagged. The formal *reduction* is not the problem; the open *mathematics* is.

## 2. The one open crux, correctly stated

`TightLooseDichotomy` = **the n=12 spectral-gap emptiness**: every primitive covering 12-set is either
- **TIGHT**: a dilated AP `c·{1,…,12}`, `M = 1/13`; or
- **LOOSE**: `∃` a real `t` with margin `≥ 2/25`.

i.e. **no 12-set has `M ∈ (1/13, 2/25)`.** (`2/25 = 2/(2n−1)` is the known second-smallest LRC value.)
This feeds LRC(14) through the compressed/spread dispatch (dominant branch already discharged, kps HYP-4087).

**The hard part is UNBOUNDED HEIGHT.** At bounded height (representatives ≤ ~52) the gap is empty by exhaustive
census (mac-mini S54/S55, opus l≥7 4.57B-leaf sweep, kps walks) — solid. The difficulty is families with **no
small representative**: kps-S11's counterexample is loose but its *witness denominator grows with height*, which
is exactly why the finite-template (Q50) refinement died.

## 3. Challenging the framing: it is ALL one phenomenon — 13-adic dilations

The borderline families — the ones that make the crux hard — are the **same object** under four names the fleet
has been circling independently:

| thread | the hard family |
|---|---|
| mac-mini MISTAKE-102 (S45/S47) | dilated APs / commensurate families that random sampling never hits |
| mac-mini CRT free-rider (S46/S47) | `c·{1..12}∪{killer}`: killer safe at a base optimum by CRT (`gcd(c,·)⊥13`) |
| opus tower-lift (S83, HYP-4126) | witnesses lift up the 13-adic tower ⟹ covers are a closed 13-adic set |
| kps CRT-lift (S11, HYP-4137) | NEEDFREE shape lifted, free moduli CRT-pinned, witness pushed to high q |

**These are one GENUS: the borderline covering families are arithmetic dilations/lifts of bounded bases, and
their loose-witness is inherited from the base.** They are NOT literally the same lemma — they are three
*complementary* dilation regimes, and it is worth being precise (this is a validity audit):
- **coprime dilation** `c ⊥ 13` (mac-mini S47): the base optimum `k/(13c)` survives because `gcd(c,killer)=1`
  spreads the killer's phase mod `c`, coprime to 13 (CRT) — the tight-AP-base free-rider case.
- **13-adic dilation** `c = 13ˡ` (opus S83): the witness lifts *up the 13-adic tower*, `13a` at level `l+1`.
- **residue-preserving lift** by `lcm(2..25)` (kps S11): the profile is preserved but the witness *moves to a
  higher free modulus* — the mechanism that killed the finite-q bound.
Together these say the unbounded-height problem is not a census run to infinity but a **reduction**: *show every
high-height loose family inherits its witness from a bounded base under one of these dilations.* The
tower-limit dichotomy is precisely the claim that **these dilations are the ONLY borderline families** — no wild
high-height family sits in the gap. Unifying the three regimes into one witness-inheritance theorem (base +
all dilations) would BE the crux. (My S47 CRT lemma is the coprime base-case, not the whole tower — do not
overstate it.)

## 4. What must be pursued (in priority order)

1. **THE TOWER-LIMIT DICHOTOMY (the crux).** The **witness-lift is already PROVED** — `speedOK13_lift`
   (opus-S84, LRCTowerLift.lean, kernel-pure, 0 sorry): a level-`l` OK-witness lifts to `13·num / 13·den` at
   level `l+1`, so **covers project DOWN the tower**. What is OPEN is the **converse dichotomy**: that the
   tower-limit covers are **ONLY 13-adic dilations** (no "wild" non-shadow high-height cover). Its first-level
   form is opus's **level-3 census mod 2197**: no non-shadow cover survives above the level-2 cover classes. If
   true, every high-height loose family is a dilation of a bounded base and inherits its witness — the loose
   branch closes at all heights. This is a **height-independent, residue-only** statement (kps-S11) — no
   MISTAKE-102 tail. *This is where fleet effort should concentrate.*
   - Concrete sub-tasks: (a) **the level-3 census mod 2197** — opus (S83) explicitly handed this to "the
     mac-mini C harness" (python was too slow at the 2028-bit class DFS); it is the witness-spectroscopy at 2197,
     the same species as the Q50 census. *This is the sharpest concrete open computation and it is in mac-mini's
     lane.* (b) fold my S47 coprime-CRT free-rider (tight-AP base case) + kps's pinned-only correction in as the
     tower base; (c) unify the three dilation regimes (§3) into one witness-inheritance theorem.
2. **`CornerLonely`** (the sub-threshold corner, killer ≤ 25B/3): the second predicate; separate, smaller,
   has its own two open items (loose-base enumeration + tables) per the S54 proof map.
3. **Formal-state consolidation** (below).

## 5. Validity assessment (the audit the owner asked for)

- **Reduction: SOUND.** `lrc14_of_dichotomy_and_corner` is 0-sorry; the citation is a proper hypothesis.
- **`native_decide` surface: LARGE (72 files).** The bounded-census and grid-table legs lean on it heavily
  (LRC14CertRoute, LRCCertTable, LRCGridValue, window packs, skeleton). native_decide is generally sound but is
  a bigger trust base than kernel reduction; the fleet correctly distinguishes "kernel-pure
  [propext, Classical.choice, Quot.sound]" results from native ones. **Recommendation:** keep the *proof of the
  two predicates* kernel-pure even if the *supporting evidence* uses native_decide; audit that no native_decide
  table silently encodes an unproven mathematical claim (spot-checks so far — klein's 144-pair AP table,
  mac-mini census — are honest data tables).
- **Bounded-height limitation: THE real gap.** Every census is bounded height. The dichotomy at *unbounded*
  height is NOT verified computationally and CANNOT be (that is item 4.1). Q50 was the attempt to make it finite
  and it FAILED. So the crux is genuinely open — do not read "census clean to height 52" as "dichotomy proved."
- **Reliability discipline (MISTAKE-101/102):** floor/spectrum claims from *random sampling* are unreliable here;
  the extremizers are arithmetic (13-adic). My high-height validation this session used **structured** CRT-lifts
  (to ~10¹⁵) — all non-AP families stayed LOOSE (≥2/25 witness), APs stayed tight — supporting the *real* surface;
  but honestly, simple L-lifts preserve the low-q witness, so they do not reproduce kps's NEEDFREE regime (kps
  already validated that one). The census's residual risk lives exactly in NEEDFREE high-height shapes = item 4.1.
- **Dead surfaces / fragmentation:** the Q50 surface is dead but *correctly marked* (kps-S11 header block). The
  broader issue is fragmentation: 23 sorries sit across ~23 files spanning the canonical node-hypotheses AND
  several older/alternative routes (Dispatch, EndgameAssembly, DominantPeel, Folding, FarPeel*, Skeleton) — some
  are live parameters, some are archaeology. **Recommendation:** a consolidation pass — one pinned top surface
  (`lrc14_of_dichotomy_and_corner`), alternative routes labeled superseded like the Q50 one, so the sorry-count
  reflects *live open obligations* not history. (I did not verify each sorry's live/dead status this session; a
  dedicated sorry-provenance audit is a good separate task.)

## 6. The honest one-paragraph state

LRC(14) is formally reduced, sorry-free, to LRC(13) + two named predicates (`TightLooseDichotomy`,
`CornerLonely`). The first predicate — the n=12 spectral gap `(1/13, 2/25)` is empty — is the crux; it is proved
at bounded height by exhaustive census but **open at unbounded height**, where the hard families are 13-adic
dilations. The finite-template shortcut (Q50) is refuted; the correct route is opus's tower-limit dichotomy
(unify with my CRT free-rider lemma), a height-independent residue statement. The formalization's best next
state = mark the dead Q50 surface, pin the canonical one, keep the two-predicate proofs kernel-pure, and
concentrate on proving the tower-limit dichotomy.
