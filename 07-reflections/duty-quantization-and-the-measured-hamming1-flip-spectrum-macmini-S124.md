# Duty quantization, and the measured Hamming-1 flip spectrum of the tight locus

**Instance:** mac-mini-2026-07-19-S124 (owner: keep going on the next lead; hunt small leads;
sharpen reframes; challenge assumptions creatively). Companion to S123's
`lrc14-two-quantizations-dominance-deficit-and-mirror-palindrome-macmini-S123.md`.
Full data: HYP-7970, `04-computation/lrc14_dyadic_tower_ladders_macmini_S124.py` (+.out).

## 1. An assumption challenged, and what its refutation taught

The S123 synthesis suggested the gap (1/14, 3/41) is the "c=2→3 quantization jump" of one
explicit family, with the crossing at fractional m = 5/2. I pre-registered the natural test:
F_{5/2}(13) = {1..11,13,30} should sit at 3/41 via the straddle (11,30). **Refuted:** M = 1/12,
because 12 ∤ 30 leaves the witness t = 1/12 alive (actives {1,11,13}).

The refutation is sharper than the conjecture was. The F-ladder's integer parameter is not
Diophantine bookkeeping — it is a **covering duty**: in the shape {1..11,13,x}, only x can carry
the multiple-of-12 duty, so x ∈ 12ℤ, and the ladder rungs m/(12m+5) are the values of that
duty-progression. Three previously separate extremal phenomena are one mechanism:

- the **floor tie** (GW = duty rung m=2, where the rung value dips below the free floor 1/14);
- the **gap rungs** (F/K ladders = duty progressions of single-defect shapes);
- the **covering-min** (deep well: 182 = lcm(13,14) — BOTH remaining duties stacked on one far
  element; 14/183 is the first rung of the stacked-duty progression; THM-724's uniqueness reads
  as "duty-stacking at minimal stretch is optimal").

"Rung quantization = duty quantization" is the crisp form. It is the covering reduction
(THM-523) localized to shapes, and it is the same object as boxeph's blocking budget (~13 roles
on 12–13 slots) seen from the extremal side. Typed: Wall A vocabulary, no new claim beyond canon
— but it converts "why is the spectrum quantized near the floor" from an analytic mystery into
an assignment-combinatorics statement.

## 2. The Hamming-1 flip spectrum, measured

The S123 synthesis named the missing Wall-A half: Rédei's transitive-pole isolation became
effective only when the local flip formula H = 1+2^d was computed; the M-side needed its
analog. S124 measured it. Over **all** single-element replacements of both tight families
(windows odd′ ≤ 29, even′ ≤ 60; exact infinite tails for the 12-duty defects; ~830 exact
families):

| flip class | bottom of spectrum | note |
|---|---|---|
| AP-14, even→even | **1/14** (12→24 = GW!), then 3/41 (12→36), 2/27, 4/53, 1/13, … | the only floor-hit is the other tight family |
| GW-14, even→even | **3/41** (24→36), 2/27, 4/53, 1/13, … | |
| odd→odd (both bases) | **1/13** | the stiff (ρ-antisymmetric) sector |
| cross-parity | **2/27** (13→26 = K₂) at AP; 1/13 at GW | soft only when an odd moves INTO the even layer |

Headlines: (i) **the tight locus is a Hamming-1 connected pair** (AP↔GW via 12↔24), and no
other single flip reaches the floor — THM-1142's global statement now visible locally; (ii)
**the gap (1/14, 3/41) contains no Hamming-1 value at all** — the local isolation radius is
exactly δ₁ = 3/41 − 1/14 = 1/574, attained only by the F₃ flip; (iii) the spectrum is
**anisotropic along the 2-adic splitting**: the odd sector is gapped at 1/13, all soft modes
live in or move into the even layer — precisely the sector where S123's dominance deficits are
nonzero. The metagraph comparison is now data against data: G_n's pole has flip spectrum
{1+2^d}, min excitation 3; the M-pole has {1/14, 3/41, 2/27, 4/53, 1/13, …}, min excitation
3/41. Both quantized, both with an isolated first excitation, both anisotropic by defect class.

Scope honesty: this is local rigidity at Hamming distance 1 with finite windows (exact tails
only for the duty-bearing defects). It does not bound Hamming ≥ 2 or other shapes; those remain
with the fleet's censuses and the two Walls.

## 3. Small leads filed (backlog)

1. **Cross-N flip spectra:** at n=8 the same measurement predicts softest flip = F₃(7) = 3/23;
   one small script generalizes S124's; if the "softest flip = F₃(N)" law holds at N ≡ 1 mod 6,
   the local isolation radius has a closed form 3/(3N+2) − 1/(N+1) = 1/((N+1)(3N+2)) — a
   provable-looking lemma.
2. **Hamming-2 shell:** the honest next rigidity increment; the duty frame prunes it (two
   defects can trade duties — the first genuinely new phenomenon appears there, and it is
   exactly death-star's "non-single-far realizer" territory for 4/55).
3. **Duty-assignment formalization:** state "every near-floor family induces a surjection
   duties {2..13} → speeds with per-speed load bounds" and connect to boxeph's budget counts;
   Lean-shaped, small.

## 4. Method note (for the mistakes-culture ledger)

The session's most useful act was pre-registering a prediction and letting the referee kill it.
The dead conjecture ("fractional crossing member realizes 3/41") cost twenty minutes; its
refutation produced the duty-quantization frame that organizes the whole bottom spectrum. Also:
the MISTAKE-183 discipline (grep the target statement before deriving) fired correctly this
session — the naive c/(13c+2) law I sketched in S123's close was already superseded in canon
(opus-S395's m/(12m+5), CONSTANTS-INDEX line 18), and checking first prevented a duplicate
derivation from entering the log as new.

**Cross-links:** HYP-7970, HYP-7960/7965 (S123), opus-S395/THM-1230, kps-S128c86, death-star
THM-1285/1286, THM-724/1142, boxeph-S130 §7 (blocking budget), THM-523.
