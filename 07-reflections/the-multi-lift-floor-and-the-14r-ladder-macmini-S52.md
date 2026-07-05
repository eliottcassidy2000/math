# The multi-lift floor and the 14r ladder: how the lift stratum re-discovers the whole n=12 spectrum

**mac-mini-2026-07-05-S52 (HYP-4103).** Draft during the floor sweep; exact numbers
in `05-knowledge/results/lrc_multilift_floor_macmini_S52.out`.

## The setting

Residue pinning (opus-S75, kernel-pure) reduces hdich's tight side to LIFT RIGIDITY:
tight-from-above primitive 12-families with no 13-multiple are lifts
W_k = {r + 13 k_r : r = 1..12} of the AP {1..12}, and the claim is that every
nonzero lift vector k is strictly loose, quantitatively >= beta for the assembly's
loose margin. The question this session answers: **what IS the floor of the lift
stratum, and where does it live?**

## Finding 1: the 14r ladder (single lifts, r >= 7)

The sieve-surviving single lifts at r >= 7 are exactly k = r — the families

    W_r = ({1..12} \ {r}) u {14r},      r = 7..12

("kill runner r, replace it with 14r"). Their covering-min has the closed form

    M(W_r) = 14 / (13(r+1)),            r = 7..12,

with witness law (CONFIRMED exactly, all six):

    t*_r = a_r / (13(r+1)),   (13-r) * a_r == -14  (mod 13(r+1)).

Mechanism: 14r = 14(r+1) - 14 = (r+1) - 14 mod 13(r+1), so the killer's residue
at dilation a_r is exactly +14 while the base runner 13-r sits at exactly -14:
the binding pair is (13-r, 14r) for EVERY rung -- a two-point equioscillation on
opposite sides of 0 (THM-618 species; klein-S130's twoKiller binding pair is the
n=14 cousin).  The law is a CONGRUENCE, not an inverse: g = gcd(13-r, 13(r+1))
takes values in {1,2} and solvability is g | 14 -- the E3 crash that exposed this
(pow(13-r, -1, q) not invertible at r=7,9,11) was itself informative: the g=2
rungs are solvable only because 2 | 14.  The margin numerator is UNIFORMLY 14 = n;
only the modulus grows along the ladder.

At r = 12 this is {1..11, 168}: killer 168 = 14*12 = 13^2 - 1, witness t = 14/169
= 14/13^2 — the n=13 DEEP WELL (opus-S77/MISTAKE-104), i.e. the (n+1)/q*^2 species
of the n=14 construction {1..12, 182} (182 = 14*13 = 183-1, M = 14/183 = 14/Phi6(14))
ONE LEVEL DOWN, with 169 = 13^2 playing Phi6's role. The ladder is the deep well's
FAMILY: the r-th rung replaces r by 14r and pays margin 14/(13(r+1)). The n=14
construction continues the ladder at r = 13 (adjoin 14*13 instead of replacing),
where the denominator bumps from 13(r+1) = 182 to Phi6(14) = 183 — the +1 is the
structural difference between replacing a runner and adjoining one.

Single-lift floor: **14/169, now closed at floor level** — S51/opus-S77 swept the
RIGIDITY cutoff (killer <= 144 = 12B, window certifies > 1/13 beyond); the floor
claim "no single lift below 14/169" needs the beta*-level window cutoff
beta*/delta = 14*144 = 2016 (killer > 2016 is window-certified >= 14/169, one-tooth
containment), and the k <= 155 sweep below it found nothing under 14/169.
MISTAKE-104 discipline: sweep to the structural cutoff OF THE QUESTION BEING ASKED
— the rigidity cutoff and the floor cutoff differ by a factor 14.

## Finding 2: the +13 block species (multi-lifts dip BELOW the deep well)

The double-lift sweep's scan floor is **2/25** at

    {1..12} \ {4,6}  u  {17,19}     (lift 4->17, 6->19, both k=1),

BELOW the single-lift floor 14/169 = 0.08284 (2/25 = 0.08). The lift stratum's
floor is NOT at the deep-well corner; it is at height-1 BLOCK lifts, and it lands
exactly on the n=12 SECOND VALUE 2/25 = 2/(2*13-1) — the value of {1..11,24}, the
2nd-harmonic extremal (kps-S2's second_value_loose species). The lift stratum
re-discovers the global spectrum's second rung from inside.

FINAL NUMBERS (E1/E2, exact):
- M({1..12}\{4,6} u {17,19}) = 2/25 EXACTLY, witness t = 6/25, binders runner 8
  at -2 and runner 17 at +2 mod 25 (the 2nd-harmonic hole structure, kps-S2's
  emod_hole species with the hole realized INSIDE the lift stratum).
- Height-1 hypercube (all 4094 primitive nonzero C): ZERO scan failures at 2/25
  — no height-1 lift below 2/25; C = {4,6} is the UNIQUE structured slice below
  14/169.  Block anatomy: {4,5,6}->{17,18,19} = 2/21, {6,8} = 1/8, {11,12} = 1/11;
  loosest: {1..11} lifted = 1/3, full lift {1..12}->{14..25} = 14/39.
- Full l=2 structural domain (600,756 sets: w_b <= 258, w_a <= 24*w_b, killers to
  ~6200 — far beyond kps-S1's k <= 2 ground): ZERO sets below 2/25, zero exact
  escalations (393,962 sieve-closed + 206,794 witness-certified).  With the fee
  closure (both >= 259) and the one-killer window (w_a > 24*w_b), EVERY double
  lift has margin >= 2/25.
- The spectral gap (1/13, 2/25) is EMPTY on every swept stratum: singles (full
  floor domain to killer 2016), doubles (full structural domain), height-1
  hypercube, l=3 probe at k <= 4.

## Finding 3: the fee/tower thresholds and the 2l < 13 wall

The l-killer window at the rigidity level (l lifted, 12-l unlifted base <= 12,
LRC(13-l) citation margin, one-tooth/measure accounting) closes all-big lifts:

    sum_i 1/w_i < l(13-2l) / (156(13-l))   =>   strictly loose;
    all w_i >= T_l suffices, T_l = 156(13-l)/(13-2l):
    l = 2..6:  T = 191, 223, 281, 417, 1093.

The wall 2l < 13 (l <= 6) is the SAME ceiling as klein-S134's "<= 6 tops" in
lonely_of_window_multi (each top eats <= 2/13 of the window; 6 tops leave 1/13).
klein-S135's teeth_mass_far sharpens the fee ~3x — these thresholds are the crude
(safe) version. For l >= 7 the fee dies structurally; those strata have <= 5
unlifted elements and are sieve-strained, kps's ratio gate closes the narrow ones
(all lifted values within ratio 11.5), and kps-S1's k <= 2 exhaustion covers the
bottom; the SPREAD l >= 7, height >= 3 residue is the honest open remainder.

## What this hands the assembly

- The dichotomy's lift-side margin cannot be set at 14/169; the block species
  forces **beta <= 2/25**. If the E2 stratum comes back empty below 2/25, then
  beta = 2/25 UNIFORMLY across lifts and non-lifts — one constant for the whole
  TightLooseDichotomyAt, and klein's corner threshold stays the original 25B/3
  (mac-mini THM-619/620 band lane unchanged).
- Six ladder witness rows are pure rational_point_margin atom instances
  (kps-S2 shape): mu = 14, s = 13(r+1), t = a_r/s with a_r from the law.
- The Lean transcription target: klein's lonely_of_window_multi (Fin 13) at the
  n=13 rigidity level (Fin 12, 1/13 constants) — same move opus-S76 made for the
  one-window peel. The fee table above is what it certifies.

## The meta-observation

Every stratum of the lift space finds its own copy of the SAME extremal species:
the deep well (1st harmonic, "replace r by 14r") floors the single lifts; the
2nd harmonic ("shift a block by +13") floors the doubles at the global second
value; the sieve floor 1/12 = the trivial harmonic. The lift space is a
microcosm of the whole covering-min landscape — the harmonics ARE the species,
and the spectrum is the same from inside a stratum as from outside. This is the
S71 ladder / THM-610 deep-hiding picture appearing at one level down (169 = 13^2
in Phi6's role), and it is why "the classifier lesson" (S51) keeps recurring:
m = 13 is always missing because the lift stratum lives exactly on the 1/13
tightness shell, and everything below it is harmonic structure.
