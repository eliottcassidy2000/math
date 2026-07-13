---
source: mac-mini-2026-07-09-S65 (cont.54, 2026-07-12)
tags: [lrc14, M-spectrum, farey, stern-brocot, ostrowski, crux, deep-well, mediant]
---
# The M-spectrum near 1/14 is a Farey tree, and the crux value is the top Ostrowski rung

klein-S266 verified the LRC(14) M-spectrum just above the wall 1/14 and found the values
1/14 < 3/41 < 2/27 < 14/183 < 1/13 < 3/37. This session recognizes that spectrum as a classical
object — a Farey / Stern-Brocot mediant tree — which pins the crux value to a named rung.

## The two structures coincide
**Ostrowski ladder (my S38):** the covering-min rungs are M_k = k/(13k+1). This gives
M_1 = 1/14 (the AP, the wall), M_2 = 2/27, ..., M_14 = 14/183 (the deep well {1..12, 182}). The
ladder's two ends are the AP (non-covering, tight) and the deep well (covering, the DC minimum).
**Farey / Stern-Brocot mediant tree** rooted at the AP value 1/14 and the peeled-12-AP value 1/13:
- 2/27 = mediant(1/14, 1/13) -- IN the spectrum
- 3/41 = mediant(1/14, 2/27) -- IN the spectrum
- and so on: every spectrum value is a mediant of the two roots' descendants.
The two coincide: the Ostrowski rungs M_k are the "spine" of the Farey tree between 1/14 and 1/13,
and the mediants fill in the near-AP families (3/41, 3/37, ...). Each spectrum value corresponds to
a specific near-AP family -- a "compression scar" where a single element is detuned off the AP.

## What this says about the crux
The LRC(14) residual (klein-S265's 5->2, sharpened cont.53 to primitive) is:
[primitive non-covering => M >= 1/14, the AP wall, sieve] + [primitive DC => M > 1/14].
The DC minimum is the DEEPEST Farey rung that stays COVERING = **14/183 = M_14 = the deep well**
{1..12, 182} (182 = 14*13 supplies both 13 and 14). So:
> **crux value = min M over primitive DC = 14/183 = the top Ostrowski rung** (> 1/14 by 1/(14*183)).
Compressed DC (max <= 13*min) bottoms at 1/13 (klein-S131/S266, extremal 2*{1..12}u{13}); non-
compressed DC peels toward the deep well 14/183 (drop 182 -> {1..12}, M=1/13, then the 182 pulls it
DOWN to 14/183 = 14/182 minus a hair). The far element LOWERS M here (14/183 < 1/13) -- the mirror
of the density side where the far element RAISES J (cont.43). Same decorrelation, opposite sign on
the two functionals: J (coverage moment) rises, M (min clearance) can dip, both staying loose.

## Why this is the right frame
The extremals being Farey mediants of the AP is the exact combinatorial shadow of "the extremals are
{k*alpha} three-gap configs" (klein's algebraically-special reflection + my cont.44): the three-gap
theorem's gap lengths ARE continued-fraction convergents, and the M-values of near-AP families ARE
their mediants. So the crux M-spectrum is not a mysterious list -- it is the Farey neighborhood of
the AP, and the crux value 14/183 is the deepest covering rung. LRC(14) = "no primitive covering
family sits below the deep-well rung," a statement about the Farey/Ostrowski tree.

-> S38 Ostrowski-ladder reflection, klein-S266 (M-spectrum), THM-527 (three-gap rigidity), cont.44
(three-gap coverage), cont.43 (far element raises J -- the mirror sign), THM-612/HYP-2621 (tight
locus {AP,GW}). File: lrc14_Mspectrum_farey_macmini_S65cont54 (+ out).


## Confirmation (cont.55 + klein-S267, independent convergence): the far-element floor IS the Ostrowski forcing
The crux value 14/183 is forced by the covering constraint, verified two independent ways this session:
- **mac-mini cont.55 (far-element floor):** over primitive COVERING {1..12, f}, the only covering f are
  f = 182m (covering needs 13|f AND 14|f since {1..12} lacks both, so 182=13*14 | f), and
  M({1..12, 182m}) = 14m/(182m+1) -- the Ostrowski ladder in m, INCREASING, minimized at m=1: 14/183.
  Dense scan f=13..400: MIN over primitive+covering = 14/183 at f=182 (deep well), nothing lower.
- **klein-S267 (same, parametrized by 13k):** {1..12, 13k} realizes rung k/(13k+1); covering <=> 14|k;
  smallest k=14 => covering-min 14/183 = n/Phi6(n), Phi6(14) = 183 = 13*14+1. Resolved three conflicting
  floor claims (kps 1/12 Vmax<=22, boxeph 1/13 Vmax<=30, deep well 14/183 needs Vmax=182) as NESTED box
  artifacts (MISTAKE-141 thrice); reconfirmed 14/183 the global covering-min (15k families + CRT-escape,
  nothing below).
So the AP (rung k=1, 1/14, tight, non-covering) and the deep well (first covering rung k=14, 14/183) are
the two ends of the Ostrowski ladder, and the covering constraint (need a 14-multiple) FORCES the family
off the tight rung up to the first covering rung -- 14/183 = n/Phi6(n) is forced, not coincidental. The
crux is: no primitive covering family beats the first covering Ostrowski rung. File:
lrc14_farelement_floor_macmini_S65cont55 (+ out).
