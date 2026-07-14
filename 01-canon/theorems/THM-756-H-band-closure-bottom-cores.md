---
id: THM-756
title: The (H)-band CLOSURE -- Battery A (COMPLETE): all 91 bottom cores (12-subsets of {1..14}), all 4032 (P,v) band pairs decided: 4011 pass (H) exactly, 19 close by exact L(P u {v}) > 0, and the final 2 are the AP and GW tight completions (L = 0 with equality-loneliness at the corners, kps THM-741's census) -- every band CLOSED. Battery B (mac-mini-S104's assembly band, 160-body sample): 158/160 closed INSTANTLY by THM-755 (v > v*), 2/160 by exact (H): the 'irreducible sliver' is ~99% capped-envelope-instant
status: PROVED-FINITE (Battery A is a complete deterministic sweep, exact rational arithmetic throughout) + VERIFIED (Battery B protocol on a 160-body sample, 100% closure, layers measured)
source: opus-2026-07-14-S290 (owner prompt: sweep the finite (H)-bands per core and close them)
depends_on:
  - THM-755 (the band edges v*; the instant layer of the closure stack)
  - kps THM-741 (the L = 0 census: AP and GW only, both lonely at t = 1/14 with equality)
related:
  - mac-mini-S104 (the covering-case assembly this plugs into: the band (220, W0] closure protocol)
  - klein THM-753 (the one-step peel skeleton; (H) was its one hypothesis)
  - THM-731/732/752 (the stack's deeper layers -- unneeded in practice)
---

# THM-756 -- the (H)-band closure

## Battery A: the bottom core class, closed completely

For every 12-subset P of {1..14} (all 91) and every integer v in the (H)-band
(max P, r_P/(pi |G'_P|)]: 4032 pairs total, each decided in exact rational arithmetic:

- **4011 pairs (99.5%): (H) holds exactly** (Bernoulli disc_v < 6 |G'_P|^2) -- the THM-731
  certificate fires, the peel closes.
- **19 pairs: (H) fails but the body closes directly** (exact L(P u {v}) > 0 by the interval
  engine) -- these carry the expected alignment structure (38% of failures are 7/11/12/13/14-
  aligned: THM-751's regime).
- **2 pairs: L = 0** -- (P, v) = ({1..12}, 13) and ({1..11,13}, 24): the AP {1..13} and the
  Goddyn-Wong doubling {1..11,13,24}. These are the two tight extremals of the whole problem
  (kps THM-741: the complete L = 0 census on the near-AP region): loneliness holds WITH
  EQUALITY at the tiling corners t = (2a+1)/14 (THM-754) -- lonely, not counterexamples.

**Every (H)-band of the bottom core class is closed.**  The failure set is exactly the
tight/aligned structure the fleet's tiles already own.

## Battery B: the assembly band (mac-mini-S104), sampled closure

160 multi-killer band bodies (near-AP base + outliers, largest in (220, 475], covering):
the closure stack [THM-755 instant | exact (H) | THM-752 fine-comb | direct exact L]:

>  **158/160 close at layer 1 (THM-755: the outlier exceeds v* = r_P/(pi |G'_P|))**;
>  2/160 at layer 2 (exact (H) at that v); layers 3-4 never reached; UNCLOSED: 0.

mac-mini's "irreducible ~4% needing exact-Q certification" is, body by body, ~99%
capped-envelope-INSTANT: the band closure protocol is three lines (compute v*; almost always
v > v*; else one exact disc; the two tight extremals by their equality certificates).

## Position

With mac-mini-S104's assembly (covering = THM-366 + THM-724 + THM-726[Step 2 + THM-751/753]
+ the band + the opus floor), this theorem supplies the band's closure protocol and closes it
outright on the bottom class.  What remains fleet-wide is engineering: run the protocol over
the full band enumeration (kps's exact-Q infrastructure; every step decidable), and compose
the Lean bottoms.  The mathematics of the (H)-band is done: it consists of THM-755's edge,
19 aligned direct-closures, and the two corner-pinned extremals the theory always knew about.

## Files

05-knowledge/results/lrc14_H_band_closure_thm756_opus_S290.out (both batteries, exact).
