---
id: THM-857
title: Scale-one Hamming-six exact component recursion
status: RESERVED / COMPUTATION IN PROGRESS — no radius-six closure claim is made by this stub
source: codex-2026-07-15-S10 H6 continuation
depends_on: [THM-815, THM-845]
related: [THM-770, THM-820, THM-856, HYP-6820]
verification:
  - 04-computation/lrc13_scale_one_hamming_six_component_recursion_codex_S10.cpp
---

# THM-857 — Scale-one Hamming-six exact component recursion

This file reserves the theorem namespace while the exhaustive exact recursion
is running.  The proved input is THM-815(C): after numerically ordering the six
proper lifts, the longest strict-safe component gives a finite bound for every
next lift.  The companion source evaluates those bounds with exact rational
intervals and no height cutoff.

Two checkpoint facts have been reproduced:

- the global first-prefix census is `194,040`, while applying each deletion
  core's own longest-component cap sharpens the traversed first layer to
  `83,881` states;
- on the missing-odd root, the recursion finds the expected covering terminal
  `14,16,18,20,22,24=2[12]`.

The safe-band operation is cached by speed.  This is lossless because the
operation bank depends only on the inserted speed; the dynamic state is the
exact residual interval union.  On the first root the cache changes runtime
from roughly `38` seconds to `6` seconds while preserving every depth counter.

The theorem claim, frozen census, independent replay, and stored certificate
remain deliberately absent until all `924` deletion cores finish.  In
particular, this stub does **not** yet exclude primitive Hamming-six covers and
does not alter the global `n=12` sporadic-branch status.

