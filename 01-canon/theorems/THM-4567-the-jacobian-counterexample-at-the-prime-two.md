---
id: THM-4567
title: "THM-1300's Keller map at the prime 2: F mod 2 is dominant and purely inseparable of degree 2 (F_2 y^2 + F_2^3 + F_1^2 F_3 = 0 mod 2); the number of 2-adic integral preimages of w depends only on w mod 8; the image of Z_2^3 has Haar measure 11/32 and a random 2-adic point has 0, 1, 2 other preimages with probabilities 7/16, 3/8, 3/16; and for every odd w there is an integer triple collision F(0, 2w, -(63w^2+1)/4) = F(1, (w-3)/2, (13-3w)/2) = F(-1, (w+3)/2, (13+3w)/2) = ((w^2-1)/4, 2w, 0)"
status: "PROVED (Hensel with det JF = -2; the mod-2 identity verified symbolically; the collision family by substitution). FINITE-EXACT: the mod-8 fibre census and the image fraction 176/512 = 11/32 (stable for k = 3..6; re-counted independently here mod 2^k, k = 3..6); the merge-partner law (two independent checks by the session's alg reader, one counting 2-adic roots of the fibre cubic). Found by the session's reader of openai/math #193 (Serre positivity), which itself gives no leverage here."
session: mac-mini-2026-10-07-oaimath3
source: 05-knowledge/results/oai3_two_orbits_twos_and_threes_20261007.md
scripts:
  - 04-computation/experiments/oai3_20261007_readers/alg/ (c1_fibre_multiplicities.py, c2_mod2_and_2adic.py, c3_check_cubic.py, c4_integer_collisions_and_mod2.py, c5_mixed_char_conservation.py, + .out)
related:
  - THM-1300 (the map: u = 1 + xy, F = (u^3 z + y^2 u (4+3xy), y + 3x u^2 z + 3x y^2 (4+3xy), 2x - 3x^2 y - x^3 z)), THM-1310 (odd p: N in {0,1,3}), THM-1345 (plane family), THM-4562 (fibre counts over F_q)
---

# THM-4567 — the Jacobian counterexample at the prime 2

## Statements (PROVED unless marked)

1. **Mod 2.** `F_2 y² + F_2³ + F_1² F_3 ≡ 0 (mod 2)` as polynomials (all coefficients even). So `F mod 2` is dominant and purely inseparable of degree 2. Compare odd `p`, where `F` is étale mod `p` (THM-1310).
2. **2-adic fibres.**
   * Since `det JF = −2`, Hensel shows that the number `N(w)` of preimages of `w` in `Z_2³` depends only on `w mod 8`. Each preimage fills exactly 2 residue classes mod 8.
   * The mod-8 census (FINITE-EXACT) has fibre sizes `{2: 112, 4: 48, 6: 16}`. So `μ(F(Z_2³)) = 176/512 = 11/32`, the same at every level `2^k`, `k = 3..6`.
   * A Haar-random `P ∈ Z_2³` has 0, 1 or 2 other integral preimages ("merge partners") with probabilities `7/16, 3/8, 3/16`. An injective map would have image measure 1/2.
3. **Conservation (reader's computation, FINITE-EXACT).**
   * On the 3 residue classes `(0,1,0), (0,1,1), (1,0,1)` the mod-2 fibre is isolated of length 2, and every point has exactly one partner in its own disc.
   * On the other 5 classes the intersection is improper and there is no in-disc partner.
   * This is the precise sense in which "merges are intersections" for this map. The Collatz map has no such finite flat structure: its fibres are infinite.
4. **Integer collisions.** For every odd `w`:

       F(0, 2w, −(63w²+1)/4) = F(1, (w−3)/2, (13−3w)/2) = F(−1, (w+3)/2, (13+3w)/2) = ((w²−1)/4, 2w, 0),

   e.g. `F(1, −1, 5) = F(−1, 2, 8) = F(0, 2, −16) = (0, 2, 0)`. This is THM-1345's plane family at `s = 2`. The smallest integer collision is `F(2, −1, 2) = F(0, −1, −4)`. THM-1300 records only the non-integral collision.
5. **Local multiplicities.** At the three points of THM-1300's triple collision each local length is 1 with no higher Tor (FINITE-EXACT). Serre's intersection multiplicity adds nothing for this étale map.
