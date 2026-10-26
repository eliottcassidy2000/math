# Pre-registration (opus-2026-10-26-S57), committed before reading the runs it names

Seen before writing this: the exact counts to 10^9 (they equal OEIS A057978), and the 10^6 profile
(`untouchable_profile_1e6` printout: residue-class shares, per-k abundancy bins, shape table).
Not yet seen: `untouchable_profile_1e9.out`, the 10-vertex tournament scan, the long V-ray run.

- **U-P2 ({2,3,5,7} = {3,5} x {2,7} split).** For even k in [10^8, 10^9], the untouchable share of each class
  r mod 210 is predicted by the product share(r mod 15) * share(r mod 14) / share(all even k). The median absolute
  relative error over the 105 even classes is <= 10%. Failure means the 3,5-part and the 2,7-part interact.
- **U-P3.** At 10^9, the mod-15 class with the largest share is k = 6 (mod 15), with share > 0.45. There both
  "easy" preimages die: 5 | k-1 and 3 | k-3. The smallest share is at k = 4 or 14 (mod 15), with share < 0.20.
- **U-P4.** The abundancy correlation persists at 10^9: the share for x = k/s(k) in [0.4, 0.6) is >= 1.5 times the
  share for x in [1.8, 2.0).
- **U-P5.** The multiplicity |s^{-1}(k)| of even k is over-dispersed at 10^9: variance/mean > 1.2.
- **T-P1.** No strong tournament on 10 vertices with a nontrivial automorphism has tc = H/|Aut| in {7, 21}.
- **R-P1.** On the V-ray from 7 (root-27 section, THM-4771(a) / root27 packet B3), the largest deficit excursion
  max_{i<j<=J}(K_j - K_i - 3(j-i)) over the first J = 10^5 returns lies in [15, 35]. The Haar model gives
  ln(J)/ln(phi) ~ 24, from the root y = phi of y^3 - 2y^2 + 1 = 0.
