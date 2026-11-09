audit K2 (2026-11-08) outputs. Produced by the auditor's programs in ../ (see ../aud_run2.sh).

o1213_lt1125.txt        aud_sc13 13 1125, merged: "E h count strongcount" (order 12) and "O h count strongcount" (order 13),
                        every self-converse tournament with H <= 1124 (counts of reduced labelled codes)
o1213.log               the four production parts (42,704 parents; 3,157,912 order-11 intermediates; 317,728 order-12
                        tournaments; 29,016 level-2 and 96,512 level-4 direct-DP checks, no failure)
aud_sc13_validation.log the validation runs: NO = 5, 7, 9, 11 with no H bound, and NO = 11 with LIM 1125
v13_5.txt, v13_7.txt, v13_9.txt   validation outputs (identical to audit K1's reduced censuses of orders 4-9; the
                        NO = 11 output equals audit K1's anti_10_red.txt and o11_full.txt and is not duplicated here)
aud_lists2.out          aud_lists2.py: comparison with sc_odd13_1124.out, statement 5's lists, delays, composite-only, f_sc,
                        the [1131, 2001) residue, and the walk-file check
aud_walk_check.out      aud_witness16.c on the packet's sc_walk.out (169 lines)
aud_witness_recheck.out audit K1's aud_witness.c on the corrected sc_witness.out (857 lines)
aud_sandwich.out        aud_sandwich.py: 711 and 855 as C3 => Z => C3^op at order 13
