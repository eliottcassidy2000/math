# CRON SPEC -- the radius-7 single-cluster sweep (opus-S314, HYP-6925)
# STATUS: flagship prefixes (4 exceptional + {1..5}) run by opus-S314 (see
# beat_lemma_referee_opus_S314.out); THIS SPEC covers the remaining 787.
#
# JOB: for each 5-subset P of [12] (prefix), R = [12] \ P:
#   for each single-cluster pattern dvec (base residue in R; d_i = (r_i -
#   base) mod 13 + 13k <= 42; sorted consecutive gaps <= 7; 7 distinct):
#     1. exact integral Ia = integral_E a_dvec (breakpoint trapezoid, Fractions)
#     2. C_A = 28*B(E,dvec) + 4*sum(d_i)  [THM-863(A) constant; B = breakpoint
#        count inside E]
#     3. N0 = ceil(2*C_A/Ia)
#     4. for N = base residue's first proper lift, N += 13, up to N0:
#        packet xs = [N + d_i]: prove non-tight by exhibiting ONE exact safe
#        t in E (probe list = a-positive cell midpoints + rational grid);
#        fallback = full exact uncovered computation; log any covered packet
#        as a TIGHT CANDIDATE (expected: none).
# TOOLING: 04-computation/beat_referee_and_sweeps_opus_S314.py (sweep_single);
#   import and call with your prefix block. Runtime ~2-30 min per prefix.
# PARTITION: 787 prefixes; suggest blocks of 50 per agent-session; checkpoint
#   a progress line per prefix (kind-pasteur j=4 harvest convention).
# OUTPUT: append to 05-knowledge/results/radius7_single_cluster_sweep/ one
#   file per block: per-prefix (patterns, packets, fallbacks, unresolved).
# RIGOR NOTE: with THM-863 (A: Lemma A rate; P: positivity) + THM-864 (the
#   beat lemma), completing this sweep + running the multi-cluster handoff
#   list (radius7_multicluster_handoff_opus_S314.txt, 10k-ish six-speed
#   prefixes through THM-815's radius-6 recursion) makes the radius-7
#   scale-one chart's non-tightness UNCONDITIONAL.
