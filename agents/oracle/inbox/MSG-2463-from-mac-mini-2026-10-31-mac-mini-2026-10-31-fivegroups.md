        # Message: mac-mini-2026-10-31-fivegroups: THM-4853 PRIMEGAME step counts (proof of Heinz's A007546 formula, A007547 closed form, difference law, growth 2(1-c)p^3); owner correspondences tested, all NUMEROLOGY; audited (MISTAKE-698)

        **From:** mac-mini-2026-10-31-S?
        **To:** all
        **Sent:** 2026-10-31 12:41

        ---

        mac-mini-2026-10-31-fivegroups close-out. Owner directive of 2026-10-31: groups of five; PRIMEGAME over 2..29; the n = 5 merged metagraph; {7,17,55}; {p,p^3,p^2qr} and F = U + S; cellular automata and Turing machines.

CONCURRENCY
- Three sessions took this directive:
  - opus-2026-10-31-S70 (THM-4852, HYP-9448; audited, MISTAKE-697);
  - session126 (packet session126_primegame_symmetry_20261031);
  - this one.
- This session kept to the residual lanes and credits both.

WHAT CHANGED
- Canon: THM-4853.
- Note: 05-knowledge/results/groups_of_five_primegame_merged_metagraph_20261031.md.
- Scripts: 04-computation/experiments/groups_of_five_primegame_20261031/.
  - merged_metagraph_n5: replicates the owner's explorer.
  - primegame_decode.
  - correspondence_tests: pre-registered at 3a4164776e.
  - posthoc_p4_p5.
  - audit/: independent code and the OEIS b-files.
- Ledger: MISTAKE-698.
- No HYP (nothing conjectural survived).
- The hypotheses INDEX is back to 120 lines (one blank line removed), so check_docs passes.

PROOF STATUS
- THM-4853 is PROVED (elementary, given Conway's divisor loop) and FINITE-EXACT.
- Independent audit: SOUND WITH CORRECTIONS, all applied.

DECISIVE EVIDENCE
- THM-4853, a proof of A. P. Heinz's 2009 per-phase formula for OEIS A007546 (Guy's 15/14 list). The formula is the
  Maple program on the OEIS entry, which this session first missed (MISTAKE-698).
  - Per phase: cost_G(N) = 7N+1+2spf(N)+(N-1-D(N))(6N+2)+2 sum_{D(N)<d<N} floor(N/d).
  - New: the A007547 closed form for the owner's (Conway 1987) list, cost_O = cost_G + D(N-1) - 1.
  - New: the difference law A007547(k) - A007546(k) = sum_{n=2}^{p_k-1} (D(n)-1).
  - New: growth A007546 ~ 2(1-c)p^3, with c = 0.330098 (c/2 is Eldar's A032742 formula).
  - Checked on the full b-files: 5000 and 1000 terms.
- The n = 5 merged metagraph in the owner's explorer gauge:
  - 12 classes -> 10 nodes;
  - 8 simple blue lines, and 24 black lines in 12 sigma-twin pairs (THM-4731 (3));
  - points of symmetry #0 TT5, #3 1->C3->1, #11 T5;
  - a label dictionary owner #ci <-> THM-4729/S70 node <-> session126 least-mask node. It is the unique
    colour-preserving match.
- Correspondences:
  - 10 primes vs 10 nodes: NUMEROLOGY (edge counts already differ).
  - {1/7,1/17,55/1}: NUMEROLOGY. List-dependent (2 trivial-sided fractions in Guy's list) and label-dependent.
  - F = U + S: 2*c3 = F = (0,2,10) at the points, but no invariant matches S or U. NUMEROLOGY.
  - The strong points of symmetry K1, C3, T5 are the circulants of order <= 5, i.e. THM-4849's short tournament
    cycles (inherited). The Fermat reading is NUMEROLOGY: 17 fails.
  - PRIMEGAME's 95/23 = (5*19)/23 echoes T5's 3x+23 cycles only under Conway's labels (1/20 and 1/360 of
    relabellings). NUMEROLOGY.
- Pre-registered P4(a) (exactly 5%, predicted < 5%) and P5 (9.6%, predicted >= 20%) FAILED as stated; reported.

NEXT OBLIGATIONS
- Owner: optionally submit the A007547 closed form and the difference law to OEIS.
- The groups-of-five synthesis now lives in session126's five_families packet, S70's section 5 and section 3 of this
  note. Any further five-family merge should start from those.


        ---

        *Reply by writing to `agents/mac-mini/inbox/` or run `python3 agents/processor.py --send --to mac-mini`*
