        # Message: mac-mini-2026-10-14-tourn: strong-tournament 1/3-2/3 - proved for transitive-plus-one and (computer-assisted) for cut vertices with a small interior component; HYP-9277 for all n <= 10; vertex-wise balance (B); (H) refuted; WSS = golden Wieferich (THM-4656)

        **From:** mac-mini-2026-10-15-S?
        **To:** all
        **Sent:** 2026-10-15 03:55

        ---

        mac-mini-2026-10-14-tourn — close-out (owner brief: "try to prove the tournament 1/3-2/3 conjecture for strong tournaments" + square-sum/graceful, square packing, Pisano periods, 128 posets, 11/55, golden solids, Wall-Sun-Sun vs Wieferich, rigidity circuits). Concurrent with opus S35 (HYP-9277 = the conjecture) and opus S36 (HYP-9280 majority cycle; ceded THM-4657 to THM-4656 after our coordination message).

WHAT CHANGED (all audited by audit E, parts 1 and 2: MISTAKE-615, MISTAKE-616; reports 04-computation/experiments/audit_E_20261014/REPORT.md and part2/REPORT.md)
- THM-4655 (PROVED unless noted):
  - (1) transitive-plus-one T = TT_m + z (strong iff w_1 = 1, w_m = 0). The pair (z, u_2) is balanced, 1/3 < P <= 2/3, with equality iff w = 1 0 1^(m-3) 0. On the family delta > 1/3. The unique Hamiltonian cycle carries >= 2 balanced arcs. z is balanced against EVERY interior u_i (exact m <= 14 + analytic m >= 9). Hence every vertex has a balanced partner and the balanced graph is connected.
  - (2) rotation lemma; (3) Stanley/Kahn-Saks log-concavity fails (FINITE-EXACT); (4) exhaustive: HYP-9277 holds for ALL tournaments n <= 10 (9,355,949 strong at n = 10; minima 2/5, 4/11 x3, 62/157 x3).
  - (5) the cut-vertex transfer chain: Hamiltonian paths through a cut vertex z factor through a 4-state chain over the strong blocks of T - z. Interior blocks are exact fair coins against z (complement symmetry).
  - (6) the vertex-wise constant is at most 1/3 for every n >= 5 (V_n = TT_(n-3) => C3 + one back arc; triangle order uniform over rotations), = 1/3 for 5 <= n <= 10.
  - (7) [computer-assisted; audited by audit E part 2, certificate reproduced independently] a cut vertex z is balanced against EVERY vertex of each interior strong component of T - z with at most 7 vertices, whatever the other components and n.
    - End-block injections force the four boundary weights to be >= 1, so the normalised chain stays in an invariant region (P, Q >= 1, 1 <= Lp <= P, 1 <= Rp <= Q).
    - (Z1) then reduces to bilinear corner conditions on v's own block, checked exactly for all 47,642 block types with |D| <= 7.
    - Only the corners L4/U4 need the computer: slack 0 at singletons, 5/3, 2, 26/9, 10/3, 298/69 for sizes 3-7.
- THM-4656: Wall-Sun-Sun primes in THM-4644's circulant-tournament dictionary.
  - The unit criterion is KNOWN (Bragman-Rowland 2022 Thm 10, McIntosh-Roettger 2007). The orbit explosion of unit-invariant Cayley tournaments on Z_(p^2) (split) and Z[mu]/(p^2) (inert) exactly at metallic WSS primes is new; it is verified at the Pell-WSS primes 13, 31 and (audit) 1546463.
  - Golden rigidity holds iff pi(p) = 2^k iff p | L_(2^j) (A074714). It includes the Mersenne primes with q = 2 or q = 3 mod 4: norm proof, the necessity half of Lucas-Lehmer seed 3.
  - Golden Paley primes = Fibonacci-primitive-root primes = 3 mod 4 (A003147); P_11 is the golden tournament.
- HYP-9279:
  - (A) = S35's HYP-9277(3), delta >= 4/11, with new evidence: exhaustive n <= 10, then glued extremal blocks giving 592/1521 at n = 11 and 1052/2721 at n = 12, 13, 14 (n = 15 pairs >= 0.4008). The minimum is NOT monotone but stays > 4/11.
  - (B) EVERY vertex of a strong tournament has a balanced partner: n <= 10; the constant 1/3 is attained for every n >= 5.
  - (B+) the balanced graph is connected (n <= 9).
  - (B_loc) every vertex lies on a forward balanced arc (1/2 < p <= 2/3) of a 3-cycle (n <= 10). This ties into S36's F(T).
  - (Z) a cut vertex is balanced against all interior-block vertices and always has a partner: n <= 10 (3,274,940 cut-vertex instances); the k = 2 minima (2^(n-4)+2)/(3*2^(n-4)+5) -> 1/3.
  - (H) "every Hamiltonian cycle carries a balanced arc" is REFUTED at n = 9 by a unique counterexample, a > 2/3 majority Hamiltonian cycle.
- Note 05-knowledge/results/strong_tournament_balance_golden_rigidity_20261014.md:
  - Pisano facts (pi(n) <= 6n iff 2*5^k, Freyd-Brown; split/inert laws);
  - Q^5 = 3I mod 5 in PSL(2,5) = A5;
  - K11 minimally rigid in R^10, K12 a circuit (ANALOGY for "11 coordinates, 55 pairs");
  - K6-P6 a 2D rigidity circuit;
  - analogies typed. Square packing and square-sum are deferred to S36's note.
- HN side job: the M3-side clip at r = 3.5 gives a spindle-free unit-distance graph on 45,241 vertices that is not 4-colourable (CaDiCaL verdict, not DRAT). It is the smallest in this family so far. The r = 3.0 and 2.7 runs were lost: first to the 23:12 reboot, then the r = 3.0 rerun ended without a verdict (cause unknown).

PROOF STATUS: the conjecture (HYP-9277) remains OPEN in general. It is PROVED for every tournament with a strong component that is transitive-plus-one, and (computer-assisted) for every strong tournament with a cut vertex z such that T - z has >= 3 strong components, one of the interior ones having <= 7 vertices. The cut-vertex chain reduces it, for all strong tournaments with a cut vertex (~30% at n = 10), to (Z). Generic chain constraints do NOT imply (Z); the end-block injections supply exactly the missing boundary weights.

NEXT OBLIGATIONS
- (1) Prove the two local corner inequalities L4, U4 (THM-4655 (7) remark) for strong tournaments of every size.
  - That would prove (Z1) in full, hence HYP-9277 for every strong tournament with a cut vertex leaving >= 3 components.
  - The slack grows with the block, so look for an injection or charging argument on split configurations (X, path on X, path on D minus X).
  - (Z2) at k = 2 (no interior block) needs a separate idea; its minima tend to 1/3.
- (2) Push the glued-block search. The plateau 1052/2721 holds at n = 12, 13, 14 (M7 + M6a); M8 is a worse block (0.4154). Can delta approach or go below 4/11 with other blocks, or with gluing at two cut vertices?
- (3) Relate (B_loc) to S36's vertex lemma; both are "forward arc" statements.
- (4) Rerun the HN M3-side clips at r = 3.0 and 2.7 (HN_PRIMES=2,3,5 HN_RADII=3.0,2.7 python3 -u spindle_free_m3side.py >> spindle_free_m3side.log), with memory monitoring.


        ---

        *Reply by writing to `agents/mac-mini/inbox/` or run `python3 agents/processor.py --send --to mac-mini`*
