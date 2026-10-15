        # Message: mac-mini-2026-10-14-wilson: close-out - K6-P6 Paley terminal = Legendre bilinear form (THM-4653) ~ Sym^2 of 32.3.c.a (HYP-9276, descent form); safe Mahler codes: bounded entry + 2-adic memory (THM-4654); Chowla 21/20-700 ledger; D19/Giuga/HN audited (MISTAKE-610..612)

        **From:** mac-mini-2026-10-14-S?
        **To:** all
        **Sent:** 2026-10-14 19:48

        ---

        mac-mini-2026-10-14-wilson close-out.

The owner gave three briefs:
1. The repair algebra at -1 and the Wilson/{2,4,8}/Fermat dictionary.
2. Double length: Guy D19, Giuga/Korselt and parity, CH/Bott/RH, and progress on Collatz/HN/JC2 where possible.
3. A specialised arithmetic evaluation of K6 - P6; stronger restrictions on the irrational Collatz sources selected by positive Mahler integers; and how 21/20 and 700 deploy in final.pdf (OpenAI, "Ordinary two-point correlations of multiplicative functions").

Every theorem below was independently audited:
- audit A for the first half (MISTAKE-606);
- audit B for D19/Giuga/HN (MISTAKE-610);
- audit C for K6-P6 (MISTAKE-611);
- audit D for THM-4654 and the PDF ledger (MISTAKE-612).

PART 1 (audited)
- THM-4640: the repair algebra at -1.
  - Collisions are increasing {2,3}-representations.
  - Primitive collisions are ladders.
  - The Gersonides carries {2,4,8} are not generators.
  - -1 switches cover almost every trailing-ones source (c(E) -> 3/8).
- THM-4641: the Wilson dictionary.
- HYP-9270: N(E) ~ 2^E E^(-0.47).

PART 2 (audited, MISTAKE-610)
- THM-4650, D19 special loci.
  - No rational-distance point lies on the sides, midlines, diagonals or circumcircle.
  - The sides, midlines and diagonals are KNOWN (Yang Ji 2021). New here: the circumcircle via Ptolemy, and the curves.
  - The midline curve is the Velu 2-isogeny image of the congruent-number curve y^2 = x^3 - x.
  - Rank 0 is unconditional (L(E,1) != 0).
  - The rectangle family and rank table; torsion Z/2 x Z/4 at 3/2, 5/6, 8/3.
- THM-4651: Giuga's sum at p^2 is the Wilson quotient (KNOWN; Volpatti's list 1277, 780887, 56151923, 11774118061).
  - The parity divide: a counterexample must be an odd Giuga number.
  - Only N(N-1) is Giuga; N(N+1) is the next primary pseudoperfect number.
- Synthesis note, with attributions corrected.
- HN.
  - Exact spindle-free rebuild: 64,513 vertices.
  - 48,121-vertex M3 u psi P1 (CaDiCaL; independently kissat404, 8,786 s).
  - Clipped pieces: 47,545 and 46,585 vertices, still 5-chromatic (CaDiCaL).
  - Still running after close-out: the |p|^2 <= 2.5 clip (45,409 vertices) and a queued M3-side clip, spindle_free_m3side.py. Their logs are in hn_records_20261014/ and still to be committed.
  - Only 72 of the 120 hinge edges are at radius 2.
  - Heule's 1,441 remains far away.

PART 3
- THM-4653 (audited CORRECT): the Paley terminal V_q = I_(K6-P6)(q) = q(q-1)(-eps M_q + 7 + 13 eps - 2 eps q) for every odd q.
  - M_q = -sum chi(a x (x-1)(x-a)) t(a) t(x) = -q^3 sum phi(x)phi(1-x) 2F1(x) 3F2(x).
  - It is evaluated by FFT in O(q log q).
- HYP-9276 (FINITE-EXACT, independently): M_q = -b(q^2) - 2q - 1.
  - g = LMFDB 32.3.c.a = eta(4z)^2 eta(8z)^2 (theta(4z)^2 + i(theta(2z)^2 - theta(4z)^2)), McCarthy-Papanikolas' f_2.
  - At Paley primes: V_q = q(q-1)(|b(q)|^2 - q^2 - 7).
  - Checked at 10,590 primes to 3e5 and 22 prime powers (Frobenius trace on Sym^2 rho_g).
  - So the degree-10 run coefficient carries a GL(3) (Sym^2 = Ad x chi_-4) coordinate beside the GL(2) Apery wheel.
  - PROVED: M_q == 2 mod 8 (q = 1 mod 4) and 16 | M_q (q = 3 mod 4).
  - Closest prior statement: McCarthy-Papanikolas section 6 (F_(p^2), numerical).
  - Four proof routes are listed.
  - New, modular-form-free descent form: q^3 sum phi(x)phi(1-x) 2F1(x)3F2(x) = H_(q^2) + chi_-4(q) q^2 + 2q + 1, with H_Q = McCarthy's 3F2*(chi4,phi,phi;eps,eps|1) over F_(q^2). Checked directly at q = 5..43.
  - McCarthy-Papanikolas Thm 1.8 extends to 20 prime powers Q = 1 mod 4, including the inert squares p^2 with p = 3 mod 4. Expected route: Katz hypergeometric sheaves over Z[i][1/2] + Chebotarev.
  - Geometric guess: the K6-P6 4-fold is Res_(Q(i)/Q) of the 3F2(1) surface.
- THM-4654 (Mahler/Collatz; audited, MISTAKE-612): safe 3/2 codes.
  - Slack automaton m -> min(3m/2 - e, 1), which is THM-4072's follower graph.
  - Block law: no 11; a gap 2 forces a gap >= 6.
  - At most floor(l/3)+1 ones.
  - The Collatz maps contract along safe codes: slopes <= 3(3^(1/3)/2)^l, offsets < 11/10. So rational shadows ENTER a bounded window.
  - Their MEMORY is 2-adic (Terras: floor(log2(2.1v))+1 future bits). The pre-audit 'memory by contraction' lemma was vacuous.
  - Complexity-approximation inequality, unconditionally linear in v and log in |u| (interface: 14Q). It is within 1 bit of tight on data; 0 violations on 4,980,199 pairs and on 2.1M audit cases.
  - Linear complexity gives a finite 2-adic exponent for z_A. That is the generic case, NOT a restriction.
  - Sturmian codes are impossible (Dubickas 2009).
- Chowla PDF ledger.
  - 1/700 is non-binding: Lemma 4.4 needs A >= 560; the paper takes A >= 1000; c = delta/(6WA).
  - 21/20 = 1 + 1/20 pays for the L/eta bins, with a margin ~220x the need; the balanced optimum is 23/20.
  - 6/5 forces the modulus 5.
  - Part II: MRT's e^(-M/20) and (log X)^(-1/700) are used only qualitatively.
  - rho = 1/20 needs 0.0008 < rho < 0.125.
  - Lemma 6.5 has a second exp(-L^(21/20)).
  - No 20 or 700 in the paper is sharp.
  - (4.3) is the characteristic polynomial of U_(3,4).
  - Links:
    - nonbacktracking operators <-> K6-P6 (skeleton only);
    - singleton centering <-> isolated-edge vanishing (generic);
    - prohibited words <-> the safe Mahler language (ANALOGY).
  - COINCIDENCES: the constants 20/21 in K6-P6; mod 5.

NEXT OBLIGATIONS
- Prove HYP-9276 via its descent form (a pure finite-field hypergeometric identity across F_q in F_(q^2)), plus the Katz/Chebotarev extension of MP Thm 1.8. Alternatives: Hijikata; Faltings-Serre after a Katz rank computation.
- HN: phase-1 random-half cores on M3 u psi P1, aiming below 1,441, plus DRAT certificates.
- Mahler: a super-linear factor-complexity bound for Mahler codes of positive integers (Dubickas gives >= 1.71 l), or a Mahler-side argument for low-complexity codes.
- Paley: the twelve-step terminals (opus S35, HYP-9275) are likely GL(2) x GL(3) or Sym^3 objects.
- Wilson completeness is expected FALSE; Giuga is open; JC2 is verified to degree 108.
- The environment is the local mac-mini, not a cloud machine.


        ---

        *Reply by writing to `agents/mac-mini/inbox/` or run `python3 agents/processor.py --send --to mac-mini`*
