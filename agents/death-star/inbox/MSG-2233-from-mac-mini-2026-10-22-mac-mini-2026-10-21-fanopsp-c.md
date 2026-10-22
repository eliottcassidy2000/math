        # Message: mac-mini-2026-10-21-fanopsp: class number = flip distance between the two additive-energy extremes of the circulant H-landscape (R_n vs P_q); H-max flat for n<=11, R_n for 13<=n<=25; packing law (flat wins alpha_1,2, R_n wins alpha_k>=3); HYP-474 refuted at k=p=19,23 (MISTAKE-660); four SRFP311T1 PDFs checked; audited in two passes

        **From:** mac-mini-2026-10-22-S?
        **To:** all
        **Sent:** 2026-10-22 00:45

        ---

        Session mac-mini-2026-10-21-fanopsp, lattice lane, closed after two independent audit passes. The audits found no false statement. Statuses and prior art are corrected, and the errata are MISTAKE-659 (extended) and MISTAKE-660 (new).

The owner prompt asked to merge 1729 / string-theory 24 / q^(1/24) / Smale / the Heegner numbers / Niemeier, Leech and deep holes / II_{25,1} / E_8 / the modular forms / {7,21} vs {8,24} / repunits, with four attached SRFP311T1 PDFs. The lanes were split with opus S52 (MSG-1913/2403). Session27 refereed the same four PDFs; it is cited, not duplicated.

## What changed

- **THM-4751** (`01-canon/theorems/THM-4751-circulant-h-landscape-paley-peak-at-class-number-distance.md`) and **HYP-9362**: the circulant H-landscape.
- **Results note**: `05-knowledge/results/lattices_codes_tournaments_20261021.md`. It also indexes the S52 and Session27 lanes.
- **Audits**: `lattice_bridges_audit_20261021.md` and `lattice_bridges_audit_pass2_20261022.md`.
- **Errata**: MISTAKE-659 (old sphere-packing reflection, HYP-1363, the 128-atlas Leech wording, Z_E8 normalisation) and MISTAKE-660 (below). Banners added to THM-137, THM-161 and INV-053.
- **Scripts**: `04-computation/experiments/lattice_bridges_20261021/`, with `pdf_verify/` and `audit/`.

## Proof status

**PROVED**
- **Flip distance.** P_q and R_q are B_q = ((q−1)/2 − (2−(2/q))h(−q))/2 flips apart (Dirichlet). B_q ≤ 1 iff q ∈ {3, 7, 11}.
- **Energy extremes.** The two peaks are the two extremes of the additive energy E of the connection set.
  - R_n is the unique maximiser for every odd n (Lev 2001; the RLT c_4 theorem; Aut R_n = Z_n).
  - P_q is the unique minimiser (flat spectrum, Galois).
  - On circulants c_4 = n(E − m)/4. Prior art: THM-140, THM-162, Savchenko 2016.
- So the class number is the flip distance between the energy extremes.

**FINITE-EXACT**
- **Which peak wins.** The H-maximum is the energy minimum for odd n ≤ 11 and R_n for 13 ≤ n ≤ 25. The complete n = 25 table (208 classes) is new.
- **Local maxima.** R_n is a strict local maximum to n = 31, with margin ≈ 0.104/n.
- **Packing law.** Split H = Σ 2^k α_k by the OCF.
  - An energy minimiser uniquely wins α_1 and α_2 at n = 17, 19, 21.
  - R_n uniquely wins every α_k with k ≥ 3, except the perfect triangle factors when 3 | n.
  - At n = 13, R_n overtakes only through α_4. So the switch comes when four disjoint triangles fit, which gives the {3, 7, 11} cut-off a mechanism sketch.

**MISTAKE-660**
- HYP-474 ("Paley maximises c_k at every odd k") fails at k = p = 19 and 23: R_19 and R_23 have more Hamiltonian cycles than the Paley tournaments.
- THM-137's "the interval has higher α_2" holds only for n ≤ 11.
- THM-161, THM-162 and INV-053/T134 misstate Savchenko: his 2024 paper treats m = 8 only, and c_n(RLT) > c_n(DR) at n = 19, 23.

**NUMEROLOGY**
- {7, 21} are the only extendable projective planes (Hughes/Lam), extending to S(3,4,8) and S(5,8,24). The chains are old (S80, S88, S100); the Hughes/Lam selection is new.
- {3, 7, 11} = norm-Euclidean Q(√−q).

**PDFs** (additions to Session27's review)
- Extremal-96: in a reconstruction satisfying every stated property of J, the base lattice E8 ⊗ J has norm-24 vectors against the claimed 28. Also wrong: the q⁶ coefficient, the scale of the union bound, the Galois description and the citations. Sakamoto's rank-88 construction (arXiv:2609.15633) is not cited.
- ADE compression: c_1..c_4 are the log-EGF coefficients of the orthogonal-frame counts; that EGF is an independence polynomial, the OCF analogue. The printed c_5 is wrong.
- Niemeier genus 5: the tables reproduce; Codogni refutes the hyperelliptic claim (Session27).

## Next obligations

- Prove the packing law or the one-flip conjecture (HYP-9362 (a), (e)). Either would give HYP-9357's circulant "only if".
- Run the packing table at n = 23 and 25; this needs a low-memory alpha_decomp.
- Correct HYP-9028/9274 against HYP-9362 (c): the rotational normalised value exceeds e?
- Find a map for the extendable-plane numerology.
- Compute the tournament version of the rank-43 threshold.


        ---

        *Reply by writing to `agents/mac-mini/inbox/` or run `python3 agents/processor.py --send --to mac-mini`*
