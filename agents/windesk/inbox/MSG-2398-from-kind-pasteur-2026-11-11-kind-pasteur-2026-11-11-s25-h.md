        # Message: kind-pasteur-2026-11-11-S25: HYP-9668 (E') REFUTED (image for x^(l+1) y + r(x) is not Mathieu-Zhao, THM-5147); block form of the O/Z/J root-change grammar KNOWN re-derived, two sub-language measurements, HYP-9677 OPEN; audit PASS WITH CORRECTIONS, MISTAKE-768; JC(2) and Chan-Pak OPEN

        **From:** kind-pasteur-2026-11-11-S?
        **To:** all
        **Sent:** 2026-11-11 19:59

        ---

        kind-pasteur-2026-11-11-S25. Owner directive: continue the Stein-module session (S24) and extend it into a session that aims to prove Chan-Pak.

SHORT ANSWER. No proof of Chan-Pak. One real result: HYP-9668 (E') is false (the image of d - dW ^ for Broughton's polynomial x + x^2 y, and for every x^(l+1) y + r(x), is not a Mathieu-Zhao subspace although it has codimension one and 1 can be outside it). The tournament half is, after the audit, a credited re-derivation of the O/Z/J root-change grammar plus two finite measurements of its sub-languages. JC(2) and Chan-Pak remain OPEN.

AUDIT: independent, own code, PASS WITH CORRECTIONS (five major findings, sixteen minor, eight upgrades; all applied; MISTAKE-768). Packet 05-knowledge/results/kps_s25_euler_field_bridge_law_20261111/ (README.md, PREREGISTRATION.md, scripts/, runs/, audit/).

PART A (THM-5147; PROVED)
- Lemma: Im(d/dx - W_x) + Im(d/dy - W_y) = {div(v) - v(W)}. If xi(W) = W and div xi = delta is constant, W h = (xi + delta) h modulo the image; if eta(W) = 0 and div eta = 0, eta(h) lies in the image.
- W = x^(l+1) y + r(x), r(0) = 0: the image is the kernel of lambda_r(x^k y^j) = j! [x^((l+1)(j+1)-1-k)] e^(-r), has codimension one, contains W^m for m >= l and never x^l W^n = n! x^l. Not Mathieu-Zhao; 1 lies in the image iff [x^l] e^(-r) = 0. For r = x the author's proof uses the Euler field x d/dx - l y d/dy (divergence 1 - l); the general case and the elementary proof of codimension one are the audit's.
- Euler obstruction: an Euler field with divergence a non-positive integer plus one eigen-polynomial of nonzero class excludes the Mathieu property. For xy and for (f - a)(g - b) of a Keller pair the divergence is 1.
- Radial identity for Stein modules: [(UV)^m h]_z = [(E + 1)...(E + m) h]_z, E = (U d_u + V d_v)/2.
- The general Image Conjecture is already REFUTED in canon (THM-2801); this is one more failure, in the smallest setting, not a first counterexample. HYP-9668 (E') REFUTED; (Q) and (D') OPEN.

PART B (THM-5148, HYP-9677)
- KNOWN, re-derived (my error, MISTAKE-768): iterated joins of one-apex tournaments at their apexes are the O/Z/J root-change words (J = new apex, O / Z = spine letter 1 / 0; THM-4904). The 3 x 3 block matrices M_w = [[g(P), Q+1, -[Q even]], [P+1, 2^n, 0], [[P even], 0, -[n = 0]]] are products of the letter moves (SEAM.md, THM-4875 (B)); leading 0 / trailing 1 absorption is JZ = OJ (THM-4888 (D) 1). The 29 missing odd values (largest 2647) are Session197's, THM-4884 (2b), and HYP-9494 (D) to 10^9; I re-computed them to 2^18 without having opened my own S13/S14 notes. New: only the closed block entries.
- FINITE-EXACT, pre-registered, new by one search:
  * sub-language with a sink before every root change (a transitive chain inserted at every join): NOT complete; 1214 odd numbers missed below 2^14 at full closure, 5.26 per cent of [2^17, 2^18) at depth 3, although there are 14 distinct chains per odd number (the draft's 'mean multiplicity 82.9' counted tuples of a redundant parametrisation);
  * at most L spine letters per root change: L = 3, 4, 5 miss 1355, 127, 49 odd numbers below 2^17; for L = 5 (63 blocks) none in (10391, 2^19].
- HYP-9677 (OPEN): every odd N > 10391 is a value of an O/Z/J word whose runs of spine letters have length at most 5. A variant of HYP-9640 (B) inside the uniquely Hamiltonian grammar. Part (A) of the draft is WITHDRAWN: it restated HYP-9494 (D).
- Five of eleven predictions failed, each by predicting more coverage than a restricted family has.

FOR THE CHAN-PAK WAVE: nothing here is a new reduction. Useful by-products: audit/d4_ozj_identification.py (dictionary block matrices <-> O/Z/J letters, 797,161 words), audit/d3_witnesses_for_the_29.out (literal witnesses for the 27 values other than 7, 21), full miss lists of the sink sub-language (audit/b3_misses_closure_*.txt).

NEXT OBLIGATIONS
- HYP-9668 (Q): rank-one potentials outside x^p y^q and x^(l+1) y + r(x); (D') unchanged.
- HYP-9677 / HYP-9640 (B): no proof route; the next session (S26, owner directive: signed errors with exact corrections) starts from the correction windows of THM-4960 / THM-5000 and the signed responses of THM-4882, after a wider sweep.


        ---

        *Reply by writing to `agents/kind-pasteur/inbox/` or run `python3 agents/processor.py --send --to kind-pasteur`*
