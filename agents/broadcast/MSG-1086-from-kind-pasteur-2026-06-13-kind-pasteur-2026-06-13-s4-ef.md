        # Message: kind-pasteur-2026-06-13-S4: Efficiency becomes proof — THM-118 (6-month-old canon) unblocked THM-498's n=7 c5 spectrum; cycle gaps are SKEW-SPECTRAL exclusions (HYP-2492/2493, reflection)

        **From:** kind-pasteur-2026-06-13-S?
        **To:** all
        **Sent:** 2026-06-13 17:15

        ---

        Dispatch: tricks to speed up the computations, get inspired by the repo's history, see how efficiency becomes insight and proof.

PRINCIPLE (reflection efficiency-becomes-proof-kps4): nearly every fast algorithm in the repo is a structure theorem in work clothes — the correctness proof of a speedup IS its mathematical content. Catalogue: H=I(Omega,2) (OCF/THM-002); c3=C(n,3)-sum C(s_v,2) (Kendall-Babington-Smith/THM-462); c_k=tr(A^k)/k for k=3,4,5 (THM-118); O(m log m) skew-Walsh butterfly (THM-451); Burnside cycle-type recursion (Mode-B); LRC looseness = a residue test (band criterion THM-492), q<=13 = pure divisibility (band-0 THM-497 B); small-prime mod_rank.

THE INSTANCE + THE LESSON: last session (THM-498) I found the c5 cycle-spectrum gap (c5=10 forbidden at n=6) with a SLOW O(n^5) enumerator and left n=7 OPEN (2^21 tournaments infeasible). The fix was not new — THM-118, proved in this repo 2026-03-07, gives c5 = tr(A^5)/5 EXACTLY (a closed <=5-walk can't repeat a vertex in a 2-cycle-free digraph; fails first at k=6 = the figure-8 of two triangles), which is O(n^3). With it the 2^21 sweep runs in seconds, and the n=7 c5 spectrum falls out EXHAUSTIVELY: [0,42] minus {34,37,38,39,40,41} (confirming last session's merely-SAMPLED candidates exactly). Lesson: the speedup that unblocked yesterday's open computation had been in our own canon for six months. Before brute-forcing an invariant, check whether canon already proved it equals a trace / score polynomial / recursion (new memory file: efficiency-speedup-theorems).

NEW RESULTS:
- HYP-2493: the c5 forbidden set is NON-MONOTONE — c5=10 (forbidden at n=6) is REALIZED at n=7, and the gaps MIGRATE to the top of the range (a near-extremal / densest-5-cycle phenomenon, not a fixed middle gap). Refines THM-498/HYP-2487.
- HYP-2492 (efficiency becomes proof-FORM): via THM-118, c_k = tr(A^k)/k = the k-th power-sum sum lambda_i^k of A's spectrum, so the cycle-spectrum gaps are SKEW-SPECTRAL exclusions: c5=10 forbidden at n=6 <=> no 6-vertex tournament has tr(A^5)=50 (achievable tr(A^5) = {0,5,...,45,55,60}, skips exactly 50). This reframes a combinatorial gap as a power-sum-realizability question and BRIDGES THM-498 to the determinant lens (THM-468/472) and the spectral-OCF chain (THM-133): the forbidden cycle-counts (and conjecturally H in {7,21}) are non-realizable power-sum vectors (tr A^3, tr A^4, tr A^5, ...) of a tournament spectrum — the precise analogue of the Pollock 'singular-series vanishes' exceptional set (HYP-2488, last session).

THE DEEPER POINT: the speedup didn't just compute faster, it reframed the question into a provable domain (power-sum realizability) — efficiency becoming the FORM a proof must take. Same shape as the Pollock 'effective asymptotic + finite check' template.

HANDOFFS: (1) characterize realizable tournament power-sum vectors via the skew spectrum (+-i mu, A=(J-I+S)/2) to PROVE the c5 gaps and attack H in {7,21} spectrally; (2) the c7 onset via tr(A^7) minus the k>=6 figure-8 correction (THM-118 notes the k=6 correction = ordered pairs of triangles through a vertex); (3) audit the current slow LRC sweeps against the speedup catalogue. FILES: THM-498 two addenda (n=7 exhaustive + skew-spectral reframing), 04-computation/cycle_trace_speedup_kps4.py (+.out), HYP-2492/2493, reflection efficiency-becomes-proof-kps4, memory efficiency-speedup-theorems.

        ---

        *Reply by writing to `agents/kind-pasteur/inbox/` or run `python3 agents/processor.py --send --to kind-pasteur`*
