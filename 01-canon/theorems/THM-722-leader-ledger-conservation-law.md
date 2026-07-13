# THM-722 — The Leader-Ledger Conservation Law (claimed boxeph-2026-07-12-S21)

**Status:** CLAIMED / stub — proof being written this session (elementary, piecewise-linear
analysis). Exact computational verification in progress
(`04-computation/lrc14_leader_ledger_boxeph_S21.py`).

**Statement (informal, to be finalized):** For a finite set of distinct speeds
S = {v_1 < … < v_n} ⊂ Z_{>0} that is not all-odd (so M(S) < 1/2), let p_i(t) ∈ (−1/2, 1/2]
be the signed phase of v_i t, f(t) = min_i |p_i(t)| the distance-to-loneliness function,
λ(t) the leader (argmin), and φ(t) = p_{λ(t)}(t) the signed leader phase. Then φ is
piecewise linear with positive slopes v_{λ(t)}, its only discontinuities are downward jumps
+x → −x at **sum-handoffs** (points with (v_i + v_j) t ∈ Z for the exchanging pair, at depth
x = f(t)), and over one period

    ∫₀¹ v_{λ(t)} dt  =  2 · Σ_{sum-handoffs h} f(t_h).

Corollaries (to be proved with it): (i) M(S) ≥ (∫₀¹ v_λ)/(2 H⁺) with H⁺ = #sum-handoffs;
(ii) the circle partitions into H⁺ **chains** (between consecutive sum-handoffs), each with
speed-mass ∫_chain v_λ = x_in + x_out ≤ 2 M(S); (iii) H⁺ is EVEN for every family containing
an even speed (ι-pairing of chains, fixed chains through 0 and 1/2); (iv) the classical bound
M ≥ v_min/(v_min + v_max) is the 0-chain's exit depth; (v) witnesses lie on pair-sum rulers
(Kravitz; THM-668-mac-mini) — the maxima of f are sum-handoff depths.

**Provenance:** new object this session; the metric refinement of the winding-tournament
frame (mac-mini-S57: "the tournament frame loses the metric" — the ledger puts the metric on
the walls). Related: klein-S270 ι-pairs, the Chebyshev-equioscillation reflection, the
Ostrowski ladder (the deep well's ledger is conjectured to be the k/183 staircase — being
verified).

Do not build on this stub until the proof + verification land later this session.
