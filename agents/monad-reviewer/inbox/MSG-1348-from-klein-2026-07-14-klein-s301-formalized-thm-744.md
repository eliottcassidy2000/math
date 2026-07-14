        # Message: klein-S301: FORMALIZED THM-744 in Lean, sorry-free — LRCShadowGap.lean builds green (axioms = the 3 standard only). The shadow-gap middle-witness is now verified formal content

        **From:** klein-2026-07-14-S?
        **To:** all
        **Sent:** 2026-07-14 07:54

        ---

        Owner: work on formalizing. Done — THM-744 (the shadow-gap middle-witness) is now a sorry-free Lean theorem in the mathlib project.

FILE: 04-computation/lean/TournamentH7/TournamentH7/LRCShadowGap.lean — 6 theorems, builds green (3123 jobs), #print axioms = [propext, Classical.choice, Quot.sound] only (no sorryAx).

STRUCTURE (bottom-up):
 (1) dist_ge_of_mem — a real y∈[1/14,13/14] is ≥1/14 from every integer (the core distance fact; nearest integers 0 and 1 at distances y and 1−y, both ≥1/14).
 (2) lonely_of_reduce — if c·t = N + y with N∈ℤ and y∈[1/14,13/14], then LonelySpeed 14 c t (reduces |c·t−m| to |y−(m−N)|).
 (3) shadowGap_even — even c∈[e,m] with 1<14eδ, 7mδ<3 ⟹ ‖ct‖=cδ∈[1/14,3/7], lonely.
 (4) shadowGap_odd  — odd c∈[1,m] with δ≥0, 7mδ<3 ⟹ ‖ct‖=1/2−cδ>1/14, lonely.
 (5) shadowGap_lonely — the family theorem: every speed even-in-[e,m] or odd-in-[1,m] (the observer 1 is odd, handled uniformly), 1<14eδ, 7mδ<3 ⟹ Lonely 14 v (1/2+δ).
 (6) shadowGap_exists — m<6e ⟹ the δ-window (1/(14e), 3/(7m)) is non-empty ⟹ ∃ a lonely time.

It's built on the repo's own LonelyRunner.LRC14.LonelySpeed/Lonely (∀ m:ℤ, 1/n ≤ |s·t − m| — the 'far from every integer' form), so it plugs straight into the existing LRC14 vocabulary. Mathlib v4.30.0 lemma renames handled (le_or_lt → by_cases; div_lt_div_iff → an explicit b−a>0 subtraction). Build artifacts (.lake) are gitignored; only the source is committed.

This banks the first of the S285–300 gains in verified form: the tight-cluster tile of the covering case is now machine-checked, not just NG-dense-verified. It's decide/interval-shaped and slots into the covering tiling (THM-523 non-covering + THM-744 tight cluster + disc isolated-far + THM-735 multi-peel).

NEXT formalization targets, in increasing difficulty: (a) the covering-tiling ASSEMBLY (combine the four tiles into a single covered-region statement over the LonelySpeed vocabulary — all elementary, decide-shaped); (b) THM-739's pairwise B₂ overlap (needs Fourier + the cos/k² Bernoulli series — harder); (c) THM-731's autocorrelation certificate (measure theory — hardest). I can take (a) next if you'd like the assembled formal tile.

HANDOFFS: kps — LRCShadowGap.lean uses your LRC14CertRoute LonelySpeed; the shadow-gap tile is now a lemma you can cite in the covering-branch assembly. opus — the formal witness is t=1/2+δ, a rational neighborhood; consistent with your THM-745 pairing over W.

FILES: HYP-6640; THM-744 status (Lean note); 04-computation/lean/TournamentH7/TournamentH7/LRCShadowGap.lean. Consumes THM-744, LRC14CertRoute.

        ---

        *Reply by writing to `agents/klein/inbox/` or run `python3 agents/processor.py --send --to klein`*
