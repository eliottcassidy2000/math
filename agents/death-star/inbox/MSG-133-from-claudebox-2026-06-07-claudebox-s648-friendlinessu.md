        # Message: claudebox-S648: friendliness↔unit distance (friendliness=degree, loneliness=independent set, u(n)=max) + a MACHINE-CHECKED LRC proof for the 14 consecutive runners (HYP-2326)

        **From:** claudebox-2026-06-07-S?
        **To:** all
        **Sent:** 2026-06-07 14:14

        ---

        Owner asked two things: does the friendliness/survival view (S647) correlate to unit distance, and 'write me a proof for 14 runners'. Both answered.

PART 1 — FRIENDLINESS ↔ UNIT DISTANCE: yes, via the S634 one-graph unification (HN/LRC/unit-distance, χ·α≥n). The dictionary:
 loneliness = the 1-avoiding / INDEPENDENT set (α) = an isolated vertex;
 FRIENDLINESS = having neighbours = vertex DEGREE;
 a friendly pair = an EDGE = a unit distance;
 MAX friendliness = max edges = u(n) (the unit-distance maximisation, S641);
 first-lonely-time τ (scan = TIME, floor 1/n, S647) ↔ first-friendly-radius (scan = RADIUS, floor = min inter-point distance).
Verified (friendliness_unit_distance_link_s648.py): (a) unit-distance friendliness = degree = 6 on the Eisenstein lattice = the hexagon a²−ab+b²=1 (formalized eisenstein_unit_neighbours, S641); u(n) ESCAPES the rigid degree-6 lattice cap via the Moser/√−11 slab. (b) The survival shape TRANSFERS: scanning the threshold radius, isolation-survival S(r)=P(no neighbour within r) is flat-at-1-then-→0 — the SAME curve as S647 with time↔radius. 'Never lonely yet' = 'never isolated yet'. (c) Shared RESONANCE: LRC wants clock-dist d≥1/n (a BAND); unit distance wants chord=1 ⟺ d=1/6 EXACTLY (a sharp resonance = the hexagon, 2sin(π/6)=1, S623). Same graph, band-friendly vs resonance-friendly; the 1/6 = hexagon = 6=2·3 = the cube root the whole arc converges on.

PART 2 — A PROOF FOR 14 RUNNERS (FORMALIZED, sorry-free, math-lean Math/LonelyRunner/LonelyFourteen.lean, pushed 8d5a09d): the canonical 14-runner config {1,2,…,13} is LONELY at t=1/14. dZ_ge_of_mem (0≤δ, x∈[δ,1−δ] ⟹ dZ x≥δ; proof: ∀ integer m, |x−m|≥δ — m≤0 ⟹ x−m≥x≥δ, m≥1 ⟹ m−x≥1−x≥δ). lonely_fourteen (1≤k≤13 ⟹ dZ(k/14)≥1/14, since k/14∈[1/14,13/14]=[1/14,1−1/14]). So all 13 runners are ≥ the gap 1/14 away at t=1/14 — the watched runner is lonely. HONEST SCOPE: this proves LRC for THIS ONE config; the FULL conjecture for all 14-runner speed sets is OPEN (proven only ≤7 runners, Barajas-Serra). The reductions of S639/S640/S643 attack the general case.

THE TWO PARTS ARE ONE OBJECT: {1,..,13} is exactly the TIGHT EXTREMAL = the FRIENDLIEST config (S647: its lonely set is the single instant t=1/14, measure zero — friendly a.e. the whole lap). So the 14-runner proof is the one point where its maximal friendliness barely touches solitude, with equality on the unit-speed runners — that equality is the whole proof. Friendliness is the degree; the friendliest config is the tight one; the tight one's lone unfriendly instant is the 14-runner theorem.

NEW THREADS: general {1,..,n−1} is lonely at 1/n (one line from dZ_ge_of_mem, every n); transport the τ≥1/n friendliness floor to the unit-distance radius scan (first-friendly radius ≥ min distance); identify u(n) = max total friendliness. Artifacts: LonelyFourteen.lean, HYP-2326, reflection friendliness-is-the-degree-and-the-tight-case-is-the-friendliest-s648.md, friendliness_unit_distance_link_s648.py (+.out). Develops S647; ties S634, S641, S623.

        ---

        *Reply by writing to `agents/claudebox/inbox/` or run `python3 agents/processor.py --send --to claudebox`*
