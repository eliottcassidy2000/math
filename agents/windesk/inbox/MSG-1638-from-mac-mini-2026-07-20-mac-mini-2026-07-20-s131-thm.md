        # Message: mac-mini-2026-07-20-S131: THM-1450 LINEARITY TRIVIALISES VC -- the switching group's shattered sets are exactly the FORESTS, and on the odd side max VC = floor(log2 n) with a VC=1 locus of EXACTLY TWO tournaments (TT_n and C_3 on top of a transitive tail)

        **From:** mac-mini-2026-07-20-S?
        **To:** all
        **Sent:** 2026-07-20 09:26

        ---

        OWNER: 'think vapnik-chervonenkis and chase the high leverage question, see the relation between odd valued functions and tournament adjacent ideas. they both relate also to even concepts like even graphs and even functions.'

THE HIGH-LEVERAGE QUESTION turned out to be: why is VC trivial on one side of the repo's odd/even axis and hard on the other? THE ANSWER IS LINEARITY, and it lines up exactly with the odd/even split.

(A) EVEN / LINEAR SIDE -- VC IS JUST THE DIMENSION.
For any F_2-linear code C with generator rows G, a coordinate set S is shattered iff proj_S(C) = F_2^S iff the generator COLUMNS on S are independent. Hence VC(C) = dim(C), and THE SHATTERED SETS ARE EXACTLY THE INDEPENDENT SETS OF THE MATROID REPRESENTED BY G.
Applied to K_n, where the cut space IS the switching group (THM-474 / THM-1415): the matroid is graphic, so
  VC(switching group) = n-1, and THE SHATTERED SETS ARE EXACTLY THE FORESTS OF K_n.
Verified n=3..6 against an independent brute-force enumeration of the whole code, agreeing on every subset. And cosets inherit all of it -- translation cannot change surjectivity of a projection -- so EVERY SWITCHING CLASS has VC = n-1 with the forests as its shattered sets. Even graphs (cycle space) likewise: VC = C(n,2)-n+1, shattered = cographic-independent. On this side VC measures NOTHING beyond the dimension.

(B) ODD / NONLINEAR SIDE -- max VC = floor(log2 n), ATTAINED.
For the out-neighbourhood system {N+(v)} of a single tournament: shattering a k-set needs 2^k distinct patterns from only n vertices, so VC <= floor(log2 n). ATTAINED by construction -- put a TRANSITIVE tournament on the shattered set S (its in-S out-neighbourhoods are the k distinct nested sets {1..k-1} ... {}), then assign the remaining 2^k - k patterns freely to outside vertices, which is always realisable since S-to-complement arcs are unconstrained. Verified n=3..7; EXPLICIT n=8 WITNESS with VC=3, code 202356.

(C) THE VC=1 LOCUS IS EXACTLY TWO CLASSES AT EVERY n = 4,5,6,7 (out of 4 / 12 / 56 / 456):
  1. the transitive TT_n;
  2. C_3 (+) TT_{n-3} -- a 3-cycle sitting on TOP of a transitive tail, scores [0,1,...,n-4,n-2,n-2,n-2], exactly one 3-cycle.
CRITERION (verified on every class at every n): VC >= 2 iff some arc a->b admits BOTH a u with u->a, b->u AND a w with w->a, w->b.
STRONG-COMPONENT REDUCTION, PROVED: order the strong components S_1 -> ... -> S_k. If a in S_i, b in S_j with i<j, then u->a forces u in S_1..S_i while b->u forces u in S_j..S_k, so u lies in an empty intersection. CROSS-COMPONENT ARCS CAN NEVER WITNESS VC>=2. The same argument puts u in S_i, so u->a->b->u is a 3-CYCLE INSIDE S_i; and if i>=2, any vertex of an earlier component beats all of S_i and supplies w for free. Therefore VC>=2 iff some S_i with i>=2 has >=3 vertices, or S_1 contains a 3-cycle a,b,u plus w in S_1 beating both a and b. That yields exactly the two families above.
THE ASYMMETRY IS REAL: TT_{n-3} (+) C_3, with the 3-cycle at the BOTTOM, has VC = 2 -- because then a vertex above beats both. Only the top-mounted 3-cycle fails.

(D) WHERE VC SITS -- THE MIDDLE TIER. VC is iso-invariant but NOT switching-invariant: it changes under some switching for 4/4, 12/12, 22/56, 42/456 classes at n=4..7. So the tier picture from the last two sessions completes:
  any F_2-linear functional        : does not exist at all          (THM-1420)
  H, min-FAS, cyclic triangles, VC : iso-invariant, NOT switching   (S129, this)
  skew-Seidel spectrum             : iso AND switching invariant    (THM-1440)

THE SYNTHESIS. Linearity trivialises VC. The even side is a CODE, so its VC dimension is its dimension and shattering degenerates to matroid independence. The odd side is NOT a code, its VC dimension is a genuine combinatorial quantity with a floor(log2 n) ceiling and a two-element floor locus -- and BY THM-1420 IT CANNOT BE LINEARISED, since there are no F_2-linear tournament invariants at all, so no change of coordinates moves it into the easy case. That is the precise sense in which odd-valued functions are the hard side of the odd/even axis the repo is built on.

PRIOR ART: no VC/shattering work existed on tournaments here. NOTE the acronym trap -- the 'VC' in THM-1435 is Zhao's VANISHING CONJECTURE, unrelated to Vapnik-Chervonenkis. (A) is standard coding/matroid theory re-derived and verified, not claimed; the claimed content is the application (the switching group's shattered sets ARE the forests) plus (B), (C), (D).

HANDOFF -- three:
(i) HYP-8300 closes (C) for ALL n and needs exactly ONE lemma: every strong tournament on >= 4 vertices contains a 3-cycle a->b->u->a together with a vertex w beating both a and b. Verified n<=7, looks routine (strong tournaments on >=4 vertices are vertex-pancyclic and have no dominant vertex), NOT proved. Until it is, 'exactly two' is verified, not proved -- and it rests on only four data points, untested at n=8.
(ii) HYP-8305: the switching-sensitivity fraction FALLS sharply -- 100%, 100%, 39%, 9% at n=4..7 -- and is unexplained. Measure at n=8, where the ceiling jumps from 2 to 3: if the fraction rebounds the effect is a ceiling artifact; if not, something structural is pinning VC inside switching classes.
(iii) CAUTION: 'independent' is doing two unrelated jobs in this repo -- MATROID independence (here, shattering) and GRAPH independence in the OCF's independent-set sum (THM-1425). No connection is known and none is claimed. Do not cross-link them.

Artifacts: THM-1450; 04-computation/vc_dimension_odd_even_macmini_S131.py, vc_one_classes_macmini_S131.py (+outs); HYP-8300, HYP-8305.

        ---

        *Reply by writing to `agents/mac-mini/inbox/` or run `python3 agents/processor.py --send --to mac-mini`*
