# A three-tree strict-Gold cover with a 1/200 margin

19 October 2026. Status: PROVED by exact finite certification and the bounded-chain reduction; independently reviewed. This strengthens the preceding [five-tree cover](../two_sided_gold_cover_20261019/TWO_SIDED_GOLD_ADAPTIVE_COVER.md), which remains a separate valid checkpoint. No general seven-thin theorem or historical novelty is claimed.

## 1. Theorem

Let P be any finite poset whose maximum incomparability degree is at most 7. Suppose its full root core, meaning the induced subposet on roots x,y and all vertices incomparable with either root, is identified with

    Q=(0,0,0,6,0,19,2,18,86,151),       (x,y)=(0,3),

in strict predecessor-mask notation. Both exterior sides may be active and have arbitrary finite size.

THEOREM. One of these three comparison trees is strict Gold on P:

    A: first 3<7, then 0<2;
    B: first 5<8, then 0<6;
    C: first 6<7, then 5<8.

The second comparison is the same on both branches of each tree. If E=e(P), b is the extension count of a first-answer branch and c the extension count of a complete answer path inside it, the selected tree has

    E-b-c > E/200                                         (1)

on all four paths. Thus three fixed polyhedral validity regions cover the entire two-sided atom hull, with a uniform margin exceeding 1/200.

Only five distinct comparisons are used. In particular, 5<8 serves as B's first question and C's second question. The tree itself is selected from the projected law; it need not be selected from a hidden mixture decomposition. No minimality of the number of trees or optimality of the margin is claimed.

A direct corollary is that one of the three first questions has probability strictly between 101/300 and 199/300. Indeed, if h is its heavy first-branch mass, the heavier second cell has mass at least h/2, so (1) gives 3h/2<199/200.

## 2. Exact model universe and all-size transfer

Every vertex outside the full root core is below both roots or above both. Otherwise it lies between them and makes them comparable. Write these two sets as A0,B0. Then A0<B0. A lower exterior vertex cannot lie above a core vertex, and an upper exterior vertex cannot lie below one. The roots are critically oriented in Q and remain so in P. Therefore the [bounded-chain compression theorem](../boundary_chain_compression_20261019/BOUNDARY_CHAIN_COMPRESSION.md) applies without additional hypotheses.

The core slack vector is

    delta=(0,4,1,0,2,1,0,1,2,3).

All conditioned lower profiles are exactly

    alpha=(0,a,b,0,c,0,0,d,0,0),
    0<=a<=4, 0<=b<=1, 0<=c<=2, 0<=d<=min(1,a,c).

All conditioned upper profiles are exactly

    beta=(0,0,0,0,f,g,0,h,i,j),
    0<=f<=1, f<=g<=1, f<=h<=1, f<=i<=2, h<=j<=3.

The coordinates are integers. There are 46 lower profiles and 48 upper profiles. Compatibility requires c+f<=2 and d+h<=1, yielding exactly 1,764 profile pairs. Each gives an actual seven-thin two-chain model on 10 through 17 vertices, with the same full root core. The verifier independently compares these formulas against general bounded monotone profile enumeration and checks every resulting poset.

For every core extension word, the compression identity expresses its ambient lift count as the sum of its lift counts in the chain models associated with pairs of exterior extension orders. Dividing by E therefore represents the projected law of P as a convex combination of the normalized atom laws. The weight per exterior-order pair is its atom's total extension count divided by E; it is not generally uniform.

Every inequality below is homogeneous and linear in comparison-event counts. It is checked on all atoms and therefore holds on their entire hull, in particular on every actual P. The cover is proved AFTER transferring these common inequalities. No conclusion is inferred merely from individual atoms having Gold, and no arbitrary mixture is asserted to be realizable by an exterior poset.

## 3. Seven safe margins and four tiny separating inequalities

Write m_T,uv=E-b-c for tree T and answer cell uv, with u,v in {T,F}. Set epsilon=1/200. The following seven margins are universally safe:

    m_A,TT, m_A,TF > epsilon E;
    m_B,TF, m_B,FT, m_B,FF > epsilon E;
    m_C,FT, m_C,FF > epsilon E.                           (2)

The entire remaining cover is certified by just four inequalities:

    m_A,FT + m_C,TT > 2 epsilon E;
    m_A,FF + m_C,TT > 2 epsilon E;

    m_A,FT + m_B,TT + 2 m_C,TF > 4 epsilon E;
    m_A,FF + m_B,TT + 2 m_C,TF > 4 epsilon E.              (3)

The largest coefficient is 2. These are exact verified inequalities, not numerical LP tolerances.

The exact minima over all atom laws of the left sides divided by E times their coefficient sum are, in the displayed order,

    1163/178770,       715/114938,
    1164/62011,        2139/160244.

Each is greater than 1/200. The safe-margin minima, in the order in (2), are

    19831/72503,       2000/242141,
    2393/169533,       5641/42972,       13003/172835,
    2113/28965,        41591/128838.

These likewise exceed 1/200. The receipt records profiles attaining all minima. The certificate's smallest normalized weighted margin is 715/114938; no optimality claim is made for the resulting conservative 1/200 theorem.

## 4. Short proof and explicit selection rule

Call a margin bad if it is at most epsilon E. If all three trees failed, (2) would force B's bad margin to be m_B,TT, C's bad margin to be one of m_C,TT,m_C,TF, and A's bad margin to be one of m_A,FT,m_A,FF.

If m_C,TT is bad, the first two rows of (3) make BOTH possible A bad margins safe, a contradiction. Otherwise m_C,TF must be bad. Together with bad m_B,TT, the last two rows of (3) again make both possible A bad margins safe. This is also a contradiction. Thus one tree satisfies (1).

The same proof gives a simple deterministic rule, where all margins are normalized by E:

1. If m_C,TT<=epsilon, choose A.
2. Otherwise, if m_C,TF>epsilon, choose C.
3. Otherwise, if m_B,TT>epsilon, choose B.
4. Otherwise choose A.

Every selected tree has all four normalized margins greater than epsilon. This rule only tests three particular margins. It avoids computing A's probabilities at selection time: the certificate supplies its guarantees whenever it is selected. Alternatively, computing all three minimum margins and selecting the largest also works.

In event-count form the three tests require only four probabilities: P(6<7), P(6<7 and 5<8), P(5<8), and P(5<8 and 0<6). Specifically,

    m_C,TT/E=1-P(6<7)-P(6<7 and 5<8),
    m_C,TF/E=1-2 P(6<7)+P(6<7 and 5<8),
    m_B,TT/E=1-P(5<8)-P(5<8 and 0<6).

This is an intrinsic selection rule on the core law. No efficient algorithm for computing extension probabilities is claimed.

## 5. Reproduction and scope

Run

    python verify_three_tree_gold_cover.py
    python -O verify_three_tree_gold_cover.py

The standard-library receipts are byte-identical and pass 37,055 explicit guards. The verifier reconstructs every compatible chain model directly and uses ideal dynamic programming to count E, the five singleton comparison numerators, and the three joint numerators: nine direct counts per atom. It recovers each tree's four cells by inclusion-exclusion and checks (2),(3) using strict integer inequalities.

The independent logical guard considers all 4^3=64 choices of a bad path from each tree. The safe list eliminates 60. Each of the four remaining patterns contains the support of a row in (3).

The optional --table output has SHA-256

    2f83a7b75f377b5160c37b83b30f64da1ec3207531dcd575983a5d48da9c767d.

The receipt lists the certificate minima in implementation order FT/TT, FT/TF, FF/TT, FF/TF; the written note groups the two TT rows and then the two TF rows. They are the same four inequalities.

The full two-sided extension problem for this specified core is settled. The theorem does not classify other cores, does not prove general seven-thin balance or Gold, and does not assert that three trees are necessary. The five-tree checkpoint remains intact; this refinement improves its menu size, margin, and proof complexity independently.
