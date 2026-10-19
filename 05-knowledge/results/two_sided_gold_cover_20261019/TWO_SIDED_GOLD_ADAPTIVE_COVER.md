# Five trees cover every two-sided exterior of the ten-vertex core

19 October 2026. Status: PROVED by the bounded-chain reduction and a finite exact certificate; independently reviewed. No general seven-thin Gold theorem, general balance theorem, or historical novelty is claimed.

## 1. The arbitrary-size two-sided case is closed for this core

Let P be a finite poset of maximum incomparability degree at most 7. Its full root core at roots x,y means the induced subposet on x,y and every vertex incomparable with at least one of them. Suppose this core, with roots identified as (0,3), has strict predecessor masks

    Q=(0,0,0,6,0,19,2,18,86,151).

These roots are critically oriented: D(0) is contained in D(3), and U(3) is contained in U(0). Both inclusions hold in every ambient poset with this full core. For an extension chosen uniformly from P, a comparison a<b is a binary question about its order in that extension.

THEOREM. One of the following five two-question trees has strict Gold on P. Each tree uses its listed second question on BOTH first branches.

    T1: first 3<7, then 6<7
    T2: first 6<7, then 3<5
    T3: first 5<8, then 3<7
    T4: first 3<7, then 5<7
    T5: first 3<5, then 2<4.

More quantitatively, put E=e(P). For every answer path in the selected tree, let b be the number of extensions in its first-answer branch and c the number in its complete two-answer cell. Then

    E-b-c > E/10000.                                      (1)

There is no size bound on P or either exterior, and both exteriors may be active. All questions involve core vertices. This closes the two-sided extension problem for this specific full core, not the classification of seven-thin cores.

The result is stronger than individual realizability requires: the same five validity regions cover the ENTIRE convex hull of the 1,764 bounded two-chain atom laws. No arbitrary mixture is identified with an ambient poset. The previously developed factorization and nonlinear realizability conditions are not needed for this particular closure.

The [one-sided theorem](../one_sided_exterior_gold_20261019/ONE_SIDED_ARBITRARY_EXTERIOR_GOLD.md) supplied one tree on each one-sided universe. The [common-first-pair obstruction](../adaptive_gold_core_selection_20261019/ADAPTIVE_GOLD_CORE_SELECTION.md) showed that one first core question cannot work for all two-sided attachments. Here the first question may be selected from four distinct core comparisons, and the full tree from the five displayed choices. No minimality of these numbers is claimed.

## 2. Complete finite model family and exact transfer

Every vertex outside the full core is comparable with both roots. It must either be below both or above both: being between them would make the roots comparable. Write the two sets as A and B. Then A<B. A vertex of A cannot be above a core vertex q, because that would make q below both roots, contrary to membership in the full core; dually for B. Thus the [bounded-chain compression theorem](../boundary_chain_compression_20261019/BOUNDARY_CHAIN_COMPRESSION.md) applies.

The core incomparability slack in label order is

    delta=(0,4,1,0,2,1,0,1,2,3).

Every conditioned lower exterior order has one of exactly 46 profiles

    alpha=(0,a,b,0,c,0,0,d,0,0),
    0<=a<=4, 0<=b<=1, 0<=c<=2, 0<=d<=min(1,a,c).

Every conditioned upper exterior order has one of exactly 48 profiles

    beta=(0,0,0,0,f,g,0,h,i,j),
    0<=f<=1, f<=g<=1, f<=h<=1, f<=i<=2, h<=j<=3.

All variables here are integers. Pair the profiles only when alpha(q)+beta(q)<=delta(q) for every q. Equivalently, the only additional joint restrictions are

    c+f<=2,       d+h<=1.

There are exactly 1,764 compatible pairs. Each gives an actual seven-thin model C_(alpha,beta) with Q, a lower chain of length max alpha and an upper chain of length max beta. The model has between 10 and 17 vertices. With lower chain c_1<...<c_s and upper chain d_1<...<d_t, its relations with Q are

    c_i<q exactly when i<=s-alpha(q),
    q<d_j exactly when j>beta(q),

and all lower-chain vertices lie below all upper-chain vertices. The verifier independently enumerates all bounded monotone profiles, checks equality with the two displayed formulas, checks every actual model's thinness and full core, and counts its extensions directly.

The compression identity says that for every core word w,

    C_P(w)=sum_(sigma in LE(A), tau in LE(B))
                  C_(C_(alpha_sigma,beta_tau))(w).          (2)

Consequently the normalized core law of P is a convex combination of the atom laws. Each occurrence of a conditioned exterior-order pair has mixture weight

    e(C_(alpha_sigma,beta_tau))/e(P).

These are extension-count weights, not uniform weights. Equation (2) preserves all joint comparison-event counts simultaneously. Thus every homogeneous linear inequality verified on all 1,764 atom count vectors transfers to every P in the theorem, as well as to the entire normalized atom hull.

## 3. The eight-inequality certificate

For tree Ti, write its four cell counts as c_i,TT, c_i,TF, c_i,FT, c_i,FF, where the first letter is the answer to its first question. Define the four Gold margins

    m_i,uv=E-(c_i,uT+c_i,uF)-c_i,uv.

Let epsilon=1/10000. The certificate first verifies the following ten inequalities universally over the atom family:

    m_1,TF > epsilon E;
    m_2,FT, m_2,FF > epsilon E;
    m_3,FT, m_3,FF > epsilon E;
    m_4,TT, m_4,TF, m_4,FT > epsilon E;
    m_5,FT, m_5,FF > epsilon E.                            (3)

The following eight weighted inequalities also hold universally. In each row the right-hand side is epsilon E times the sum of the coefficients on the left.

    I1: m_1,TT + m_4,FF > 2 epsilon E
    I2: m_1,FT + m_2,TF > 2 epsilon E
    I3: m_1,FF + m_2,TF > 2 epsilon E
    I4: m_1,FT + m_3,TT > 2 epsilon E
    I5: m_1,FF + m_2,TT + m_3,TT > 3 epsilon E
    I6: 3 m_1,FF + 4 m_2,TT + m_3,TF > 8 epsilon E
    I7: 6 m_1,FT + m_3,TF + 3 m_5,TT > 10 epsilon E
    I8: 64 m_1,FT + 58 m_2,TT + 26 m_3,TF + 21 m_5,TF
                                                   > 169 epsilon E. (4)

All coefficients are nonnegative integers. The largest is 64. A weighted average of allegedly bad margins cannot exceed epsilon E, so each row prohibits its whole support from being bad simultaneously.

For auditability, the exact minimum over all atom laws of each left side divided by E times its coefficient sum is, in order I1 through I8,

    2542/38111,
    540/30163,
    2205/55193,
    2481/163474,
    362/30461,
    391/109746,
    872/380105,
    3205/27793233.

Every fraction is strictly greater than 1/10000. These are minima of the particular displayed certificates, not a claim of the optimal five-tree uniform Gold margin. The supplied receipt also lists all ten exact safe-path minima and explicit profile indices attaining every recorded minimum.

## 4. A short elimination proof of the five-region cover

Call a margin bad if it is at most epsilon E. Suppose for contradiction that every tree has at least one bad margin.

1. By (3), T4's bad margin must be m_4,FF. Inequality I1 then prevents m_1,TT from being bad. Thus T1 must have a bad FT or FF margin.
2. If m_2,TF were bad, I2 and I3 would prevent both remaining possible T1 bad margins. Hence T2's bad margin must be m_2,TT.
3. If m_3,TT were bad, I4 would prevent T1's FT margin from being bad, and I5, together with the already bad m_2,TT, would prevent T1's FF margin from being bad. Hence T3's bad margin must be m_3,TF.
4. With m_2,TT and m_3,TF bad, I6 prevents m_1,FF from being bad. Therefore m_1,FT is bad.
5. Inequality I7 now prevents m_5,TT from being bad. Inequality I8 prevents m_5,TF from being bad. The other two T5 margins are safe by (3). This contradicts the supposed failure of T5.

Therefore one tree has all four margins greater than epsilon E. This proves (1) and the theorem. Notice that the conclusion uses the linear inequalities AFTER transferring them to P. It never commits the invalid inference that individually Gold atoms have a Gold mixture.

Equivalently, write R_i for the four-inequality region m_i,uv>epsilon E in the normalized atom hull. The proof shows that R_1 through R_5 cover that entire hull. An observable selection rule is simply to compute the five minimum margins and choose a tree attaining their maximum. This requires six distinct one-comparison probabilities and five joint probabilities. It does not require knowing a hidden exterior-profile decomposition, and no claim of efficient extension counting is made.

## 5. Independent exact reproduction

Run from this directory:

    python verify_two_sided_gold_cover.py
    python -O verify_two_sided_gold_cover.py

Normal and optimized receipts agree on 52,951 explicit guards. Every guard is an explicit runtime check; none relies on Python assertions. The verifier uses only the standard library and literal core, tree and eight-row certificate data. It does not load the numerical discovery arrays, an LP solver, or another research packet's executable code.

For each of the 1,764 actual models, the verifier computes 12 counts directly by ideal dynamic programming: the total E, six distinct comparison numerators, and the five needed joint numerators. The four cells for each tree are recovered by inclusion-exclusion. The DP recurrence, after adding the requested comparison edges to the predecessor constraints, is

    F(S)=sum_(q available after S) F(S union {q}),
    F(all vertices)=1.

Adding comparison edges need not trigger a separate transitive-closure computation: legal topological orders already respect every implied relation; a cycle leaves no complete extension and hence gives count zero. This is independent of the gap-count formula used during discovery.

The verifier then checks (3) and (4) with strict INTEGER comparisons. Finally it independently enumerates the 4^5=1,024 possible choices of one bad path per tree. Exactly 1,000 are eliminated by the safe-path list, and each of the remaining 24 contains the complete support of one of the eight weighted inequalities. This finite logical check supplements the five-step proof above.

The optional command

    python verify_two_sided_gold_cover.py --table

emits the complete exact atom count table as JSON lines. Its SHA-256 is

    d4caae3e8354652af27e52ed8c33e213ff12e462168933719abeeb44c6f65547.

The table is reproducible output, not required input. The small frozen payload retains the proof, verifier, exact receipts, and their hashes; numerical mixed-integer and LP discovery files remain outside it.

## 6. Scope and next structural question

This theorem settles arbitrary two-sided exteriors of the specified core, including all its fixed-interface realizable profile laws. It supplies no classification of other cores. The 1/10000 margin and five-tree menu are convenient exact choices, not proved optimal. The prior obstruction still rules out a universally valid single first core question, but does not establish a lower bound of five trees or four distinct first questions.

The key reusable instrument is the cover certificate: identify margins that are universally safe, then rule out simultaneous failure of the remaining trees by a few common positive linear combinations. This is different from demanding one universal Gold tree, and it can succeed even when first-question balanced-pair lists are disjoint. For another core the whole atom hull may genuinely fail; only then would its actual factorization and nonlinear interface restrictions be needed to repair the relaxation.
