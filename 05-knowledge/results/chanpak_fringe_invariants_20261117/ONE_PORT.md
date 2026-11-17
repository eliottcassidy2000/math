# Exact Hamiltonian substitutability with one exposed port

A core consists of D and one distinguished vertex v. Each external vertex must be uniform to D, but the direction of its arc with v is arbitrary. The exterior can itself be an arbitrary tournament. Define the following sequences for r ordered nonempty paths covering D:

- R_r: unrestricted.
- E_r: the first path ends at a vertex beating v.
- S_r: the first path starts at a vertex beaten by v.
- B_r: the first path ends at a vertex beating v and the second starts at a vertex beaten by v. This is used only for r>=2.

These are exact and complete: two cores give equal H in every common exterior of the indicated kind iff the four sequences coincide. For two cores with at most n vertices in D it is sufficient to use the explicit tests below, requiring at most 2n−2 ordinary exterior vertices if n>=2 (and two for n=1). This is scalar Hamilton tomography with one exposed port, not adjacency reconstruction.

## Sufficiency

Contract each maximal consecutive run in D in a full path to the symbol D. Classify the resulting exterior word, with v appearing once. Every contact between an ordinary exterior vertex and a D run depends only on the specified uniform direction. The only arc tests on vertices inside runs are immediately before or after v, and give exactly R, E, S or B. If two runs touch v, they are distinct runs. Permuting the labels of ordered nonempty run slots lets the affected slots always be called first (or first and second). Each exterior word therefore contributes one of the four counts.

## Necessity: transverse probe algebra

An ordinary probe pair consists of transitive sets A and F, with D→A→F→D. Set alpha=1 if v beats all of F; set beta=1 if all of A beats v. Their complement signs specify all opposite arcs. The four logical sign settings are four physical tournaments. A first difference subtracts sign 0 from sign 1; a mixed difference subtracts twice.

In the contraction, each regular separator between consecutive D runs uses at least one A followed by one F. A separator with v alone contributes a B statistic. A two-vertex separator containing v is either (v,F) or (A,v), contributing E or S. A separator containing v but touching neither adjacent D run uses length at least three (A,v,F at minimum). These simple length/side constraints give the following triangular decoder.

1. **Recover R at every r.** A mixed difference in alpha,beta kills every path with any D-v contact (only one remaining v edge could depend on the variables), and also every path with v terminal. For r>=2 use |A|=|F|=r−1. The surviving paths have at most r runs: v does not itself save an ordinary separator vertex. At maximum all exterior vertices are used in separators, one of which is A,v,F and the other r−2 are A,F. The opposite order F,v,A cannot give a shortest separator because D loses to F and beats A. The mixed coefficient is (r−1)((r−1)!)² R_r plus already known lower R. To obtain the base R_1, instead use two exterior vertices X,Y both losing to D, and vary v→X and v→Y separately. The mixed difference forces the words D,X,v,Y or D,Y,v,X, each of coefficient −1, giving −2R_1.

2. **Recover E_1,S_1.** Use X→Y, with both exterior vertices losing to D. Let x=[v→X], y=[v→Y]. The coefficients of E_1 and S_1 are respectively x and 1+(1−x)+2(1−y). For E_1 the single D run starts the path, then v and the full X→Y. For S_1 the exterior portion before v can be empty or any increasing prefix-subset whose last vertex beats v; the unused vertices appear after the D run. The only two-run marked case has v between the D runs and all ordinary exterior at the end, independent of the two signs. Terms with no D-v contact are known from R. Thus vary y first to read S_1, then vary x to read E_1.

3. **Recover E_r,S_r for r>=2.** Use r−2 A's and r−1 F's (thus 2r−3 ordinary exterior), and take the alpha difference. Every B term drops because both v edges go into D. Among surviving paths there are at most r D runs, because all r−1 gaps now need at least two vertices each and the total of ordinary exterior plus v is 2r−2. At r runs all gaps have exactly two vertices, so the special gap must be v,F: otherwise the supply of A's is insufficient. The leading coefficient is (r−1)(r−2)!(r−1)! E_r. All lesser E,S terms and all R terms have already been recovered. Exchange set sizes and take the beta difference to read S_r instead; the special gap is A,v with the same positive coefficient.

4. **Recover B_r for r>=2.** Use balanced |A|=|F|=r−2. There are at most r runs; at r all gaps are regular A,F except for one consisting solely of v, hence the leading coefficient is (r−1)((r−2)!)² B_r. Subtract known smaller B terms and all already determined R,E,S terms.

Run steps 3 and 4 in increasing r. Coefficients of lower statistics, if needed algorithmically, are computed by literal enumeration of exterior words with nonadjacent repeated D symbols followed by the sign tests and linear differences. All arguments use signed sums only as observations of actual tournaments; no signed quantity is asserted to be an attainable Hamilton count. The ordinary unmarked uniform module theorem reads just R; the new point is a finite complete decoder when one port is exposed nonuniformly.

## Active, genuinely nonisomorphic interchangeable cores

Take D with outgoing masks (14,28,8,16,5), reading zero-based vertex bits. At the exposed vertex v=5 use outgoing mask 10 in the first core and 5 in the other. Each contact mask has both directions, and D has no proper nontrivial homogeneous module. The complete core rows are respectively

    A=(46,28,40,16,37,10),   B=(14,60,8,48,37,5).

Both have H=29. They are not even isomorphic after global reversal: the unique outdegree-one vertex of A beats a vertex of outdegree three, whereas that vertex of B and of its converse beats a vertex of outdegree two. Their common exact sequences, with r=1,...,5, are

    R=(9,66,186,240,120)
    E=(5,38,110,144,72)
    S=(4,28,76,96,48)
    B=(—,20,56,72,36).

These finite profile entries can be computed by cutting permutations of the five core vertices at all choices of run boundaries, checking every uncut edge and the indicated boundary directions. Thus these two distinct six-vertex tournaments are interchangeable in every one-port exterior, not just a bank of tested exteriors.

The scope boundary is physical: append vertex 6 beating only old vertex 0. This violates uniformity at just that one old vertex (besides the first port that is allowed to vary). The Hamilton counts become 115 and 75. So exposing a second port can genuinely distinguish the two, even with a single-vertex test.

A tiny unit-digit read is also physical. Append a vertex with outgoing mask 51, beating old vertices 0,1,4,5 and losing to 2,3. The two interchangeable precursors now read 135 and 133. This exposes two additional D exceptions rather than just one. The unit separation is a finite exact seed, not a repeatable or growing digit theorem.
