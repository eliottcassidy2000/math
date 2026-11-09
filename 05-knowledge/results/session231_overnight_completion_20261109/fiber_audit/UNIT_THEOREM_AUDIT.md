# Audited positive theorem: arbitrarily long consecutive-odd fibers

## Verdict and precise scope

PASS. There is an explicitly certified family of strong tournaments whose Hamilton-path counts fill arbitrarily long finite runs of consecutive odd integers, with all members of a run having a common full deletion interface. The family is **saturated**: its common deletion-interface gcd is 37, not 1. Its baselines are large and not freely prescribed. This does not prove all-target Chan–Pak coverage, eventual coverage of every odd integer, a bounded-start interval, or an optimal vertex bound.

The all-k conclusion is proved by a boundary-cancellation identity, not inferred from checking finitely many k. `unit_audit.py` is an independent standalone stdlib-only checker with an owned copy of the input, UNIT_WITNESS.json. It imports none of the producer's modules. Normal and optimized runs match exactly. The literal graph-level count check is for the 12-vertex k=0 case; longer cases are justified by exact composition and the all-k proof below, not by purported generic enumeration on 1,200 or more vertices.

## Literal core and actual directed loops

Let C be the seven-vertex tournament with outgoing adjacency masks

    [8,9,3,52,7,87,31], cut mask30.

Literal permutation enumeration of Hamilton paths and ordered covering pairs of paths gives

    (h,q,A,G,T,S,E) = (37,740,269,422,132,12,21).

Its row-major upper-triangle code is1147012. The checker verifies this encoding and all tournament arcs independently.

At the augmented carry state (A,G,p)=(269,422,1), the following two 108-step directed closed walks have gains +34 and −90. Their columns are singleton cut-bit words:

    X+ = 011000101101101101000100111100111011110100011011100001011100001110110101000110001010001010111010000111100000
    Y+ = 011100111111000011110100010101101110000111010011011011000101111011100111101001111000011000011011101001001001

    X− = 011000101101101101000100111100111100001101000011000100001011100100011110100011110001001010111010000111100000
    Y− = 011100111111000011110100010101110101110011101010010001111000101100100100111010000100011000011011101001001001

Define X=(X+)^8(X−)^3 and Y=(Y+)^8(Y−)^3. Both have length d=1188. Each of the eleven segments is literally checked to close at the same augmented root, using only forward edges

    A' = 2A+hy−qx,
    2G' = G+qy−hx,
    p' = y.

Integrality and −h<A',G'<q are checked for every edge. No inverse edge is assumed and no path is reversed. The local weight

    w = Ay+hp y−xG'−hx(G' mod2)

sums to +34 on each of the first eight segments and −90 on each of the last three. Thus the concatenated gain is8·34−3·90=2. Direct whole-word evaluations independently certify both constituent gains and the final gain, so the result does not rest solely on the cocycle bookkeeping.

The endpoint bits are first(X)=0, last(X)=0, first(Y)=0, last(Y)=1. In particular the correct augmented root has p=1, not0. The boundary bits used below are u=last(X)=0 and v=first(Y)=0.

## Global boundary identity

For a singleton word Z write its state as

    (1,2^|Z|,A_Z,G_Z,T_Z,S_Z,E_Z).

These quantities are independently computed without calling the general composition function: A_Z is its forward binary integer, G_Z its reverse binary integer, S_Z and E_Z its first and last bits, and

    T_Z = Σ_{i<j} z_i z_j 2^(j−i−1) + Σ_i z_i z_(i+1).

The closed carry equations imply C Y and X C have equal h,q,A,G. Their unflanked T difference is

    δT = A G_Y − A_X G + h(T_Y−T_X) + E S_Y − E_X S.

Their S difference is S−hS_X and their E difference is hE_Y−E. Adding singleton boundary bits u=E_X and v=S_Y therefore gives

    Δ = δT + u(S−hS_X) + v(hE_Y−E)
      = A G_Y − A_X G + h(T_Y−T_X) − huS_X + hvE_Y.

The physical core's S,E terms cancel exactly. For the explicit words above, exact evaluation gives Δ=2. Since u=v=0 here, the final two terms also vanish. Consequently s0 C Y s0 and s0 X C s0 have identical full interfaces and T values differing by exactly2, with the left side larger.

The full-interface equality and the actual difference are checked separately by composing the complete core and word states. This is an independent global control on the directed-loop calculation.

## All-k interval theorem

For every integer k>=0 and j in {0,...,k}, form the ordinal deletion graph

    L_(k,j) = s1 ⇒ s0 ⇒ X^j ⇒ C ⇒ Y^(k−j) ⇒ s0.

It has 1188k+10 vertices. All singleton words describe literal transitive ordinal blocks with selected vertices given by their bits; C has the fixed physical graph and cut above. Every listed construction is therefore an actual tournament/cut object, not a bookkeeping state.

Successive positions differ by replacing C Y with X C. That pair has equal h,q,A,G, so those coordinates remain equal under the common ordinal contexts. The outside singleton boundaries give S=h=37 and E=0 for every j.

For every adjacent replacement, the left common context ends with0: either the explicit leading s0 or the last bit of X. The right common context starts with0: either the explicit trailing s0 or the first bit of Y. Both contexts have Hamilton count1. Thus the same global boundary calculation applies at every position and for every k, with no multiplier or changing endpoint correction:

    T(L_(k,j+1)) = T(L_(k,j)) − 2.

Their common full interface has

    h=37,
    q=740·2^(1188k+3),
    S=37, E=0,

and common A,G. In particular gcd(h,q)=37 for every k. Saturation is genuine and must not be described otherwise.

Apply star closure by adding vertices v,r with every L vertex beating v, v beating r, r beating exactly the selected cut vertices, and the unselected L vertices beating r. The count identity is

    H(star(L)) = 2A+h−T.

Hence the k+1 star counts are exactly

    N_k, N_k+2, ..., N_k+2k,

on 1188k+12 vertices each. Strongness holds for every k because r reaches the first selected singleton, that singleton reaches all later L vertices, all L vertices reach v, and v reaches r.

The baseline is positive: any Hamilton path through C, preceded and followed by the ordinal singleton vertices and then v,r, is a Hamilton path of the star. Thus there are at least h=37 such paths independently of the cut. Oddness can also be seen directly from the baseline formula below, without invoking a general tournament parity theorem.

This proves arbitrarily long consecutive-odd runs in actual strong-tournament Hamilton counts and arbitrarily long consecutive step-two full-interface T fibers. It does not turn those finite runs into an all-target theorem.

## Explicit baseline and its limitation

Let W=Y^k. For k>0 compute A_W,G_W,T_W by the singleton-word formulas above; for k=0 set these three scalars to0. Then exact composition gives

    N_k = 6996·2^(1188k) +148A_W−1749G_W−37T_W−939.

The first bit of every nonempty W is0, so G_W is even. T_W is always even: appending bit b adds b(A_prefix+last(prefix)), whose two summands have the same parity. The other terms show N_k is odd. This also covers k=0, where N_0=6057.

The checker directly enumerates Hamilton paths in the actual 12-vertex k=0 star and obtains6057, and checks its strong connectivity from every vertex. Exact arithmetic for k=0,...,8 confirms the common interface, every count N_k+2j, the baseline formula, and the gcd37.

Already k=1 has a 1201-bit baseline, while supplying only a two-value run on1200 vertices. For k=8 the baseline has9517 bits while the interval has only9 values. These examples make the central caveat concrete: the construction guarantees interval length, but does not permit choosing the interval's starting point or filling the gaps between these baselines.

## Files and reproduction

    python3 -B unit_audit.py
    python3 -O -B unit_audit.py

Run from the audit directory with UNIT_WITNESS.json next to the script. Results are UNIT_AUDIT.json and matching normal/optimized logs. The receipt records every rooted-segment gain, exact root and endpoint bits, physical core, checks through k=8, literal base-star adjacency/count, and a hash of the full checked forward-edge trace. No producer files were changed and no commit was made by this audit worker.
