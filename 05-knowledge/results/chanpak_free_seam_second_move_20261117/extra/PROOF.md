# Farey endpoint passports for a nontransitive second move

The reflected center is D=(54,32,2,119,6,20,55) with apex v outgoing mark 26. Its signed profile is H(D)=5, optional two-slot count P_D=256, weighted starts and ends L_D=R_D=142, signed values L_D^−=−2, R_D^−=0, and κ_D=91. All bits are zero-based outgoing masks.

## Signed counts are ordinary opposite endpoint counts

For a tournament X and a chosen starting vertex z define M(X,z) as the signed sum of ordered covers by a nonempty path starting at z and an optional second path, sign `(-1)^length(first)`. Let E(X,z) count full paths ending z. Then **M(X,z)=−E(X,z)**. The singleton first path gives −H(X−z). Otherwise use the second vertex t, necessarily z→t: deleting z gives precisely the cover counted by M(X−z,t), with flipped sign. Induction gives

    M(X,z)=−H(X−z)+sum_(z→t) E(X−z,t)=−E(X,z),

where the last equality partitions full paths on X−z by their last vertex and appends z on exactly the missing part. At order one the result is immediate. Summing over marked vertices gives L_X^−=−e_X, where e_X is the full path-end count in the v-outgoing mark set. The terminal version follows by converse.

## Two independent marked strong interiors are physical

For a marked X, let X* mean its converse with complemented mark set, equivalently the v-contact reversal. Form deletion X→D→Y* or Y→D→X*, with the corresponding marks, then adjoin v. Write h_X=H(X), e_X as above, and d_X=H(X+v)−h_X; analogous notation holds for Y. The exact odd-cycle interval transfer gives the first H as

    h_X h_Y κ_D + h_D(d_X h_Y+d_Y h_X)
    + 142(L_X h_Y+L_Y h_X) + 256 L_X L_Y
    + R_D^− L_X^− h_Y + L_D^− L_Y^− h_X.

Every term except the last two is symmetric in X,Y. The reflected minus original difference is exactly

    2(h_Y e_X−h_X e_Y).

This is an exact determinant rather than a scalar or modular approximation. A transitive prefix with arbitrary v-beats word can be placed before either inside block X,Y: it leaves h and e unchanged (all full paths end in the terminal block, and the prefix has one path). Thus the two arbitrary prefix words survive the switch while changing the symmetric baseline. In the actual tournament these arbitrary strong flank interiors are not ordinal additions to the original host. Positive v contact to the first flank and negative v contact from the last make the whole graph strong; a strong first/last inside block reaches its contacts even without a prefix.

Let R_m be transitive on 0,...,m+1 except that m+1→0, m≥1. It is strong. It has h_m=2^m+1; its full path-end counts are

    (1,1,2,4,...,2^(m−1),1).

An ordinary increasing path ends at the last vertex; a path using the back edge chooses freely which smaller interior vertices occur after it and before its prescribed terminal vertex. This gives the displayed powers, including the terminal zero path and the natural path. Subset sums of these counts realize any integer e in [0,h_m] by a literal tag set. For any coprime h_m,h_n, Bézout produces a determinant of magnitude 1 within this interval: choose e_m≡h_n^(-1) mod h_m with 1≤e_m<h_m, then e_n=(h_n e_m−1)/h_m lies in [0,h_n), and signs may be exchanged. An especially simple infinite subfamily uses m,n=m+1, e_m=1,e_n=2 since 2h_m−h_(m+1)=1. For instance use tag 0 in R_m and tags 0,last in R_(m+1). The original reflection gap is exactly two at every such height and under both arbitrary prefix words. This proves a nontransitive variable-height second construction, not exact placement for all targets.

## One freely chosen strong interior and one ordinary flank

Alternatively keep Y transitive with positive binary word w: h_Y=1,e_Y=[w odd]. Choose R_m inside X with mark {0} if w is even, or the complement of {0} if w is odd. Then e_X=1 or h_m−1, so the determinant is respectively +1 or −1; no synchronization between m,w is required. Prepend any transitive binary prefix u to X, including the empty prefix u=0. Let q=2^m, h=q+1 and P=4q+2(3^m−1). The last identity comes from expanding odd packs in R_m: select an odd subset of the m interior vertices for the unique possible nontrivial cycle, plus the two outer vertices. If w is even, the interior R_m passport for v has L=2q−1, R=P−2h, R^−=−q, and d=2q−2. Indeed a subset path starting at 0 is increasing, with arbitrary subsequent vertices except that choosing only the far end fails; a subset path ending 0 uses the back edge or is singleton, giving weighted end total h. The full path starting at 0 is unique. Marked-start cycles through v then select only the even subcount, giving d=L−1. If w is odd, complementary marks replace L by P−h−L, R by h, R^− by −1, and d by q.

The full left flank has weighted starts ℓ=Pu+L and
d_X=d+h(F(u)−1)+uR−[u odd] R^− by the adjacent signed interval law. The two heights are read by the symmetric formula above and the determinant. For even w and u=0 the lower count reduces to

    2^m(380+5F(w)+654w)−66+5F(w)−114w.

Examples of literal twin-prime counts in this new family are 42461/42463 at (m,u,w)=(1,0,21), 26699/26701 at (2,0,10), and 1525031/1525033 at (3,2,25). These are finite witnesses, not twin-prime infinitude. The held-out targets 8191, 65537, 458759 do not land in this specialized one-interior family under the exact bounded decoder in `try_target.py`; failure here says nothing about the full spectrum or the two-strong-interior family.

## Reflection-covariant guards

Let B=X→D→Y*, B'=Y→D→X*. Because L_D=R_D and reflected marks, the total weighted end count of B equals the weighted start count of B', and their signed counterparts coincide. For example the ordinary equality reads

    R(B) = R_X h_D h_Y + P_X R_D h_Y + P_X P_D L_Y
         = L(B').

In the signed equality, all but the first block's terms vanish for R by parity balance of an optional earlier slot, and similarly all but the last for L. A transitive guard with arbitrary v-beats word z may be prepended to B, and its marked converse appended to B'. Both counts then shift by

    h_B(F(z)−1)+z R(B)−[z odd] R^−(B).

Thus the guard preserves every Farey determinant repair. The move is reflection-covariant (different physical side on the two hosts), not one fixed labeled edit. If there is already an arbitrary transitive prefix on that same side, the guard can be absorbed into it; it is not a third intrinsic degree of freedom in that case. Neither this guard nor strong-flank variation has been shown to fill the absolute target gap.

## A scoped no-go for a single extra apex reading old contact bits

Keep the original seven-vertex reflected center and free transitive flanks (before the Farey change). An added vertex s uses one of four fixed Boolean functions of [v→x] to decide [s→x] on each flank independently: constant 0, constant 1, copy, or complement. Its seven center contacts and s-v arc are arbitrary but fixed. Thus a labelled routine has 4²·2⁸=4096 possibilities, and we allow a different routine on each reflected host.

The exact C++ certificate `apex_no_go.cpp` groups routines by their fourteen shift values on the printed fixed test addresses. It leaves only nine ordered candidate pairs, all of which use constant rules on both flanks of both hosts. In particular no routine reading old apex bits on either flank can commute for **every** admissible pair. The certificate does not assert that the nine remaining constant-routine pairs all commute beyond the tests.

The count for each routine is exact: if s→i and j→s, a newly formed odd cycle containing s has a nonempty even-length path from i to j in the old host. Its contribution is two times the number of such paths, weighted by the full Hamilton count on the leftover vertices, by the odd-pack identity. For each i,j this coefficient is computed by ordinary subset path DP. Additional independent literal DP computes every base count and checks sample extended instances. This is a finite separating witness for a universal claim in a finite action space, not a nonexistence theorem for arbitrary selective contacts or for multiple next vertices.
