# Positive defect transport without a zero reset

2026-11-04. **PROVED exact identities and an unbounded named-state transport; FINITE-EXACT controls.** No closed two-digit realization grammar or full Hamilton spectrum theorem is asserted.

## 1. A positive formula on the deleted physical root

Let B be a nonempty tournament with physical mark w, U=B−w, and I=N_B^−(w). All names and incidences are retained. Work in the squarefree commutative ring R_U=Z[z_v:v∈U]/(z_v²). Let G_ij enumerate all simple directed i-to-j paths, including singleton paths. Write

    P=Σ_ij G_ij, e_I=Σ_{j∈I,i∈U}G_ij, s_I=Σ_{i∈I,j∈U}G_ij,
    s_O=Σ_{i∈U\I,j∈U}G_ij.

Products impose disjointness; [U] extracts the monomial on every vertex of U. Empty slots have weight1. Then

    H(B)=[U](1+e_I)(1+s_O),
    E_w(B,2)=[U](1+e_I)(1+P),
    b_w(B):=E_w(B,2)−H(B)=[U](1+e_I)s_I.             (1)

Proof: split a Hamilton path at w. Its left slot is empty or ends at an actual in-neighbor, while its right slot is empty or starts at an actual out-neighbor. For E_w, remove w from the end of its first slot; the second slot is unrestricted. Subtract and use P=s_I+s_O. This also works for B={w}, when U is empty and b=0.

Thus defect is a nonnegative integer counting an explicit family of bad joins. In particular it is not an observer phase, and cannot be set to zero by changing the name of the state. If I={a}, the two nonempty-slot product e_I s_I vanishes because both paths contain a. Therefore b_w(B) is exactly the number of Hamilton paths of B−w starting at a. For larger I, the cross-port joint terms remain essential.

## 2. Exact fixed-mark birth transport

Add a physically named vertex x, keep w fixed, and on U set a_i=[i→x], q_j=[x→j], ε=[x→w]. Define

    X=1ᵀGa, Y=qᵀG1, C=qᵀG1_I, D=1_IᵀGa.

The new in-neighbor set is I'=I∪{x} if ε=1, otherwise I'=I. The full endpoint update gives

    e'_I'=e_I+z_x(1+X)(C+ε),
    s'_I'=s_I+z_x(D+ε)(1+Y).

Substitution in (1), using z_x²=0, proves the all-order law

    b_w(B+x)=[U]{(1+X)(C+ε)s_I+(1+e_I)(D+ε)(1+Y)}. (2)

Every factor uses the actual old endpoint matrix. In particular D and C cannot in general be reconstructed from H, b or separate one-port marginals. Together with the prior named endpoint birth/reversal/deletion rules, (2) is an exact unbounded algebraic transport system: arbitrary finite operation words can be evaluated without scalar or phase resets. It is not a theorem that an operation word realizing either desired ternary digit exists, and it does not supply native complement legality certificates.

## 3. An arbitrary-length defect-preserving construction

Add x with w→x and x→every vertex of U. Then a=0, q=1, ε=0, so (2) specializes to

    b_w(B+x)=b_w(B)−S_I(U),                        (3)

where S_I(U) is the number of Hamilton paths of U starting in I. This is a literal permanent birth, preserving every old edge and name. If B is strong and has at least two vertices, B+x is strong: x reaches all old vertices except w directly, reaches w through the old strong B, and all old vertices reach x through w.

Repeat this birth t≥1 times. The newcomers form a transitive block TT_t above U, all beaten by w. Put h=H(B), f=E_w(B,2), A=H(U) (A=1 when U is empty), and S=S_I(U). Then

    H(B_t)=h+S+(2^t−1)(f−A),
    E_w(B_t,2)=A+2^t(f−A),
    b_w(B_t)=b_w(B)−S.                             (4)

If U is empty, B_t is transitive with w as its source, so h=f=A=1 and S=b=0 verify all three lines directly. For the remaining proof assume U is nonempty. Proof of the first line: split a path at w. If the right slot contains a newcomer and the left old slot is empty, all newcomers must lie right, contributing A. If the right slot contains a nonempty subset of newcomers and the left old slot is nonempty, there are2^t−1 assignments of newcomers and the old-slot contribution is [U]e_I(1+P)=f−A. If the right slot has no newcomer, all newcomers lie left; the old contribution is h−S_O(U). Since A−S_O(U)=S, the formula follows. Equation(3) gives the constant defect because no Hamilton path of TT_t→U starts in I. Adding the two gives the second line.

This is a genuinely unbounded named construction, but its numerical growth is binary and it does not close both ternary digits. In particular, once a nonempty transitive prefix is already present, further such births preserve the actual defect exactly, rather than repairing it.

## 4. Multiport lift: the missing third coordinate

For strong suffix C, a nonempty physical port set W⊆C, and s≥1, form T_s by TT_s→C, a root r which beats TT_s and C\W, and W→r. On C define

    h=H(C), F=Σ_{w∈W}E_w(C,2),
    J=[C]e_W s_W.

Here F is a sum of single-end responses, while J contains the directed, disjoint, jointly constrained end/start responses. Then

    E_r(T_s,2)=h+2^sF,
    H(T_s)=h−J+2^sF,
    b_r(T_s)=J.                                    (5)

For defect, (1) applies to U=TT_s→C. A Hamilton path of U cannot start in W. In the remaining two-slot product all prefix vertices must occur in the left slot because the right starts in C. Thus the coefficient is precisely J. For E_r, split the two ordered slots across TT_s→C; assigning the prefix vertices contributes2^s to the suffix end-marked response, and the case with root first and all U in the other slot contributes h.

For two such states, with α=h−J, the necessary and sufficient condition for a shared-clock branch H'_s=3H_s+d at every s≥1 is

    F'=3F, α'=3α+d.                                (6)

It does not force J'=0 or J'=3J. Root defect is a third coordinate to transport separately. The sufficiency is substitution; necessity follows by subtracting the identities at two clocks. Equality of scalar count at one clock establishes neither condition.

These formulas explain exactly why a genuine multiport+4 endpoint may exist while the earlier single-port zero-ready search excludes it. Neither class can be substituted for the other without retaining J and the underlying physical relations.

## 5. Verification and scope

`verify.py` constructs path polynomials by direct enumeration and compares extracted coefficients against independent endpoint subset DP. It exhausts all labelled tournaments of orders1 through4, every physical mark, and every birth cut:285 defect identities,4306 birth identities,855 arbitrary-prefix checks at lengths1,2,3. It also verifies strongness preservation whenever the source is nontrivial strong. Checks use explicit exceptions and remain enabled under `python -O`. Normal and optimized output are recorded separately.

A supplemental bounded census `search.cpp` exhausts8192 ordered two-birth extensions of each of the two named Session180 prototype outputs, keeping every old edge and mark0. It finds six raw(53,63) hits from(17,21), no(55,63); and two raw(59,63) hits from(19,21), no(61,63). The last two have a global sink and are not strong. These raw count tests are not an additional closure theorem; the grammar worker's separate treatment classifies the continuation scope.

The prior F49/F343 coordinates and named81/27 hostile matrices are not altered here. Equation(2) reads the full named endpoint matrix before taking coefficients, so no equality of field phase, automorphism order or H is used to identify different physical states. Current candidate phases must still be recomputed from current endpoint fields; no finite phase label repairs a missing J or full joint interface.

## 6. The actual512425 endpoint: carry3690, do not erase it

The sibling multiport packet supplies a literal named17-vertex candidate, retaining physical root39, prefix50→51 and in-port setW={4,21}. `candidate_transport.py` independently recomputes the suffix endpoint DP, both directed cross-port coefficients, full H and E_root(2), and strongness from its actual rows. The receipt gives

    h=18583, F=124383,
    J(4,21)=1206, J(21,4)=2484, J=3690,
    H=512425=3·170807+4, E_39(2)=516115.

The diagonal joint terms vanish because both paths would contain the same physical port. Equation(5) now gives, at every clock s≥1,

    H_s=14893+124383·2^s,
    E_39(T_s,2)=18583+124383·2^s,
    b_39(T_s)=3690.

This is an exact defect-preserving infinite extension of the actual +4 candidate, not an appeal to a fresh zero-ready representative. Its shared-clock offset from the original source family2195+42153·2^s is

    H_s−3(2195+42153·2^s)=8308−2076·2^s.

For offset4 this forces2^s=4, so s=2 is the unique positive clock. Offset2 would require2^s=8306/2076, not an integer. Thus this new multiport return is a genuine +4 construction at its displayed clock but is not a uniform ternary branch under prefix cloning.

The next concrete algebraic target for a clock-uniform continuation from this actual candidate is a realizable marked multiport suffix with

    F'=373149,
    h'−J'=44681 (+2) or44683 (+4),

together with its actual J' and all named endpoint relations. A claim of zero-ready output would additionally need J'=0; neither the scalar target nor the field codecs imply this. The physical transformation, native complement guards when used, and closure of its output class still need proof. This target permits a nonzero transported defect and is strictly more flexible than silently requiring or assuming zero-readiness.
