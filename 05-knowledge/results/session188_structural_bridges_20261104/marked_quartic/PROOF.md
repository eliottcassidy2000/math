# Endpoint-marked module profiles factor named shared-label quartics

**PROVED:** an exact marked-profile factorization for every quotient tournament, including prime quotients and arbitrary pair-free children; and an all-order support-localization rule for endpoint-balanced children. **EXACT CONTROLS:** a pair-free nine-name example and a pair-free fifteen-name prime-five substitution, retaining actual physical query counts. **NOT CLAIMED:** universal separation of arbitrary prime substitutions, a full-response collision, or permission to omit the remaining cubic or quartic coordinates.

## 1. Marked ordered covers

For a nonempty child A and a physical vertex u in it, define E_u(r) to be the number of ordered spanning covers by r nonempty directed paths whose first path ends at u. Define S_u(r) similarly with the first path starting at u. Set these sequences to zero at r=0 and beyond |A|. The choice of first slot is conventional: permuting ordered cover slots gives the same count for any specified slot.

Their normalization is

    sum_u E_u(r)=sum_u S_u(r)=q_r(A).

For n=|A|>=2, the next-to-last cover coordinate has the elementary score formula

    E_u(n-1)-S_u(n-1)=(n-2)!(d_in(u)-d_out(u)).          (1)

Such a cover contains exactly one two-vertex path. When u is a singleton, its contribution is the same on both sides. Otherwise that path is an incoming or outgoing edge at u, respectively, and the remaining slots have (n-2)! orders. This proves (1).

Call A endpoint-balanced when E_u(r)=S_u(r) for every u,r. A vertex-transitive tournament is endpoint-balanced: automorphisms distribute designated endpoints uniformly, so both quantities equal q_r(A)/|A|. Formula (1) also shows every nontrivial endpoint-balanced tournament is regular and therefore has odd order. Deleting any one vertex then leaves an even-order card; for every remaining vertex, its marked endpoint sequences differ at the card's next-to-last coordinate, since d_in-d_out is odd. This identifies a guaranteed nonzero marked-card input, but does not by itself prove that every outer kernel detects it.

## 2. The arbitrary-quotient run kernel

Let T=Q[A_1,...,A_s] be a substitution on fixed physical child sets, and let U=Q^op[A_1,...,A_s] retain all child interiors. Let B be the adjacency matrix of Q, let Z=diag(z_1,...,z_s), and put

    G=(I-ZB)^(-1)Z.

G_ij is the formal generating series of nonempty child-run words starting in color i and ending in color j. Every transition follows Q and visits a different child; each occurrence contributes its child's variable. Reversing a run word gives G_(Q^op)=G^t. All functionals below truncate beyond each finite child size, so these are finite coefficient calculations, not analytic infinite sums.

Write

    c_i=sum_h G_hi,       r_i=sum_h G_ih,
    K_ij=c_i r_j-c_j r_i.                                (2)

For sequences alpha and beta assigned to different children i and j, define L_(alpha,beta) on a monomial z_1^t1...z_s^ts by multiplying alpha(t_i), beta(t_j), and q_(t_h)(A_h) for all other children. No nonempty child has a q_0 contribution.

Let C_uv count ordered spanning two-path covers whose first path ends at u and second starts at v. For u in A_i and v in A_j, i!=j,

    [C_uv(T)+C_vu(T)]-[C_uv(U)+C_vu(U)]
       = L_(E_u,S_v)(K_ij)-L_(S_u,E_v)(K_ij).            (3)

Proof: the generating function of the two path slots for C_uv is c_i r_j. Each slot is nonempty. Read runs in the first slot and then the second slot; the run ending at u and the run starting at v occupy fixed positions in their respective children's ordered covers. Their exact numbers are E_u and S_v, independent of those positions by slot permutation. Because the marked vertices are in different children, no joint same-child marked profile is needed. This proves C_uv(T)=L_(E_u,S_v)(c_i r_j). Quotient converse transposes G, replacing c_i r_j by c_j r_i. The reversed physical endpoint order uses the swapped marked profiles. Subtracting gives (3).

The kernel is an antisymmetric outer factor. Equation (3) pairs it with the antisymmetric marked-profile tensor E_u tensor S_v minus S_u tensor E_v. It is not a product of uncontrolled endpoint marginals from the full ambient tournament.

## 3. Exact named shared-label quartic

Let w,u,v be distinct physical old names, with u and v in different children. Suppose deleting w leaves every child nonempty. Use the shared-label coefficient identity from `response/PROOF.md`:

    [a_w a_u c_w c_v] F_T = C_uv(T-w)+C_vu(T-w).          (4)

After deletion, the same quotient Q remains, with the one affected child replaced by its physical vertex-deleted card. Form all q, E and S profiles on these cards. Applying (3) to them gives the exact coefficient

    [a_w a_u c_w c_v](F_T-F_U)
       = L_(E_u,S_v)(K_ij)-L_(S_u,E_v)(K_ij),            (5)

where every object on the right now belongs to the deletion cards. This is an integer identity for the full Boolean response polynomial. It does not assume modulo-four or modulo-eight separation and does not discard any cubic or quartic interaction.

If both marked children remain endpoint-balanced after deleting w, their antisymmetric tensor vanishes. In particular, when all original children are endpoint-balanced and w lies in neither marked child, the coefficient in (5) is exactly zero, for every quotient Q and arbitrary profiles in all other children. Thus cross-child shared quartics can be supported only when the repeated name w lies in one of the two marked children.

This rule applies, in particular, to arbitrary vertex-transitive children, which are pair-free when nontrivial: a homogeneous pair would have two different outdegrees, contradicting regularity. The rule is not asserted for shared quartics whose two singly occurring names lie in the same child; those require joint two-marked profiles. Disjoint-label quartics and every cubic remain part of the full response and are not removed by this reduction.

## 4. Exact pair-free examples

For Q=C3 and three directed-triangle children, with names 0,1,2 in the first child and 3,4,5 in the second,

    [a_0 a_1 c_0 c_3](F_T-F_U)=-24.

The retained complete nine-name response table has 216 nonzero shared-label quartics and 324 nonzero disjoint-label quartics. All 216 shared terms have the repeated name in one of the two different marked children, as predicted. The same marked positions with repeated name6 in the third child give zero.

For the regular prime five-vertex quotient with i->i+1,i+2 modulo five and five directed-triangle children, the same physical coefficient is

    [a_0 a_1 c_0 c_3](F_T-F_U)=-170964.

This is a fifteen-name strongly connected, pair-free example. The coefficient is computed both by the marked run kernel and by full-card endpoint DP. The independent literal checker retains all sixteen physical interpolation counts on each side for every displayed coefficient, including the zero controls. No finite classification or general reconstruction conclusion follows from these examples.

## 5. Reproduction and residual

`verify.py` computes marked profiles by subset Hamilton DP and ordered-cover convolution; it checks all 1,099 labelled child tournaments through order five. For arbitrary Q it enumerates color-word counts by a bounded count-vector recurrence, then evaluates (2)–(5). Full-card endpoint counts supply a separate verification independent of the child-run decomposition. The complete earlier nine-name coefficient table is checked rather than replaced by a summary.

Compile `literal.cpp` with a standard C++17 compiler and run it to reproduce `literal.json` (the local build uses `g++ -std=c++17 -O3 -Wall -Wextra -Werror literal.cpp -o literal`). Every extension has at most seventeen names, so counts are at most 17! and the signed sixteen-term interpolation fits safely in 64-bit integers. It evaluates the actual port counts by integer subset DP and takes ordinary finite differences. Then run `python verify.py > receipt.json` and `python -O verify.py > receipt.optimized.json`.

The open next step is to establish a non-annihilation theorem for the outer kernels on these marked-card tensors, or to exhibit an exact collision showing such a theorem cannot hold. A nonzero marked profile alone is not treated as a separator. The present packet is a structural reduction for genuinely named, nonuniform observables on arbitrary prime quotients, complementary to the coarse-response formulas.
