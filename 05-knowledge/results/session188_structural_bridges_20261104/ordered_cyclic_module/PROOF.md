# A cyclic triangle over an ordered exterior module

## 1. Physical statement

Let T have a homogeneous directed triangle A={p,q,r}, with p->q->r->p. Put R=T-A. Suppose R has a module M=X union Y with nonempty sides X->Y; thus every name in E=R-M sees all of M uniformly. The triangle sees the sides cyclically:

    A->X->Y->A.

Choose physical names i in X and j in Y. Their complete one-marked ordered-cover profiles must be endpoint-balanced: E_i^X(h)=S_i^X(h)=b_X(h), and E_j^Y(h)=S_j^Y(h)=b_Y(h), with b_X(1),b_Y(1)>0. Covers have h ordered nonempty directed paths, and a marked endpoint belongs to one prescribed slot. Vertex-transitive X,Y satisfy these hypotheses. No condition is imposed on E or on its incidences toward A.

Let T^A reverse only the triangle. Then the full named two-port responses satisfy the strict exact disjoint-quartic inequality

    [a_p a_i c_q c_j](F_T-F_(T^A))>0.                 (1)

If A, X and Y are children of a fixed quotient and the whole quotient is conversed, with all child interiors held fixed, the same coefficient is obtained. More generally the same restoration holds whenever the target has the same unmarked child q profiles, the same balanced profiles at i,j, and the same orientation of A after quotient converse. These are statements about the exact physical coefficient, not a blockwise or quadratic truncation.

The proof below yields a nonnegative rational generating function. In particular it proves positivity at all orders, without an assumption that the exterior is transitive or that any child has two vertices.

## 2. Eliminate the arbitrary exterior

Work in the squarefree algebra on variables z_v for v in E. Let G_uv sum directed simple paths from u to v in E, including singleton paths. Let E_v=sum_u G_uv and S_u=sum_v G_uv. Empty paths are represented only by the explicit constants below.

Write I_A={v in E:v->A}, O_A={v in E:A->v}, and I_M={v in E:v->M}, O_M={v in E:M->v}. Define the nonnegative E-path polynomials

    alpha = 1+sum_(v in I_A) E_v,    beta = 1+sum_(u in O_A) S_u,
    a = 1+sum_(v in I_M) E_v,        b = 1+sum_(u in O_M) S_u,
    sigma = sum_(u in O_M,v in I_M) G_uv,
    mu = sum_(u in O_A,v in I_M) G_uv,
    kappa = sum_(u in O_M,v in I_A) G_uv,
    tau = sum_(u in O_A,v in I_A) G_uv.

Introduce commuting variables x,y counting successive X- and Y-runs. These variables are not squarefree. Put

    h=x+y+xy,       D=1-sigma*h,
    U=mu*h+x(1+y), V=kappa*h+y(1+x).

After eliminating E, the X,Y transition matrix is [[0,1],[0,0]]+sigma*J. Its run-path matrix N is

    N_xx=x(1-sigma*y)/D,  N_yy=y(1-sigma*x)/D,
    N_xy=xy(1+sigma)/D,   N_yx=xy*sigma/D.

Its column sums, row sums and total are

    C_x=x/D, C_y=y(1+x)/D,
    R_x=x(1+y)/D, R_y=y/D, H=h/D.

The same Schur/run decomposition as in cyclic_child_disjoint/PROOF.md gives the quantities entering the triangle kernel on the full exterior R:

    alpha_R=alpha+a(C_y+kappa*H),
    beta_R=beta+b(R_x+mu*H),
    tau_R=tau+N_xy+mu*C_y+kappa*R_x+mu*kappa*H,
    E_x=a*C_x, E_y=a*C_y, S_x=b*R_x, S_y=b*R_y,
    rho_x=N_xx+mu*C_x, rho_y=N_xy+mu*C_y,
    ell_x=N_xy+kappa*R_x, ell_y=N_yy+kappa*R_y.

Every decomposition includes empty exterior segments and arbitrary later returns to X,Y through E. The inverse D is a formal series; sigma has positive E-degree, so it is nilpotent in the squarefree E-algebra. These formulas are exact, not a bounded walk approximation.

## 3. A coefficientwise positive identity

The exact disjoint triangle kernel is

    K=2 alpha_R beta_R (N_xy-N_yx)
      -(1+2 tau_R)(E_x S_y-E_y S_x)
      +alpha_R(rho_x S_y-rho_y S_x)
      +beta_R(E_x ell_y-E_y ell_x).

Direct substitution and cancellation give the identity

    K/(xy) = (2 alpha beta+alpha b+a beta)/D
           +[alpha b U+a beta V
             +ab((kappa+mu+2 tau)h+2x+2y+xy)]/D^2
           +2abUV/D^3.                                  (2)

For completeness, this can be checked by grouping the numerator over D^3 into the four independent products alpha beta, alpha b, a beta and ab. The first is 2D^2. The next two are D(D+U) and D(D+V). The last is

    D((kappa+mu+2 tau)h+2x+2y+xy)+2UV.

All four groups reproduce the displayed expression. This check is an identity over the integer polynomial ring in the eight interface quantities and x,y; it does not rely on tournament inequalities.

The series expansions of D^(-1), D^(-2), D^(-3) have nonnegative integer coefficients. Every interface quantity is nonnegative, hence (2) is coefficientwise nonnegative in x,y and the E variables.

Apply the functional x^s y^t -> b_X(s)b_Y(t), with b(0)=0, and extract the full E-support. This gives exactly (1). Indeed each skeleton has a unique prescribed X-run and Y-run carrying its respective physical marked endpoint. Permuting the other run slots gives the ordered-cover counts, without an additional factor s or t. Endpoint balance permits the same b whether the marked run starts or ends at its physical name. No unmarked covering coordinate or jointly marked response has been discarded.

The coefficient of xy in K is 2 alpha beta+alpha b+a beta. Moreover [E](alpha beta)=H(E with one extra vertex having A's incidences), which is positive by existence of a Hamilton path in every tournament. Therefore

    (1) >= 2 b_X(1)b_Y(1) H(E+A_singleton)>0.           (3)

This proof also covers E empty, singletons X or Y, and arbitrarily large or nontransitive sides. For E empty the exact polynomial is

    K=xy(4+3x+3y+5xy+2x^2 y+2xy^2+2x^2 y^2).

The theorem fails to cover a transitive attachment X->A->Y by this same positive argument; that different interface is not silently included.

## 4. Prime-quotient original-module root families

Let Q be a prime tournament of order m>=3. Suppose it has three distinct names z,u,v with z->u->v->z and {u,v} a module of Q-z. Inflate z by C3, u and v by arbitrary nonempty vertex-transitive children, and all other names by arbitrary nonempty children. Fix all physical names.

For every XOR-family of modules of this original substitution that includes the whole root, the source and target full named responses differ. One fixed battery of at most m+20 actual queries suffices against every such target.

Proof: proper original modules lie inside individual children. If scalar H increases under the family, the scalar query separates. Otherwise module_decomposition/PROOF.md gives equality of all corresponding child q profiles. vertex_transitive_extremes/PROOF.md preserves every physical one-marked profile in the two marked children under their internal original-module families, including their whole converse. If the triangle child is selected, it is restored after the root converse; the exact marked-profile restoration in cyclic_child_disjoint/PROOF.md and (1) give a strictly positive disjoint quartic. Its sixteen finite-difference queries include (0,0), so no extra scalar query is needed. If the triangle child is unselected, queries uniform on each other child restore the target to the whole-converse response, and the exposed-anchor battery in returns_next/UNIFORM.md detects it. There are m+2 compressed blocks, giving at most m+4 anchor queries. Their union therefore has at most m+20 queries.

This is root-family detection relative to the original source. For arbitrary vertex-transitive marked children it is not a claim of full orbit injectivity or of a common reconstruction battery.

## 5. All critical all-C3 families

The elementary concrete T,U,W families in the earlier packets all have the cyclic ordered-module witness above:

- T_(2r+1), arcs forward distance 1 through r: take z=0,u=r,v=r+1. The pair is homogeneous after deleting 0.
- U_(2r+1), reverse chain A_0,...,A_(r-1), forward chain B_0,...,B_r, with A_i->B_j iff j<=i: choose z=A_i,u=B_i,v=B_(i+1), for any 0<=i<r.
- W_(2r+1), forward chain 0,...,2r-1 and a hub beating even names and losing to odd ones: take z=hub,u=2s,v=2s+1, for any 0<=s<r.

These witnesses include r=1. The separate critical_quotient_synthesis/PROOF.md identifies the three concrete families with the classical critical-prime classification; that sourced classification is not re-proved or enlarged here.

Consequently every T,U,W quotient inflated entirely by directed triangles has all 2^(m+1) original-module XOR states pairwise distinct in its full named two-port response. Its only proper nontrivial modules are the triangle children. A difference family omitting the root is separated by an isolated selected triangle. A difference family including it is separated by Section4. The same argument works from every state: a quotient converse reverses the ordered witness pair, and swapping u,v restores its cyclic convention; child converse preserves its marked profiles. No relabeling of physical responses is used to identify two states.

The result is an all-order pair-free original-module orbit theorem. It is not arbitrary-prime quotient closure or unrestricted named reconstruction.

## 6. Scope

The exact positivity theorem allows arbitrary exterior E and any nonempty balanced marked sides. The root-family theorem additionally needs prime quotient and vertex-transitive marked children to justify preservation under their internal original-module XORs. The full orbit corollary uses all-C3 children.

The local triangle-converse assertion includes arbitrary exterior E, so it applies in any larger tournament where its two module hypotheses remain true. The quotient-converse and root-family assertions are whole-core statements, with separately justified ordered strongly-connected-component padding allowed. No arbitrary strong outer-context lift of those latter comparisons is claimed: a positive exact quartic need not be nonzero modulo eight, and the even-module residue transfer does not preserve its sign automatically. No global closure is asserted.
