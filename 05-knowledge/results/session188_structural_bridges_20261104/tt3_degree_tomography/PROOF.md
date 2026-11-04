# Transitive-triple degree tomography and an equal-profile exterior skew

**Research packet: independent review pending.** The identities below are derived from the frozen full quartic normal form and exact ordered-module Schur decomposition. They do not replace those sources. No unrestricted reconstruction or arbitrary-prime closure is claimed.

## 1. Three marks inside a physical triple

Let A={p,q,r} be a module of a tournament T, with arbitrary exterior R. The names p,q,r are distinct. For i in R, consider the complete physical quartic

    L_r,i(T)=[a_p a_q c_r c_i] F_T.

On R let G_uv be its squarefree simple-path polynomial, including singletons, E_v its column sum and S_u its row sum. Put

    alpha=1+sum_(v->A) E_v, beta=1+sum_(A->u) S_u,
    tau=sum_(A->u,v->A) G_uv,
    ell_i=sum_(v->A)G_iv, rho_i=sum_(A->u)G_ui.

Then, writing d_A(r) for r's indegree inside A,

    L_r,i(T)=[R]{alpha S_i(d_A(r)+2tau)
                  +beta E_i(2-d_A(r)+2tau)
                  +2alpha beta(ell_i+rho_i)}.             (1)

This is valid for all eight orientations of the three physical internal edges, not only transitive triples. Consequently, if T' changes only the orientation inside A,

    L_r,i(T)-L_r,i(T')=(d_A(r)-d_A'(r)) K_i,
    K_i=[R](alpha S_i-beta E_i).                           (2)

Proof: after eliminating R, the effective A-path matrix is
M=(I-Z(B_A+tau J))^(-1)Z. Write C,R for its column and row sums, and H for its total. The eight-term three-cover normal form gives

    alpha(C_p M_qr+C_q M_pr)(S_i+beta H ell_i)
    +alpha beta(C_p R_q+C_q R_p)R_r rho_i
    +alpha beta C_r(C_p R_q+C_q R_p)ell_i
    +beta(E_i+alpha H rho_i)(M_rp R_q+M_rq R_p).

Extract the full three-name support. Every displayed H correction has degree at least four and vanishes. The first bracket has coefficient d_A(r)+2tau; the last has coefficient 2-d_A(r)+2tau. Each triple product has coefficient two. This proves (1), with all exterior support still extracted jointly.

Adjoin a virtual singleton z to R with A's original incidences. Then K_i is exactly C_zi(R+z)-C_iz(R+z), where C_uv counts ordered spanning two-path covers with first path ending at u and second starting at v. Thus (2) is degree tomography through an exterior two-cover skew, not a product of independent marginals.

## 2. A binomial positivity lemma

Let q_0=0 and q_1,...,q_m be positive, with q_k=0 outside 1<=k<=m. Assume q is log-concave: q_k^2>=q_(k-1)q_(k+1). Define L(x^s y^t)=q_s q_t. For every r>=0,

    L(-(x-y)^2(x+y+xy)^r)>=0.                            (3)

Indeed expand (x+y+xy)^r as sum_c binom(r,c)(xy)^c(x+y)^(r-c). Write N=r-c+2. The corresponding summand is

    sum_(j=0)^N binom(N,j)
      [N-(2j-N)^2]/[N(N-1)] q_(c+j)q_(c+N-j).           (4)

Log-concavity and interval support imply that g(j)=q_(c+j)q_(c+N-j) is symmetric and nondecreasing toward the center. The same holds for f(j)=N-(2j-N)^2. Under the binomial weights binom(N,j)/2^N, the mean of f is zero. Their covariance is nonnegative: pair two independent indices and average (f(j)-f(k))(g(j)-g(k)), which is nonnegative because both depend monotonically on distance from the center. Therefore the weighted mean of fg is nonnegative, proving (3). Zeros at the support boundary cause no exception. For r=0 the value is exactly 2q_1^2.

## 3. Strict skew for equal log-concave witnesses

Suppose R has a module X union Y, with X->Y, A->X and Y->A. Put E=R-(X union Y). Assume the COMPLETE ordered nonempty path-cover profiles of X and Y are identical to a log-concave q as in Section2. No condition is imposed on E or its incidences toward A.

Set K_X=sum_(i in X)K_i and K_Y=sum_(i in Y)K_i. Then

    K_X-K_Y >= 2 q_1^2 H(E+z_M)>0,                     (5)

where z_M is a singleton with the common incidence of X union Y toward E.

Use the eight nonnegative E-interface polynomials alpha,beta,a,b,sigma,mu,kappa,tau defined in ordered_cyclic_module/PROOF.md, Section2. Here alpha,beta refer to A; a,b to X union Y. Put h=x+y+xy, D=1-sigma h, U=mu h+x(1+y), V=kappa h+y(1+x). Exact elimination gives

    K_X=L[E]{x[b alpha(1+y)-a beta]/D
                 +ab x[V(1+y)-U]/D^2},
    K_Y=L[E]{y[b alpha-a beta(1+x)]/D
                 +ab y[V-U(1+x)]/D^2}.                (6)

The aggregate endpoint sums use q_s(X),q_t(Y), with no extra slot multiplier. The functional is symmetric because the two complete profiles agree. Hence the terms proportional to (x-y) and (x-y)h vanish, leaving

    K_X-K_Y=L[E]{xy(alpha b+a beta)/D
      +ab[-(x-y)^2+xy(x+y+2xy)+(mu+kappa)xyh]/D^2}.     (7)

All terms other than -(x-y)^2 are coefficientwise nonnegative. Expand D^-2=sum_(r>=0)(r+1)sigma^r h^r; sigma is nilpotent in the E-squarefree algebra. Section2 makes each remaining functional nonnegative. Its r=0 contribution is 2q_1^2[E](ab). The exact path decomposition at z_M gives [E](ab)=H(E+z_M)>0. This proves (5) at every order.

No equal-profile assumption is silently dropped. Unequal TT sides admit K_X=K_Y=806 and even K_X-K_Y=-18 in exact nine-name controls. Neither implies that both skews vanish. Individual K_X can also be negative or zero.

Transitive tournaments of every order have log-concave q_k=k!S(n,k). An elementary proof uses q_(n,k)=k(q_(n-1,k)+q_(n-1,k-1)): convolution with (1,1) preserves a nonnegative interval-supported log-concave sequence, and multiplication by k preserves log-concavity. For convolution, expand its log-concavity difference into the two adjacent old differences and q_k q_(k-1)-q_(k-2)q_(k+1), all nonnegative by the old ratio ordering. Directed triangles also qualify, with q=(3,6,6).

## 4. Proposed root-family consequence, awaiting independent scope review

Let Q be prime, with distinct z,u,v such that z->u->v->z and {u,v} is a module of Q-z. Inflate z by a physical TT3, and u,v by arbitrary nonempty children with the same log-concave complete profile. Other children are arbitrary nonempty tournaments. Fix any ambient root context. Consider an XOR family of ORIGINAL source modules inside this core, containing the whole core.

If scalar H strictly increases, the scalar physical query detects it. Otherwise hereditary original-module monotonicity gives equal corresponding child profiles. In particular the target distinguished triple is transitive, though its physical order may have changed.

The all-cover top-quartic converse lemma and uniform restoration of every other child identify a quotient-converse quartic comparison with a local comparison between source A and the converse of target A. This applies also to the saturated 3+1 quartics: identify c_i within X or within Y, while the other two a marks and one c mark remain three distinct names of A. Bidegree at most(2,2) excludes extra alias routes. Enlarge E by the ambient names; X union Y remains a module after deleting A.

By (2) and (5), equality of all six quartics (three choices of r, two marked sides) forces equality of all three indegrees between A and the converse of target A. A transitive triple is uniquely determined by its named indegree vector. Thus target A must be precisely A^op. Queries nonuniform only on A then restore all other children to the whole-converse comparison; the established uniform exposed-anchor separator detects it, with its physical ambient lift.

This would prove source-relative root-family detection, not arbitrary pairwise reconstruction of outputs with larger child profiles. A conservative battery is the 96 finite-difference queries for the six quartics, a scalar query, and at most m+4 compressed anchor queries, totaling at most m+101; duplicate queries need not be removed. The concrete critical T/U/W witnesses from critical_quotient_synthesis make all-TT3 sources an unbounded family. No completeness at arbitrary prime quotients is asserted.

This section remains explicitly pending review of the restoration/anchor integration. Sections1--3 should also be independently checked before being cited as a new theorem. The frozen TT3 circulation and marked-path-duality packets are unaffected.
