# Optional rate extension: almost square-root normalized error

**PROVED.** Under the exact same hypotheses as PROOF.md, for every eta>0,

    H(T_n)/mu_n - 1 = o(n^(-1/2+eta)) almost surely.

It suffices to consider 0<eta<1/2. Continue to write kappa-1=Q_n+R_n.
Let B be any oriented incidence matrix of the complete graph K_n: its columns
have a +1 and a -1. If x lists the independent arc signs, the score vector is
Bx. Thus Q_n=x^T C x, where

    C=(2I-B^T B)/[n(n-1)].

C has zero diagonal and trace zero. Because BB^T=nI-J, its eigenvalues are
(2-n)/[n(n-1)] with multiplicity n-1, and 2/[n(n-1)] with multiplicity
(n-1)(n-2)/2. In particular ||C||_op<=1/n for n>=3, and
v:=2 tr(C^2)=Var(Q_n)<=2/n.

For any t with 2|t| ||C||_op<1, Gaussian replacement gives

    E exp(t Q_n) <= E exp(t g^T C g).

Indeed, conditional on all other variables, a zero-diagonal quadratic form
is A x+B. Replacing a sign by a standard Gaussian increases its exponential
moment because cosh(tA)<=exp(t^2 A^2/2). Iterating proves the inequality.
The bound on t ensures finite Gaussian moments at every stage. Diagonalizing
and expanding -log(1-u), using tr(C)=0, gives

    log E exp(t Q_n)
      <= v t^2/[2(1-2|t| ||C||_op)].

Chernoff with t=delta/[v+2||C||_op delta], and separately with -t, yields

    P(|Q_n|>=delta)
      <=2 exp(-delta^2/[2(v+2||C||_op delta)])
      <=2 exp(-n delta^2/[4(1+delta)]).

Together with E R_n^2<=2/[n(n-1)], this proves for epsilon>0

    P(|kappa-1|>epsilon)
      <=2 exp(-n epsilon^2/[16(1+epsilon/2)])
        +8/[epsilon^2 n(n-1)].

Take epsilon=a n^(-1/2+eta) for any fixed a>0. The first bound is summable
because it is eventually bounded by 2 exp(-c_a n^(2eta)); the second is
O_a(n^(-1-2eta)), also summable. Apply Borel-Cantelli for every positive
reciprocal integer a. This proves the claimed little-o rate. As in the main
report, there is no required relation between the tournaments at different
orders, and no independence across n is used.

The complete-convergence result in PROOF.md does not depend on this optional
extension. No law of the iterated logarithm or uniform-tail asymptotic for the
entire path polynomial is claimed here.
