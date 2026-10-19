# A finite odd-part reduction for a cofactor deformation

19 October 2026. Proved reduction; global convergence for the full deformation family is a research hypothesis, not established here.

For a fixed positive integer t define, on positive even n=cp with c=2^a exact and p odd,

    F_t(n)=(c-1)p+2tc-1.

The t=1 case is the globally classified cofactor map. For t>1 this is an explicit arithmetic deformation, not the aliquot sum. Its fixed points are exactly c(2tc-1). All individual core branches have slope (c-1)/c<1, but their intercepts grow with the core, so branchwise contraction alone does not prove convergence.

**Two-step theorem.** Suppose F_t(n)>n and the odd part k of p+1 satisfies k>=2t+1. Then F_t²(n)<n.

**Proof.** Put q=p+1=dk with d=2^v2(q). Growth says q<2tc. Since k>=2t+1, it follows that d<c; both are powers of two, so 2d divides c. The first image is

    n_1=c(p+2t)-p-1=c(q+2t-1)-q,

whose exact dyadic core is d: n_1/d=(c/d)(q+2t-1)-k is odd. Direct calculation gives

    d(F_t²(n)-n)
      =(c+d)(2t(d-1)-p)+2(t-1)d+p+1.

Since p=dk-1>=(2t+1)d-1, the coefficient in parentheses is at most -d-2t+1. Since p+1<2tc,

    d(F_t²(n)-n)
      < -(c+d)(d+2t-1)+2(t-1)d+2tc
       = c(1-d)-d²-d <0.

This proves the theorem for every t and every allowed input, without an exponent bound.

**Finite residual family.** Outside fixed points and immediate descent, every input not covered by this theorem has

    p+1=2^b k,       k odd,       1<=k<=2t-1.

There are only t odd coefficients k to handle for a fixed t, although the two dyadic exponents remain unbounded. At t=1 there is only k=1, precisely the Mersenne-odd-part lane that the global six-step proof resolves. At t=2 the surviving lanes are k=1 and k=3.

This is a structural reduction, not a finite-state proof: the remaining exponent pairs can still carry unbounded information, and transitions between their regions have to be proved. A useful next target is a bounded-block descent theorem outside finitely many explicitly periodic starts for t=2. The broad hypothesis that every F_t orbit is eventually periodic should remain labeled conjectural until those exceptional transitions are controlled. No finite experiment in this note establishes it.
