# Nonlinear matching bootstrap: sharp weighted stability for every n >= 140

Status: PROVED, independently audited and accepted. This new packet
uses the independently accepted tensor matching theorem in
`../tensor_matching/PROOF.md`; it does not modify the frozen n>=748 or
n>=526 proofs. Its improvement is a new nonlinear input, not a refinement
of the old c<1/4 constant alone.

Write R=floor(n/2) and

    kappa_(2R)   = (4/R)[1-(3/4)^R-(1/4)^R],
    kappa_(2R+1) = (4/R)[1-(3/4)^R].

For every n>=140, K_n=kappa_n. Positive-variance equality consists exactly
of ordinal fair pairs, with one singleton in odd order. This proof leaves
orders10..139 unresolved unless separately treated.

## 1. Accepted inputs and the retained nonlinear matching theorem

As before, q_e=p_e(1-p_e), V=sum q_e, Delta=H-det(I+A)-abs(det A),
B=E Delta(T), and W is the sum of q_e q_f over incident unordered edge pairs.
The Session138 inputs give

    B>=W>=2V^2/n-V/4,
    Delta>=B+g(V),       g(V)=(1-exp(-4V))/4.

A descending-weight greedy matching has selected weights w_i, assigned
 discarded masses S_i, total matching mass M, and D=sum S_i=V-M. Its groups
have disjoint counted wedge pairs and satisfy

    W_i>=w_i S_i,
    W_i>=w_i S_i/2+S_i^2/4,
    W>=D^2/(2n).

The retained matching bounds are Delta>=B+kappa_n M, and
Delta>=B+a_r M if its size r has n>2r, where

    a_r=(4/r)[1-(3/4)^r].

The NEW tensor theorem, after first rounding all edges outside a retained
matching, gives

    Delta>=B+1-product_i(1-q_i)-1_{n=2r}product_i q_i
         >=B+1-exp(-M)-4^(-R).                          (1)

The last uniform loss is valid because the product subtraction is present
only if r=R and n=2R, when product_i q_i<=4^(-R). Otherwise that subtraction
is zero. The matching probabilities are unchanged when its complement is
rounded, and the terminal binary expectation is retained exactly once.

## 2. Candidate configurations and a bootstrap rule

Assume V>0 and Delta<=kappa V, where kappa=kappa_n. Since R>=70,

    kappa R<=4,       kappa<=4/R<=2/35,
    kappa<=8/(n-1).

The accepted exponential estimate exp(4)>50 gives g(1)>49/200>2/35.
Since g(V)/V decreases, V<=1 is impossible. Thus V>1 and

    V>49R/800>=343/80>4.

Since exp(16)>50^4>250000, we have

    kappa V>=B+g(V)>B+c_0,       c_0=249999/1000000.       (2)

Whenever kappa V>B+c for c>0, use B>=D^2/(2n) and complete a square:

    M=V-D>c/kappa+D^2/(2n kappa)-D
      >=c/kappa-n kappa/2
      >=cR/4-4-2/R=:m_c(R),                             (3)

because n<=2R+1. The function m_c(R) increases with R. Combining(1),(3)
gives the bootstrap rule

    kappa V > B+1-exp[-m_c(R)]-4^(-R).                   (4)

Starting from c_0, the following three advances hold uniformly for R>=70:

    c_0 -> c_1=7/25 -> c_2=57/100 -> c_3=99/100.

At R=70 the corresponding mass lower bounds m_c(70) are

    969951/2800000,        61/70,        333/56.

To check an advance c->c', take respectively degrees2,3,4 in the positive
exponential Taylor sum T_k(m)=sum_(j=0)^k m^j/j!. Exact rational arithmetic
verifies

    T_k(m_c(70)) [1-c'-4^(-70)] >1.

Since exp(m)>T_k(m), this implies
1-exp[-m_c(70)]-4^(-70)>c'. Monotonicity of m_c and decrease of4^(-R)
extend each step to every R>=70. Therefore every candidate satisfies

    kappa V>B+c,           c=99/100.                    (5)

This is the qualitative strengthening: the binary baseline is retained
while the scalar gap rises from about1/4 to99/100.

## 3. Retained mass, light matching edges, and their deficit budget

The variance-wedge bound gives

    V<=n/8+n kappa/2<=n/8+4n/(n-1).                     (6)

The full greedy matching gives

    B<=kappa D.                                         (7)

Fix

    d=19/20,       t=d/8=19/160,
    lambda=d/64=19/1280,       h=d/16=19/320.

Suppose L>=1 selected greedy edges have w_i<=t, with total weight U.
Delete them, leaving a matching of size r_h and mass M_h. By(3),(5),

    r_h>=4M_h=4(V-D-U)
        >4[c/kappa-n kappa/2-tn/2]
        >=c(n-1)/2-16n/(n-1)-dn/4.                      (8)

The last expression has positive derivative c/2-d/4+16/(n-1)^2.
At n=140 it exceeds19 by12229/27800. Thus integer r_h>=20.
Set m=R-L>=r_h, so m>=20 and n>2m. Since a_r decreases with r,

    kappa V>=a_(r_h) M_h>=a_m(V-D-U).

With a=3/4 and epsilon=a^(R-L), the same exact rearrangement as in the
n>=526 proof, now using kappa V>c, yields

    cL-4U<4epsilon V+4D.                                (9)

## 4. The geometric error is paid per deleted edge

For every R>=70 and 1<=L<=R-20,

    4V a^(R-L)/L < c-d=1/25.                            (10)

Indeed(6) and monotonicity of n/2+16n/(n-1) for n>=7 give

    4V<=Q_R:=R+33/2+8/R.

The function log[a^(R-L)/L] is convex in L, with second derivative1/L^2.
The left side of(10) is at most

    max(Q_R a^(R-1), [Q_R/(R-20)]a^20).

Both expressions decrease with integer R>=70. For the first,
Q_(R+1)<Q_R+1 and Q_R>R imply successive ratio<(3/4)(1+1/R)<1.
For the second,

    Q_R/(R-20)=1+[73/2+8/R]/(R-20),

whose positive excess decreases. Both endpoint quantities at R=70 are
strictly below1/25, by exact rational checks. Combining(9),(10) gives

    sum_light(d-4w_i)=dL-4U<4D.                         (11)

Each summand is at least d/2>0, so D>0 in this branch.

## 5. The existing local dual certificate closes the argument

For each light group w<=t, the stronger wedge invoice and the identity

    [wS/2+S^2/4]-tS+lambda(d-4w)
       =(S/2-t+w/2)^2+w(2t-w)/4>=0

give W_i>=tS_i-lambda(d-4w_i). Heavy groups give W_i>=tS_i. Summing and
using(11),

    W>=tD-lambda sum_light(d-4w_i)
      >(t-4lambda)D=hD.

The exact endpoint margin is

    h-4/70=1/448>0.

Thus W>hD>kappa D>=B>=W, contradicting(7). Hence L=0.
All selected weights then exceed t>kappa. If D>0, summing W_i>=w_iS_i
again gives W>kappa D, a contradiction. Therefore D=B=0, and the uncertainty
graph is a matching. Its accepted sharp matching-support bound proves
Delta>=kappa V; ordinal fair pairs attain it.

For equality with V>0 the same conclusion D=B=0 holds. The accepted binary
equality classification and H=1 scalar minimization give exactly maximal
ordinal fair pairs, with the singleton in odd order.

## Verification and scope

`verify_n140.py` checks every rational bootstrap, endpoint, and polynomial
identity with explicit exceptions retained under Python -O. Its deliberately
lower threshold139 fails the first bootstrap step; that rejects this
certificate at139, not the stability theorem itself. Uniform extension is
proved by the displayed monotonic arguments. No new graph census or
floating-point optimization is used.
