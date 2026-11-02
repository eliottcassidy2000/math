# Unused matching capacity proves sharp stability from order 42

**PROVED, with a same-agent audit and exact rational checks.** This extends
the incoming [order-48 proof](../stability_748_guarded_encodings_20261102/STABILITY_48.md)
to orders42..47. Together with the accepted3..9 and48+ results, the remaining
range is10..41. No independent-agent or formal audit of this extension is claimed.

For a weighted tournament A, put V=sum p_e(1-p_e),
Delta=H(A)-det(I+A)-abs(det A), R=floor(n/2), and

    kappa=(4/R)[1-(3/4)^R-1_(n even)(1/4)^R].

Then Delta>=kappa V for every n>=42. Positive-variance equality consists
exactly of ordinal maximal fair pairs, with one singleton at odd order;
at V=0 equality is transitivity. The coefficient is optimal.

## Inputs and bookkeeping

Use the descending-weight greedy matching. Its selected weights w_i lie
in(0,1/4]; assign each discarded edge to the first selected edge incident
to it. Let S_i be its group's discarded mass, M=sum w_i, D=sum S_i,
V=M+D, B the expected binary Delta under independent rounding, and W the
sum q_e q_f over incident unordered distinct edge pairs. Pad to R groups
with (w,S)=(0,0) if needed. The accepted bounds are

    B>=W>=2D^2/n,
    Delta>=B+(1-exp(-4V))/4,
    Delta>=B+kappa M,
    Delta>=B+1-product_i(1-w_i)-1_(n even)product_i w_i.

The last product subtraction may be inserted even when the matching is
not perfect: zero padding makes it zero then. These are the retained
binary baseline, degree cancellation, and independently audited
[tensor matching theorem](../session139_incidence_moduli_and_stability_20261102/tensor_matching/PROOF.md).
The group contributes at least

    G(w,S)=max(wS,wS/2+S^2/4).

The two expressions are alternatives, never added. Distinct groups have
disjoint counted wedges, as in the order48 proof.

## 1. A capacity potential valid on the full continuum

For0<k<=1/5,0<=w<=1/4,S>=0,

    G(w,S)>=k S-k^2(1-4w).                         (1)

If w>=k, use wS>=kS. If w<k, minimizing G-kS over S yields

    -c_k(w), where
    c_k(w)=(k-w/2)^2           for0<=w<=2k/3,
            2w(k-w)           for2k/3<=w<=k.

In the first interval,

    c_k(w)-k^2(1-4w)
      =w[k(4k-1)+w/4]<=wk(4k-5/6)<=0.

In the second interval completing the square gives

    c_k(w)-k^2(1-4w)
      <=k^2[(1+2k)^2/2-1]<=-k^2/50<0.

Thus(1) is proved without a grid argument. Summing once per group gives

    B>=kD-k^2(R-4M).                              (2)

The potential R-4M is exactly unused matching capacity. It keeps information
that a light-edge count loses. The condition k<=1/5 is sufficient;
the same proof reaches k<=(sqrt(2)-1)/2, but that improvement is not needed.

## 2. A finite rational bootstrap for R=21,22,23

Assume V>0 and Delta<=kappa V. Write

    K_R=(4/R)(1-(3/4)^R), E_R=4^(-R), N_R=2R+1.

For these three R, kappa<=K_R<1/5. Since
(1-exp(-5))/5>K_R and the function(1-exp(-4V))/(4V) decreases,
V>5/4. Therefore kappa V>B+c0, c0=2483/10000.
Both strict exponential comparisons follow from the positive degree12
Taylor lower bound for exp(5), checked rationally in verify.py.

Whenever kappa V>B+c, degree cancellation gives

    M>c/kappa+2D^2/(n kappa)-D
      >=c/kappa-n kappa/8
      >=c/K_R-N_R K_R/8=:m_R(c).                  (3)

AM-GM, with all R slots including zeros, strengthens the tensor reserve to

    Delta>=B+1-(1-M/R)^R-E_R.                    (4)

On0<=M<=R/4 its rightmost term increases. Starting with c0, choose

    c_next=floor(10^6[1-(1-m_R(c)/R)^R-E_R])/10^6.

Each listed step strictly increases c and is justified by(3),(4), even
though the floor may equal the unrounded expression: M>m makes the
reserve comparison strict. Exact integer arithmetic reaches c>=99/100
after8,5,4 steps respectively for R=21,22,23. The complete fractions,
not floating values, appear in verification.json. Thus M>m_R(99/100).

## 3. Concavity closes the entire compatible mass interval

By(2), using k=kappa, a counterexample would require

    F_n(M):=1-(1-M/R)^R-1_(n even)4^(-R)
              -kappa^2 R-kappa(1-4kappa)M <=0.    (5)

This follows by combining Delta<=kappa(M+D) with(2),(4), using the
parity-appropriate subtraction in(4). F_n is strictly concave on[0,R/4].
Its endpoint at R/4 is exactly zero by the definition of kappa.
For each n=42..47, the exact rational check at m_R(99/100) is positive
(greater than1/25). Concavity therefore forces F_n(M)>0 for
m_R(99/100)<=M<R/4, a contradiction. The only possible equality mass
is M=R/4, which forces all R selected weights to equal1/4.

Finally B<=kappa D, whereas W>=sum w_i S_i=D/4 and kappa<1/4.
Consequently D=B=0. Apply the inherited matching equality classification.
The ordinal fair-pair examples attain the coefficient, proving sharpness.
The incoming48+ theorem completes the stated all-order range.

## 4. Failure boundary and next test

This certificate does not reach order41. At R20 its initial bootstrap
value decreases; verification.json records the exact failed comparison.
That is a limitation of this proof, not a counterexample tournament.
An explicit rational scalar configuration at n41 is recorded separately
to test whether these lower bounds alone can settle the next order.
The next useful information is actual degree/group compatibility or
directed-cycle baseline beyond W, not another unsupported averaging step.

The general contribution is(1)--(2), not the constants42 or5. Tensor parity
supplies a nonlinear product; unused capacity turns local losses into one
linear term; AM-GM and concavity then make two endpoint checks decisive.
