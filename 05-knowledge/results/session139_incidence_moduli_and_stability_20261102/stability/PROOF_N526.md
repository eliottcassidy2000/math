# Light-mass dual certificate: exact weighted stability for n >= 526

Status: PROVED, pending independent review of this packet. This improves
the earlier new shared-budget n>=748 proof, without enumerating graphs.
The remaining unproved orders are10..525 (unless separately resolved).

For R=floor(n/2), write

    kappa_(2R)   = (4/R)[1-(3/4)^R-(1/4)^R],
    kappa_(2R+1) = (4/R)[1-(3/4)^R].

For every n>=526 the theorem is

    H(A)-det(I+A)-abs(det A) >= kappa_n sum_e p_e(1-p_e).

The constant is sharp. Equality with positive variance consists exactly
of ordinal fair pairs, with one singleton in odd order.

## 1. Inputs and notation

Put q_e=p_e(1-p_e), V=sum q_e, Delta=H-det(I+A)-abs(det A),
B=E Delta(T), and W=sum q_e q_f over unordered incident edge pairs.
The accepted Session138 inputs are:

    B>=W,             W>=2V^2/n-V/4,
    Delta>=B+g(V),     g(V)=(1-exp(-4V))/4.

For every retained matching of r positive edges and mass M,

    Delta>=B+kappa_n M;
    Delta>=B+a_r M if n>2r,
    a_r=(4/r)[1-(3/4)^r].

The second coefficient decreases with r: it is the average of the first
r terms of the decreasing sequence1,(3/4),(3/4)^2,... .

Use descending-weight greedy matching. Let w_i be its selected weights,
S_i the total mass of edges discarded at selection i, M=sum w_i, and
D=sum S_i=V-M. All assigned discarded weights are at most w_i. Count
within group i the wedges containing its selected edge, and the wedges
within each of its two endpoint stars. Different groups use disjoint
edge sets, so their counted wedge pairs do not overlap. If the star masses
are X_i,Y_i, then

    W_i>=w_i S_i,
    W_i>=w_i S_i+(X_i^2+Y_i^2-sum_assigned q_e^2)/2
        >=w_i S_i/2+S_i^2/4.                            (1)

Here sum_assigned q_e^2<=w_i S_i and X_i^2+Y_i^2>=S_i^2/2.
In particular W>=D^2/(2n), since there are at most n/2 groups.
These are the accepted greedy-group bounds, with the stronger invoice
retained rather than weakened to S_i^2/4.

## 2. Candidate configurations and retained-baseline coupling

Assume V>0 and Delta<=kappa V, where kappa=kappa_n. Since R>=263,

    kappa R<=4,
    kappa<=8/(n-1),
    kappa<=4/263.

The partial exponential sum sum_(j=0)^7 4^j/j!=16319/315>50 gives
c0=g(1)>49/200>4/263. Since g(V)/V decreases, V<=1 is impossible.
Thus V>1 and kappa V>=g(V)>49/200, whence

    V>49(n-1)/1600>4.

Since exp(16)>50^4>250000,

    g(V)>c,                  c=249999/1000000.           (2)

The wedge lower bound gives

    V<=n/8+n kappa/2<=n/8+4n/(n-1).                     (3)

Retaining the full greedy matching gives

    B<=kappa D.                                         (4)

Retaining B in the nonlinear inequality and completing a square gives

    kappa V>=B+g(V)>D^2/(2n)+c,
    V-D>c/kappa+D^2/(2n kappa)-D
       >=c/kappa-n kappa/2.                             (5)

This coupling is stronger than independent upper and lower estimates on D,V.

## 3. Delete light selected edges, keeping their exact mass

Fix

    d=487/2000,       t=d/8=487/16000,
    lambda=d/64=487/128000,       h=d/16=487/32000.

Suppose L>=1 selected edges have w_i<=t. Let U be their total weight,
and let r_h,M_h be the size and mass of the retained heavy matching.
Then U<=tL<=tn/2. By (5),

    r_h>=4M_h=4(V-D-U)
        >4[c/kappa-n kappa/2-tn/2]
        >=c(n-1)/2-16n/(n-1)-dn/4.                      (6)

The last expression is increasing for n>=526: its derivative is
c/2-d/4+16/(n-1)^2>0. At526 it exceeds17 by964339/1680000.
Since r_h is an integer,

    r_h>=18.                                            (7)

Set m=R-L. Then m>=r_h>=18 and n>2m. Monotonicity of a_r and the
retained matching inequality give, with a=3/4 and epsilon=a^m,

    kappa V>=a_(r_h) M_h>=a_m(V-D-U),
    kappa V(R-L)>=4(1-epsilon)(V-D-U).

Rearranging, then using kappa R<=4 and kappa V>c,

    kappa V L-4(1-epsilon)U
       <=[kappa R-4+4epsilon]V+4(1-epsilon)D
       <=4epsilon V+4D,
    cL-4U<4epsilon V+4D.                                (8)

Unlike a count-only argument, this retains U exactly.

## 4. The finite geometric tail is paid per deleted edge

For every R>=263 and 1<=L<=R-18, we prove

    4V a^(R-L)/L < c-d.                                 (9)

From (3), and because n/2+16n/(n-1) increases for n>=7,

    4V<=n/2+16n/(n-1)<=Q_R:=R+33/2+8/R.

The function log[a^(R-L)/L] has second derivative1/L^2>0. Its maximum
on the real interval[1,R-18] is therefore at an endpoint. Hence the
left side of (9) is at most

    max(Q_R a^(R-1), [Q_R/(R-18)]a^18).                 (10)

Both expressions decrease with integer R>=263. For the first,
Q_(R+1)<Q_R+1 and Q_R>R, so its successive ratio is less than
(3/4)(1+1/R)<1. For the second,

    Q_R/(R-18)=1+[69/2+8/R]/(R-18),

whose positive excess decreases. Thus only the exact endpoints R=263
need checking. Both are strictly less than c-d=6499/1000000. The tighter
one has exact slack

    (c-d)-[Q_263/245](3/4)^18
      =1846127022469099/27674621771776000000>0.

Combining (8),(9) yields the exact light-mass deficit budget

    sum_light(d-4w_i)=dL-4U<4D.                         (11)

In particular D>0 because every light w_i<=d/8 makes d-4w_i>=d/2>0.

## 5. A pointwise dual certificate spends that budget once

For every light group w<=t and S>=0, (1) gives W_i>=wS/2+S^2/4. The
following identity is the new local dual certificate:

    [wS/2+S^2/4]-tS+lambda(d-4w)
       =(S/2-t+w/2)^2+w(2t-w)/4 >=0.                   (12)

Thus W_i>=tS_i-lambda(d-4w_i). For heavy groups, W_i>=w_iS_i>=tS_i.
Summing all groups and using (11),

    W>=tD-lambda sum_light(d-4w_i)
      >(t-4lambda)D=hD.                                (13)

But

    h-4/263 =81/8416000>0,

so W>hD>kappa D>=B>=W, contradicting (4). Therefore L=0.

## 6. Closure and equality

All selected weights now satisfy w_i>t>4/263>=kappa. If D>0, summing
W_i>=w_iS_i gives W>kappa D, again contradicting (4). Hence D=0 and
B=0. The uncertainty graph is a matching, and the accepted matching-support
bound gives Delta>=kappa V. Ordinal fair pairs provide the matching upper
bound on K_n, proving sharpness.

For equality with V>0 the same reasoning gives D=B=0. By the accepted
binary equality theorem (n>=4), every supported completion is transitive.
The accepted H=1 scalar minimization then yields exactly maximal fair-pair
ordinal configurations, with one singleton for odd n.

## Scope and verification

No graph census, floating-point optimization, unproved Pfaffian unit claim,
or Riemann-hypothesis analogy enters this proof. `verify_n526.py` checks
all numerical endpoints and the polynomial identity exactly, with active
checks retained under Python -O. Uniform extension follows from the stated
monotonic arguments, not from a finite sample. This proof alone leaves
orders10..525 open.
