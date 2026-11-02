# Exact weighted tournament stability from order 320

Status: **PROVED**, combining accepted Session138 inputs with the local
dual identity of concurrent Session139 and the retained baseline developed
in this packet. Same-agent proof audit and exact rational verification;
no independent-agent or formal audit of this refinement is claimed.

Later refinement: [STABILITY_238.md](STABILITY_238.md) combines this flexible
dual with the incoming stronger discarded-degree bound. The open intervals
below describe this320 proof; the current unresolved interval is10..237.

## Theorem and credit

For every real weighted tournament A of order n>=320, put

    Delta=H(A)-det(I+A)-|det A|,
    V=sum_{i<j} A_ij A_ji, R=floor(n/2),
    kappa_n=(4/R)[1-(3/4)^R-1_{n even}(1/4)^R].

Then **Delta>=kappa_n V**, with equality at V>0 precisely for ordinal
sums of R fair pairs and, when n is odd, one singleton. Hence K_n=kappa_n.
At V=0 equality is precisely transitivity. With the inherited orders3..9,
the remaining unresolved interval is **10..319**.

This improves the sufficient threshold of [our first proof](STABILITY_466.md).
While that checkpoint was being integrated, concurrent Session139 supplied
an [independently audited n>=526 proof](../session139_incidence_moduli_and_stability_20261102/stability/PROOF_N526.md)
at revision9c7eb0b13e. Its new ingredient is a pointwise light-group dual
certificate retaining the exact light mass U. The improvement here combines
that ingredient with the -BR baseline term and frees the dual parameters.
Both contributions are needed for the displayed certificate. The previous
thresholds remain valid sufficient results, not retracted theorems.

## Inputs and putative extremizers

Use the notation and accepted inputs of [STABILITY_466.md](STABILITY_466.md):
B=E Delta(T) under independent binary rounding, W the incident-edge
variance-product sum, and D the mass discarded by descending-weight greedy
matching. Keep its disjoint group contributions W_i, selected weights w_i,
and assigned masses S_i. The stronger accepted group inequality is

    W_i>=w_i S_i,
    W_i>=w_i S_i/2+S_i^2/4.                              (1)

For completeness, split assigned edges between the selected edge's two
endpoint stars, of total weights X_i,Y_i. Their within-star wedges plus
the selected/assigned wedges give

    W_i>=w_i S_i+(X_i^2+Y_i^2-sum_assigned q_e^2)/2.

Greediness gives sum q_e^2<=w_i S_i, and X_i^2+Y_i^2>=S_i^2/2,
proving (1). Contributions from distinct groups use disjoint edge pairs.

Assume V>0 and Delta<=kappa V. Now R>=160 and kappa<=1/40.
The same exp(8) partial-sum argument as in the first proof gives, with
c=999/4000,

    V>2,   kappa V>B+c,   B<=kappa D,
    W>=D^2/(2n),         V<=n/8+n kappa/2.                (2)

Indeed f(2)>c and kappa<c/2 exclude V<=2. Retaining the full greedy
matching gives B<=kappa D. The wedge lower bound gives the last estimate.

## Keep the exact light mass and its baseline

Choose the rational parameters

    t=25001/1000000, d=667/5000, lambda=t^2/d.

Let L be the number of selected weights w_i<=t and U their sum. Suppose
L>=1. Delete those edges, retaining h edges of mass M_h=V-D-U.
Since U<=tR, (2) and completion of a square give

    M_h>c/kappa + D^2/(2n kappa)-D-tR
       >=c/kappa-n kappa/2-tR.

As weights are at most1/4,

    h>=4M_h>(c-4t)R-16-8/R=:g(R).                        (3)

The function g is strictly increasing for R>=160, and
g(160)=98867/12500>7. Thus h>=8. Put m=R-L>=h and
epsilon=(3/4)^m. In particular n>2h and m>=8.
The matching coefficient a_h decreases with h, so

    (kappa V-B)(R-L)>=4(1-epsilon)(V-D-U).

Expand while retaining B:

    L(kappa V-B)-4(1-epsilon)U
      <=[kappa R-4+4epsilon]V+4(1-epsilon)D-BR
      <=4epsilon V+4D-BR.

By (2), this implies

    cL-4U<4epsilon V+4D-BR.                              (4)

## Uniformly pay the geometric tail per removed edge

From (2), n<=2R+1, and kappa<=4/R,

    4V<=n/2+2n kappa<=Q_R:=R+33/2+8/R.

For every integer 1<=L<=R-8,

    4V(3/4)^(R-L)/L
       <=max{Q_R(3/4)^(R-1), Q_R(3/4)^8/(R-8)}<c-d.    (5)

The first bound follows because the successive ratio of
(3/4)^(-L)/L is (4/3)L/(L+1), an increasing function of L;
a sequence with increasing ratios has its interval maximum at an endpoint.
Both endpoint bounds decrease with R>=160. For the first,
Q_(R+1)<Q_R+1 and Q_R>R, so its next/previous ratio is less than
(3/4)(1+1/R)<1. For the second,

    Q_R/(R-8)=1+(49/2+8/R)/(R-8),

whose positive excess decreases. At R=160 both inequalities are strict
rational comparisons. The tighter slack is

    (c-d)-[Q_160/152](3/4)^8
       =1681793/24903680000>0.

Combining (4),(5) yields

    sum_light(d-4w_i)=dL-4U<4D-BR.                       (6)

Here d-4t>0, so D>0.

## A flexible local dual closes the argument

For 0<=w<=t and S>=0, the exact identity is

    wS/2+S^2/4-tS+lambda(d-4w)
       =(S/2-t+w/2)^2+w(t-4lambda-w/4).                 (7)

We used lambda d=t^2. The right side is nonnegative because

    3t/4-4lambda=25001/2900000000>0.

Thus each light group has W_i>=tS_i-lambda(d-4w_i).
Each heavy group has W_i>=w_i S_i>=tS_i. Summing disjoint contributions
once and using (6),

    B>=W>(t-4lambda)D+lambda BR.                         (8)

The coefficient t-4lambda=9075363/1450000000 is positive. If
1-lambda R<=0, (8) is immediately impossible. Otherwise use B<=kappa D
and D>0 to obtain

    kappa>t-lambda(4-kappa R).                           (9)

But lambda<1 and

    0<4-kappa R<=8(3/4)^160<1/100000000.

The last inequality is the integer check 800000000*3^160<4^160.
Since t-1/40=1/1000000,

    t-lambda(4-kappa R)>t-1/100000000>1/40>=kappa,

contradicting (9). Hence no light selected edge exists.

All selected w_i>t>kappa. If D>0, then B>=W>=sum_i w_i S_i>kappa D,
contradicting (2). Therefore D=B=0, uncertainty is a matching, and the
accepted matching-support bound proves Delta>=kappa V. The accepted H=1
classification gives the stated equality cases and attainment, exactly
as in the first proof.

## Boundary, recurrence and next obligation

This proof is uniform for all n>=320. The constants fail to certify the
next orders; that is not a counterexample to the conjectured sharp formula.
The coupled proof now has two negative terms that must survive any further
optimization: the baseline -BR in (6) and the exact light mass -4U.

For all R on the proved large-order range, set
o_R=R K_(2R+1)/4, e_R=R K_(2R)/4. The exact formula gives

    o_R=(3/4)o_(R-1)+1/4,
    e_R=e_(R-1)-(3/16)e_(R-2)+3/16,
    K_(2R+1)-K_(2R)=4^(1-R)/R.

The recurrences apply when every indexed global coefficient is in a proved
range; on the H=1 locus their fair-pair formulas hold at every admissible R.
The even correction is a full pair-cover determinant term, not a theorem
identifying the weighted locus with the merged metagraph's complement edges.

`verify_stability_320.py` independently checks (7) as a polynomial identity,
every rational endpoint and the all-order monotonicity conditions, with
finite replay at orders320..1499. The next substantive problem is to retain
more of the actual (w_i,S_i) distribution or directed-triangle incidence,
rather than merely retune the same coarse mass bounds at orders10..319.
