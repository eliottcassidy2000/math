# Exact weighted tournament stability from order 238

Status: **PROVED**, conditional only on the accepted Session138 inequalities
listed below. Same-agent proof review with exact rational certificates and
independent native graph tests; no independent-agent or formal audit of this
final combination is claimed.

## Statement and combined provenance

For every real weighted tournament A of order n>=238, let

    Delta=H(A)-det(I+A)-|det A|,
    V=sum_{i<j} A_ij A_ji, R=floor(n/2),
    kappa_n=(4/R)[1-(3/4)^R-1_{n even}(1/4)^R].

Then Delta>=kappa_n V, sharply. Equality with V>0 occurs precisely at
ordinal sums of R fair pairs and, at odd order, one singleton. At V=0,
equality is precisely transitivity. Thus K_n=kappa_n for every n>=238;
together with the inherited orders3..9 the unresolved interval is10..237.

The source chain is explicit:

- Session138 supplies binary rounding, matching support, and equality.
- Concurrent Session139 supplies the exact light-mass dual and per-edge
  geometric-tail argument in its526 proof.
- This packet's [466 proof](STABILITY_466.md) keeps the -BR baseline;
  its [320 proof](STABILITY_320.md) frees the local dual parameters.
- Incoming [discarded-degree cancellation](../stability_priority_moduli_20261102/STABILITY_N500.md),
  used in that packet's276 proof, improves W>=D^2/(2n) to W>=2D^2/n.

Combining the last two mechanisms yields238. This is a sufficient cutoff,
not a proof that237 fails. The earlier320,314,276 and other cutoffs remain
valid historical bounds. We make no external novelty claim.

## 1. The two graph inequalities retained

Independently round A to a binary tournament T and let B=E Delta(T).
For q_e=A_ij A_ji, let W=sum q_e q_f over unordered incident edge pairs.
Select a descending-weight greedy matching. Its discarded edges have total
weight D and weighted degrees d_v. Every discarded edge f was deleted by
a selected edge e with q_e>=q_f; the selected/discarded wedge pays q_f^2.
These wedges are distinct. The discarded/discarded wedges are exactly
(1/2)sum d_v^2-sum_discarded q_f^2. Therefore

    W>=(1/2)sum_v d_v^2>=2D^2/n.                         (1)

This is the incoming cancellation lemma, repeated here to expose every
dependency. It holds for all nonnegative weighted simple graphs.

For a selected edge of weight w, let S be the mass assigned to it and
W_i its within-group wedges. With endpoint-star masses X,Y, greediness
gives sum_assigned q_f^2<=wS, so

    W_i>=wS,
    W_i>=wS+(X^2+Y^2-sum_assigned q_f^2)/2
        >=wS/2+S^2/4.                                  (2)

Different groups use disjoint edge pairs. Equations(1),(2) are separate
valid bounds on the same W; their right sides are not added together.

The accepted tournament inputs are

    B>=W>=2V^2/n-V/4,
    Delta>=B+f(V), f(V)=(1-exp(-4V))/4,
    Delta>=B+kappa_n M for any retained matching mass M,
    Delta>=B+a_h M if it has h edges and n>2h,
    a_h=(4/h)[1-(3/4)^h].

See the named Session138 dependency links in [the first proof](STABILITY_466.md).

## 2. At least ten heavy selected edges

Assume V>0 and Delta<=kappa V, where kappa=kappa_n. Now R>=119.
Set

    c=999/4000, t=4/119, d=1793/10000, lambda=t^2/d.

Since exp(8)>47259/35>1000, f(2)>c. The decreasing ratio f(V)/V and
kappa<=4/119<c/2 exclude V<=2. Thus

    kappa V>B+c, B<=kappa D, V<=n/8+n kappa/2.            (3)

The middle inequality retains the full greedy matching; the last follows
from the global wedge lower bound.

Suppose L>=1 selected weights satisfy w_i<=t, and let U be their exact
total weight. Deleting them retains h<=R-L edges and mass M=V-D-U.
By (1),(3), U<=tR and completion of a square,

    M>c/kappa+2D^2/(n kappa)-D-tR
     >=c/kappa-n kappa/8-tR.

Since each weight is at most1/4,

    h>=4M>(c-4t)R-4-2/R=:g(R).                          (4)

Here we used 1/kappa>=R/4 and n kappa/2<=4+2/R. The function g increases
for R>=119 and

    g(119)=4618839/476000>9.

Consequently h>=10 and m:=R-L>=h. In particular n>2h and M>0.

## 3. Pay the geometric tail without discarding B or U

The coefficients a_h decrease with h. For epsilon=(3/4)^m,

    (kappa V-B)(R-L)>=4(1-epsilon)(V-D-U).

Expanding as in the320 proof, and using kappa R<=4,

    cL-4U<4epsilon V+4D-BR.                              (5)

From (3),

    4V<=Q_R:=R+33/2+8/R.

For 1<=L<=R-10, the sequence (3/4)^(-L)/L has increasing successive
ratios (4/3)L/(L+1). Its maximum on this interval lies at an endpoint.
Therefore

    4epsilon V/L
       <=max{Q_R(3/4)^(R-1), Q_R(3/4)^10/(R-10)}<c-d.   (6)

Both endpoint expressions decrease with R>=119. The first has successive
ratio less than(3/4)(1+1/R)<1. For the second,

    Q_R/(R-10)=1+(53/2+8/R)/(R-10),

whose positive excess decreases. At R=119 both strict comparisons are
rational; the tighter slack is

    (c-d)-[Q_119/109](3/4)^10
       =6985054879/17001349120000>0.

Equations(5),(6) give the baseline-preserving light budget

    sum_light(d-4w_i)=dL-4U<4D-BR.                       (7)

Since d-4t>0, D>0.

## 4. The exact dual contradiction

For 0<=w<=t and S>=0,

    wS/2+S^2/4-tS+lambda(d-4w)
      =(S/2-t+w/2)^2+w(t-4lambda-w/4)>=0.               (8)

The identity uses lambda d=t^2; nonnegativity follows from

    3t/4-4lambda=101/25390673>0.

Use (8) on light groups and W_i>=w_i S_i>=tS_i on heavy groups.
Their wedge contributions are disjoint. With (7),

    B>=W>(t-4lambda)D+lambda BR,
    t-4lambda=213468/25390673>0.                         (9)

If1-lambda R<=0, (9) is immediately impossible. Otherwise (3) and
kappa<=4/R give

    (t-4lambda)D < B(1-lambda R)
                  <=(4/R-4lambda)D
                  <=(t-4lambda)D,

because R>=119=4/t. This is a contradiction even at the endpoint; it
does not require a small-tail estimate for4-kappa R.

Thus L=0. All selected weights exceed t>=4/R>=kappa. If D>0, their
selected/assigned wedges give B>=W>kappa D, contradicting (3). Hence
D=B=0. The uncertainty graph is a matching; the accepted matching theorem
gives the bound. Zero binary gap at n>=4 forces supported transitive
completions; the accepted H=1 scalar classification gives exactly the
maximal fair-pair equality locus and attainment.

## Audit and next boundary

`verify_stability_238.py` checks every rational endpoint, the dual identity,
and the stronger graph bound on all unweighted graphs through order5,
all729 four-vertex graphs with edge weights0,1,2, and balanced complete
bipartite controls. It also replays the scalar certificate at238..1499.
The all-order conclusion is the argument above, not that finite replay.

The coarse scalar proof still does not decide10..237. The next useful
object is actual directed incidence or a tighter joint distribution of
selected weights and discarded masses; native H, full determinants, and
the same source tournament must survive that next compression.
