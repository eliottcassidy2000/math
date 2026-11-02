# A two-regime local invoice proves sharp stability for every n>=48

**Status: PROVED; same-agent proof audit with exact rational certificates and direct disjoint-group checks. No independent-agent or formal audit of this new theorem is claimed.** Accepted inputs are
the tensor reserve, stronger discarded-degree estimate and matching/equality
theorems already reviewed in [STABILITY_66.md](STABILITY_66.md). The
concurrent [n56 proof](../stability_priority_spectrum_20261102/STABILITY_56.md)
supplies the flexible-dual framing. Here exact matching slots and a stronger
local dual, using the maximum of two valid invoices, lower the threshold.

For a real weighted tournament A, write V=sum_e p_e(1-p_e),
Delta=H-det(I+A)-abs(det A), R=floor(n/2), and

    kappa_n=(4/R)[1-(3/4)^R-1_(n even)(1/4)^R].

**Theorem.** For every n>=48, Delta>=kappa_n V with optimal coefficient.
Positive-variance equality is exactly maximal ordinal fair pairs and, for
odd n, one singleton. At V=0 equality is exactly transitivity. Together
with inherited n=3..9, the unresolved interval is10..47.

## 1. The local improvement and its exact scalar boundary

Each greedy group has selected weight w and discarded mass S. Two valid
invoices are W_i>=wS and W_i>=wS/2+S^2/4. Do NOT add them. Use their
maximum

    G(w,S)=max(wS,wS/2+S^2/4).

For every t>0, put d=5t and lambda=t/5. For0<=w<=t and S>=0,

    G(w,S)>=tS-lambda(d-4w).                           (1)

**Proof.** Set u=w/t in[0,1], v=S/t>=0. Dividing the desired difference
by t^2 gives

    max(uv,uv/2+v^2/4)-v+1-(4/5)u.

The two invoices switch at v=2u. Minimizing first over v gives

    u(1/5-u/4)                      if0<=u<=2/3;
    2(u-7/10)^2+1/50                if2/3<=u<=1.

The first is at least u/30, and the second is positive. This proves(1)
on the full continuum, not merely a rational grid.

More generally if lambda=t^2/d, the same calculation with beta=4t/d
shows that this scalar dual is valid iff beta<=2sqrt(2)-2, equivalently
d/t>=2+2sqrt(2). Indeed the two minimum polynomials become u(1-beta-u/4) and
2u^2-(2+beta)u+1. For0<beta<=2sqrt(2)-2 the first is
nonnegative on[0,2/3] and the second is nonnegative everywhere.
Above this boundary, beta<=2 gives a negative value at
u=(2+beta)/4 in[2/3,1]; beta>2 gives a negative value at u=1.
The rational choice5 is safely above that exact threshold.
For d/t=24/5 the scalar values u=17/24,v=17/12 give deficit-1/288.
This hostile is about the two-invoice relaxation, not an actual tournament.

## 2. Candidate bootstrap, uniformly for R>=24

Let B=E Delta(T) for the common independent binary rounding, W the incident
variance-pair sum, M the greedy matching mass and D=V-M. Accepted inputs:

    B>=W>=2D^2/n,     W>=2V^2/n-V/4,
    Delta>=B+(1-exp(-4V))/4,
    Delta>=B+kappa_n M,
    Delta>=B+1-exp(-M)-4^(-R).                         (2)

The last bound is the independently audited tensor matching reserve;
the perfect-matching subtraction is handled exactly as in STABILITY_66.md.
Assume V>0 and Delta<=kappa V. Since kappa<=4/R<=1/6 and
g(1)>49/200>1/6, monotonicity of g(V)/V excludes V<=1. Hence
kappa V>B+49/200. Whenever kappa V>B+c, completion of the square gives

    M>c/kappa-n kappa/8>=cR/4-1-1/(2R)=:m_c(R).        (3)

The following uniform bootstrap follows from(2),(3):

    49/200 -> 7/20 -> 13/20 -> 47/50 -> 99/100.

At R24 the four mass bounds are539/1200,259/240,691/240,5543/1200.
For each advance c->c_next the degree12 positive Taylor sum satisfies

    T12(m_c(24))[1-c_next-4^(-24)]>1.

The verifier supplies exact rational slacks. Increasing m_c(R) and
decreasing4^(-R) extend the comparisons to all R>=24. Thus all candidates
satisfy kappa V>B+c with c=99/100. Also B<=kappa D.

## 3. Exact slots give twelve heavy edges

Choose

    t=1/6, d=5/6, lambda=1/30, h=t-4lambda=1/30.

Suppose L>=1 selected weights are at most t, with total U. Retain r_h
heavy edges with mass M_h=M-U. Since U<=tL<=t(R-r_h) and r_h>=4M_h,

    r_h(1-4t)>(c-4t)R-4-2/R=:g(R).                    (4)

The positive coefficient is1-4t=1/3. At R24,

    g(24)-11/3=1/100>0.

As g increases, r_h>=12 uniformly. Put m=R-L, so12<=r_h<=m<=R-1 and
n>2m. Applying the accepted nonperfect matching coefficient a_(r_h),
decreasing it to a_m and keeping B gives the exact rearrangement

    cL-4U <4D-RB+4V(3/4)^m.                           (5)

This retains the same negative baseline as the concurrent56 proof. The
stronger retained-edge count comes from using the exact slot constraint.

## 4. Pay the tail, then switch local invoices

From B>=2V^2/n-V/4 and B<=kappa V,
4V<=n/2+2n kappa<=Q_R=R+33/2+8/R. Convexity in L and decreasing endpoints give

    4V(3/4)^m/L
      <=max{Q_R(3/4)^(R-1),[Q_R/(R-12)](3/4)^12}
      <c-d=47/300.                                   (6)

For the second monotonicity use
Q_R/(R-12)=1+(57/2+8/R)/(R-12); the first uses
Q_(R+1)<Q_R+1. At R24 the two values are about0.05463 and0.10779,
both strictly below47/300; the actual checks use exact fractions.
Equations(5),(6) imply

    dL-4U<4D-RB.

Because d-4t=1/6>0, this forces D>0. Apply(1) to each light group and
W_i>=w_i S_i>=tS_i to each heavy group. The counted wedges are disjoint: assign a discarded edge to the first
greedy selected edge incident to it. Within each group, wS counts the
selected/discarded pairs. If the two arm masses are a,b, the same group
also has at least wS+(a^2+b^2-wS)/2>=wS/2+S^2/4, since
every discarded weight is at most w. Taking their maximum spends the
group once; no wedge can belong to two different groups.
Therefore

    B>=W>=tD-lambda(dL-4U)>(t-4lambda)D+lambda RB,
    (1-lambda R)B>hD.                                (7)

If lambda R>=1 this is impossible. Otherwise B<=kappa D<=4D/R yields

    h<(1-lambda R)kappa<=4/R-4lambda
      <=1/6-4/30=1/30=h,

also impossible, because R>=24=4/t. Thus L=0. Then w_i>t>=kappa for
every selected edge, so D>0 would imply W>=sum w_iS_i>kappa D>=B.
Hence D=B=0. Matching-support stability and the accepted H=1 equality
classification give the theorem and its exact equality cases.

## 5. Scope and the next real obstruction

This newly resolves48..55 beyond the incoming56 theorem. The chosen
certificate fails at47, where R23 no longer satisfies4/R<=t and the
same slot calculation does not force twelve heavy edges. This does not
refute the desired inequality there. Further progress can retain the full
matching product, the actual light mass, and the unspent group squares;
any mixture of triangle and wedge payments must still count overlaps.

The radical2+2sqrt(2) here is an exact boundary of a two-regime quadratic
optimization. The rational5 is a sufficient surrogate. That is a genuine
numerical mechanism, not an identification with all other occurrences of5,
11 or the golden ratio. The raw-H spectrum and all other major conjectures
remain at their separately stated frontiers.
