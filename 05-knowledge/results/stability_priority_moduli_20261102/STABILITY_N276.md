# Sharp weighted tournament stability for every n >= 276

Status: **PROVED** by the uniform argument below, with exact arithmetic
verification. Not yet independently reviewed by another agent or person.
This is the current result of this packet; `STABILITY_N500.md` preserves
its earlier intermediate proof. The remaining interval is **10..275 OPEN**.

## Inheritance and the theorem

For every real weighted tournament A of order n>=276, put

    Delta=H(A)-det(I+A)-abs(det A),
    V=sum_(i<j) A_ij A_ji,  R=floor(n/2),
    kappa_n=(4/R)[1-(3/4)^R-1_(n even)(1/4)^R].

Then Delta>=kappa_n V, with sharp coefficient. For V>0 equality holds
exactly at ordinal sums of R fair pairs and, for odd n, one singleton.
At V=0 equality is the inherited deterministic transitive case.

The dependencies and definitions of B=E Delta(T), W, greedy matching,
and assigned discarded masses are specified in `STABILITY_N500.md`.
In particular the accepted Session138 inputs are

    B>=W>=2V^2/n-V/4,  Delta>=B+g(V),
    g(V)=(1-exp(-4V))/4,
    Delta>=B+kappa_n M for any matching of mass M,
    Delta>=B+a_m M for a nonperfect matching of size m,
    a_m=(4/m)[1-(3/4)^m].

This packet's discarded-degree cancellation proves W>=2D^2/n, where
D is the mass discarded by greedy matching. Its disjoint group budget gives

    B>=eta D-L eta^2                                   (1)

for any eta>0, where L selected edges have weight at most eta.
These two graph lemmas are proved in `STABILITY_N500.md`, section 1;
they hold beyond tournament variance weights as well.

**Incoming work incorporated before publication:** commit `9c7eb0b13e`,
`../session139_incidence_moduli_and_stability_20261102/stability/PROOF_N526.md`,
sections 2 and 4, supplies the retained-baseline square completion and the
convex endpoint bound for a geometric tail *per removed edge*. Its n>=526
theorem uses the earlier weaker W>=D^2/(2n). Here the stronger factor-four
bound, retention of B in both matching estimates, and the same per-edge
tail mechanism give n>=276. The prior mechanism is credited, not presented
as newly discovered by this packet. The original Session138 and Session139
proofs and frozen audits are unchanged.

## 1. Candidate configurations and coupled matching mass

Assume V>0 and Delta<=kappa V, with kappa=kappa_n. Set

    c=24999/100000, eta=2/69, a=7/60, delta=1/59.

The exact Taylor check exp(12)>25000 gives g(3)>c. Since g(V)/V decreases,
V<=3 would imply Delta/V>c/3>8/275>=kappa. Therefore V>3 and

    kappa V-B >= g(V)>c.                               (2)

The wedge estimate gives V<=n/8+n kappa/2<=n/8+4n/(n-1).
Retaining the full greedy matching gives B<=kappa D. The new degree bound,
combined with (2), also yields

    V-D > c/kappa + 2D^2/(n kappa)-D
        = c/kappa - n kappa/8
          + [2/(n kappa)](D-n kappa/4)^2
        >= c/kappa - n kappa/8.                        (3)

The square completion retains dependence between V and D. Bounding them
independently loses three units of matching mass at the asymptotic scale.

## 2. At least fifteen retained heavy edges

Suppose L>=1. Remove selected edges of weight <=eta. The remaining matching
has size r_h<=R-L and mass M_h>=V-D-eta L. From (3),

    4M_h > c(n-1)/2 - 4n/(n-1) - 2eta n.               (4)

The right side is increasing, since its derivative is
c/2-2eta+4/(n-1)^2>0. At n=276 it exceeds 14 by 31599/88000.
Every edge has variance <=1/4, so r_h>=4M_h>14. Thus

    r_h>=15,  m:=R-L>=15,  n>2m.

The coefficient a_j is the mean of the first j terms of the decreasing
sequence 1,3/4,(3/4)^2,..., hence decreases. Keeping B and replacing r_h
by the larger m therefore gives, with epsilon=(3/4)^m,

    (kappa V-B)(R-L) >= 4(1-epsilon)(V-D-eta L).

The retained mass lower bound is positive. Rearrangement and kappa R<=4 give

    L[kappa V-B-4eta(1-epsilon)]
       <= 4epsilon V +4D-RB.                           (5)

## 3. The tail is small per removed edge, uniformly in n

For either parity n<=2R+1 and

    4V<=n/2+16n/(n-1)<=Q_R:=R+33/2+8/R.

For 1<=L<=R-15, convexity of
log((3/4)^(R-L)/L) bounds its maximum by its two endpoint values. Hence

    4epsilon V/L <= max( Q_R(3/4)^(R-1),
                         Q_R(3/4)^15/(R-15) ) < delta. (6)

Both expressions decrease for integer R>=138. For the first,
Q_(R+1)<Q_R+1 and Q_R>R imply a successive ratio less than
(3/4)(1+1/R)<1. For the second,

    Q_R/(R-15)=1+(63/2+8/R)/(R-15),

whose positive excess decreases. The endpoint inequalities at R=138
are exact rational checks; the tighter slack is

    1/59-Q_138(3/4)^15/123
        =18763479311/119479547723776>0.

By (2), the bracket in (5), after charging delta per removed edge, exceeds
c-4eta-delta. Its excess over a is

    c-4eta-delta-a=58643/135700000>0.

It follows that

    aL<4D-RB.                                         (7)

In particular D>0.

## 4. One shared budget forces the uncertainty graph to be a matching

Combine (1) and (7):

    D eta(a-4eta)<B(a-R eta^2).                         (8)

Here a-4eta=1/1380>0. If a-R eta^2<=0, (8) is impossible.
Otherwise B<=kappa D<=4D/R gives

    B(a-R eta^2) <= D(4a/R-4eta^2)
                 <= D(eta a-4eta^2),

because R>=138=4/eta. This contradicts the *strict* inequality (8),
including the endpoint R=138. Therefore L=0.

Every selected weight exceeds eta>=4/R>=kappa. If D>0, assigned wedges
give W>=sum_i w_i S_i>eta D>=kappa D>=B>=W, impossible. Hence D=B=0.
The matching-support inequality proves Delta>=kappa V. Equality forces
all positive-probability binary completions to be transitive; H=1 and the
inherited scalar theorem give the exact fair-pair equality locus. The
ordinal fair-pair construction attains kappa, proving sharpness.

## Boundaries

This is an all-real-weight, all-n>=276 proof, not a census or an analytic
continuation from small n. It neither resolves n=10..275 nor asserts the
threshold is optimal. The coefficient 2 in W>=2D^2/n is asymptotically sharp
for general weighted graphs; further improvement must use other structure
or a sharper combined budget. The current endpoint estimates stop short
near n=274: with eta=4/137 and a fifteen-edge floor, even the idealized
c=1/4 leaves 1/4-32/137-Q_137(3/4)^15/122<0. This is a limit of this
particular scalar test, not a mathematical obstruction to stability there.
