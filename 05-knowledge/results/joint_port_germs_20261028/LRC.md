# Positive boundary germs and a height-adapted Gibbs certificate

**PROVED elementary deductions / CITED inherited Gibbs and contact mechanisms /
FINITE-EXACT controls. Same-agent audit.** Local labels L1-L6 belong here.
These are detection and continuation results; LRC(14) remains open.

## L1. A positive first germ replaces reflection cancellation

Let V be a nonempty finite set of positive integer speeds, M=max(V), and
0<delta<=1/2. Write

    f_V(t)=min_(v in V) ||vt||,
    G_delta(V)={t in R/Z : f_V(t)>=delta},
    L_V(delta)=length(G_delta(V)).

The set is a finite union of closed intervals, including possible isolated
points, and excludes0. For a component [a,b] inside(0,1), define

    P(a)={v : frac(va)=delta},
    N(b)={v : frac(vb)=1-delta}.

Both are nonempty. For an isolated point a=b both sets occur at that point.
At delta=1/2, an active speed belongs to both sets.

**Theorem.** For sufficiently small epsilon>0,

    L_V(delta-epsilon)=L_V(delta)+epsilon D_V(delta),
    D_V(delta)=sum_components (1/max P(a) + 1/max N(b)).

In particular D>=0 and

    G_delta nonempty  iff  D>0  iff  (L,D)!=(0,0).

**Proof.** The allowed interval endpoints for a speed v are
(k+delta)/v and (k+1-delta)/v. Only finitely many endpoint-order changes
are possible as delta varies. Near the chosen threshold, all inactive
constraints have slack. At a left endpoint the new lower boundary is

    max_(v in P(a)) (a-epsilon/v) = a-epsilon/max P(a).

At a right endpoint it is b+epsilon/max N(b). The same calculation expands
an isolated point into an interval. Distinct components remain separated
for small epsilon. Summing their changes proves the formula. Every nonempty
component has two positive contributions, and an empty compact good set
stays empty for a small threshold change. QED.

The minimum rather than maximum of the reciprocal rates is essential:
all active inequalities must hold simultaneously. This is a conjunction,
not a sum of independent possible expansions.

The pair (L,D) distinguishes the previous reflection hostiles:

| V | delta | L | D |
|---|---:|---:|---:|
| {1,2} | 1/3 | 0 | 3 |
| {1,2} | 1/4 | 1/4 | 3 |
| {1,2} | 2/5 | 0 | 0 |

The natural reflection has trace zero in all three cases. For AP13 at
delta=1/14, there are six isolated witnesses and

    (L,D)=(0,23324/6435).

Uniform dilation preserves this germ: L_(cV)(delta)=L_V(delta), by the
measure-preserving circle covering t->ct, so D also agrees. This improves
on unnormalized component counts, which repeat under dilation. It does
not give a positive lower bound independent of the velocity height.

## L2. An explicit epsilon and a quantitative measure gap

Suppose delta=a/q in lowest terms. Then

    epsilon=1/(4qM)

is always small enough for the exact identity in L1.

**Proof.** Endpoint order changes occur only when two affine endpoints
coincide, or at0 or1/2. Their threshold values have denominators at most2M:
the relevant denominators divide v+w or |v-w|, with the same-speed
exception delta=1/2. A distinct such level b/c differs from a/q by at
least1/(qc)>=1/(2qM)=2epsilon. Thus there is no lower critical level
between delta-epsilon and delta. The volume is affine there and continuous
at delta, giving the identity throughout this interval. Continuity follows
also because f_V has no nonempty interval on which it is constant:
its finitely many linear pieces have nonzero slopes.

If G_delta is nonempty, pick t_0 in it. Since f_V is M-Lipschitz on the
circle, the arc of radius epsilon/M about t_0 lies in G_(delta-epsilon).
Its length is2epsilon/M. Therefore

    G_delta empty       => L_V(delta-epsilon)=0,
    G_delta nonempty    => L_V(delta-epsilon)>=1/(2qM^2).

This is a height-dependent separation, not a uniform small-period claim.
Normalizing by gcd(V) first can improve the height used in the bound.

The rational contact mechanism is classical and is also recovered in the
concurrent [Session91 LRC note](../session91_packing_golden_carriers_20261029/geometry/LRC_REFLECTION_AND_GOLDEN_GAP.md).
The present use keeps all endpoint-order changes, not just maxima, to
certify an exact first germ and its measure gap. No external priority is
asserted for this elementary arrangement argument.

## L3. One rational Gibbs parameter now detects weak witnesses

For h=delta-epsilon define the danger depth and its probability generating
function by

    H_h(t)=sum_(v in V) 1_{||vt||<h},
    Z_(V,h)(z)=integral_(R/Z) z^H_h(t) dt.

Take

    z=1/(4qM^2).

Then

    G_delta(V) nonempty  iff  Z_(V,delta-epsilon)(z)>z.

More quantitatively, the two possibilities satisfy Z<=z and Z>=2z.

**Proof.** In the empty case, L2 gives zero good-time measure at h, so
H_h>=1 almost everywhere and z^H_h<=z. In the nonempty case, the zero-depth
mass is at least1/(2qM^2)=2z, and Z is at least that mass. QED.

THM-2143
[Gibbs zero-count certificate and tight-boundary limitation](../../../01-canon/theorems/THM-2143-artanh-riesz-gibbs-certificate-soundness-boundary.md)
already proves the strict-branch criterion with an existential rational
fugacity. L2 supplies an effective height-dependent threshold relaxation
and mass gap, making a single explicit fugacity sound and complete for
the original weak threshold. This does not refute its warning: the
threshold has changed and the necessary height information is retained.

For LRC(14), substitute q=14 after primitive normalization. A proof of the
resulting positive Gibbs inequality for every thirteen-speed row would
prove the conjecture. No such uniform inequality is proved here.
The formula uses the full correlated depth distribution; replacing it by
independent runner probabilities or a mean depth is not justified.

The quantitative gap permits a certified numerical decision if the absolute
error in Z is less than z/2, by comparing with3z/2. The required accuracy
depends on M; no polynomial bit-complexity or fixed-resolution result follows.

## L4. Directional germs recover owner sets by Mobius inversion

Relax the individual thresholds to delta-epsilon lambda_v, lambda_v>=0.
The same endpoint argument gives the directional coefficient

    D_V(lambda)=sum_components (
        min_(v in P(a)) lambda_v/v
        + min_(v in N(b)) lambda_v/v ).

For a subset S of the labelled speeds, choose lambda_v=v for v in S
and0 otherwise. Then D_V(lambda) is the number of endpoint sides whose
entire active owner set is contained in S.

If C(S) is that count and m(A) is the multiplicity of the owner set A,
then C(S)=sum_(A subset S) m(A). Ordinary Boolean-lattice Mobius inversion
recovers every m(A). This is an exact recovery of the side-owner multiset,
not of the pairing of sides into components or their physical locations.
The directed tests in the checker use an epsilon bounded by1/(4qM^2)
so distinct affine endpoints cannot change order under these rate choices.

This is a genuine link to the repository's relation/character vocabulary:
the response operator decides which owner sets survive a chosen perturbation.
It does not identify those owner hyperedges with graphic cycles, Fano lines,
or tournament cyclic triples.

## L5. Even the complete threshold-volume germ forgets phase

Let

    V={1,3,6},    W={2,3,6},    delta=1/3.

Their good sets are respectively

    {4/9,5/9},    {2/9,7/9}.

Each has one point with positive owner3 and negative owner6 and one with
the roles reversed. The first speed is inactive at both surviving points.
Identify the first coordinate (speed1 in V, speed2 in W). For sufficiently
small signed threshold relaxations h=(h_*,h_3,h_6),
both volume functions are exactly

    2 max(0, h_3/3+h_6/6).

Thus their whole local multivariate volume germs agree, not merely L and D.
The inactive speed stays slack in a neighborhood; the other potential
contacts remain excluded by a positive gap. The two active linear endpoint
inequalities give the displayed formula.

Yet appending speed7 at the same threshold kills both V points:
their new distances are1/9. It preserves both W points, whose new distances
are4/9. This is a phase-sensitive future test with a genuinely new speed
on both sides. No threshold-volume germ with the phases forgotten is
closed under this operation.

This is a small explicit continuation analogue of THM-2050
[period14 top germs and global loneliness](../../../01-canon/theorems/THM-2050-period14-top-germs-do-not-determine-global-loneliness.md).
That theorem already warns that complete local data need a global sidecar.
Our germ in L1 measures the whole good set, avoiding restriction to just one
old clock layer, but it still cannot discard the absolute phase when the
operation is addition of a new runner.

A second hostile has positive bulk: {1,3} and {1,6} have
(L,D)=(1/9,2/3) at1/3, while adding4 kills only the first.
These examples concern a fixed threshold. After adding a speed, we do not
silently replace it by the different LRC threshold associated with the new
number of runners.

## L6. A closed owner carrier on the isolated-witness stratum

If L_V(delta)=0, G_delta is a finite set. Retain for each t in it

    (t, P_t, N_t),
    P_t=max{v:frac(vt)=delta},
    N_t=max{v:frac(vt)=1-delta}.

To append a new speed w at the same threshold:

- Delete the record if ||wt||<delta.
- Keep it unchanged if ||wt||>delta.
- At equality update P_t to max(P_t,w) when frac(wt)=delta, and update N_t
  analogously when frac(wt)=1-delta. At delta=1/2 both updates apply.

This is closed under every finite sequence of such additions. The new good
set is a subset of the old finite set, and L1's uniform germ at each surviving
point uses exactly these two maxima. No new witness can appear upon adding
a constraint. Full directional perturbations instead require the complete
active owner sets, not just their maxima.

Reflection sends (t,P,N) to(1-t,N,P). Hence one may store a canonical
representative t<=1/2 with multiplicity2, except at1/2 where it is1.
For AP13 the three representatives are

    (1/14,1,13), (3/14,5,9), (5/14,3,11),

each of multiplicity2. Their weights sum to23324/6435.
This is an exact symmetry compression; it does not force the carrier to
be nonempty. Further minimization for a restricted future-speed menu is
not claimed.

The carrier does not close deletion of a speed, movement of the threshold,
or a change to a lower runner-dependent threshold: these can recreate bulk
intervals. L5 proves why eliminating phase is invalid even in this narrow
closed stratum.
