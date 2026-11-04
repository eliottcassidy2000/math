# Direct two-query homogeneous-pair separation

**PROVED: P1-P3, at every finite order.** The exclusion itself is inherited
from **C2**, [two_port_eight_core_20261103/algebra/PROOF.md](../two_port_eight_core_20261103/algebra/PROOF.md).
We reproduce its swap identities and replace the cross-cover measurement
by two actual port assignments. No literature priority is claimed.

## Definitions and inherited mechanism

All old vertices are named. F_T(a,c) counts directed Hamilton paths after
adding x,y with x->y, with a_v=1 iff v->x and c_v=1 iff v->y.
The two masks are independent, including at the same old name.
S_v and E_v count Hamilton paths starting and ending at v.
C_uv counts ordered spanning covers by two **nonempty** directed paths,
the first ending at u and the second starting at v; C_vv=0.
Put R_v=E_v+sum_w C_vw.

A homogeneous pair is a two-vertex module: both names have identical
incidences to every exterior vertex. For a->b let Q contract the pair to z.
Contraction keeps every other name fixed.

The closest mechanism is incoming **C2** above, together with
**THM-4145**, [rooted-homogeneous-pair-expansion-two-defect-formula](../../../01-canon/theorems/THM-4145-rooted-homogeneous-pair-expansion-two-defect-formula.md),
and the linear response tomography in [Session27, section7](../session27_repunits_lattices_20261022/session27-response/proof.md).
The least-used sidecar here is the *row sum* of C, including its endpoint
completion E, rather than one cross-cover contrast. It has a direct
positive partition of the quotient Hamilton paths.

## P1. An exact calibrated query

For every tournament T and named vertex v,

    F_T(e_v,0)=H(T)+2R_v(T).                         (1)

**Proof.** The port y beats every old vertex and has only x as an
in-neighbor. The port x beats every old vertex except v. A Hamilton word
has exactly one of these three forms:

* x,y,W, with W a Hamilton word of T;
* y,P,x,R, with P ending at v and R possibly empty;
* P,x,y,R, with P ending at v and R possibly empty.

The latter two families each count E_v+sum_w C_vw: R empty contributes
E_v, otherwise its first vertex is the named w. All required port edges
are automatic and the decomposition is unique. Thus (1). QED.

This is an equality of integer query values, not a parity condition.
In particular c=0 does not make y a universal source because x->y.

## P2. Strict named order for a homogeneous pair

For every homogeneous pair a->b,

    S_a(T)-S_b(T)=S_z(Q),
    E_b(T)-E_a(T)=E_z(Q),                           (2)
    C_bw(T)-C_aw(T)=C_zw(Q),
    C_wa(T)-C_wb(T)=C_wz(Q)       (w exterior),     (3)
    C_ab(T)=C_ba(T)=H(Q).                          (4)

**Proof.** These are the identities in incoming C2, with (4) also in
THM-4145. Swapping a,b cancels words or covers where they are not
consecutive. This includes covers placing them in different paths.
Uniform exterior incidences preserve every affected arc. The unmatched
objects have consecutive block a,b, which contracts to z. For example
in C_bw-C_aw the unmatched block ends the first path, leaving C_zw(Q).
The endpoint identities and the other cross-cover identity are the
start/end versions of this argument. For (4), concatenate the two paths
at their designated a,b or b,a boundary and contract that boundary.
Both are bijections to Hamilton words of Q; no actual arc at the cut
is required. QED.

Consequently (including the canceling internal entries in (4)),

    R_b(T)-R_a(T)=E_z(Q)+sum_w C_zw(Q).             (5)

Every Hamilton word of Q either ends at z or has one particular
successor w with z->w. Cutting just after z gives

    H(Q)=E_z(Q)+sum_(z->w) C_zw(Q).                (6)

All C entries are nonnegative, so (5) is at least H(Q)>0. The final
positivity needs only existence of a tournament Hamilton path:
insert a new vertex before the first vertex it beats in an existing
Hamilton word, or append it if there is none.

Combining (1) and (5) proves the direct certificate

    F_T(e_b,0)-F_T(e_a,0)
      =2(E_z(Q)+sum_w C_zw(Q)) >= 2H(Q)>0.         (7)

The bound is attained, for example, for the order-two pair, where Q is
a singleton and the difference is 2. More generally equality in the
bound is equivalent to sum_(w->z) C_zw(Q)=0. No strongness hypothesis
is used. If the same named pair is homogeneous in U with b->a, the
same difference for U is strictly negative. This remains true if its
contracted quotient has completely different arcs.

Thus the two query vectors cannot agree or be positive scalar multiples.
This gives the same projective exclusion as incoming C2 through a direct
two-query measurement.

## P3. Ambient module converse and the remaining scope

Let M be a module in an arbitrary tournament T. If a,b form a homogeneous
pair in M, they also do in T: exterior vertices see all of M uniformly.
Reversing every internal arc of M reverses a->b and preserves pair
homogeneity, including its incidences to the other vertices of M.
Apply (7) in the two *actual ambient* tournaments. Their full named
two-port responses differ. There is no restriction on ambient order,
quotient strongness, or other arcs, and no assumption that exact core
response equality descends through an unknown scalar.

This rules out the entire homogeneous-pair branch. Combined with
[T5 and the finite structural reduction](../two_port_twelve_principles_20261103/STRUCTURE.md),
it leaves a strong, nonselfconverse, **prime** twelve-vertex core.
This agrees with incoming **C4**, which already establishes that conclusion.
A prime core remains **OPEN**; nothing here proves universal reconstruction.

For the whole-core comparison Delta=F_T-F_(T^op), **T4** in
[the full cubic proof](../two_port_twelve_principles_20261103/PROOF.md)
shows that Delta has degree at most three, with all named mixed cubic
coefficients retained, including repeated old labels across the ports.
Equation (7) is a separating linear-coordinate measurement of that full
function, not a replacement of it by low-order sidecars. A partial module
reversal inside an asymmetric ambient tournament need not have a cubic
contrast; (7) remains valid directly without such a degree claim.

## Hostile controls and reproduction

A transitive four-vertex tournament with outgoing rows (10,8,11,0) has
pair 0->1 in a middle ordered component. Reverse only that pair, giving
(8,9,11,0). Both have H=1, S=(0,0,1,0), E=(0,0,0,1).
Thus endpoint equality alone does not detect the reversal. The two actual
query vectors at names0,1 are respectively (5,9) and (9,5).

Run `python check.py --output check.json`, and the same command with
`python -O`. The exhaustive universe is all 1,099 labelled quotients
of orders1..5, each distinguished root: 5,405 expansions. It checks
42,600 cross-cover identities, 2,618 actual port counts, 293 independent
literal permutation-and-cut sidecars, 58 literal extended Hamilton counts,
and 121 prescribed nested ambient-module comparisons. No filter is applied
to the rooted quotient universe. All gates raise errors under -O.

The quotient-side DP and literal word/cut enumeration are distinct
implementations. They are finite controls for the bijective proof,
not a substitute for its arbitrary-order quantifiers. This session does
not claim a separate-agent proof audit.

