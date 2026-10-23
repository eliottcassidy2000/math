# Paired-charge circuits, binary cuts, and the complete dyadic carry compiler

23 October 2026. **PROVED:** the general paired-charge criterion, rowwise
shared-coordinate circuit law, binary-cut/carry compiler and conditional
injective realization rule. **FINITE-EXACT:** the compressed original
order-eight labels, their 1,695 fibres, and the stated positive/hostile tests.
No unconditional all-order native supply theorem is asserted.

## 1. Inheritance and the genuinely different coordinate

The closest mechanism is the all-modulus row-image criterion in
`session44_guarded_native_moves_20261023/electrical-duality/UNIVERSAL_ELECTRICAL_GUARDS.md`.
The nearest obstruction is the fixed-source tariff in
`session44_guarded_native_moves_20261023/rectangle-bank/THEOREM.md`: at the
original sources 1238 and 2246, a full-bank character makes every net target
coefficient even, but the selected target ports have capacity one. The
canonical hostile is the weighted C4 scalar-phase/current distinction, not
merely the earlier one-character failure. The least-used sidecar is the
complete terminating dyadic filtration and its actual marked harmonic lifts.

The concept board is: paired-charge image / target-owned rows / shared
source-target supports / binary cut cycles / higher-bit carry / authentic
injective path resources. Anchor: multi-source cancellation. Niche: compress
the full original dyadic bank, not just one primitive character. Wildcard:
charge-balanced composite source divisors.

The new coordinate is the finite group **2K_(2)**, the image of doubling in
the original 2-primary critical group. It is not the odd-primary critical
group, a native path space, or a list of independent per-rectangle tariffs.

## 2. All-size paired-charge theorem

Let L be a connected finite graph's integral symmetric Laplacian, let
D^0 be its integral zero-sum divisors, K=D^0/LZ^V, and K_(2) its 2-primary
part. Write pi_2 for primary projection. Fix a root r and define

    Q_b = 2^b K_(2),
    ell_b(v) = 2^b pi_2([e_v-e_r]) in Q_b.                 (1)

Let F be any integral symmetric matrix with zero row sums; positive
conductances and native feasibility will be imposed separately. Then

    ker(L mod 2^k) subset ker((L+2^b F) mod 2^k)
    for every k >= 1

if and only if every physical row i satisfies

    sum_j F_ij ell_b(j) = 0 in Q_b.                       (2)

In particular, vertices which occur as both sources and targets are NOT
split into independent copies. Their contributions are added in the same
row before (2) is evaluated.

Proof. For each k the finite perfect dot-product pairing gives
(ker L)^perp=im L over Z/2^k. A zero-sum row d annihilates the original
harmonic kernel iff [d] lies in 2^k K. The intersection of 2^k K over all k
is precisely the odd-primary subgroup, since multiplication by two is an
automorphism there and eventually kills K_(2). Apply this to every row of
2^b F. Equation (1) rewrites its 2-primary class as (2). This proves both
directions, with no inference of equality of harmonic kernels. QED.

If the exponent of K_(2) divides 2^a, the condition at k=a already implies
all k. Odd moduli are NOT covered. In abstract invariant-factor coordinates,
if K_(2)=direct-sum Z/2^{a_i}, then

    Q_b = direct-sum Z/2^{max(a_i-b,0)}.                   (3)

An owner-marked map, not just these factor orders, is necessary to use it.

For authenticated rectangles let u_t=e_{S_t}-e_{R_t},
v_t=e_{U_t}-e_{V_t}, and

    F = sum_t a_t (u_t v_t^T + v_t u_t^T).                (4)

The exact debt at a physical vertex j is

    sum_t a_t (u_t[j] ell_b(v_t) + v_t[j] ell_b(u_t)).     (5)

This is an additive finite-group circulation law with a capacity sidecar.
There is no rule requiring each summand separately to vanish. If source
and target role sets are disjoint and M=sum_t a_t u_t v_t^T, (5) becomes
M ell_b(targets)=0 and M^T ell_b(sources)=0. With overlapping roles only
(5), not these separate equations, is valid.

## 3. Primitive composite charges evade a single-pair tariff

For any two disjointly supported zero-sum integer divisors u,v, write
F=u v^T+v u^T and g_u=gcd_j |u_j|, g_v=gcd_j |v_j|, with gcd(0)=0.
Equation (2) is equivalent to

    g_u ell_b(v)=0 and g_v ell_b(u)=0.                    (6)

Indeed source-support rows are u_j ell_b(v); their joint vanishing is the
gcd condition by integral Bezout. Target rows give the other condition.
For primitive u,v, the condition is exactly ell_b(u)=ell_b(v)=0.

A useful nontrivial circuit has ell_b(U)=ell_b(V) and four distinct source
vertices A,B,C,D outside those targets such that

    ell_b(A)+ell_b(B)=ell_b(C)+ell_b(D).                   (7)

Let u1=e_A-e_C, u2=e_B-e_D, v=e_U-e_V. The sum of the two amplitude-one
rectangles has F=(u1+u2)v^T+v(u1+u2)^T and is safe by (6). If all four
source labels are distinct, each individual rectangle is unsafe. Their
nonzero debts cancel at the SAME two target rows, not in an unowned scalar
sum. This can evade a single-source-pair coordinate tariff without even
amplitude. Native realization still requires four independently available
source transitions and certified nonrecycled literal resources.

A dual construction exchanges the source and target roles. More general
menus are governed by (5), not merely this sufficient rank-two pattern.

## 4. Exact original order-eight compression

The pinned original filtration in
`metagraph_order8_electrical_20261023/README.md` and its receipt proves

    K_(2) = (Z/2)^3358 + (Z/4)^4 + Z/256.

The image Q_1 therefore has order 2^11=2048 and type

    Q_1 = (Z/2)^4 + Z/128.                               (8)

Use Y0,...,Y4 from `radical_lifts_mod4`, Z from
`eighth_radical_lifts_mod8[0]`, and X from the last primitive higher lift,
all rooted at 3527. Literal multiplication by all 438,216 original edge
rows verifies LY_i=0 mod4, LZ=0 mod8, LX=0 mod256. The five parity vectors
Y_i mod2 are independent. All 3,528 coordinates verify

    X = Y0+Y4 mod2,
    X = Z+2Y0+2Y3 mod4.                                 (9)

The complete marked label can thus be taken as

    ell(j) = (Y0[j],Y1[j],Y2[j],Y3[j] mod2; X[j] mod128). (10)

This is an isomorphism-coordinate version of (1), not a claim that the
coordinates in the receipt are a Smith transformation.

Here is a proof of completeness, avoiding a selected-bank assumption.
The original grounded H256 has order 2^3374. The kernel of coordinatewise
reduction H256 -> (Z/128)^V consists of 128 times the original grounded
binary kernel, of order 2^3363. Its image therefore has order 2048.
The vectors 64Y0,...,64Y4 generate a subgroup of order 32 in that image.
X mod128 has order 128 because some coordinate is odd, and
64X=64(Y0+Y4). Its intersection with the even subgroup has order exactly
two: the even subgroup has exponent two, and 64X is nonzero. Thus these
vectors generate 32*128/2=2048 elements and are complete. Relation (9)
allows Y4 to be removed. Consequently, for every integer symmetric
zero-row-sum F,

    2F retains ALL original dyadic banks
      iff F Y_i=0 mod2 (i=0,1,2,3), and F X=0 mod128.     (11)

The all-dyadic conclusion uses the terminating exponent-256 computation,
not an extrapolation beyond the selected primitive character.

The same argument gives a precise tower. Modulus two is automatic; for
k>=2 preservation of the full original H_(2^k) is equivalent to

    F Y_i=0 mod2 (i=0,1,2,3),
    F X=0 mod 2^min(k-1,7).                              (12)

For 2<=k<=8, scale Y_i by 2^{k-2} and reduce X modulo 2^k. The reduced
harmonic image has order 2^{k+3}; those displayed generators attain it.
For k>8 use saturation in Section 2.

The original 3,528 vertices occupy exactly 1,695 labels. Fibre sizes are
1:642, 2:545, 3:310, 4:139, 5:44, 6:15. Full lists and the complete marked
coordinate table are saved in `node_labels.json`. The four finite-order
coordinates are additive residues, not opaque class names.

## 5. Binary cut circuits and higher carry, with every row retained

For any integer symmetric zero-row-sum F and binary potential y, let
C={v:y_v=1}. Then

    (F y)_j = sum_(v across the cut from j) F_jv mod2.    (13)

Thus F y=0 mod2 iff the mod-two graph of odd F_jv crossing that cut has
even degree at EVERY vertex. Equivalently its crossing edges form an
Eulerian subgraph, hence a union of cycles. A lone crossing edge, or any
odd-degree crossing row, gives an explicit local obstruction. The signs
are retained in F for higher layers even though parity forgets them.
Apply (13) to the four Y parity cuts and to X mod2. This yields a common
binary feasibility gate: a cancellation must share a physical row and a
cut debt; two defects at different target rows cannot cancel.

Write X mod128=sum_(t=0)^6 2^t x_t with x_t binary. Define c_0=0. Then
F X=0 mod128 iff every one of the following seven tests succeeds:

    F x_t+c_t is even at every row,
    c_(t+1)=(F x_t+c_t)/2.                               (14)

Inductively, after t successful tests,
sum_(s<t) 2^s F x_s = 2^t c_t. Adding the next term proves (14) and the
final divisibility. Negative entries are ordinary integers; no nonunit
is cancelled modulo a ring. Replacing a failed odd numerator by a floor
would be unsound. The carry vector is part of the certificate.

The compiler aggregates (4) first, retaining a single row for every
physical vertex, evaluates (10), and produces the four parity-cut ledgers
plus seven carry ledgers. A failure includes its row, character coordinate
and nonzero residue. A success can additionally be checked by literal
multiplication on the 2,048 generated reduced harmonic vectors. This is a
small marked-module/carry certificate, not one-hot SAT or a Graver restatement.

## 6. Exact capacity boundary and a sufficient native wiring rule

The electrical theorem alone does not select paths. An authenticated menu
must specify, for every actually used copy of a rectangle, its literal
before/after Hamiltonian paths at each source, actual complement owners,
full shared-coordinate record, kernel/automorphism lifts, converse records,
and paired-black resource orbit IDs. A coefficient two needs TWO supplied
copies; a large aggregate edge conductance does not create the second copy.

Within a fixed authenticated menu, the exact certificate has three parts:

1. The union of the supplied before records and the union of the after
   records each define a single value at every physical shared-coordinate
   key, and retain all prescribed pins and required shared equalities.
2. Actual source slots and literal resources are injective in each endpoint
   packet, with the declared unchanged complement explicitly reserving its
   used resources. If the menu requires all selected old/new witnesses to
   be mutually distinct, enforce that stronger union condition as well.
3. Its compiled F passes (5), or equivalently (11) for the original atlas,
   and the actual colored multiplicity difference is exactly 2F with the
   sign convention of (4). In particular no loop or converse factor is lost.

These conditions are necessary and sufficient for the stated simultaneous
endpoint replacement within that authenticated finite menu to retain all
requested dyadic original harmonic vectors and the stipulated resource and
shared-coordinate requirements. Necessity reads the actual endpoints.
For sufficiency select exactly the supplied paths at the named slots;
parts 1 and 2 make the union well-defined and injective, the unaffected
complement is retained, and part 3 proves the harmonic inclusion. This
asserts a simultaneous move, not a sequence of legal one-slot intermediates.
If the supplied endpoints are only selected observer packets, it certifies
only that packet, not a full native graph or global bijection.

For a charge-balanced circuit (7), a particularly transparent sufficient
resource hypothesis is four sources with literal U<->V alternatives, all
four source owners distinct, both endpoint records preserving the same
entire shared record, and all selected F/R resource orbits disjoint from
one another and the fixed complement. Then the signed source partition
chooses the four alternatives directly, and no port is recycled.

The usual capacity inequalities are necessary but not sufficient without
these explicit resources. Nor do separately compatible local shared records
imply compatibility of their union. These are the exact boundaries of the
all-size wiring rule.

## 7. Controls and non-consequences

The checker keeps the old fixed pair 1238/2246 and targets 881/2762: the
primitive X guard passes modulo eight but a Y parity row fails. It tests
the new algebraically safe 913/2701 by 398/528 rectangle, but does not
claim its paths are authenticated by this packet; the native producer owns
that evidence. It constructs a four-source relation from original marked
labels and checks that two individually unsafe rectangles cancel on their
shared targets. That latter certificate is explicitly algebraic unless a
native path menu is attached.

A source/target-overlap hostile uses u=v=(1,-1), modulo eight and x=(0,4):
the two row contributions cancel although the separate source-gap tariff
would forbid it. The compiler retains both terms. The weighted C4 with
weights01=03=1,12=23=2 has grounded H8 equal to 0000,0040,0424,0464.
Targets1 and3 have identical scalar evaluations, but the current of
character0424 from0 lands as 4400 versus4004. Thus scalar phase equality
never replaces the owned row-current criterion.

Kernel preservation means INCLUSION of the original harmonic bank in the
new bank. No equality of kernels, equality of critical groups, odd-primary
preservation, universal existence of native options, or unconditional
all-size tournament wiring follows.

## 8. Reproduction and source scope

Run from this packet with the repository supplied explicitly:

    python3 -B check.py --repo /path/to/math --out normal.json
    python3 -O -B check.py --repo /path/to/math --out optimized.json
    cmp normal.json optimized.json

`recheck_filtration.py` reruns the original pinned order-eight calculation
and compares its whole result with the original receipt. It deliberately
does not rerun the earlier order<=7 control file, which is absent from this
checkout. Direct invocation of the old checker's main program hit that
missing historical control input; no source was edited. The main checker
here independently multiplies the full original edge matrix by every used
harmonic lift and verifies the marked generator/compression relations.
Source pins and a final manifest keep these dependencies explicit.
