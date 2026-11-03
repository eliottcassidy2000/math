# Incoming source returns and the polarity boundary

Read after clean rebases throughd9d9f3e483 and25c9e96557. These connections retain
their source attribution and do not alter O3's all-order H35 obstruction.

## I1. The incoming zero return extends both new ladders indefinitely

[Session163 FIRST_STEP](../session163_ternary_triples_20261103/strong_returns/FIRST_STEP.md)
and [MOVING_GUARD](../session163_ternary_triples_20261103/strong_returns/MOVING_GUARD.md)
give an exact two-ear return. For a guard q with unique in-neighbour r,
adding x,y with q->x->y->q and x,y beating every other old vertex gives

    H_new=3H+2*#paths of T-q starting at r.

The new guard x has unique in-neighbour q, and deleting x leaves y a
source. Its next correction is therefore zero. The construction repeats,
preserves strongness, and triples the fixed-port observer in the query
polarity specified there. This is the same literal pair of births as
our reset, specialized to a unique-in-neighbour input; it is credited
to the incoming proof rather than rediscovered here.

Apply it **after our first reset and before the final converse**. That
first reset always supplies exactly the required new guard x, even
when the original selected root had two in-neighbours. Keep our spare
old port separately. For every k>=1 and j>=0 this yields actual strong
tournaments with counts

    H_minus(k,j)=3^j(9*2^k-1),
    H_plus (k,j)=3^j(9*2^k+1).                         (I1)

The source program is exact: after the first pair of commands, repeat
`(0,2,1),(0,2,2)` j times, then take global converse and select slot0.
Every repetition retains the spare and moves the guard. The same joint
polynomial operator applies at all sizes. This extends an entire
two-parameter family, rather than a finite list of examples.

The important boundary is polarity. Before final converse, the incoming
three-coordinate observer at the spare is multiplied by3 each time.
After converse, its dual observer generally is not. For example our
minus k1 return has v=(17,16,17), while one zero continuation has
v=(51,115,116), not3v. Converse exchanges the start and end queries;
the incoming commuting-query hypothesis does not transfer for free.
Our full joint program transports the correct dual observation and
F49 phase. This is a concrete reason to retain both polarized responses.

`return_program.py` checks the exact incoming3x observer before converse,
the full joint computation afterwards, independent subset counts on
small cases, and both k20,j20 strong sources of order65. No claim is
made that this zero-branch closure supplies nonzero-branch closure.

## I2. Finite observation and source realization remain separate

The incoming [MRB/square-sum finite-horizon theorem](../mrb_square_probe_spectrum_20261103/HORIZON.md)
recovers every k-birth scalar response by signed tensor interpolation
from cuts of size at most two, on the same actual named source. Its
local cause is degree two at each newborn's path incidences. Our capped
atlas uses a different positive quantity: bounded ordered joint response
forces a logarithmic bound on a terminal word. One is exact observation
at a chosen horizon; the other is source classification at a chosen
response budget. Combining them does not turn signed interpolation into
an attainable source, and it cannot create a ready source at H35.

The incoming [MRB-star family](../session162_mrb_square_20261103/mrb_star/NOTE.md)
also omits H35 at every size despite a literal strong H35 witness. Its
family obstruction is distinct from our rejection of every possible
1/2-ready presentation. We keep both scopes. The supplied witness is
not needed to establish ours, since RETURN independently constructs H35.

The [H497 packet](../rooted_carry_closure_20261103/README.md), initially a
finite checkpoint atd9d9f3e483, was promoted with its proof write-up in
fffa1865ba and extended in25c9e96557. It now supplies actual balanced
primitive returns, shared phase legalizers, and, after adjoining the
same moving-guard reset, the complete2401-element M2(F7) residue monoid.
These are proved within their stated languages. Its balanced-pair
phase is different from the golden Q4 observer used here; the packet
itself gives literal zero-pair phases with nonzero golden phases.

The most productive next combination is now specific: use the incoming
fixed-horizon probes to search for nonzero-offset operations on the
graded return states, then certify a returned **actual source**. The
polarity counterexample above is a cheap hostile for any proposed
compression that keeps only one of the two endpoint observers.

## I3. A seven-to-one obstruction to identifying the two field charts

The incoming [PC3 field pairing](../rooted_carry_closure_20261103/PHASE_49_50_51.md)
uses shared F7 matrices on two F49 coordinates. Every such linear
action commutes with componentwise Frobenius. Its42 nonrational
directions therefore form21 Frobenius pairs. The later
[complete monoid](../rooted_carry_closure_20261103/SINGULAR_RETURN.md)
retains this commutation property, even when matrices are singular.

Our plus controller's21-state action has a different exact law:
`sigma G sigma G^(-1)=tau`, where tau is a nonzero trace-zero
translation of order7. Suppose a proposed map from our phases into
the incoming paired chart intertwines G with one shared F7-matrix
action and sigma with componentwise Frobenius. Commutation in the
target forces the map to identify p with tau(p). Therefore it must
collapse each seven-element translation orbit, and has at most three
distinct images on our21 states. In particular it cannot be faithful.
This proof only uses equivariance; on the image the G action is
invertible because G^3 is the identity there.

Thus the incoming actual residue monoid does not by itself realize our
synthetic phase controller with both field operations identified.
The strongest survivor is a bridge retaining an additional source
coordinate or implementing a non-Frobenius-fixed affine translation.
In homogeneous field coordinates our G has an off-diagonal constant
`-1+2theta`, outside F7, which pinpoints the missing operation. This is
a concrete transport obstruction, not a disagreement between the two
correct appearances of21,42 and F49.
