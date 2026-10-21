# External interaction owners for the 94 residual intervals

2026-10-20. **PROVED:** guarded native inverse formula, exact ledger
transfer and scale-ray credit persistence. **FINITE-EXACT:** at N=4096,
94 distinct old/new interaction addresses pay all94 whole intervals left
by the frozen [flow bank](../annulus_residual_flow_20261020k/README.md).
This enlarges the ownership certificate; it does not reduce physical
energy. The uniform153N estimate remains **OPEN**.

## Native source and clock guards

Retain the inherited627-point grounded frontier V and first-arrival
canonical primitive map P, not literal odd Collatz U. If a is1 or2 and
j>=0, put

    s = [2^a(3u+1)4^j-(2^a+3)]/9.                 (1)

Require u primitive, integral positive s, s primitive and its first
branch a=2 when s=1mod8, otherwise a=1. Then P(s)=u with edge clock
a+2j. Indeed (3s+1)/2^a=B^j(u), B(z)=4z+1. Each nonfinal
B iterate is5mod8 and strips exactly once; the primitive u stops the
stripping. This is the inherited inverse-P mechanism, with all guards
retained, not a new convergence assertion.

The credit compiler additionally requires u in V, s outside V,
u<4N<=s<256N and

    beta(A_(64N)(u),(-1)^(a+2j) A_(64N)(s))=-3.   (2)

Here Q(x,y)=x²-xy+y² and beta is its polarization. The old source u has
depth/clock0; s has depth1 and the displayed clock. Thus this actual
old/new pair contributes exactly -3*1_[1,infinity)(d). Neither source
norm is spent. Conditions(1)-(2), especially the band/charge test, are
essential; inverse integrality alone supplies no negative interaction.

## Finite completion and exact transfer

The bounded scan stops at new source144369 after finding94 owners. For
example (u,s)=(87,59619), a=1,j=5 has terminal87 and clock11. All94
unordered pairs are distinct, old/new, and therefore outside every
all-new sibling-pair address. Different pairs may share a source: these
are owned interaction terms, not withdrawals of source norm. Every row
retains mod15, mod35, actual terminal, clock, depth and unsigned charge.

Every residual positive interval starts at depth>=1. Assign one distinct
credit to each interval. The permanent negative step contains its whole
interval. Write C_d=-282 for d>=1 and C_0=0, U_d for the old residual
positive debit, J_d for the full original sibling-bank debit, and R_d
for every other original interaction. Then

    R'_d=R_d-C_d,  J'_d=J_d+C_d,
    I_d=R_d+J_d=R'_d+J'_d,
    J'_d <= U_d+C_d <= 0.                         (3)

All257 physical horizons are checked, including the late cut115 and the
tail after the last sibling interval. Retaining the transferred slack
gives exactly the same bound R'+U+C=R+U. Dropping that nonpositive
slack gives the weaker bound I<=R'; nevertheless this finite control
still has max_d(R'_d-M_d)=-201392. No improvement has been obtained by
recycling R. Reusing these94 pair terms in a second external bank would
be double counting. Scalar payments may cross terminals; no vector is
moved between terminals.

## Infinite source-marked scale rays

One three-strip lift is s->64s+7(2^a+3), retaining parent and adding6
to the clock; its charge need not stay fixed. After THREE such lifts,
put F=64^(3k), N_k=FN and

    s_k=F s+(2^a+3)(F-1)/9,   j_k=j+9k.          (4)

The additive term lies in[0,F), so every integer annulus-band boundary
is preserved. Also F=1mod27, hence s_k=smod3. The six-band native charge
law therefore retains the new charge, while the inherited64-scale law
retains each fixed old charge. The clock increases by18k, preserving its
sign. Each of the94 unique pairs (u,s_k) still supplies -3 from depth1.
Indeed max V=63665 and the least selected s is33265: for k>=1,
s_k>=262144*33265>max V; k=0 was separately guarded. Thus no lifted
source is already a frontier terminal. Same-branch injectivity and the
disjoint first-branch residue classes preserve the94 distinct owners.
This is an infinite credit family; no assertion says it covers the new
scale's growing residual bank or controls its other interactions.

The modular reader computes the numerator of(1) modulo9m before exact
division by9. It retains the actual joint mod15/mod35 registers. Controls
use k=0,1,2,7 and moduli15,35,105, without materializing any long orbit.

## Validation

Run `python -B check.py --verify` and `python -B -O check.py --verify`
from this directory. The source receipt hash is authenticated; every
selected route and charge is replayed, with an independent charge
formula. Type, source-branch, premature-credit and duplicate-address
controls are retained. The discovery scan is limited to primitive new
sources below256N and stops once94 owners are found; production never
assumes an unresolved route terminates.

Generated receipt:4,356 checks; verification adds one equality gate.
LF SHA256: `d702d85fd51e0f9778ac1fb9a784a1426b3af76ea4aa0cb1d5a32392711eef82`.

Independent peer proof/API and optimized saved-receipt audit: **PASS**.
The phase-lift first-arrival boundary and scalar-transfer direction were
checked explicitly; no physical-energy improvement is inferred.
