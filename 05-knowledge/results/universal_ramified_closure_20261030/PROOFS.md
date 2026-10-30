# The universal ramified symmetry group and its selector geometry

2026-10-30. **PROVED scoped, same-agent review; FINITE-EXACT checks
separate. No literature-priority claim.** This is the group common to
every bank, not the full group of each individual bank.

## U1. Source, clocks, and the exact universal group

Fix p>=5 prime, k>=1, and D with v_p(D)=1. Work in

    R=(Z/p^k)[delta]/(delta^2-D), q=p^k, N=q^2,
    rho(x+y delta)=x mod p,
    g=(1+delta)/(1-delta), C=<g>,
    A={z -> g^t z+b : t mod q, b in R}.

The denominator1-delta is a unit. The norm is Q(x+y delta)=x^2-Dy^2.
For nonzero t modulo q,

    v_delta(g^t-1)=1+2v_p(t).                         (U1.1)

Indeed g-1=2delta/(1-delta) has valuation1. For an element w of positive
delta-valuation, the linear term p*w in (1+w)^p-1 has strictly smaller
valuation than every higher term when p>=5. Raising to p therefore
adds2 to this valuation until truncation. A power coprime to p does
not change it. Consequently g has order q, has norm1, and is1 modulo
delta. Conversely, for each y mod q, x^2=1+Dy^2 has exactly one solution
x=1 mod p: lift one digit at a time using the unit derivative2x.
Thus C is the ENTIRE norm-one unit subgroup with residue1.

This includes the previous Fibonacci source: for delta=2mu-a and
D=a^2+4 simply ramified, mu^4 generates this same subgroup C. The
new theorem also applies to p=7,11,... via the displayed Cayley clock.
The inherited200 sqrt(pN) construction is NOT asserted for these
additional sources without checking its separate balance hypotheses.

Color every ordered pair(u,v) by the literal orbit(v-u)C, including
the diagonal color. Since |C| is odd, no nonzero orbit contains its
negative. Choosing one orbit in each opposite pair gives an actual
tournament; let B be this complete independent-sign bank.

Set eta=g^(p^(k-1))-1. Its valuation is2k-1, eta^2=p eta=0, and
eta*z=eta*rho(z). For h:F_p -> F_p put

    sigma_h(z)=z+eta*h(rho(z)),
    F={sigma_h : h:F_p -> F_p}, H=F A.

The earlier [S3b](../ramified_symmetry_spectral_20261030/PROOFS.md)
proves F=C_p^p preserves every ordered orbit: equal-residue differences
are unchanged; a unit difference v+eta*d equals g^(s p^(k-1))*v
for s=d/rho(v). A normalizes F, and A intersect F consists of the
linear-affine functions h(z)=u+vz. Hence |H|=p^(3k+p-2).

**Universal closure theorem.** At every k>=1,

    intersection_(T in B) Aut(T)
      = Aut(the fully colored difference relation)
      = A^(2) = H.                                   (U1.2)

Here A^(2), the standard two-closure, means permutations preserving
every orbit of A on ordered pairs. There is no further universal
symmetry hidden above the fibre shifts. Individual tournaments may
have larger groups because they merge colors.

### Proof of the upper bound

Preservation of every bank is equivalent to preservation of every
ordered orbit color: distinct orbit pairs can have independently chosen
signs, and opposite orbits always have opposite signs. The colors are
precisely A's ordered-pair orbits. It remains to show their automorphism
group is contained in H.

Translate an automorphism f so f(0)=0. Its value at1 lies in C, so
multiply by an inverse clock power to arrange f(1)=1 as well.
Every clock orbit preserves Q, so Q(f(u)-f(v))=Q(u-v). Polarizing,
using that2 is a unit, gives preservation of the bilinear form

    B((x,y),(x',y'))=xx'-Dyy'.

Pairing with1 shows that the first coordinate of f(x,y) is exactly x.
Write f(delta)=(0,c). Norm preservation gives c^2=1 mod p^(k-1),
and its directed orbit gives c=1 mod p. For k>=2 the two roots of
unity over the odd prime-power ring imply c=1 mod p^(k-1).
Pairing a general f(x,y) with f(delta) therefore fixes its second
coordinate modulo p^(k-1). Thus f(z)-z lies in the socle line eta F_p.
At k=1 this same conclusion follows already from the first coordinate.

Write f(z)=z+eta*h_z. If v=z-w has delta-valuation1, a nonzero change
eta*(h_z-h_w) has valuation2k-1. But every nonzero clock displacement
(g^t-1)v has EVEN valuation2+2v_p(t), and cannot equal it. Hence h_z=h_w.
Any two points of the same residue fibre connect by at most two such
valuation-one differences: choose a third point whose delta coefficient
mod p differs from both of theirs. There is such a value since p>=5.
Consequently h_z depends only on rho(z), proving f is in F. Undo the
normalizations to obtain H. This proves(U1.2). QED.

The quadratic norm forces all freedom into its bilinear radical; the
ODD/EVEN valuation boundary of directed orbit displacements couples
that radical across each residue fibre. Norm preservation alone does
not supply this coupling. The square-zero layer is not discarded.

## U2. Every finite named-selector stabilizer

Let S be a nonempty finite set of named vertices, and let j be the
number of distinct residues rho(S). Pointwise stabilizers in H obey:

    j>=2: |H_(S)|=p^(p-j).                           (U2.1)

If j=1, choose s0 in S and put
tau=min_(s in S, s!=s0) v_delta(s-s0), with tau=2k for a singleton.
Then

    |H_(S)|=p^(p-2+floor((tau+1)/2)).                 (U2.2)

The floor equals k for tau=2k. These formulas include repeated named
vertices by replacing a tuple by its support. The empty stabilizer
constraint leaves the entire group H.

**Proof.** Translate s0 to0. A root-fixing map can be written as
g^t*z+eta*h(rho(z)), with h(0)=0; this representation has multiplicity p,
because a top clock power is itself a linear fibre shift.
If S contains a unit v, fixing it forces (g^t-1)v into eta F_p, hence
t is divisible by p^(k-1). The whole map is then a fibre shift and
fixing S kills exactly j independent function values. This proves(U2.1).

If all pins have residue0, the h term vanishes at every pin. The clock
must satisfy g^t*v=v for every difference v, equivalently
1+2v_p(t)+tau>=2k. There are p^floor((tau+1)/2) clock choices and
p^(p-1) functions with h(0)=0; divide by the p-fold duplication to
obtain(U2.2). QED.

The groups themselves have types C_p^(p-j) in(U2.1), and
C_(p^ell) times C_p^(p-2) in(U2.2), ell=floor((tau+1)/2).
For the latter, the clock commutes with fibre shifts and a complement
to their shared top clock is given by h(0)=h(1)=0. These data determine
the stabilizer type, NOT the full configuration orbit. For example,
the ordered pairs(0,1) and(0,2) have the same occupancy and stabilizer
order, but norms1 and4 distinguish their directed difference colors.

Here tau is literally the lowest-common-ancestor depth in the nested
cosets modulo delta^t. No arbitrary rooted tree is substituted for the
ring. Once two residue fibres are met, only the occupied-fibre count
remains. Within one fibre, the common ancestor depth is indispensable.

**Sharp base size.** Exactly p named vertices are necessary and
sufficient to remove every universal automorphism. Any unoccupied
residue fibre supports a nonzero shift fixing all pins. A transversal
with one vertex in each fibre has trivial stabilizer by(U2.1).
This is a point-selector statement, not a lower bound on binary ports.

## U3. Free auxiliary points, but not free direct signatures

Take the minimal base B0={0,1,...,p-1}. Knowing images of these p named
vertices under a color-preserving map determines the images of
B1={delta,1+delta,...,p-1+delta} uniquely using their pair colors to B0.
Then pair colors to B0 union B1 determine every vertex uniquely.

**Proof in the fixed-base frame.** Pair norms to0 and1 fix the real
coordinate of j+delta to j. Its color relative to j is delta C.
If g^t*delta has zero real coordinate, its delta coefficient in g^t-1
must be divisible by p^(k-1), forcing t divisible by p^(k-1); that
power fixes delta. Thus j+delta itself is forced. Alternatively use
the norm root c=1 mod p^(k-1) and the valuation-one obstruction in U1.
The same argument works after undoing a known universal transport.

For a general vertex, pair norms to0 and1 fix x, and pair norms to
delta fix y modulo p^(k-1). Only an eta displacement can remain.
Let j=x mod p. If y is nonzero mod p, the difference to j has valuation1;
otherwise the difference to j+delta does. Its directed color excludes
any nonzero eta displacement. Hence the vertex is determined. QED.

The first-stage pair signature need NOT itself be injective. At p=5,
k=2,D=5, the five points(5,0),(5,5),(5,10),(5,15),(5,20) have identical
colors to all five integer landmarks. Yet their pointwise stabilizer
is trivial. The auxiliary delta landmarks resolve this collision.
Thus a group-theoretic base is not automatically a one-round resolving
set; closure under the actual relation supplies the extra information.
The second stage needs no additional chosen selector.

## U4. Exact joint-port enumeration and minimal recurrences

Let O_r be the number of H-orbits of ordered r-tuples of vertices,
allowing repetition. For r>=1,

    |H| O_r = p^(p-2) sum_(m=2..k) L_m p^((2k-2m+1)r)
               +(N/p)^r sum_(z=0..p) binom(p,z)(p-1)^(p-z) z^r,
    L_m=(p-1)p^(3m-2).                               (U4.1)

For k=1 the first sum is empty. Burnside on ordered tuples uses
Fix(f)^r, not the total number of cycles. In the inherited affine
inventory, all fixed-point types with m>=2 survive with multiplicity
p^(p-2). The identity and the m=1 fixed type are precisely the fixed
elements of A intersect F, and are replaced by F: a function with z
zeros fixes zN/p vertices and has multiplicity binom(p,z)(p-1)^(p-z).
All other elements have no fixed points. This also proves(U4.1) at
k=1, where the FULL cycle formula from the previous packet fails but
its fixed-point inventory still survives.

For fixed p,k, O_r for r>=1 has minimal rational recurrence order

    p+k-1, with characteristic polynomial
    product_(z=1..p)(X-zN/p)
       * product_(m=2..k)(X-p^(2k-2m+1)).             (U4.2)

All displayed roots are distinct positive integers and all their
coefficients in(U4.1) are positive, so a Vandermonde argument proves
minimality. The arbitrary r=0 convention O_0=1 is separate. Ordered
tuples of DISTINCT pins are counted by replacing each Fix(f)^r by
its falling factorial (Fix(f))_r in the same fixed-point inventory.

An independent derivation partitions all tuples by U2. For j>=2
residue fibres their labelled mass is

    (p)_j S(r,j) (N/p)^r,

where S(r,j) is a Stirling number of the second kind. Multiply by
p^(p-j) and divide by |H|. For tuples in one fibre of exact common
depth tau=1,...,2k-1 the labelled mass is

    N[p^((2k-tau)(r-1))-p^((2k-tau-1)(r-1))].

Multiply by(U2.2) and divide by |H|. The N fully diagonal tuples
make one final orbit. These strata reproduce(U4.1); at r=2 the result
is2p^k-1, exactly the number of directed difference colors.

For named binary ATTACHMENTS, use the complete Q-color cycle sum in
the preceding packet S3b, Q=2^r; it is a different action from r point
selectors. U1 now certifies that formula for the FULL universal group,
without changing it into the full group of an arbitrary individual bank.
For the pinned625 source, whose individual group is already H, both
counting laws are full graph-group laws. Its prime-child count from
the previous packet remains valid, with no new A000568 term claimed.

## U5. A logarithmic tournament packet realizes the universal action

There are M=p^k-1 pairs of opposite nonzero clock orbits. This follows
from(U1.1): each valuation stratum gives its number of elements divided
by its clock period, totaling2(p^k-1) directed orbits.
Assign each opposite pair a different binary word of length
r=1+ceil(log2 M), choosing one word from each antipodal pair of words.
Use the r bits as orientations in r separate named tournaments.
Opposite differences receive complementary words. Every ordered color
is recoverable from the joint word, so the common automorphism group
of this ordered tournament packet is EXACTLY H.

One chosen bank may be prescribed as the first layer; use its direction
as the first bit and encode the M pair names in the remaining bits.
The other layers are actual tournaments but are not necessarily prime.
The bit bound is minimal for INJECTIVELY CODING all2M directed colors,
because r bits have at most2^r words. It is not a minimal number of
banks needed merely to realize H as an automorphism group: one pinned
625 bank already does that.

## U6. The universal kernel and the bank action compose exactly

Let Aff(R^+) be the affine group of the additive Z/p^k-module, and
let M be the normalizer of A inside it. Then

    H intersect Aff(R^+)=A,
    H is normal in H M,
    1 -> H -> H M -> M/A -> 1.                       (U6.1)

The quotient M/A is exactly M's faithful action on the bank B.
For every bank state T,

    |(H M)_T|=|H| |(M/A)_T|,
    |(H M)_T|/|M_T|=p^(p-2).                        (U6.2)

Thus adjoining the nonlinear universal group changes stabilizer
weights but does NOT change bank-state orbits. No splitting of the
extension and no equality(H M)_T=Aut(T) is assumed.

**Proof.** If a fibre shift is affine, h(z)-h(0) is additive on F_p,
hence linear; these are exactly the shifts already in A. Removing
the affine factor from any element of H proves the intersection.
Normalizing A also normalizes its ordered-pair orbits and therefore
their full automorphism group H=A^(2). The kernel of the M-action on
every bank is M intersect H=A, by U1. The isomorphism and stabilizer
formulas follow. The same statements hold for any subgroup of M
containing A. QED.

This composes with the incoming [Session111 bank recurrence](../session111_scale_connections_20261030/BANK_RECURRENCE.md),
whose normalizer acts on orientation coordinates rather than on vertex
selectors. It neither replaces that recurrence nor counts its group
again. At p=5,k=2, its odd effective subgroup E=C5^2 acts on four fixed
bits and two independent ten-bit blocks, each comprising two parallel
five-cycles. This actual representation is inherited from
[ACTION_CARRIERS.md](../session111_scale_connections_20261030/ACTION_CARRIERS.md).

For binary bank states, a block is fixed by its nonidentity coordinate
shift in4 of its1024 assignments. Consequently the exact strata are:

| stratum | effective stabilizer order | labelled bank states | E-orbits | known lifted vertex stabilizer |
|---|---:|---:|---:|---:|
|neither block fixed|1|16646400|665856|5^9|
|one block fixed|5|130560|26112|5^10|
|both blocks fixed|25|256|256|5^11|

The orbit sum692224 is inherited; the new interpretation retains its
universal vertex kernel. All2^24 banks are included, without a prime
or balance filter. For Q unsigned colors the labelled stratum with
j fixed blocks is Q^4 binom(2,j)(Q^2)^j(Q^10-Q^2)^(2-j), and its
orbit count is that integer divided by5^(2-j). Q=2^r can be interpreted
as r jointly named bank assignments in the inherited fixed-origin gauge.
These strata give certified subgroup stabilizers, not an all-bank
classification of full tournament automorphism groups.

## U7. Full unsigned closure: a two-bit phase and its sharp depth boundary

Now color a difference only by its opposite pair v*(C union -C),
and let W preserve every such unsigned color. Let c be ring conjugation
delta -> -delta. At EVERY k>=2,

    W=H semidirect <-1,c>,  W/H=C2^2, |W|=4|H|.     (U7.1)

This extends the incoming [SIGN_QUOTIENT](../session111_scale_connections_20261030/SIGN_QUOTIENT.md)
from the additive normalizer to ALL permutations preserving these
unsigned colors. Its two orientation probes are inherited, not newly
claimed for that narrower normalizer. W is an even-order symmetry
group of an unsigned colored relation, not of an individual tournament.

**Proof.** An unsigned-color automorphism preserves Q. Translate0 and
use a multiplier in C union -C to fix1. Polarization fixes the first
coordinate as in U1. The image of delta has second coordinate congruent
to either1 or-1 modulo p^(k-1); use c if necessary to choose1. The
remaining map has form z -> z+eta*h_z. For a difference v of valuation1,
membership in a NEGATIVE clock orbit would change it by
-(g^t+1)v, of valuation1. This cannot equal a nonzero socle displacement
of valuation2k-1>=3. Only the positive clock orbit remains, so U1's
odd/even argument forces h to depend solely on the residue. This proves
the upper bound. Both -1 and c preserve unsigned colors and normalize H;
their four-element group intersects the odd group H trivially. QED.

Conjugation sends delta^d*u into(-1)^d times its directed clock orbit:
the quotient of the conjugated unit by u has norm1 and residue1. Thus
-1 reverses every directed color, and c reverses just the odd-valuation
colors. The directions of1 and delta distinguish the four cosets in
W/H. Exactly two bits recover this GLOBAL orientation phase after
factoring out H; they do not identify an arbitrary vertex permutation.

At depth k=1 the conclusion is false. Unit colors record only the
nonzero residue difference, while within each residue fibre the colors
are the signed difference pairs of a p-cycle. Therefore

    W=(D_(2p)) wr D_(2p),  |W|=(2p)^(p+1),          (U7.2)

in its imprimitive action on p fibres of size p; D_(2p) has order2p.
Each fibre can reflect independently. H has index2^(p+1) here and is
not in general normal in W; a quotient C2^(p+1) is NOT asserted.
The p-cycle color already forces each component's dihedral group, so
the displayed wreath group is an upper as well as a lower bound.

**Free phase by selector placement.** For k>=2, W still has minimal
point-base size p. Use the transversal

    {0,1,2+delta,3,...,p-1}.

After fixing0 and1 a W element is a fibre shift followed possibly by c.
The pin2+delta excludes c, because its displacement2delta cannot be
cancelled in the top socle when k>=2. The occupied fibres then kill F.
The earlier unoccupied-fibre lower bound proves minimality. Thus tilting
one of the existing p pins removes the global phase without another pin.
The integer transversal instead leaves c. At k=1, each fibre's dihedral
action requires two pins, so the sharp unsigned base size is2p.

The finite checks exhaust W's two-point stabilizer at25 and625. Their
orders are4000 and250, giving full unsigned orders1000000 and7812500.
At625,3750 of4000 initially possible local reflection choices are
rejected by joint landmark colors: only the two global phases remain.
This is a literal alignment constraint, not independent sign averaging.

## U8. Generality, literature, and stopping boundaries

The standard ambient language is a cyclotomic association scheme,
Cyc(C,R), whose color-preserving group is the two-closure of A.
[Evdokimov--Ponomarenko](https://arxiv.org/abs/math/0608469), Sections1
and3.2, already study these groups and explain nonlinear symmetries
from an impure unit subgroup containing1+I for a nonzero ideal I.
Here I=eta F_p. Their necessary normality obstruction applies;
their Galois-ring classification is not applied to this ramified ring.
U1 supplies its own exact upper bound for this quadratic source.
This is an application-specific proof and enumeration, with no claim
that two-closure or the impurity mechanism was invented here.

The restriction p>=5 in the displayed Cayley-clock construction is
material. At p=3,k=2,D=6, g=(4,5) mod9 has order3, not9. The norm-one
kernel is not generated by that clock; U1's parameter assumptions fail.
Also Q alone cannot replace directed orbit colors, and the entire bank
cannot be replaced by one arbitrarily selected tournament.

The81/27 control remains outside this ramified source but is replayed
as a hostile to recovering full groups from coarse triple observations.
The incoming Session111 joint-code and normalizer-bank actions are
distinct: binary port correlations and the action on bank coordinates
are not the vertex point-selector action counted in U4.

Actual complements can leave these cyclotomic schemes. Their full
target relation and selected path remain necessary. Minimal general
complement/re-marking carriers, classification of individual ramified
bank groups at every depth, and unrestricted A000568 remain OPEN.
