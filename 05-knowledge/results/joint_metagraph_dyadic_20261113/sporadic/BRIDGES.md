# Signed pair lifts and eta coefficients beyond fixed finite-matrix recurrences

**Status: PROVED IN NOTES, Propositions1–4 including3a, after independent audit.**
Literature inputs remain **CITED**; specified computations **FINITE-EXACT**.
See the root [audit](../AUDIT.md) and [promotion record](../PROMOTIONS.json).
No tournament spectrum theorem, new moonshine theorem, or group action
inferred from matching prime factors is claimed.

## Inheritance and the concept board

Closest mechanisms are [THM-024, involution anti-automorphisms](../../../../01-canon/theorems/THM-024-involution-anti-aut.md),
[THM-4857, tournament automata and QR codes](../../../../01-canon/theorems/THM-4857-the-gf2-automaton-of-a-tournament-classifies-the-merged-metagraph-exactly-through-order-five-circulants-are-additive-cellular-automata-and-paley-tournaments-give-quadratic-residue-codes.md),
and the [existing prime-level eta recurrence](../../weighted_crt_clocks_20261028/ETA_RECURRENCE.md).
The [current A5 frame packet](../../session314_metagraph_constraint_algebra_20261113/HISTORICAL_BRIDGES.md)
already connects K5, ten K3,3 frames, and two distinct Petersen constructions.
We retain its distinction between a canonical graph and a marked switching representative.
The [cyclic Frobenius kernel sidecar](../../session97_overnight_arithmetic_carriers_20261029/supersingular/CYCLIC_KERNEL_RECURRENCE.md)
already constructs the Hamming and Golay extensions and rejects their use as a characterization of Ogg primes.

Canonical hostile: a reflection on five letters has one fixed point, but
its induced action on ten unordered pairs has two. Corrected near miss:
prime support is not an action, and a Petersen minor does not transport
orientations or boundary responses. Least-used relevant sidecar: the
composite eta levels 4,9,25, with the previously omitted exponent parity.

| Live concept | Representation / predicate | Operation and missing coordinate | Decisive probe |
|---|---|---|---|
| Tournament recurrence | finite matrix output, rational generating function | increase block length; retain actual boundary state | exponential modular-cusp singularity |
| A5 / Petersen | five letters and their ten pairs | pair lift; retain stabilizer action on a pair | signed reflection obstruction |
| K3,3 | two triples, C3 x C3 action | complete all missing pairs; retain two internal triangles | eight orientations, each H=9 |
| 23 to 24 | Paley row code to parity extension | binary span then one parity coordinate | complete 4096-word enumerator |
| Literal set 4,2,5,3,7,13 | integer eta ratios | retain level and exponent parity | add 9 and25 as odd-exponent controls |
| Ogg and sporadic primes | supersingular Frobenius versus group-order support | quotient by Fricke / forget level | 11 versus4; 29 versus Baby Monster |

Anchor: exact recurrence interfaces. Niche: signed actions under the
vertex-to-pair functor. Wildcard: the literal unusual set as an eta
exponent-parity filter. These are different mathematical objects; the
connection contracts below say exactly what transfers.

## 1. Signed orbital criterion and the pair-lift obstruction

Let a finite group G act on a finite set V and let epsilon:G->{+1,-1}
be a homomorphism. A **signed invariant tournament** means that g preserves
every arc when epsilon(g)=+1 and reverses every arc when epsilon(g)=-1.
This is a global reversal sign, not an independently chosen edge sign.

**Proposition 1, PROVED IN NOTES (standard orbital construction, stated explicitly).**
Such a tournament exists iff, for every unordered pair P={u,v}, the
permutation sign of the action of its setwise stabilizer G_P on its two
points equals epsilon restricted to G_P. If the condition holds, the number
of labelled tournaments is exactly 2^r, where r is the number of G-orbits
on unordered pairs.

*Proof.* Choose a direction on one pair in each orbit. Transport it using g,
reversing once if epsilon(g)=-1. Two transporters agree precisely when their
quotient in G_P acts on the two endpoints with sign epsilon. This proves
well-definedness, necessity, and the two independent choices per orbit. QED.

For epsilon=+1 this is equivalent to the permutation image of G having odd
order. An even-order image has an involution, which swaps some pair. An
odd-order image has no element swapping a pair, because such an element
would have even permutation order. This is the familiar tournament
automorphism obstruction; it is inherited, not a new theorem.

**Proposition 2, PROVED IN NOTES (the additional connection here).** Let m>=3 be odd.
The dihedral group of order2m acts on Z/m by x->a+x and x->a-x. Give
rotations sign+ and reflections sign-. Then:

1. On the original m vertices there are exactly 2^((m-1)/2) signed invariant
   tournaments, the circulant tournaments.
2. On the binomial(m,2) unordered-pair vertices, the rotation subgroup alone
   admits exactly 2^((m-1)(m-2)(m+1)/8) invariant tournaments.
3. The full induced signed dihedral action admits none when m>=5. At m=3
   it admits two.

*Proof.* On Z/m, nonidentity rotations fix no vertex; a reflection fixes
exactly one vertex and exchanges every other point with a distinct point.
An unordered-pair stabilizer is either trivial or generated by the
reflection exchanging its endpoints, with the required negative sign.
There are (m-1)/2 distance classes.

No nonidentity rotation fixes an unordered pair of original vertices: a
setwise fixation either fixes a vertex or exchanges two, both impossible
for a nonidentity odd-order rotation. Thus the cyclic action is free on
the new vertices, and on their unordered pairs as well. If v=m(m-1)/2,
the number of pair orbits is binomial(v,2)/m, the claimed exponent.

Reflection x->-x fixes the (m-1)/2 new vertices {a,-a}. For m>=5 it fixes
two distinct new vertices pointwise, so it would have to reverse the arc
between them while leaving its endpoints fixed. This violates Proposition1.
At m=3 the pair action is equivalent to the original action. QED.

At m=5 this gives **4** signed dihedral tournaments on K5's vertex set,
**512** C5-invariant tournaments on Petersen's ten pair vertices, and
**0** with that signed dihedral lift. The choice of edges of the Petersen
graph is not used to hide missing tournament pairs: the target tournament
has every one of the45 pairs oriented.

Every homomorphism A5 into the automorphism/anti-automorphism group of a
tournament is trivial: its composite sign map to C2 is trivial, since A5
is simple nonabelian, and then a nontrivial image would embed A5 in an
odd-order group. Therefore allowing one global converse sign does not
repair the full A5 obstruction. C5 is a precise surviving action.

For K3,3, retain its two parts and the C3 x C3 action rotating each. There
are three unordered-pair orbits: the edges of each internal triangle and
all nine crossing pairs. The eight invariant tournament completions are
ordinal joins of two directed triangles, in either order. Every one has
H=3*3=9 and is nonstrong. Thus the same original undirected graph can
support a useful odd action while its completion has very restricted
Hamilton behavior.

**Connection contract.** Source: a dihedral action on a cyclic tournament.
Target: its induced action on pairs, including Petersen at m=5. Map:
g.{u,v}={gu,gv}. Preserved: the group action, cyclic orbits, and disjointness.
Destroyed: signed orientability, because fixed-point counts change.
Sidecar: epsilon and the action of every two-point stabilizer. Cheapest
test: count two fixed pair-vertices of one negative involution. This is
useful before imposing a symmetry reduction on any tournament recurrence.

## 2. The literal set, composite eta levels, and a recurrence boundary

Let L={N>1:24/(N-1) is a positive integer}. Elementary divisor arithmetic gives

    L={2,3,4,5,7,9,13,25}.

The user's literal set S={4,2,5,3,7,13} is exactly the subfamily where the
exponent r_N=24/(N-1) is **even**. The two omitted levels9,25 have exponents3,1.
The equality S={d+1:d|12} is already in [THM-4838](../../../../01-canon/theorems/THM-4838-von-staudt-clausen-as-reciprocal-sum-thresholds-orbifold-euler-characteristics-triangle-groups-2-3-p-and-the-twelve-behind-the-owners-sets.md).
This is an exact property of the supplied set, not a claim that this was its
original reason for occurring or that it classifies tournament counts.

For every N in L, set

    f_N(tau)=(eta(tau)/eta(N tau))^r_N=q^-1 F_N(q),
    F_N(q)=product_(d>=1,N does not divide d)(1-q^d)^r_N,
    F_N(q)=sum_(n>=0) a_N(n) q^n,   a_N(0)=1.

**CITED.** These eight functions are principal moduli for X_0(N), as listed
in [Aricheta, AppendixB](https://arxiv.org/html/1809.07421).
This is X_0(N), not its Fricke quotient. A second primary route is
[Savitt, Theorem1.1 and Remark1.2](https://math.jhu.edu/~savitt/papers/pdfs/newman.pdf):
the two eta congruences here are -24 and24; weight is0 and N^r_N is a square,
so the character is trivial. The only pole is the simple pole at infinity.
For prime N all other cusps have positive order; for N=4,9,25 the intermediate
denominator sqrt(N) gives order0. Hence the function has degree1 to P^1.

**Proposition 3, PROVED IN NOTES (prime case inherited; composite extension explicit).**
For all n>=1,

    n a_N(n) = -r_N sum_(k=1..n)
       [sigma_1(k)-N*1_(N|k)*sigma_1(k/N)] a_N(n-k).

*Proof.* Take q times the logarithmic derivative of the displayed product.
The coefficient of q^k is -r_N times the sum of divisors of k not divisible
by N. Such excluded divisors are exactly N times divisors of k/N. Equate
coefficients of q F'_N = F_N (q F'_N/F_N). QED. This argument works for
composite N too; it does not assume the complement of N-multiples is the
set of integers coprime to N. For instance level4 includes d=2.

**Proposition 3a, PROVED IN NOTES (integral roots and exact dyadic depth).** More
generally let N>1 and r>=1 be any integers, and put
F(q)=product_(N does not divide d)(1-q^d)^r. For every integer s>=0,

    F has a 2^s-th root in 1+q Z[[q]]  iff  2^s divides r.

Moreover the largest 2-power by which all nonzero exponents of F mod2 are
divisible is exactly 2^v2(r). Its first nonconstant odd coefficient is at
q^(2^v2(r)). If F=sum a(n)q^n then

    v2(a(n)) >= max(0,v2(r)-v2(n)),   n>=1,

with v2(0)=infinity; equality holds at n=1.

*Proof.* The q coefficient of F is -r. If F=G^(2^s) with G=1+bq+..., it
equals2^s b, proving necessity; dividing every Euler exponent by2^s gives
an integral root and proves sufficiency. Write r=2^v u with u odd and
H=product_(N does not divide d)(1-q^d). In F2[[q]],
F=H(q^(2^v))^u. Since H=1-q+..., the first nonconstant term is exactly
q^(2^v), proving both support claims. The logarithmic derivative gives
n a(n)=-r times an integer; divide the resulting2-adic valuation inequality
by n. Finally a(1)=-r is nonzero and gives equality. QED.

For our eight Hauptmoduls the maximal integral-root depths in level order
2,3,4,5,7,9,13,25 are **3,2,3,1,2,0,1,0**. Thus the literal set S is
exactly the integral-square family. Equivalently, among these eight all
odd-index coefficients are even precisely at the six supplied levels.
Levels9 and25 fail already at q. This is a precise dyadic correspondence
between exponent division, formal roots, Frobenius support, and a coefficient
valuation bound. It does not equate eta coefficients with tournament counts.
The prime-level recurrence is inherited; this formal-root interpretation was
not found in the bounded search of the cited eta packet and its nearby
Frobenius/power sidecars. No historical novelty for the elementary power
argument is claimed.

**Proposition 4, PROVED IN NOTES (recurrence-type obstruction).** F_N(q) is not a
rational function of q. In particular its coefficients cannot, even
eventually, have the form u^T M^n v for one fixed finite matrix M over C,
or satisfy a homogeneous constant-coefficient linear recurrence.

*Proof.* Put tau=1/N+iy with y>0 and gamma=[[1,0],[N,1]] in Gamma_0(N).
Then z=gamma^-1 tau=-1/N+i/(N^2 y). Modular invariance and the simple
q^-1 pole give

    |f_N(tau)|=|f_N(z)| ~ exp(2*pi/(N^2*y))  as y->0+.

Meanwhile q=exp(2*pi*i*tau) tends radially to exp(2*pi*i/N), and multiplying
by q does not change this growth. A rational function has at most a
finite-order pole at that point and hence grows at most polynomially in
1/y, a contradiction. A fixed matrix output has generating function
u^T(I-qM)^-1v; an eventual recurrence only adds a polynomial initial part.
Both are rational. QED.

This is not an obstruction to nonlinear recurrences, growing state,
nonconstant transition matrices, or coefficient extraction from a larger
object. It identifies a necessary difference between a fixed finite-port
tournament transfer and the eta recursion. Adding or removing finitely
many initial coefficients does not repair the mismatch.

**Connection contract.** Source: eta product with retained N. Target: exact
coefficient recursion with a growing convolution history. Map: logarithmic
derivative. Preserved: every coefficient, divisor filter, and level.
Destroyed if one keeps only r: the forbidden-divisor schedule. Needed
sidecar: N and all earlier coefficients (or another proven sufficient
state). Cheapest hostile to a fixed-state transfer: the modular cusp above.

## 3. What 23/24 and the sporadic prime layers really retain

The 23->24 bridge is already exact in [THM-4857](../../../../01-canon/theorems/THM-4857-the-gf2-automaton-of-a-tournament-classifies-the-merged-metagraph-exactly-through-order-five-circulants-are-additive-cellular-automata-and-paley-tournaments-give-quadratic-residue-codes.md):
take the binary row span of the Paley23 adjacency matrix and append a parity
coordinate. Our complete word enumeration independently recovers dimension12
and extended weights

    1 + 759 z^8 + 2576 z^12 + 759 z^16 + z^24.

This retains a code; it forgets the adjacency matrix and generally many
tournament distinctions. No M24 action on the original tournament follows.
Indeed any faithful action of an even-order group on a tournament is
excluded above. The code's automorphisms and the tournament's automorphisms
are different groups. Levels23,24,25,26 also cannot be merged as a uniform
eta-ratio family: only25 belongs to L. Other legitimate modular roles
require their own group/level specification.

**CITED**, with current primary sources checked: the prime support of the
[Monster](https://brauer.maths.qmul.ac.uk/Atlas/spor/M/) is
{2,3,5,7,11,13,17,19,23,29,31,41,47,59,71}; the
[Baby Monster](https://brauer.maths.qmul.ac.uk/Atlas/spor/B/) omits exactly
{29,41,59,71}. Ogg's condition says that **all** supersingular j-invariants
in characteristic p belong to F_p, not merely F_(p^2), and is equivalent
to genus0 of X_0(p)^+. [Duncan--Ono, Theorem1.1](https://arxiv.org/pdf/1411.5354)
gives a genuine map from suitable order-p McKay--Thompson series, via U(p)
and reduction mod p, to rational functions whose poles are the supersingular
j-invariants other than0 and1728. It does not identify tournament path
counts with those coefficients or construct a sporadic action on prime labels.

The first25 primes are the primes below100; their intersection with Ogg's
list has15 elements, leaving {37,43,53,61,67,73,79,83,89,97}. The detailed
square-sum graph and its hostile null models are inherited from
[THM-4860](../../../../01-canon/theorems/THM-4860-the-prime-square-sum-graph-is-bipartite-by-residue-mod-4-two-primes-1-mod-3-are-never-adjacent-chebyshev-bias-obstruction-to-prime-square-sum-paths-none-below-615731-and-the-25-primes-below-100-by-ogg-layer.md),
not recomputed or promoted here. In particular4 is a composite eta level,
while11 is an Ogg prime of positive X_0(11) genus. Substituting one for the
other erases precisely the distinction this packet needs.

### The recovered 23/24/25/26 chain has an actual lattice map

This chain is **inherited PROVED/CITED**, not a new result of this session.
The [audited Golay-to-Leech construction](../../session97_overnight_arithmetic_carriers_20261029/leech/GOLAY_LEECH_SUM_GUARD.md)
supplies the missing map after the23->24 parity extension. For its extended
code C, take integer24-tuples x with common coordinate parity m in {0,1},
residue support {i:x_i=m+2 mod4} in C, and sum x_i=4m mod8; divide by sqrt8.
The result is the Leech lattice Lambda. The sum guard and odd-coordinate glue
are essential: ordinary binary Construction A instead has48 norm-two roots
and gives the A1^24 Niemeier lattice. Thus a code length24 is not by itself
a Leech construction.

Now add the rank-two hyperbolic plane U with Gram matrix [[0,-1],[-1,0]].
The orthogonal sum Lambda plus U has **rank26 and signature(25,1)** and is
the even unimodular lattice II_(25,1). Conversely, in its standard diagonal
coordinates with norm sum_(i=0..24) x_i^2-t^2, put

    w=(0,1,...,24;70),       z=(1,0,2,...,24;70).

The coordinate sums are300, the parity check is300-70=230 even, and
sum_(j=0..24) j^2=4900=70^2. Hence w and z are null and w.z=-1.
They span a unimodular hyperbolic plane; each is primitive because pairing
with the other is a unit. The complementary projection is explicitly

    x -> x+(x.z)w+(x.w)z.

It lands in their common orthogonal complement, and integrally splits off U.
The null quotient (w^perp intersect II_(25,1))/Z w is the Leech lattice.
This identification is the classical construction in
[Conway--Sloane, *Lorentzian forms for the Leech lattice*, p.216](https://neilsloane.com/doc/Me91.pdf);
the direct-sum form and sign convention are also explicit in
[Borcherds, *Automorphic forms and Lie algebras*, section1](https://math.berkeley.edu/~reb/papers/cdm/cdm.pdf).

[THM-4752, tournament scores in II_(25,1)](../../../../01-canon/theorems/THM-4752-tournament-scores-in-ii-25-1-backward-length-and-the-a24-deep-hole.md)
already gives the actual tournament interface: on25 vertices the score lift
v_T=(s;70) has norm -2c3(T), so the lift is null exactly for transitive T.
For general cyclic T, v_T is not in v_T^perp, so v_T^perp/Zv_T is not a
defined quotient. The always-defined A24-cusp projection in that theorem is
a different map, retaining scores. Its canonical hostile has identical
25-vertex score lifts but H=657 versus643. MISTAKE-661 and the
[subsequent geometry audit](../../session28_modular_spectral_20261022/session28-lorentzian/GEOMETRY.md)
also distinguish a cusp pairing from hyperbolic distance; a null cusp lies
at infinite hyperbolic distance. The older norm sign and primitive-vector
conventions used here agree with the repaired canon.

The roles are consequently explicit: **23** is the punctured code length
and cyclic coordinate action; **24** is the extended code length and Leech
rank; **25** is the positive index and the number of score coordinates;
**26** is total Lorentzian rank. The eta level N=25 is a separate parameter.
The chain preserves a named coordinate action through gluing, then a lattice
and its bilinear form; it loses tournament arcs before reaching the lattice.
No Hamilton-spectrum or Monster action follows from those ranks.

## 4. Exact universe, controls, and remaining frontier

[check.py](check.py) uses only the Python standard library and explicit
exceptions, so `-O` does not remove validation. It compares signed orbital
transport with an independent GF(2) constraint elimination for every stated
action at odd m=3,5,7,9,11,13. For m=5 it enumerates all512 labelled
C5-invariant tournaments and counts actual Hamilton paths by endpoint DP;
these are not512 isomorphism classes. The m=3 pair lift is a positive control
for the m>=5 obstruction; A5 and the reflected ten-point action are hostiles.
K3 x K3 is not used: the six-vertex test is explicitly the two-part C3 x C3 action.
The eta test compares all coefficients through degree128 by logarithmic
recursion and direct finite-product multiplication at every N in L. It
also checks each integral dyadic root by independent polynomial powering,
the first forbidden root by its q coefficient, exact mod2 support depth,
and the coefficient valuation inequality. It includes composite4,9,25 and
the odd-exponent boundary. Complete Paley7/23
codeword enumeration supplies the coding controls.

Reproduce from this directory:

    python3 check.py > check.json
    python3 -O check.py > check.optimized.json
    cmp check.json check.optimized.json

**OPEN:** transporting a sporadic/modular object to an exact tournament
path-spectrum statement. The operational advance is narrower: symmetry
reductions must pass signed two-point stabilizers, and a modular coefficient
recurrence cannot silently be replaced by a fixed finite-matrix transfer.
These decide two proposed connections before expensive enumeration.
