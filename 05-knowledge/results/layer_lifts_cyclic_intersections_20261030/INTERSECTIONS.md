# Prime cyclic banks: point intersections after legal native moves

**PROVED**, using the **CITED** classical prime-degree theorem of Burnside.
**FINITE-EXACT** controls are separately listed below. These are bounds for
a specified family of labelled affine charts. No bound on the number of
charts needed for unrestricted A000568 is claimed.

## 1. The native objects and the retained information

Fix an odd number n of named vertices. A tournament is the vector of its
binom(n,2) edge bits, where the smaller endpoint pointing to the larger is1.
Relabeling acts affinely, with endpoint-order sign flips. Its linear part
permutes unordered edges; it fixes the all-ones vector J. Converse is
translation by J. Ties do not occur: this is the intrinsic tournament relation.

For a regular cyclic subgroup H of S_n, let F_H be the affine space of
tournaments fixed by H, and U_H its homogeneous direction space. On an odd
cycle, U_H consists of functions constant on the (n-1)/2 distance-pair edge
orbits. Thus dim U_H=(n-1)/2. F_H has this direction and is nonempty.

A **native program chart** starts in F_H, follows a specified finite
sequence of actual Hamilton-path complements and relabelings, and keeps
exactly the starting points for which every requested path is legal.
A path complement pins the path's directed edges and flips every other
edge. It is a coordinate constraint followed by a fixed bit translation.
The program chart is therefore an affine flat (or empty), with direction
contained in U_K, where K is the transported cyclic parent. If at least one
native move occurred, its direction does not contain J: the last path pins
an edge, and subsequent relabelings preserve J. Source labels, chosen paths,
endpoint ports, and the transported parent are part of the chart record.
We take exact unions only when forgetting construction history.

This is the carrier from
[the preceding packet](../anchor_content_flat_quotients_20261030/TOURNAMENTS.md).
Arbitrary substitutions or atom attachments can leave this single-parent
family; they require a separate direction-space analysis.

## 2. The prime-order point-intersection theorem

Let p be an odd prime and H,K distinct regular cyclic subgroups of S_p. Then

    U_H intersection U_K = span(J).                    (1)

Consequently two nonempty native program charts with distinct transported
parents, each having undergone a legal path complement, intersect in at most
one labelled tournament. In fact one chart excluding J already suffices.
This holds at every program depth and for different legal paths and ports.

**Proof.** The group G=<H,K> is transitive of prime degree. The classical
Burnside theorem says G is doubly transitive or has a normal Sylow p-subgroup;
see Dobson–Witte, *Transitive permutation groups of prime-squared degree*,
[Theorem1](https://arxiv.org/pdf/math/0012192), which states that input.
H and K are distinct Sylow p-subgroups of G, since a Sylow p-subgroup of S_p
has order p. Hence the normal alternative is impossible. Double transitivity
makes G transitive on unordered edges. The common homogeneous invariant
functions are therefore just the constants, proving(1). An intersection of
affine flats is empty or a translate of the intersection of their directions.
The direction of a native chart excludes J, so the latter intersection is0.

Without a path pin, arbitrary translates of U_H and U_K can meet in two
points, a converse pair; taking the translates through the same point attains
this bound. The original unshifted banks F_H and F_K are disjoint: a tournament
cannot be invariant under a doubly transitive group, which can reverse an arc.
The distinctions between original banks, arbitrary translates, and legal
native charts are essential.

The classical group input is also used in
[THM-4644, prime circulant rigidity](../../../01-canon/theorems/THM-4644-fermat-primes-rigid-circulant-tournaments-wieferich-orbit-explosion.md).
The new step here applies it to **differences** of native-chart points, not
to automorphisms of a resulting tournament after complement.

## 3. A bound on the intersection representation

Let M charts on p vertices be grouped by transported cyclic parent H, and let
V_H be the union of the charts with that parent. Let X be the set of points
belonging to at least two different V_H, with owner multiplicity
m(x)=#{H:x in V_H}. There is an exact integer indicator identity

    1_union = sum_H 1_(V_H) - sum_(x in X) (m(x)-1)1_{x}.       (2)

Its cross-parent correction has support at most
sum_(H<K) M_H M_K <= binom(M,2), by the singleton-intersection theorem.
Correct by owner multiplicity, **not** by the number of overlapping chart
pairs: several pairs can describe the same point, even within one parent.

Inside a fixed parent, separate the disjoint cosets of U_H. Choose one origin
in each. Every chart in that coset is an axis-aligned subcube in the
r=(p-1)/2 distance-orbit coordinates. Indeed each original edge pin fixes
one such coordinate, and translations change only its required value.
Intersections of these subcubes remain subcubes. There are at most3^r
nonempty coordinate descriptions, since each coordinate is free, zero, or
one. At most M_H cosets are occupied. Exact inclusion-exclusion can therefore
be coalesced into at most3^r M_H internal flat terms. Combining with(2)
gives an indicator representation with at most

    3^((p-1)/2) M + binom(M,2)                         (3)

terms before further coincidences are merged. Coefficients can be large;
(3) bounds the number of distinct retained supports, not their bit cost.

This improves a generic2^M subset-intersection bound **for this carrier**.
It does not bound M, remove the exponential dependence on p, or prove that
all tournaments lie in a small such union. Enumerating the points in each
small chart can also be cheaper than retaining its full arrangement; no
optimality claim is made for(3).

For an S_p-invariant union, apply the exact indicator to every permutation's
affine fixed set, then Burnside. For marked ports use their intended group.
One must first saturate under that group if the union is not invariant.
The incoming Session112
[cumulative-cut ranking](../session112_counting_interfaces_20261030/counting/README.md)
can use the same representation: prefix cuts only add coordinate pins,
relabeling transports parents, and the singleton bound survives. This removes
one possible intersection explosion in that oracle; the number of saturated
charts and the full-group sum can still be enormous. The separate
[shared-residue bank ranking algorithm](../session112_counting_interfaces_20261030/counting/BANK_RANKING.md)
has a different, explicitly retained bank action and does not erase this cost.

## 4. Composite orders: an all-size hostile with identical ports

Let n=3m with m odd and m>=3. On Z_n put

    a(x)=x+1,
    v(x)=x+3 if x=0 modulo3, and v(x)=x otherwise,
    b=v a v^-1.

The regular cyclic groups <a>,<b> are distinct. Their fixed tournament banks
share the family C3[R], where R is any cyclic tournament on m vertices,
with fibres given by residue modulo3. All have the common Hamilton path
P=(0,1,...,n-1). After the same legal native complement along P, the two
charts intersect in exactly

    2^((m-1)/2) labelled tournaments.                  (4)

The named endpoint ports are the same0 and n-1 throughout.

**Proof.** Write x=3t+i, i in{0,1,2}. The permutation b a^-1 translates the
three fibres by(1,-1,0). Conjugating by a rotates those components. They
generate all translations(t_0,t_1,t_2) with t_0+t_1+t_2=0 modulo m; the
projection to any two fibres is surjective. Thus all cross-fibre unordered
edges form one orbit. Within a fibre, the difference up to sign is invariant,
and translations together with a connect all three copies. There are
(m-1)/2 within-fibre edge orbits, hence(m+1)/2 edge orbits altogether.

An invariant orientation chooses a direction for the outer triangle and one
for each internal distance pair; it is exactly C3[R] or its reversed outer
triangle. The path P uses only cross-fibre edges and fixes that outer
orientation. Every internal choice remains free. Translation by the same
complement mask preserves intersection size, proving(4). Distinctness of
the cyclic groups can also be seen directly: b sends0 to n-2 and1 to2, so
it is not any constant translation a^j when n>=9.

At n=9,15,21,27 the intersection sizes are2,4,8,16. In particular the
twenty-seven-vertex example is not an isolated defect of a small search.
Its block structure causes the extra common directions at all sizes.

Do not identify this cyclic C_27 construction with the Paley tournament on
the additive group of F_27, whose regular subgroup is C3^3. The root-lattice
map q=27 -> C=7 uses the doubly regular Paley object, not this cyclic hostile.

## 5. Finite controls and a multiplicity warning

[cyclic_intersections.py](cyclic_intersections.py) enumerates every regular
cyclic subgroup at p=3,5,7, all its constant-step path charts, and every
legal Hamilton-path chart at p=5. It checks group edge orbits by union-find,
affine equations by exact F2 elimination, intersections by explicit bitsets,
and composite invariance by a direct adjacency implementation.

The all-path p=5 control has6 parents and240 chart records. The sum of parent
union sizes is360 but the labelled union has280 points. Exactly40 points
have three parent owners each. Thus the correction is40*(3-1)=80.
There are120 nonempty cross-chart pairs; subtracting120 would be wrong.
The resulting280 agrees with the preceding packet's first native layer.

The controls total35,483 explicit checks. The theorem for all prime p uses
the group argument, not extrapolation from those orders. The full composite
family uses the fibre-translation proof, not extrapolation from four examples.
