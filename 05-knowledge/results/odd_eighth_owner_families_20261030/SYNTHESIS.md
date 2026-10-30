# Odd eighths: exact maps, distinct operations, and retained owners

**PROVED elementary identities / CITED classical and repository mechanisms /
OPEN unconstructed bridges.** No literature-priority claim. The productive
connection is a shared, explicitly labelled phase carrier together with its
consumer. A matching fraction alone supplies no transport theorem.

## 1. A four-character dictionary

Let U_8={1,3,5,7}. Multiplication modulo eight makes it the Klein four group.
Define chi_-1(a)=(-1)^((a-1)/2), chi_2(a)=(-1)^((a*a-1)/8), and
chi_-2=chi_-1*chi_2. These four rows are an orthogonal Walsh basis:

| Character | 1 | 3 | 5 | 7 | Sum chi(a) zeta_8^a |
|---|---:|---:|---:|---:|---|
| 1 | 1 | 1 | 1 | 1 | 0 |
| chi_-1 | 1 | -1 | 1 | -1 | 0 |
| chi_2 | 1 | -1 | -1 | 1 | sqrt(8) |
| chi_-2 | 1 | 1 | -1 | -1 | i sqrt(8) |

This classical character/Gauss-sum table is [Keith Conrad, Example 3.4,
printed page 5](https://kconrad.math.uconn.edu/blurbs/gradnumthy/Gauss-Jacobi-sums.pdf).
The last two characters are primitive modulo eight; the first two factor
through modulo four. Our checker recomputes the table with exact integer
coordinates in (1/sqrt(2),i/sqrt(2)), without numerical complex arithmetic.

Two useful decoders, for odd v and a, are

    ||va/8|| = 1/4 - chi_2(v)*chi_2(a)/8,
    zeta_8^a = [chi_2(a)+i*chi_-2(a)]/sqrt(2).

The first retains just the reflection-even bit. The second retains all four
oriented phases. Reflection has the two orbits {1,7} and {3,5}; consequently
every reflection-even scalar observable on this set is alpha+beta*chi_2.
This last classification also appears explicitly in the incoming
[PHASE_BRIDGE, B2](../odd_eighths_uniform_families_20261030/PHASE_BRIDGE.md),
which is credited rather than reclaimed here.

The field is Q(zeta_8)=Q(i,sqrt(2)): i=zeta_8^2 and sqrt(2)=zeta_8+zeta_8^-1,
while zeta_8=(1+i)/sqrt(2). Its Galois action is a->u*a for u in U_8,
with independent signs on i and sqrt(2). These elementary equalities
explain the two quadratic bits; no elliptic-curve or Sha conclusion follows.

## 2. Three different clocks on the same four labels

| Operation | Map | Structure retained or lost |
|---|---|---|
| Geometric quarter-turn | a -> a+2 modulo 8 | A four-cycle; retains the circular ordering |
| Cyclotomic automorphism | a -> u*a modulo 8 | U_8 is V_4; every nonidentity element has order two |
| Squaring the complex root | zeta_8^a -> zeta_4^a | Two-to-one, leaving two primitive fourth roots |
| Squaring the unit label | a -> a*a modulo 8 | Four-to-one, constant at 1 |

Thus rotation, Galois action and the two squaring maps must not be silently
identified. In particular neither characteristic-zero squaring map is a
field Frobenius. This distinction is a cheap hostile for a proposed
doubling/squaring or additive/multiplicative analogy.

## 3. The tournament connection is an actual matrix identity

[PROOFS, T1](PROOFS.md) proves that the Cayley transform of the transitive
tournament matrix is the signed cyclic shift W, with W^n=(-1)^(n-1)I.
At order four its eigenphases are the four odd eighths, exactly. At every
dyadic order n=2^(a-1) they are the primitive 2^a-grid, exactly the maximizing
time set of the AP speed family {1,...,2^a-1}. The map identifies phase labels,
not vertices with runners or arbitrary attachment operations.

This extends the explicit order-four CM calculation already present in
[Session37, ARITHMETIC_CUT_COMPLETION, sections 4 and 4B](../session37_statistics_parity_cuts_20261023/arithmetic/ARITHMETIC_CUT_COMPLETION.md):
K=-S^2 has eigenvalues 3+2sqrt(2),3-2sqrt(2), each twice, and the Cayley field
is Q(zeta_8). That file also supplies the finite order-five census. The
classical classification and count of switching-transitive local orders
are already in [THM-468, determinant floor](../../../01-canon/theorems/THM-468-tournament-determinant-floor.md).
Here the explicit signed-cycle criterion explains their recurrence
L_(n+1)=2n L_n and supplies a literal operator model at every size.

The inverse Cayley map retains the entire labelled matrix. Replacing it
by eigenvalues loses the switching gauge, labels and generally much more.
Even the transitive and cyclic triangles have the same skew spectrum but
Hamiltonian-path counts 1 and 3; see also
[THM-1415, canonical star quotient](../../../01-canon/theorems/THM-1415-switching-is-the-canonical-star-quotient.md).
The native 81/27 automorphism control remains in our certificate. We do not
claim it is cospectral; its hostile predicate is equality of labelled cyclic
triples with unequal native automorphism groups.

## 4. Boundary geometry needs integer heights as well as residues

For W(8,E), the top four points are always 1/8,3/8,5/8,7/8. They do not
determine its top windows: adjoining speed 9 replaces physical owner 1 by 9.
For m=5 the empty E gives M=5/41 and E={9} gives M=5/43 after the same
tail 40. Both original cores had M=1/8 and the identical maximizing grid.

The source-to-target map is therefore

    (physical point a/8, actual directional owner c_u)
       -> (primitive residue a, reciprocal slope 1/c_u).

The residue records the phase. The owner records how far a physical time
can move before losing safety. Unit multiplication can permute the first
without preserving the second. O1-O5 classify exactly when the integer
heights can be forgotten for a query, and which future operations expose
them again. The first-tail query recovers every second-band bit, including
nonunits invisible to the whole upper chamber.

The incoming [Session115 harmonic theorem](../session115_odd_eighth_classifications_20261030/lrc/HARMONIC_COMPLETION.md)
already gives the physically weighted odd-eighth decoder

    w_8(j)=1/j+1/(8-j)=88/105+(32/105)*chi_2(j),

on j in U_8 (here j is its own inverse). Its q=7 harmonic completion and
the incoming [source-energy bridge](../odd_eighths_uniform_families_20261030/PHASE_BRIDGE.md)
must retain multiplicity, displacement address and clock parity. We do not
reclaim those formulas. P2 extends the harmonic mechanism to every second-band
core containing q+1, using the actual owner list from O1.

## 5. Niche recoveries and non-identifications

| Recovered source | What the eighths mean | Useful mechanism; boundary |
|---|---|---|
| [THM-3435, dyadic fibre grids](../../../01-canon/theorems/THM-3435-dyadic-fibre-grid-decomposition-for-literal-half-twists.md) | Half-twisted finite grid coordinates | Counts obey a three-state mod-7 carry, but literal coverage needs the affine starting coset. The order-three orbit of 2 is not a ternary branching tree. |
| [THM-080, odd Walsh transfer](../../../01-canon/theorems/THM-080-transfer-matrix-odd-walsh.md) | Exact Fourier amplitudes such as 1/8 | A named root determines traversal sign; each unrooted even component contributes its missing factor 2. Amplitudes are not circle positions. |
| [Switching pairing boundary](../switching_branch_integrality_20261020/PAIRING_BOUNDARY.md) | Odd-element proportions and complementary branch weights | An odd number of unadorned types cannot be paired. Decorated centralizer weights, not a visual four-point orbit, are needed. |
| [THM-2598, quartic V4 torsor](../../../01-canon/theorems/THM-2598-quartic-v4-resolvent-torsor-and-universal-cusp-boundary.md) | Four solution sheets with a V4 action | Choosing a sheet gives coordinates; forgetting it is a genuine quotient. No identified map to our owner integers or a Keller consequence. |
| [Unit-distance holonomy](../diagonal_palette_holonomy_20261011g.md) | Signed incidence and reflection seams | Actual fixed points must survive a quotient. This is a shared proof discipline, not a bijection with LRC windows. |
| [Reflection-parity controls](../odd_sign_spectra_trijections_20261028/LRC_UNIT_DISTANCE.md) | Trace of reflection on components | A reflection trace can miss the distinction between isolated points and intervals. It cannot replace literal closed survivors. |

The August Boolean-quartet work uses four disjoint support pairs on an
eight-sheet fibre. Its full coordinate is a choice among those support pairs,
not automatically a label in U_8. The prior [dyadic-fibre reflection](../../../07-reflections/lrc-dyadic-fibre-grid-ternary-carry-and-boolean-quartet-codex-20260815.md)
is provenance; the current theorem above and its correction lineage control
the mathematical statement.

Incoming Session115 also gives two different complete classifications:
[odd eighths in coherent Fano gauges](../session115_odd_eighth_classifications_20261030/quadratic/ODD_EIGHTH_PHASES.md)
are probabilities of Boolean character tests; [affine-cube group masses](../session115_odd_eighth_classifications_20261030/switching/ODD_EIGHTH_GROUPS.md)
come from dimensions of fixed spaces and weighted centralizers. Both are
valuable because they specify their predicates and clocks. No map has
been constructed from those probability values to our physical maximizing
points, so they are not extra instances of O1 or T1.

A later [Session115 synthesis](../session115_odd_eighth_classifications_20261030/SYNTHESIS.md)
does construct a different exact map between those neighboring fields:
for odd-order h, the Boolean event 1[N_h v!=0] is exactly the even-order
event of the affine transformation (v,h). Normalizing a Boolean gauge
requires moving its named character observer. This is now a proved
group/Boolean bridge; it still does not locate an LRC maximizing time.

That source also computes U_(3->5)chi_2=-i*chi_-2 for its particular
word-address unitary on odd eighths. The positive weight
88/105+(32/105)chi_2 therefore becomes nonreal. Our exact character code
replays its four values, while crediting the incoming operator construction.
This is stronger than a caution about a shared denominator: even on the
same named four phases, preservation of a quadratic norm need not preserve
the positive physical weight. The full clock and operation remain essential.

## 6. The shared principle and its decisive tests

For this session the exact synthesis is operational: a residue or spectral
quotient is usable when its consumer factors through it, and it is closed
only when the future operations respect the same quotient. Here we can prove
the exact boundary, rather than leave it as a caution:

* Five eventual scalar answers lift to sixteen states under future unit-owner
  additions; admitting the first-tail query forces all 128 source states.
* First-arrival functions suffice for one shot, but the {4}/{4,28} hostile
  separates after the common addition 32.
* The full Cayley matrix has a rational inverse; its eigenphase multiset
  retains the local-order predicate but not Hamiltonian counts.

The next useful theorem must identify an operation-closed carrier or prove
the necessity of its missing coordinate. Averaging, topology, a shared field
or a repeated fraction cannot supply that coordinate without such a map.
LRC(14), general tournament enumeration improvements, Collatz convergence
and Sha finiteness remain OPEN.
