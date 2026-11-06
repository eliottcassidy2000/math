# A base-eleven scalar witness bound and exact forward stabilization

Status: PROVED relative to the CITED Busch bound, with independently verified
finite witnesses and exact integer checks. Depends on the cited Busch
strong minimum and finitely many explicit small tournament witnesses. This
replaces a witness by a different tournament with the same H, and need not
preserve its strong components, isomorphism class, or metagraph position.

## Theorem

Every attained integer h>1 has a witness U of order m satisfying

    3^5 * 11^m <= 11^3 * h^5.                         (1)

Consequently the minimum witness size obeys

    m(h) <= floor(3 + 5 log_11(h/3)).                  (2)

The separate case h=1 has m(1)=1. This theorem is conditional only on h being
attained; it does not settle which positive odd h are attained. The universal
constant in (2) is not claimed optimal.

## Proof: normalize all but four scalar factors

Delete singleton SCCs from any witness and write its H as the product of
nontrivial strong-component counts. There are no strong components of order 2.
Call a block of size s and count d good when

    d^5 >= 11^s.                                     (3)

Every strong component of order s>=9 is good. Indeed Busch's exact minimum
f(s) has f(9)=75, f(10)=125, f(11)=225 and f(s+3)=5f(s). Check the three
integer inequalities f(s)^5>=11^s, and the step multiplier 5^5>=11^3.
Notably 75^5=2373046875 >= 11^9=2357947691.

For orders 3 through 8 use Busch's minima 3,5,9,15,25,45 and oddness of H.
Order3 has count 3 and order 4 count5. At order 5, counts>=11 are good.
At order 6, counts>=19 are good; replace count 15 by a 5-vertex count 15 witness.
At order 7, counts>=29 are good; replace counts 25 and27 by 6-vertex witnesses.
At order 8, counts>=47 are good; replace count 45 by a 6-vertex witness.
All four replacement blocks are good. Each comparison is an explicit integer
fifth-power inequality recorded in the checker. The only remaining nongood
scalar factors are

    3, 5, 9, 17,

with witness sizes 3,4,5,6. Pool these as a single scalar product 3^a5^b17^c,
counting a factor 9 as two factors 3. This step remembers H but not the original
component partition.

## Proof: group factors and discharge twelve residual cases

For each complete triple in a, use a 6-vertex count 27 block. For each pair
in b, use a 6-vertex count 25 block. For each pair in c, use an 8-vertex count 289
block. All three are good. The residual scalar is

    D=3^r5^t17^u, 0<=r<=2, 0<=t,u<=1.

For these twelve possibilities use the following witness orders:

    D:     1  3  5  9  15  17  45  51  85  153  255  765
    cost:  0  3  4  5   5   6   6   7   7    7    8   11

The empty residual has cost 0 and is omitted. The residual 765 is a dominance
join of count 3 on 3 vertices and count 255 on 8 vertices. All others have direct
adjacency-code witnesses in witnesses.json. Each row satisfies exactly

    3^5 * 11^cost <= 11^3 * D^5.                     (4)

The checker uses integer powers, and verifies the physical witnesses by
permutation enumeration. Concatenate all good blocks and the residual witness
by dominance joins. Path counts multiply, and block orders add. Multiplying
(3) for all good blocks and (4) for the residual proves (1). The resulting
nonempty tournament has count h, as required. QED.

## Exact forward consequence and decision complexity

For every n>=1, every globally attained h>1 satisfying

    11^3 * h^5 < 3^5 * 11^(n+1)                      (5)

already occurs at order n. Formula (1) forces a witness order less than n+1;
pad with universal sources to order n. Thus S_n and the global spectrum agree
through this explicit strict window. This is exact stabilization, not an
assertion that all odd integers in the window occur.

For any proposed h>=3, enumerate tournaments through the largest integer B
with 3^5*11^B<=11^3*h^5. The first occurrence gives m(h); if none occurs, h is
absent globally. This cutoff is valid without knowing attainability in advance.
It remains a 2^O(log^2 h) algorithm, but shrinks the exponent's constant versus
the inherited SCC-preserving cutoff floor(3log_3 h). The strong-atom divisor
algorithm may still be preferable on individual inputs.

## Provenance and boundary

Busch's minimum and the inherited SCC-preserving 3log_3 h cutoff appear in
05-knowledge/results/spectrum_basins_20261020e/spectrum.md sections2–3.
The small witnesses are extracted from the independently audited
session207/size-spectrum/census.json, then checked again here by brute force.
The first scalar-compression argument of this session gave
3^(2m-1)<=h^5. Formula (1) strengthens that bound pointwise for h>=3, since
5/log11 < (5/2)/log3 and both real-valued bounds equal3 at h=3.
The key additional operation is scalar replacement: the original reduced
core may still have 3 log_3 h vertices (a chain of triangles), so it would be
false to impose (1) on that unchanged core. No novelty claim for Busch or the
small count spectra is made.
