# Independent review: mixed births at fixed physical ports

**PASS**, bound to the six source hashes in SOURCE_SHA256SUMS. No correction was needed. The source files were copied before execution and were not changed.

## Exact scope

The interface retains one fixed set of original physical ports. Births may be left or right with arbitrary masks, but newborns are never promoted and ports are never forgotten or replaced. The all-order obstruction concerns repetition of one fixed nonempty word, from a nonempty actual initial tournament: a nonconstant eventual law H_(n+1)=3H_n+c forces one direction to occur exactly once and c=2 modulo4. It does not cover adaptive words, changing ports, reordering to a new operation model, or rebuilding graphs. Strongness is not assumed.

## Mathematical checks

1. A retained cover satisfies c-P<=1. Deleting a newborn from an internal position restores at least one port incidence and raises c by one; endpoint deletion restores incidences without increasing c. Thus every output predecessor is retained. Unordered components and repeated (*,*) components have exactly the multiplicities implemented by the indexed predecessor/successor choices. Shapes encode everything needed for both insertion directions.
2. For right births, 2k-P plus the number of nonport starts is nondecreasing, at most2k+1. Equal-potential transitions only append at a nonport end and leave the shape fixed. Their coefficient is the number d of such ends, at mostk+1, independent of mask. The left statement is its reversal. This supplies uniform finite annihilators including zero eigenvalues and Jordan blocks.
3. Opposite-side births commute on actual graphs: their new cross arc is left-to-right in either temporal order, all old incidences agree, and swapping the two new labels fixes every old vertex and port. At the cover level the same bijection applies to two insertions into each component configuration; deleting either new label gives a retained intermediate cover. Thus the formal matrices commute as well, including their formal coordinates. Grouping a period's two side-subwords does not reorder masks within either side.
4. Simultaneous triangularization gives period eigenvalues d^r e^l. A nonconstant eventual ternary affine sequence has a genuine3 mode, so at least one nonempty side length must be1. The absent-side case is treated without0^0.
5. With one left birth per period, the clone root loses only to at mostk original ports. Each internal root-to-root fragment ends at a distinct such port, giving the uniform bound j<=k+1. This is crucial: it excludes a growing-index tail in the Stirling expansion.
6. The triangular Stirling inversion uses fixed positive clone multiplicities and uniform right-only recurrences. Positivity kills j>=4 eventually; the only bounded surviving nonzero recurrence base is1, so a3 is eventually constant and even. The growth bound on a2 eliminates every base>=2. For r>=2 the a1 bound eliminates every base>=4. For r=1, a1 is itself a pure clone sequence, whose leading ternary coefficient is even. Polynomial factors and nilpotent transients cannot evade these steps. Odd Hamilton counts then force the constant term beta odd, hence c=-2beta=2 modulo4.
7. The explicit mixed +2 example is correct. After a nonempty period it has a universal transitive source prefix, and therefore is not strong. That is consistent with the theorem's arbitrary-tournament scope; it is not presented as a strong returning clock.

## Independent computation

The frozen verifier ran both normally and with Python -O. Both outputs equal the frozen receipt byte for byte.

The separate independent.py uses endpoint subset DP and unordered set partitions, not the source's path-permutation/cut enumeration. It reproduces all2,764 source profiles and33,194 literal birth vectors through order4 with up to3 ordered ports. It also checks115,711 physical two-birth commutations directly on every arc. A fifth-port stress check, beyond the source's k<=4 potential tests, verifies all9,147 shapes and3,586,272 nonzero transition entries against the potential and mask-independent diagonal laws. See INDEPENDENT.json.

The separate witness.py counts the seed's repeated-root skeletons literally, obtaining a1=179,a2=180,a3=60 and all higher coefficients zero. Thus its clone count is179+180(2^m-2)+60(3^m-3*2^m+3)=60*3^m-1. Independent forward subset DP gives179,539,1619,4859,14579 for zero through four mixed periods, and every arc matches the transitive-prefix/clone model. See WITNESS.json.

The finite checks validate the implementation and named witness. The all-k, all-period and eventual assertions are justified by the proof, not extrapolation. No claim of minimal interface dimension or unrestricted spectrum completeness is approved.
