# Infinite birth families, inverse divisor sums, and quadratic rays

**11 November 2026. PROVED scoped family formulas; FINITE-EXACT complete
coefficients in the stated band; CITED cutoff; INHERITED seed rays.**
Chan-Pak, arbitrary odd-coefficient seeds, and complete all-order
blue/black metagraph classification remain **OPEN**.

The main result is a complete classification of **one-vertex tournament
additions with target Hamiltonian-path count below 125, at every graph
size**. There are 946 nonzero source-count/target-count families, generated
by 872 primitive strong-child kernels. Each family has an explicit rational
generating function and an eventual polynomial counting its marked births.
The primitive records retain actual source adjacency and attachment masks.
The cutoff comes from Busch's strong-order bound, not from extrapolating
a small-order census.

More generally, the proofs classify every birth by its absorbed interval
of strong components, and give rational families for any fixed pair of
path counts once its finite primitive bank is known. This classifies the
vertex-addition relation, distinct from literal fixed-order blue/black lines.

[Complete birth proofs](BIRTH_THEOREMS.md) ·
[Arithmetic and prime-geometry proofs](ARITHMETIC.md) ·
[Inheritance, connection contracts, and open obligations](BOARD.md)

## Concrete results

- **All sizes, target H<125:** 946 nonzero families, with eventual degrees
  1..5 occurring 744,161,34,6,1 times. A separate child-first calculation
  verifies all 687 ordinary strong cores and 5,086 pointed strong cores.
- **33 in the growth relation:** marked births 33->101 number 16n-90
  for source order n>=8; births 33->103 number 8n-42 for n>=7.
- **Actual common sources:** all 200 middle-core types for coprime targets
  below 125 are retained. Exactly 3n-19 sources of H33 support BOTH 101
  and 103 at order n>=8, including eleven at order ten. The H11 sources
  supporting both 35 and 37 number n-4 for n>=6.
- **Inherited Chan-Pak restriction, now counted:** a one-vertex step h->3h+2 or h->3h+4 must
  absorb every nontrivial source component. Its existence is exactly
  positivity of a finite primitive kernel. Uniform positivity above the
  inherited threshold 95 is still open.
- **Size sequences through 64:** the entire absorbed-component-width
  triangle follows from the ordinary and subset-marked Burnside series.
  This avoids generating the underlying tournaments for these statistics.
  It is not a benchmark improvement over the existing A000568 engines.
- **A135244 through target 10,000:** all 9,999 reference entries reproduced
  using a one-million-source divisor sieve plus explicit prime-square and
  semiprime strata. The generic quadratic source bound is 99,980,001.
- **An exact arithmetic bridge:** k^2+1 is prime iff
  A135244(k^2+2)=(k^2+1)^2. In particular A135244(33)=13*19=247,
  explained by the closest prime pair summing to 32.
- **A002496 through input k=100,000:** all 6,656 prime values reproduced.
  The prior seeds construct all quadratic candidate values with
  k=d*2^a, odd d<=63 and a>=4. This proves neither primality along those
  rays nor infinitude of A002496.

## Sequence sample

Source order is n; a birth adds one vertex. All counts are ordinary rooted
isomorphism classes. The last three columns partition the second column.

| n | all births Q_n | singleton component | one component absorbed | multiple components absorbed |
|---:|---:|---:|---:|---:|
| 0 | 1 | 1 | 0 | 0 |
| 1 | 2 | 2 | 0 | 0 |
| 2 | 4 | 3 | 0 | 1 |
| 3 | 12 | 6 | 2 | 4 |
| 4 | 48 | 13 | 18 | 17 |
| 5 | 296 | 36 | 174 | 86 |
| 6 | 3,040 | 148 | 2,296 | 596 |
| 7 | 54,256 | 1,064 | 46,868 | 6,324 |
| 8 | 1,716,608 | 14,848 | 1,589,744 | 112,016 |
| 9 | 97,213,472 | 398,064 | 93,320,536 | 3,494,872 |

The n<=8 rows agree with every attachment orbit of every ordinary source
class; larger rows are exact evaluations of the proved Burnside identities.
The complete formulas, not a claim of global sequence novelty, define the
newly tracked statistics.

## Reproduce and query

From the repository root, with Python 3 and g++ available:

    python3 -B 05-knowledge/results/infinite_birth_families_20261111/reproduce.py

The reproduction uses the two pinned inherited census files and the saved
primary numerical controls. It compiles its C++ helper in a temporary
directory. No additional service or agent is needed.

For an exact count at any nonnegative source order, with odd
1<=source-h<=target-h<125:

    python3 -B 05-knowledge/results/infinite_birth_families_20261111/families.py --source-h 33 --target-h 101 --order 100

This returns 1510 and a primitive adjacency witness padded by 94 uniform
source vertices. The opposite ternary target 103 returns 758 at order 100.
A zero answer means no such one-vertex birth at that order; an absent
family means none at any order. It is not a prohibition of multi-step
constructions or realization of the target by a different operation.

## Evidence and limitations

The source-first census uses all 7,412 ordinary classes through order eight,
with 1,774,266 attachment orbits. There are 57,658 direct child-SCC controls,
3,402 independent literal child-H checks, and 15,624 fixed-H coefficient
checks. The child-first audit visits all 191,536 order-nine classes, counts
paths in every strong one, and independently computes the root orbits and
deletion counts of every retained core. It verifies all 872 primitive
kernel polynomials. A further 17,577 coefficients test the compiled
all-size formulas against the source-first census.

The common-source atlas receives 1,930 literal normalized or joint witness
checks, and its entire target-availability set agrees on all 1,454 eligible
source classes through eight vertices. Its scalar facets summarize common
existence; actual cuts and canonical core rows remain attached.

The 32 quadratic-ray bases receive 64 literal tail checks. The arithmetic
sieve agrees with primary OEIS files, trial-division inverse fibres, a
separate quadratic primality check, and 99 instances of the inverse-equality
bridge. Saved hashes use LF-normalized text and literal compressed bytes.
These are independent computational paths within this session, not an
external peer review or a claim of globally novel mathematics.

Full metagraph wiring, arbitrary primitive H kernels, all odd seeds, and
the required uniform Chan-Pak construction remain open. Source identities
must survive when composing the counted operations; scalar positivity
alone is insufficient.
