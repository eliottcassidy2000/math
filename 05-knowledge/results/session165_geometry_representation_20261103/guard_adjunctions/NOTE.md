# Guard adjunctions, dihedral chambers, and actual blue tilings

**Status: PROVED general identities and constructions; FINITE-EXACT controls.**
The new object is a precise connection between arithmetic guard observations,
rank-two reflection arrangements, and a complement-stable bank of genuinely
blue fixed-path tilings. It is not a FRACTRAN simulator, a theorem about the
full Hamilton-path spectrum, or a claim of literature priority.

## Inheritance and the missing coordinate

[S160's priority scanner](https://github.com/eliottcassidy2000/math/blob/17b4acad2ee77e9bcec299f9fcba816a7b2a0cf7/05-knowledge/results/session160_ordered_programs_20261103/priority_tilings/NOTE.md)
already proves the lcm feasibility test, the 2,592/512 counts for its thirteen
guards, and an observation collision. We retain those facts as inherited.
[S151's incidence](https://github.com/eliottcassidy2000/math/blob/main/05-knowledge/results/session151_spectrum_metagraph_20261103/incidence/INCIDENCE.md)
already supplies path-orbit and double-coset counts: these are not new here.
[S151's line convention](https://github.com/eliottcassidy2000/math/blob/main/05-knowledge/results/session151_spectrum_metagraph_20261103/metagraph/STRUCTURE.md)
distinguishes chord complement F from reflection S and requires S-fixed
presentations for blue lines. The construction below explicitly enforces this.

The closest mechanism is the lcm criterion; the hostile is equal observations
with different future observations. The corrected near miss is treating
complementation as arithmetic dynamics. The underused sidecar is the actual
prime-exponent vector, rather than only the currently enabled rules.

## 1. An exact adjunction, and its maximal complement-stable part

Let b_1,...,b_g be the denominators of any finite reduced FRACTRAN program.
The numerators and priority order play no role in the following static object.
For S contained in [g] and N a positive integer define

    lambda(S) = lcm{b_i : i in S},       lambda(empty)=1,
    rho(N) = {i : b_i divides N}.

View subsets as a category ordered by inclusion, and positive integers as a
category ordered by divisibility. Then lambda is left adjoint to rho:

    lambda(S) divides N  iff  S is contained in rho(N).       (1)

This is immediate from the definition of lcm. Hence c=rho lambda is a closure
operator. More explicitly it is extensive, monotone, and idempotent, since
lambda(c(S))=lambda(S). The realizable observations are exactly

    C = {S : c(S)=S}.                                        (2)

Indeed rho(N) is closed and lambda(S) itself realizes every closed S. Thus
lambda rho chooses a canonical least integer with the same observation;
that integer divides the original N. This is an exact quotient of observations,
not a quotient of ordered-program dynamics.

Put

    B = {S in C : [g] minus S also belongs to C}.              (3)

B is the largest complement-invariant subset of C. Its members will be called
**biclosed guard sets**, with respect to this particular lcm closure. If any
b_i=1 then B is empty. If B is nonempty and g>0, complementation acts freely.
In general B is not a Boolean algebra. For denominators (4,6,9), the biclosed
sets {4,6} and {6,9} have intersection {6}, which is closed but not biclosed:
its complement {4,9} forces 6 by lcm. Denominators (6,10,15) give only empty
and full as biclosed sets, despite all three denominators being pairwise
incomparable. Pairwise divisibility alone does not capture lcm entailment.

## 2. Two-prime antichains are exactly dihedral chamber banks

Fix distinct primes p,q, m>=2, and integer exponents

    a_1 > a_2 > ... > a_m >= 0,
    0 <= c_1 < c_2 < ... < c_m,
    b_i = p^(a_i) q^(c_i).

None of these denominators is 1. For a nonempty index set S with least and
greatest indices l,r,

    lambda(S) = p^(a_l) q^(c_r),
    c(S) = {l,l+1,...,r}.                                   (4)

The first identity follows from coordinatewise maxima. An index i divides
this lcm precisely when a_i<=a_l and c_i<=c_r, namely l<=i<=r.
The empty set is closed. Therefore

    |C| = 1 + m(m+1)/2,                                     (5)

and its members are exactly empty and the intervals. An interval and its
complement are both intervals exactly when it is a prefix or suffix, including
empty and full. Consequently

    |B| = 2m.                                               (6)

Join two B-states when they differ in exactly one guard bit. The result is the
cycle C_(2m), with the following order, where [j]={1,...,j}:

    empty, [1], [2], ..., [m],
    [m] minus [1], ..., [m] minus [m-1], empty.                (7)

There are no extra edges: an interior prefix has only its next larger and
smaller prefixes as one-bit neighbors, and similarly for suffixes. Complement
is the antipodal map j -> j+m on the cycle labels 0,...,2m-1.

This is also an exact rank-two root-geometry model. Set

    alpha_i = (cos((i-1) pi/m), sin((i-1) pi/m)).

The roots lie in an angular interval shorter than pi. The positive cone of a
nonempty collection contains exactly the roots between its extreme indices.
Thus b_i -> alpha_i intertwines the lcm closure with positive-cone closure.
It identifies B with the biclosed subsets of the positive roots of the real
dihedral root system I_2(m). Equivalently these are the negative-root sets of
the 2m chambers cut out by the lines alpha_i perpendicular: a generic linear
functional can change sign at most once along this angular interval. The
chamber adjacency graph is exactly (7), and reversing the functional is
complementation. This establishes the identification directly, without an
appeal to a classification theorem. The standard root-system terminology and
relation between finite biclosed sets and inversion sets are recorded in
[Barkley--Speyer](https://arxiv.org/abs/2207.05998) and
[Speyer's Coxeter notes](https://dept.math.lsa.umich.edu/~speyer/665_2019/BruhatNotes.pdf).

For (4,6,9), this is the A_2 hexagon:

    000 -> 100 -> 110 -> 111 -> 011 -> 001 -> 000.

The omitted observations are 010, whose complement is not realizable, and
101, which is itself not realizable. Their least-integer representatives are
1,4,12,36,18,9 in the displayed cycle order. Those integers are not an orbit of
a FRACTRAN program merely because their observation graph is a cycle.

## 3. Prime separation gives tensor products

Suppose the guards are partitioned into blocks, and no prime occurs in the
denominators of more than one block. Observation factors over the associated
prime parts of N. Conversely independent realizing integers can be multiplied.
Therefore

    C = product C_j,          B = product B_j.               (8)

If every block is a two-prime antichain of size m_j>=2, the one-bit-flip graph
is the Cartesian product

    Gamma = C_(2m_1) square ... square C_(2m_s),              (9)

and simultaneous guard complementation is the product of the antipodes.
This realizes a product of rank-two reflection-arrangement chamber banks.
It does not assert a root-system model for arbitrary overlapping denominator
families.

There is an elementary representation-theoretic refinement. Let V be the
complex vector space of functions on Gamma, A its adjacency operator, and F
the complement operator. For k_j in {0,...,2m_j-1}, tensor Fourier characters
are simultaneous eigenvectors, with

    A:  sum_j 2 cos(pi k_j/m_j),
    F:  (-1)^(sum_j k_j).                                   (10)

Thus V is the tensor product of the regular cyclic permutation modules,
and the F-even and F-odd subspaces each have half the dimension. The cyclic
action here rotates the enumerated observation cycle; it is not claimed to
act on integers by rational multiplication or on unmarked tournaments by
isomorphisms.

For one block an exact character/closed-walk identity is

    tr(F A^ell)
      = sum_(k=0)^(2m-1) (-1)^k (2 cos(pi k/m))^ell
      = 2m sum_(r=0)^ell binom(ell,r)
                     1[2r-ell = m mod 2m].                  (11)

The trace counts length-ell walks from a state to its complement. The second
formula counts choices of clockwise versus counterclockwise steps. In
particular it is zero for ell<m and equals 4m at ell=m. In the product bank,
put D=sum_j m_j. The distance from any state to its complement is D, and the
number of shortest such walks from each state is

    2^s D! / product_j(m_j!).                               (12)

There are two directions around each cycle, and the multinomial coefficient
interleaves the steps in the different coordinates. For (4,6,9,5,7), Gamma is
C_6 square C_4 with 24 states and 40 shortest antipodal walks from each state.

The old S160 denominator list factors into four comparable two-guard blocks
and five one-guard blocks. Formula (8) recovers its inherited counts
3^4 2^5=2592 and 2^4 2^5=512. These comparable blocks are not the antichain
blocks in (9): their biclosed sets 00 and 11 have no one-bit edge. We do not
apply the cycle-graph conclusion to them.

## 4. A literal embedding into blue marked-path tilings

Assume g≥1 in this section. On vertices 0,...,n-1 keep the path 0->1->...->n-1 forward. A tiling is its
bit vector on the nonconsecutive pairs; a 1 reverses that chord from the
transitive orientation. Let S permute chord positions by

    (i,j) -> (n-1-j,n-1-i),

and let F complement every chord bit. In the S151 convention a line
{t,F(t)} is blue exactly when S(t)=t.

The chord positions have d=floor((n-1)^2/4) S-orbits. Choose n with d>=g and
assign every such orbit to one of the g guards, using every guard at least
once. Define t(e) to put the enabled bit of that guard on the entire orbit.
Then

    S(t(e))=t(e),
    F(t(e))=t(complement e),                                 (13)

and e -> t(e) is injective. Hence every nonempty biclosed bank B embeds into
actual blue tilings, and gives exactly |B|/2 distinct individual blue lines.
One may take n=2 ceil(sqrt(g))+1. The marked path, the assigned chord-orbit
labels, and the line identities must be retained. Isomorphism classes can
identify states; they do not identify distinct individual lines in the
specified line multigraph.

**Do not confuse two graphs.** The chamber graph Gamma joins states differing
in one guard. The actual blue metagraph lines pair complementary states,
which are antipodes of Gamma. A chamber edge is not claimed to be a black
line, a blue line, or a FRACTRAN transition.

For g=3 use n=5. In lexicographic orbit order the four reflection orbits are

    {(0,2),(2,4)}, {(0,3),(1,4)}, {(0,4)}, {(1,3)}.

Assign them guards 1,2,3,1. Independent Hamilton-path DP and permutation
counting give:

| Enabled bits | Least integer | H(t(e)) | Canonical tournament code |
|---|---:|---:|---:|
| 000 | 1 | 1 | 0 |
| 100 | 4 | 11 | 12 |
| 110 | 12 | 11 | 12 |
| 111 | 36 | 9 | 10 |
| 011 | 18 | 15 | 76 |
| 001 | 9 | 9 | 8 |

The canonical code is the minimum reversed-edge bitmask over all vertex
permutations, with pairs in lexicographic order. The three blue lines are
000--111, 100--011, and 110--001. The arithmetic states 100 and 110 already
become the same unmarked tournament class, so the chamber model cannot be
recovered from that quotient. H also identifies additional states.

## 5. The exact dynamical obstruction survives the geometry

Take the ordered reduced FRACTRAN program (1/4,1/6,1/9). Integers 4 and 16 have
the same biclosed observation 100. Their first applicable rule is 1/4, but
their successors 1 and 4 have observations 000 and 100. Both remain in B.
Therefore even this dihedral blue bank has no well-defined deterministic
transition induced by the program.

The failure cannot be repaired by a finite enrichment that tracks this
program exactly on all powers of four and still outputs the guard bits.
Starting at 4^k produces k successive 100 observations, then halts at
observation 000. Equivalently, make halt absorbing and repeat 000 forever. Different k have different future observation streams. In any
deterministic quotient, two inputs in the same state must have identical
future output streams. Thus infinitely many quotient states are necessary.
This is a direct distinguishable-futures proof, not an appeal to FRACTRAN
universality. A bounded horizon admits finite quotients, and a retained
unbounded exponent register restores the distinction.

## 6. Loss ledger and reproducibility

- **Source:** positive integers with explicit denominators. **Target:** closed
  guard sets. **Map:** rho. **Preserved:** every current divisibility guard and
  therefore current first-applicable index. **Lost:** exponents above the
  thresholds, foreign primes, and generally next-state information. **Sidecar:**
  sufficiently detailed exponent registers or N itself. **Hostile:** 4 versus 16.
- **Source:** a two-prime antichain's closed sets. **Target:** rank-two root
  closure. **Map:** b_i -> alpha_i. **Preserved:** closure, biclosure, Hamming
  adjacency, complement. **Lost:** primes, exponent sizes, fractions, priority.
- **Source:** biclosed guards. **Target:** blue marked-path tilings. **Map:**
  orbit-constant chord assignment. **Preserved:** every guard, complement and
  individual line identity. **Lost upon further quotient:** station semantics,
  and possibly distinct guard states, as 100/110 demonstrate.

Run `python3 check.py` and `python3 -O check.py` in this directory; the outputs
are byte-identical. The test computes closure directly by divisibility and
checks complete arithmetic residue periods when the lcm is below 10,000.
It checks every guard subset for rank-two modules m=2,...,8, every induced
cycle edge and antipode, and the integer form of (11) against independent
walk DP through ell=2m+1. It checks all 24 product states and their 40 shortest
antipodal walks. The literal n=5 blue embedding is independently checked for
path validity, reflection, complementation, isomorphism collisions and path
counts. These finite controls corroborate the general proofs above; they are
not a census of arbitrary FRACTRAN programs or tournaments.

## Universal extension

[UNIVERSAL_CODES.md](UNIVERSAL_CODES.md) proves that every finite closure system
has an exact denominator realization, with minimum prime count equal to the
width of its proper meet-irreducible poset. It also constructs type-A banks
whose complement-realizable observations form the permutohedron, retaining
the distinction between guard adjacency and actual blue lines.
