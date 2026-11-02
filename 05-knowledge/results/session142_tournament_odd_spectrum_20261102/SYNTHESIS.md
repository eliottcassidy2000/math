# What is closed, and what is not

## Established or independently rechecked

1. **PROVED, recovered exclusion:** no tournament has ordinary Hamiltonian
   path count7 or21. Cited Busch strong-minimum/overlap results reduce the
   exclusion to small orders; the selected272 overlap presentations use two
   exact counters. Independent full labelled enumeration through order6
   corroborates the small spectra. The all-size argument is the written
   proof, not an extrapolation from that enumeration.
2. **PROVED arithmetic refinement:** for the candidate monoid A of positive
   odd integers excluding7,21, the multiplicative atoms are odd primes
   p!=7, the numbers7p for odd primes p!=3, and63,343. This is a theorem
   about A, not a claim that A is already the attained spectrum.
3. **PROVED factor-length bound:** writing a=v7(h), the maximum number of
   nonunit A-factors is Omega(h)-ceil(a/2), except h=7*3^b, b>=2, where it
   is b-1; L_A(1)=0. Every actual tournament with H=h has at most L_A(h)
   non-singleton strong components. Uniform attainment is not asserted.
4. **FINITE-EXACT carry certificates:** strong tournament adjacency rows
   attain49,63,343, with independent DP and permutation counts. Combined
   with SCC multiplication these recover the exact two-prime-lane reduction.
5. **PROVED bridge obstruction:** H_n=1+2Q_n coefficientwise over the integers.
   A square integer-polynomial response map with k raw-H coordinates has
   Jacobian determinant divisible by2^k. Normalizing the count coordinate
   removes that parity factor but not the source-domain constraints.
6. **PROVED legal-domain obstruction:** weighted orientation interpolation
   fills the complete interval between minimum and maximum counts. Actual
   one-edge templates also give normalized responses2-t and9-t, whose
   derivatives are-1, yet forbidden targets3,10 have inverse t=-1 outside
   the permitted{0,1}. Even an integral proper global inverse need not be
   a legal tournament constructor.
7. **PROVED sufficient algebraic inverse, limited scope:** an integer
   polynomial map homogeneous for matching positive input/output weights
   and with Jacobian determinant+/-1 has an integer polynomial inverse,
   by block recursion in increasing weight. This supplies no Boolean or
   nonnegative inverse guarantee; the shear(x+y^2,y) is a hostile control.

## Explicit correction, preserved lineage

The main positive-weight theorem in the inherited Jacobian
`THM-1370-elliptic-jacobian-conjecture-all-dimensions.md` survives. Its broad
sentence excluding every identity-plus-ordinary-homogeneous map from positive
weighted homogeneity does not: the shear(x+y^2,y) with weights(2,1) is a
counterexample. Claims that failure requires an indefinite grading must
also be restricted to maps that admit the asserted grading. The new bridge
proof records these corrections; the canon file now contains the independently
reviewed narrow repair and explicit positive-target/origin argument. Its
historical source snapshot is unchanged.
This is separate from the other theorem sharing the number1370, whose full
slug is `THM-1370-h-spectrum-omits-7-21-all-n.md`.

## The remaining theorem

Construct a strong tournament with H=p for every odd prime p!=7 and one
with H=7p for every odd prime p!=3. Equivalently, supply a legal gadget
atlas with exact normalized response laws and a well-founded inverse descent
covering every q>=0 except3,10, terminating in certified base tournaments.
Lattice membership, rank, a continuous image, a determinant, local inverse,
finite count bank, or arithmetic factorization alone does not prove that.

Royle's reported bounded prefix is safely inclusive80555; the selected local
THM-4104 prefix is80405. These inherited prefixes are cited, not rerun here.
Later repository constructive prefixes, when used, retain their own source
certificates and do not promote the universal claim. No unrelated current
order-eleven claim, LRC(14), planar-Jacobian conjecture or earlier session's
approval is altered or promoted by this release.

## Checkpoint2: the new structural restrictions

[Prime-lane proof](minimal_gap/PROOF.md) proves that **every** actual tournament
with H=p prime, p>3, is asymmetric, including those padded by arbitrary
singleton strong components. The free automorphism action on Hamiltonian
vertex words gives |Aut(T)| dividing H; the cited strong minimum makes the
core smaller than p, excluding an order-p permutation.

For H=7p and prime p>=11, only group orders1 or7 remain. The complete small
cyclic-orbit tables plus an injective merge of disjoint path words give the
sharp first symmetric prime quotient: every7p witness for odd prime p<593
is asymmetric, while the actual ten-vertex rows
`[142,156,184,240,225,195,135,512,255,383]` have H=4151=7*593 and exactly
seven automorphisms. These are all-witness restrictions, not existence for
every prescribed p or7p.

[Constructor-language proof](constructions/BINARY_NO_GO.md) addresses the
actual inherited append/re-mark operations S0,R1 on positive marked states
(C,L,g,t). Both choices have the same next count2F+I, where F=2C+g+t and
I=2L+g-t. Eligibility for a binary step is exactly I=+/-1. An R cannot
start three binary steps at F>=5, and at depth k>=4 all first k-1 operations
must be S. Thus finite seed banks have only O(log X) reachable binary-growth
outputs below X, despite permitting arbitrary operation depth and state
choices. This is an all-depth proof, not a finite-search extrapolation.

The two literal H=5 seeds recover rays2^(k+2)+1 and6*2^k-1. Both rays are
inherited. The new conclusion is the obstruction for this exact operation
alphabet, including arbitrary named-neutral S1/reset interleavings: those
enter an invariant state set from which binary growth cannot resume. Other
neutral gadgets, extra vertices, nonuniform attachments,
old-edge edits, extra marks or older-vertex re-marking are outside it.
A useful next construction must retain a legal marked state while adding
an operation able to renew the missing I=+/-1 branching resource.

### Current finite-prefix provenance

The cited Royle/local order11 prefixes above are not the best later repository
construction bound. The inherited [October31 constructive packet](../square_sum_spectrum_overnight_20261031d/README.md)
records all allowed odd values through99,999,999. That100million scan is
**inherited and not rerun here**, not a Session142 numerical advance and not
an infinite completeness theorem. Source-preserving individual recovered
witnesses are a separate later checkpoint.

## Checkpoint3: size-dependent module amplification

[Path-cover proof](path_cover_amplification/PROOF.md) counts incidences of a
Hamiltonian path and k-1 cuts. Forgetting the order of the k resulting
nonempty directed paths loses at most k! choices, so
`k! p_k(T) >= H(T) binom(n-1,k-1)`. In the inherited exact homogeneous
substitution kernel the coefficient is k!W_k, not W_k alone. A cyclic port
has W_1>=3 and W_2>=1, giving `H(Q[T]) >= (n+2)H(T)` and an extra1 for even n.
The context may depend on T; a finite-template restriction is not needed.
Inactive singleton-SCC ports remain exact odd multiplication. This does not
apply to the nonhomogeneous marked S0/R1 operations or arbitrary multi-port
surgery, and does not prove any additional unattainable integer.

## Checkpoint3: recover actual counts, not a new numerical frontier

[Recovered witnesses](prime80557/README.md) preserve the prior executable
seven-parent constructor and22 source snapshots. Explicit adjacency matrices
have H80557 at order19 and H100000007 at order29. Both are strong and prime-
count, hence asymmetric. Neither is asserted smallest or newly discovered.
80557 was already covered by an earlier million-prefix construction;
100000007 is only an isolated positive control beyond the inherited100million
scan, not a new interval or a global realization theorem.

The small witness is checked by ordinary full subset DP at O0/O2. The larger
one has a transitive parent after deleting its apex, and is counted directly
from the adjacency by partitioning paths at that apex and summing over the
maximum-before/minimum-after pair. A separate two-block-side enumeration
also works. No order29 full subset-DP run is claimed.80498 inherited inverse
controls and598 small graph controls are rerun; the billion-pair finite
prefix scans are not. Source paths are pinned to a published commit.

These positive controls help distinguish a *failed restricted language* from
a missing actual count. The remaining obligation is still uniform legal
source-image coverage of the infinite prime and7p lanes, not another local
rank argument or a larger list of isolated witnesses.

## CITED incoming full-state inverse: a distinct solved subproblem

The concurrent [compressed legal descent](../spectrum_descent_jacobian_20261102/DESCENT.md),
read at published commit0040dca1285542b81b305d96aa97181dfe8d8376, provides a
minimum-cost complete inverse for its declared singleton/D **full-state**
grammar. Its scalar-target entry theorem remains OPEN. In its ordered
coordinates(F,q,g,h), positivity-admissible(9,7,3,5) is not a reachable state,
even though the scalar Hamiltonian count9 is attained. Thus selecting a legal
state above a scalar target remains a separate obligation.

The accompanying [Boolean realization handoff](../boolean_realization_71_17_20261102/HANDOFF.md)
correctly separates that Bellman search over legal predecessors from the
refuted local greedy energy rule. These incoming results are **CITED here**,
not newly replayed, independently promoted or claimed as Session142 discoveries.
They complement rather than duplicate the recovered fixed-source inverse.
