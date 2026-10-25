# Square-gap prime support: exact count, valuation filter, and descent

Status: PROVED by the self-contained arguments below; finite checks are separately FINITE-EXACT. No novelty claim. The one-prime obstruction is already Proposition 11 of Voronsbekher's 2026 preprint. This is not a proof of nonexistence of general magic hourglasses or magic squares.

## Inheritance and scope

The inherited session58 hourglass work in `/tmp/session58-hourglasses/check.py` retains the seven-cell collision guard and rejects the Cain three-AP-only example as a non-hourglass. The inherited mechanism is the centered three-gap Schur relation; the canonical hostile is a three-AP packet whose top-row sum is wrong. The least-used sidecar here is the Gaussian allocation index at each center prime. The corrected near miss is that counting three available gaps is necessary, not sufficient: center root 125 has three gaps but no relation.

Read-only repository guardrails were inspected. All new files are confined to this assigned directory; publication and repository coordination belong to the parent. Live concepts: centered Schur relations; Gaussian allocation; valuation minima; common-root scaling; squarefree support. None is identified with an LRC carrier.

## 1. Definitions and exact hourglass predicate

For a positive integer e let D_e be the positive integers d for which e²−d and e²+d are integer squares. Necessarily 0<d<e²: equality would force 2e² to be a square, impossible for e>0. Therefore allowing zero as a square adds no endpoint.

There is a bijection between D_e and unordered positive integer pairs {u,v}, u≠v, satisfying u²+v²=e², given by d=2uv. Indeed (u±v)²=e²±2uv. Conversely, from a²=e²−d and b²=e²+d with 0<a<b, parity gives u=(b+a)/2, v=(b−a)/2. The pair is unique because u²+v² and uv determine the two squares.

A seven-distinct-square hourglass with center e² exists if and only if D_e contains three distinct positive elements x,y,z with x+y=z. To see this, the three opposite pairs are e²±d_i. Each full outer row sums to 3e² precisely when a signed sum of the three d_i vanishes; with positive d_i this is equivalent to one being the sum of the others. Three distinct positive gaps make all seven entries distinct. Conversely a repeated absolute gap repeats an entry. Center e=0 is excluded throughout.

## 2. Exact hereditary scaling

For every positive divisor t of e,

    D_e ∩ t² Z = t² D_(e/t).

One direction is multiplication. For the other, t² divides e² and d, hence divides both endpoint squares. Prime by prime, t²|a² implies t|a, so division gives the smaller center. This works for arbitrary t, not merely primes absent from split support, and preserves positivity, distinctness, and Schur relations.

Write e=m n where m contains exactly the powers of 2 and primes congruent to 3 modulo 4 dividing e, and n contains the prime powers congruent to 1 modulo 4. Then D_e=m²D_n. In a representation u²+v²=e², every inert prime q^b||e divides both u,v to exponent at least b: otherwise remove their common power of q and reduce modulo q to contradict that −1 is a nonsquare. If 2|e then u²+v²≡0 (mod 4) forces u,v even; iterate. The converse is scaling. In particular D_e is empty if n=1.

## 3. Gaussian index parametrization and count

Write n=∏p_i^{a_i} with all p_i≡1 (mod 4), and choose Gaussian primes π_i with π_i conjugate(π_i)=p_i. Unique factorization parametrizes every norm-n² Gaussian integer, up to a unit, by

    z(j)=∏ π_i^(a_i+j_i) conjugate(π_i)^(a_i−j_i),
    −a_i≤j_i≤a_i.

The positive gap is |2 Re(z(j)) Im(z(j))|. Multiplication by a unit and conjugation leave it unchanged. Conjugation sends j to −j. The zero vector is exactly the zero-gap class: a real or purely imaginary Gaussian integer of norm n² is a unit times n, whose conjugate prime allocations are balanced. A diagonal representation would give n²=2u², impossible. Distinct nonzero vectors modulo j~−j give distinct gaps: equal product magnitude and fixed norm determine the absolute coordinate pair, and its sign/swap symmetries are exactly units and conjugation.

Consequently

    |D_e| = (∏(2a_i+1)−1)/2.

For one split prime to exponent a, there are exactly a gaps. For three distinct split primes to exponent one, there are (3³−1)/2=13 gaps, not three. Center root 1105=5·13·17 is an explicit control.

## 4. Exact valuation where the allocation is unbalanced

At a fixed p=p_i and a=a_i, if j_i≠0 then

    v_p(d(j))=2(a−|j_i|).

If j_i=0, only the inequality v_p(d(j))≥2a is universally valid.

Proof: factor p^(a−|j_i|) out of z. When j_i≠0, the residual Gaussian integer w contains exactly one of π_i and its conjugate, not both. If p divided Re(w), then w+conjugate(w) would vanish modulo π_i while exactly one summand does, a contradiction. The same argument with w−conjugate(w) treats Im(w). Denominators 2 and 2i are units modulo the odd prime π_i. Thus neither residual coordinate is divisible by p, proving the formula. In the balanced case the factor p^a alone gives the stated lower bound.

Hostile to strengthening the balanced case to equality: e=65 and d=3000 have endpoints 35² and 85². The underlying pair is (60,25), the p=5 allocation is balanced, but v_5(3000)=3>2v_5(65)=2.

## 5. One-prime obstruction and the general tropical filter

If n=p^a, positive gaps correspond to |j|=1,...,a and have the distinct valuations 0,2,...,2a−2. A signed sum of three distinct gaps cannot vanish: its uniquely least p-adic valuation cannot cancel. This excludes every e=m p^a after arbitrary inert/2 scaling, not only pure prime powers. For a<3 cardinality already excludes a relation; a≥3 is where valuation adds information.

More generally, if three distinct Gaussian-index gaps form a signed relation, then at every split prime p_i,

    max(|j_i^(1)|, |j_i^(2)|, |j_i^(3)|)

is attained at least twice. If the maximum is positive and unique, the exact valuation formula makes the corresponding gap the unique minimum. A balanced index has valuation at least 2a and therefore cannot interfere. If the maximum is zero the condition is automatic. This is a necessary filter, NOT a sufficient relation certificate.

For center 65 all four unordered triples of gaps pass this filter, but none is a Schur triple. For center 1105, 150 of 286 triples pass, but none is a Schur triple. Thus the filter is not a solution even at two-prime support.

## 6. Primitive support reduction and incidence consequence

For a triple of gaps at stripped center n, put M_i=max_r |j_i^(r)| and t=∏p_i^(a_i−M_i). Then all three gaps are divisible by t², so the hereditary identity transports the triple to center n/t=∏p_i^M_i. The new indices are the same j_i; their ranges shrink to [−M_i,M_i]. Every prime with M_i=0 disappears entirely. Any Schur relation and all seven distinctness/positivity guards survive.

This t is exactly the common divisor of n and the six noncentral positive square roots. At a prime with M_i>0, some gap has exact valuation 2(a_i−M_i), preventing further simultaneous root division. At a prime with M_i=0, all available center powers divide every root. This also follows directly from the hereditary identity and the valuation formula.

After dividing by the common divisor of all seven roots, every surviving center prime has M_i=a_i. The maximum-twice condition then forces at least two of the three gaps to be p_i-units. Thus, in any hypothetical primitive hourglass, each center prime divides at most one of the three positive gaps:

    gcd(e, gcd(d_r,d_s))=1 for every r≠s.

The center has no 2 or inert prime, and cannot have only one split prime. These are necessary restrictions, not an existence construction. No claim here excludes arbitrary two-prime-power centers.

## 7. Independent exact checks

Run `python check.py` and `python -O check.py` in this directory. No validation depends on Python assertions. The independent endpoint-square enumeration and Gaussian-index generation agree for every e=1,...,3000, covering 3,172 positive gaps. The exact count, stripping identity, every local valuation formula, and hereditary slice for every prime divisor of every tested e all pass. There are 1,487 tested centers with exactly one split prime; their valuation lists are exactly the predicted sequences. No Schur triple occurs in this bounded universe. This is only a bounded check, not an all-center nonexistence result.

Files: check.py, check.out, check-optimized.out, certificate.json. Arithmetic uses only integers and integer square roots.

## 8. Literature provenance and limits

Primary source: Valentyn Voronsbekher, *On magic hourglasses of squares having no centre of prime-power or semiprime type*, Zenodo preprint, 13 August 2026, DOI https://doi.org/10.5281/zenodo.21910218; record https://zenodo.org/records/21910218.

The primary metadata and PDF were retrieved directly through the public Zenodo API after the web reader failed. The source PDF and extracted text were consulted locally but are deliberately excluded from this publication; retrieve the primary source through the DOI or record above. The PDF's Lemma 7 contains inert/2 stripping, Lemma 10 Gaussian prime allocations, and Proposition 11 the prime-power obstruction. It also states and proves a squarefree-semiprime exclusion; that separate proof is not independently audited or used here. The displayed positive-gap formulas in its Lemma 10 should be read with absolute values when coordinate products are negative; this does not affect its valuation argument. This note uses absolute values throughout.

The author-research page https://github.com/uwe-schwarz/magic-squares-of-squares/blob/main/research/prime-support/hourglass-sunit.md was used to locate the primary source, not to establish novelty or inherit an unverified theorem. In particular the valuation exclusion is a recovered known mechanism. The preprint is not represented here as refereed literature, and no global open-problem status is inferred from a local search.

## 9. The squarefree projective incidence lens

Assume the stripped center e is the product of k≥1 distinct split primes. Its Gaussian allocation vectors belong to {−1,0,1}^k minus zero. Reducing coordinates modulo 3 identifies this set with F_3^k minus zero; simultaneous sign identifies the two elements of each one-dimensional F_3-subspace. Thus the gap indexing is in bijection with PG(k−1,3). For k=0 the set is empty rather than a nonempty projective space. Choices of prime order or conjugate orientation induce coordinate permutations or sign changes, hence preserve projective incidence. Once those choices are fixed the correspondence is explicit; it is not an assertion that all projective symmetries preserve numerical gap values.

This gives the exact point counts 1,4,13 for k=1,2,3. For k=3, the 13 index points carry the 13 lines of PG(2,3), each containing four points, and each pair belongs to one line. These lines are index-incidence objects, not additive gap equations.

There is an arithmetic interpretation of coordinate subspaces. For a subset S of prime coordinates, j_p=0 for every p in S if and only if ∏_(p in S) p² divides the gap. The forward implication follows from the balanced-index valuation inequality; the reverse uses v_p(d)=0 when j_p=±1. Put t=∏_(p in S)p. The coordinate subspace maps under d↦d/t² bijectively onto D_(e/t), preserving the remaining indices and every actual Schur relation. In projective terms this is the natural deletion isomorphism from the coordinate subspace PG(k−|S|−1,3). When all coordinates are deleted, both sides are empty. Arbitrary projective hyperplanes do not receive this arithmetic interpretation from the present proof.

A cheapest hostile to conflating projective dependence with an additive relation occurs already at e=65. The index vectors (1,0),(0,1),(1,1) satisfy the displayed F_3 linear dependence, but their corresponding gaps are respectively 4056,3000,3696. None is the sum of the other two. Indeed the entire four-point projective line D_65={2016,3000,3696,4056} contains no Schur triple. At e=1105 the coordinate line j_17=0 is exactly

    17² D_65 = {582624,867000,1068144,1172184},

again with no Schur triple. The independent test enumerates all 13 projective lines at center 1105 and verifies four-point sizes and unique pair incidence; no line contains a Schur triple. This refutes line/dependence-implies-Schur, not the converse implication, whose global truth is not examined here.

Loss ledger: source is the Gaussian allocation set with its exact integer exponents, target the projective F_3 index set; the map preserves index identity modulo simultaneous conjugation and all coordinate-zero incidence, but finite-field addition is not Gaussian multiplication or numerical gap addition. In normalized Gaussian phase coordinates z(j)/e=∏(π_p/conjugate(π_p))^j_p, exponent addition takes place over the integers, not modulo 3. Wrapping j+j' modulo 3 discards a generally nontrivial phase factor; absolute coordinate products discard further phase information. Exact gap values are the sidecar needed to test the actual additive relation. The count 13 and line size 4 therefore have a precise, sharply limited indexing meaning, not numerological force.
