# Split-sum prime generation: the quadratic trap and exact binary balance

19 October 2026. **PROVED:** the reductions and elementary classification
below. **CITED:** Booker–Pomerance's complete squarefree-smooth residue
theorem and the repository's circulant-tournament theorem. **FINITE-EXACT:**
the companion controls and explicit sources. **OPEN:** completeness of both
the original Euclid–Mullin sequence and the greedy split-sum sequence.
No literature-priority claim is made.

## 1. The useful advance and the inherited boundary

There is a literal tournament/prime-generation bridge here, but it belongs
to the **split-sum** variant of Euclid's construction. Once the smallest
missing prime is fixed, failure to produce it is exactly confinement of
the subsequent prime outputs to a quadratic-residue subgroup. At primes
congruent to 3 modulo 4, that subgroup is the out-neighborhood of zero in
the Paley tournament. A single output outside that neighborhood forces the
missing prime at the next step.

A second result classifies exact multiplicity balance of divisor residues.
The same cyclic 2-group condition that characterizes the rigid-circulant
primes characterizes the possibility of a perfectly uniform binary
subset-product distribution. Yet a minimum-size balanced source cannot
generate its Fermat prime: it has the wrong quadratic phase. One extra
nonresidue factor is necessary and sufficient in the extremal case.

The main hostile is especially concrete. At the Fermat prime 65537, the
divisors of the primorial below 65537 have full residue support, by the
published theorem, but cannot be exactly balanced in signed residue
classes. The missing binary layer consists of the two composite residues
256 and 65281 = 97·673. The smallest prime repairing that layer is 327941.
Full support, exact balance, and the target's quadratic phase are three
different data.

Inheritance used:

- `number_tournament_bridges_20261014.md` and
  `THM-4644-fermat-primes-rigid-circulant-tournaments-wieferich-orbit-explosion.md`:
  cyclic 2-groups explain Fermat-prime circulant rigidity.
- `tournament_path_defects_20261014e.md`: a binary address is useful only
  when its decoder and its actual arithmetic predicate remain specified.
- `prime_pairs_local_symmetry_20261019/forbidden-pattern-classification.md`:
  correlation reconstruction has an odd-order hypothesis and loses
  information at even-order boundaries.
- `scaling_residues_correlations_20261019/SCALING_RESIDUES_AND_CORRELATIONS.md`:
  a complete-looking normalized kernel need not retain the required phase.
- `wilson_defect_trees_20261014c.md`: a freely selected branch is not the
  supplied unchanged source. Here branching existence is not greedy
  completeness. No Wilson-prime finiteness premise is used.

## 2. A lossless candidate-prime test

Let P be a finite set of distinct primes and N = ∏_{p∈P} p. A split source
is an unordered pair {d,N/d}, with d dividing N. Each split sum d+N/d is
coprime to N: a prime in P divides exactly one summand. Thus every prime
factor of a split sum is new to P.

For an odd prime q not in P, put

    H_q = F_q^*/{±1},          h = |H_q| = (q−1)/2,
    D_{P,q}([x]) = #{d | N : d ≡ x or −x (mod q)},
    Z_q(P) = #{d | N : q divides d+N/d}.

D counts actual divisors with their multiplicities, not merely the
attained residue classes. Every divisor is counted once in H_q. Consequently
∑_x D_{P,q}([x]) = 2^{|P|}.

### Theorem 1. Exact phase-and-root-class formula. PROVED.

If (−N/q) = −1, then Z_q(P) = 0. If (−N/q) = 1 and s² ≡ −N modulo q,
then

    Z_q(P) = D_{P,q}([s]).                              (1)

For N>1, Z_q(P) is even, and Z_q(P)/2 counts unordered split sources.
In particular, if D_{P,q} is positive on every signed residue class, then

    q divides some split sum  iff  (−N/q) = 1.          (2)

Proof. Since d is invertible modulo q,
q | d+N/d iff d² ≡ −N. A nonzero square has exactly the roots s and −s,
which together form one class in H_q. This proves (1) and (2).
Complementation d↔N/d preserves divisibility of the split sum. It has no
fixed divisor when N>1 is squarefree, proving the evenness claim. ∎

The quadratic character alone is insufficient without the support
hypothesis: with P={2} and q=11, −2 is a square modulo 11, but neither
divisor 1 nor 2 belongs to the root class {3,8}. There is no candidate 11.

### Fourier and the retained marker

The characters of H_q are precisely the multiplicative characters χ of
F_q^* satisfying χ(−1)=1. Finite character orthogonality gives the exact
identity

    D_{P,q}([x]) = (1/h) ∑_{χ(−1)=1}
                    conjugate(χ(x)) ∏_{p∈P}(1+χ(p)).   (3)

The product is the binary subset-product transform. Formula (3) measures
the divisor distribution. Formula (1) also needs N modulo q to identify
whether the root exists and which class must be evaluated. Discarding
that marker loses the prime-generation predicate even if every value of
D is retained. The [kernel-recovery companion](KERNEL_PHASE_RECOVERY.md) refines this loss statement: reflection centers recover the phase in the stated nonsaturated Fermat regimes, while root-class ambiguity can remain.

Two exact hostile pairs make this loss visible:

| q | First prime source | Second prime source | Common complete D | Z values |
|---|---|---|---|---|
| 5 | {2,3,7}, N=42 | {2,3,11}, N=66 | (4,4) | 0 and 4 |
| 7 | {2,3,5,11}, N=330 | {2,3,5,17}, N=510 | (6,5,5) | 0 and 6 |

Classes in the last column of data are ordered by representatives
1,…,(q−1)/2. In both rows q is the smallest missing prime for both sources;
the source sizes and complete signed-divisor kernels agree. In the second
row even the reason for equality is immediate: 11 ≡ −17 modulo 7.
These are genuine arithmetic sources, not arbitrary artificial arrays.

## 3. The smallest-missing-prime reduction

**Published input.** Booker and Pomerance, *Squarefree smooth numbers and
Euclidean prime generators*, Theorem 1, prove that for prime q other than
5 and 7 every nonzero residue is represented by a squarefree integer whose
prime factors are less than q. Their Section 4 records that at q=5 and 7
the only missing nonzero residue is 4. Primary text:
https://arxiv.org/html/1607.01557v4 .

For the two exceptions, passing to signed residue classes removes the
exception: −4 is 1 modulo 5 and 3 modulo 7, and both are attained.
Thus, for **every odd prime q**, the divisors of the product of all primes
less than q meet every element of H_q. This is a direct corollary of the
published theorem, not a new proof of that theorem.

Define the greedy split-sum sequence by adjoining the least prime dividing
at least one of the numbers d+N/d, d|N. This is Chua's sequence, discussed
by Booker in https://arxiv.org/html/1605.08929 . It is not the original
Euclid–Mullin rule, which only uses N+1.

**Attribution.** The eligibility test and nonresidue phase-flip mechanism
already occur in the proof of Booker–Pomerance's Corollary 8; their
Section 6 also removes the 5,7 exceptions using signed representatives.
The next theorem and corollary are the greedy-process consequences and
Paley formulation of those published ingredients. The later binary-balance
classification and explicit obstructions are separately proved here.

### Theorem 2. Quadratic trapping is the exact omission obstruction. PROVED.

Let q be the smallest prime absent from the current finite source P.

1. If q=2, the next selected prime is 2.
2. If q is odd, q is the next selected prime iff (−N/q)=1.
3. Suppose q is odd and currently blocked, so (−N/q)=−1. If the next
   selected prime is r, then q is selected on the immediately following
   step iff (r/q)=−1.
4. Consequently, if q remains omitted forever from this stage, all
   subsequent selected primes lie in the quadratic-residue subgroup
   modulo q. Conversely, a forever-residue tail from a blocked stage
   keeps q omitted.

Proof. If 2 is absent, N is odd and N+1 is even. For odd q, all primes
less than q belong to P, so the published input and the preceding signed
exception check imply full signed divisor support. Theorem 1 gives
candidate eligibility. Every prime below q already divides N and cannot
divide a split sum, so eligibility is equivalent to selection by the
greedy rule. If r≠q is selected, q remains the smallest missing prime,
and full support persists because the old divisors remain divisors of Nr.
Finally

    (−Nr/q) = (−N/q)(r/q) = −(r/q).

Apply the same selection criterion at the next step. This proves all four
claims. ∎

In particular, **a failure of completeness cannot omit just finitely many
primes**. If q is the least omitted prime, all but finitely many prime
nonresidues modulo q must also be omitted. There are infinitely many such
prime nonresidues by the elementary argument in Corollary 3 below. This
consequence is specific to the proved split-sum reduction; it is not
asserted for the original Euclid–Mullin sequence.

No factorization of an enormous split-sum product is needed to test the
smallest candidate after its source-support theorem has been established.
This does not determine which larger prime is selected when q is blocked.
That selection remains the difficult dynamical input.

### Corollary 3. An equivalent, sharply specified proof target. PROVED.

For the greedy split-sum sequence from any finite prime seed, the following
are equivalent:

    (a) Every prime eventually occurs.
    (b) For every odd prime q, infinitely many selected primes r satisfy
        (r/q)=−1.
    (c) No odd prime q has an eventual selected-prime tail entirely in
        its quadratic-residue subgroup.

Proof. If a prime is omitted, choose the smallest omitted q. The prime 2
cannot be omitted, and all primes smaller than q have occurred after a
finite time. Theorem 2 then violates (b) and (c). Conversely, there are
infinitely many prime nonresidues modulo each odd q, by an elementary
Euclidean argument: if only finitely many were present, CRT would give
a positive nonresidue modulo q coprime to their product. Its prime
factorization must contain another nonresidue prime. Therefore a sequence
containing every prime satisfies (b). Equivalence of (b) and (c) follows
because selected primes are distinct, so r=q occurs at most once. ∎

This is an equivalence for this particular greedy process, not a proof of
any of the three statements. It requires neither a limiting density nor
equidistribution among all q−1 nonzero residues. The next exact obligation
is to rule out eventual confinement to a single index-two subgroup.

### The literal Paley bridge

For q≡3 modulo 4, orient a→b exactly when (b−a/q)=1. This is the Paley
tournament: negating a nonzero difference reverses the sign. Its
out-neighbors of 0 are exactly the quadratic residues. Theorem 2 says
that omission of the smallest missing q forces the residue labels of the
entire subsequent prime tail into this one out-neighborhood. A first
escape forces q at the next step.

For q≡1 modulo 4 the same relation is a Paley graph, not a tournament;
the corresponding statement uses its neighborhood. No ties are oriented
artificially, and no Hamiltonian-path count is asserted to govern the
process. The target being excluded is q itself, the zero residue, which
is kept outside the graph-neighborhood condition until it is selected.

Example: after the first three greedy outputs 2,3,5, the missing prime
7 is blocked. The next outputs are 11,37,13, with residues 4,2,6 modulo 7.
The first two are out-neighbors of 0; 6 is an in-neighbor. The next output
is therefore 7. The checker obtains these terms from the actual split rule.

**Original-Mullin boundary.** For the N+1 rule, even infinitely many
nonresidue multipliers do not by themselves force the running product to
hit −1. Modulo 7, alternating multipliers 3 and 5 are both nonresidues,
while their cumulative products alternate 3,1 and avoid 6. This is a
residue-level countermodel, not a claimed Euclid–Mullin trajectory. It
identifies the precise implication that fails without split-sum saturation.
Booker–Simon (2026), https://arxiv.org/html/2601.21901v1 , still formulates
the original least-factor question as open; no resolution is claimed here.

## 4. Exact binary balance singles out the same Fermat primes

Consider any list t_1,…,t_k in a cyclic group Z/hZ, and let C(x) count
the binary subsets I with ∑_{i∈I}t_i=x. Repetition in the list is allowed;
addresses are distinct even if their group values coincide.

### Theorem 4. Complete cyclic binary-balance classification. PROVED.

C is constant iff both of the following hold:

1. h is a power of two, say h=2^s, including h=1;
2. for every v=0,…,s−1, at least one t_i has 2-adic valuation v
   modulo 2^s. A zero residue supplies none of these layers.

When these conditions hold, C(x)=2^k/h for all x, and k≥s.

Proof. A uniform integer count makes h divide 2^k, hence h is a power
of two. The Fourier transform is

    F(j)=∏_i (1+exp(2πi j t_i/2^s)).

Uniformity is equivalent to F(j)=0 for every j≠0. Write j=2^a u with
u odd. A factor vanishes iff jt_i≡2^{s−1} modulo 2^s, which is equivalent
to v_2(t_i)=s−a−1. Thus all characters with a fixed a vanish precisely
when that layer is represented, and every nontrivial character vanishes
precisely under condition 2. The constant value and minimum length follow
from the total 2^k. For h=1 the layer condition is empty. ∎

Apply this to H_q using a primitive root g modulo q and
t_p=log_g(p) modulo h. Formula (3) is exactly the Fourier transform in
the proof. We obtain two statements with different quantifiers:

- At the abstract residue-list level, some exactly balanced binary
  distribution exists on H_q iff q−1 is a power of two, equivalently q
  is a Fermat prime.
- For a **supplied squarefree prime source N**, exact balance additionally
  requires its prime factors to meet every layer in Theorem 4.

Changing primitive root multiplies all discrete logarithms by an odd unit,
so these layer conditions are independent of that gauge.

The first equivalence agrees exactly with the repository's established
classification of odd primes for which every circulant tournament has
only translation automorphisms. That tournament statement is inherited
from THM-4644, via Akin/Burnside/Alspach, not proved anew here. The common
carrier is the cyclic unit group with no odd-order factor. This does not
identify a tournament path with a divisor or force a prime value.

### Theorem 5. Balance plus generation costs one additional binary choice.
PROVED.

Let q>3 be a Fermat prime and h=(q−1)/2=2^s. Suppose the actual divisors
of squarefree N with k distinct prime factors are exactly balanced on H_q.
Then

    q divides a split sum iff an even number of the prime factors
    of N are quadratic nonresidues modulo q.             (4)

If q divides a split sum, necessarily k≥s+1. Equality in this
generating-source bound occurs precisely when layer 0 occurs twice
and each layer 1,…,s−1 occurs once.
In this extremal case there are exactly two ordered divisor witnesses,
forming one unordered split source.

Proof. Since q≡1 modulo 4, (−N/q)=(N/q). The latter is (−1)^{m_0},
where m_0 is the number of prime factors with odd discrete logarithm,
namely those in layer 0. Apply Theorem 1 and exact balance. Every layer
is present, so m_0≥1; an even m_0 requires m_0≥2, giving k≥s+1. At
equality the listed layer multiplicities are forced and sufficient.
Then D([x])=2^{s+1}/2^s=2, and complementation pairs the two witnesses. ∎

In particular a minimal balanced list of s prime factors has one factor
in each layer and cannot generate q. This is an exact obstruction, not
a distributional heuristic. It is also visible from complement parity:
one divisor per signed class would give an impossible odd number of
witnesses if the candidate root class existed.

The boundary q=3 has h=1, so every divisor distribution is balanced.
Generation instead requires (N/3)=−1; the one-factor source N=2 suffices.
Non-Fermat primes can certainly be generated. Exact balance is a much
stronger property than the support condition needed in Theorem 1.

## 5. The 65537 obstruction and an extremal arithmetic certificate

With primitive root 3 modulo 65537, h=32768=2^15. The primes below
65537 occupy the valuation layers 0,…,14 in the exact numbers

    3260, 1631, 817, 444, 203, 93, 49, 26,
    9, 4, 2, 1, 1, 2, 0.

Layer 14 consists of the two residue classes whose square is −1:

    256 = 2^8,           65281 = 97·673.

Neither is a prime below 65537. Thus the primewise divisor distribution
of the full primorial cannot be exactly uniform, although its support is
complete. This distinction is not remedied by merely having 6542 binary
prime choices. An extra independent factor can annihilate a Fourier mode
only if its own residue supplies the missing layer.

The least prime in these two residue classes is 327941. The companion
certifies this by a complete sieve up to that value. Adjoining it to the
full primorial supplies the sole missing layer and makes every signed
residue occur exactly 2^6528 times. The theorem supplies this exact huge
multiplicity; no enumeration of the primorial's divisors is performed.

In fact this gives an exact optimization statement: **the minimum possible
largest prime factor of any squarefree source exactly balanced on H_65537
is 327941**. Necessity follows because every such source needs a prime in
layer 14. Sufficiency is supplied by the 15-prime source below. The same
minimum largest-prime bound is attained by its 16-prime generating
extension. This conclusion uses the finite prime check only for the first
prime in a precisely specified pair of residue classes, not to infer an
unbounded pattern from a list.

A much smaller source demonstrates the sharp generating threshold. Take
one prime from each layer, in layer order:

    3, 19, 17, 13, 233, 839, 431, 157, 6373, 10163,
    257, 2, 64513, 61441, 327941.

Its 32768 divisors give one representative per signed class, but it
cannot produce 65537. Adjoin the additional nonresidue prime 5. The
resulting squarefree source is

    N = 36054720350329695914767176044436363972270.

Its 65536 divisors give exactly two representatives per signed class.
Exactly one unordered split produces a multiple of 65537:

    d = 135625235,
    N/d = 265840795411928288380603919649882,
    (d+N/d)/65537 = 4056346726458768152045471341.

This is a certificate for an allowed split from the displayed source.
It does not say 65537 is the least candidate, nor that this source is
visited by either greedy prime sequence.

## 6. Transfer contract, tests, and next work

The precise maps are:

1. A prime subset address goes to its divisor residue in H_q. Multiplicity
   and all Fourier factors are retained; the exact divisor and full unit
   residue are lost. Retaining the prime-address subset restores the
   divisor. Retaining N modulo q restores the candidate root test.
2. At a smallest missing q, the published support theorem removes the
   support uncertainty, leaving the one-bit quadratic phase. Appending
   a new prime multiplies that phase by its Legendre symbol. This is the
   exact operation behind the trapping theorem.
3. At q≡3 modulo 4, the Legendre symbol is an intrinsic Paley arc from
   zero. At q≡1 modulo 4 it is an undirected adjacency. Neither map
   imports Hamiltonian-path information or establishes recurrence.
4. Binary balancing and circulant rigidity both descend to the 2-group
   classification. This shared classification does not identify their
   source objects or remove the quadratic-phase marker.

The useful next proof target is Corollary 3(b)/(c), or the weaker
trajectory-specific assertion that every currently blocked smallest
missing q eventually receives a nonresidue output. The relevant analytic
quantity is a character of the **selected next prime**, not a character
of an arbitrary available split sum or of a randomly chosen divisor.
Booker's freely chosen branches do not supply that greedy selection law.
No bound on the waiting time is proposed from 80 initial terms.

Run `python3 -B verify_prime_branching.py` and
`python3 -B -O verify_prime_branching.py` in this directory. The universe is:

- all 74,612 cyclic residue multisets of length at most 5, with group
  orders 1 through 16, for the full balancing iff;
- 1,440 prime-source/modulus cases, using every subset of
  {2,3,5,7,11,13} after excluding the modulus, at odd primes through 101;
  literal divisors independently check quotient dynamic programming,
  candidate counts, complement parity, and the support criterion;
- signed primorial support at every odd prime through 499;
- the first 80 actual greedy split-sum steps, including all rejected
  smaller unused prime candidates and explicit divisor witnesses;
- all primewise 2-adic layers at 3,5,17,257,65537; the minimum repair prime
  through 327941; literal enumeration of the extremal 16-prime source;
- six malformed interfaces, including Boolean/composite moduli,
  repeated/nonprime factors, and a factor equal to the modulus.

All 79,609 always-active checks pass in normal and optimized Python.
The canonical JSON receipt SHA256 is
`596864b2b827961bd4cc37a531b81928d240f90ee8cfa674114d659874cadda6`.
Finite checks support unrestricted elementary proofs and the stated
finite examples; they do not upgrade an open recurrence conjecture.

The primary sources were checked directly. The old 2016 coverage
conjecture was not treated as open after finding its 2017 proof. The
largest-factor and least-factor Euclid–Mullin sequences, the split-sum
variant, and branching existence have been kept distinct throughout.

