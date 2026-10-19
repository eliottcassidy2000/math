# When the divisor kernel already determines its quadratic phase

19 October 2026. **PROVED, with independent review PASS.** This complements the [split-sum prime-generation packet](PRIME_BRANCHING_PHASE_AND_FERMAT_BALANCE.md); no completeness claim is made.

## 1. A recovery theorem for every finite cyclic group

Let t_1,...,t_k lie in G=Z/hZ, and let D(x) count binary subset sums. Write T=sum_i t_i and let

    Per(D)={a in G : D(x+a)=D(x) for every x}.

Complementing a subset proves D(x)=D(T-x). Consequently the complete set of reflection centers

    C(D)={c in G : D(x)=D(c-x) for every x}

is exactly the coset T+Per(D). Indeed composing two such reflections gives translation by their difference, and composing the known reflection with any period gives another reflection.

The integer |Per(D)| divides 2^k: D is constant on every translation orbit, each of which has size |Per(D)|. Thus the period subgroup has power-of-two order, even when h has odd factors. This is a useful restriction that does not follow merely from reflection symmetry.

If h=2^s m with m odd and s>=1, the parity of T is determined by the full array D whenever |Per(D)|<2^s. If |Per(D)|=2^s, the array's reflection centers have both parities. This last assertion concerns reflection information; it alone does not assert that every such center is realized by another binary list of the same length.

Proof: the unique subgroup of order 2^r has generator h/2^r=2^{s-r}m. Its elements are all even exactly when r<s. The period group cannot have larger order than the Sylow 2-subgroup by the preceding divisibility.

## 2. Exact recovery at Fermat primes

Let q>3 be a Fermat prime, h=(q-1)/2=2^s, and identify F_q^*/{+1,-1} with Z/hZ using any primitive root. For a squarefree prime source N not divisible by q, let D be its signed divisor-count kernel and k its number of prime factors.

The quadratic character (N/q) is (-1)^T. The following recovery statements are therefore exact:

1. If D is not uniform, it determines (N/q). Calculate C(D); all its elements have the same parity, which is the desired parity of T.
2. If D is uniform with multiplicity one, then (N/q)=-1. In particular q cannot divide any split sum d+N/d.
3. If D is uniform with even multiplicity, both phases occur among abstract binary lists of the same length k. Thus the kernel and list length cannot in general determine the phase. This does not assert both realizations in a prescribed prime set or a greedy trajectory.

For item 1, a proper subgroup of a cyclic 2-group contains only even residues. For item 2, uniform multiplicity one means k=s. The binary-balance classification requires exactly one logarithm in each valuation layer 0,...,s-1; there is exactly one odd logarithm, so T is odd. An alternative proof uses complementation: if T were even, its two fixed residue classes under x->T-x would each contain an odd number of subset addresses, impossible because complementation is fixed-point-free when k>0.

For item 3, k>=s+1. A balanced list with odd total is the standard layer list

    1,2,4,...,2^(s-1)

followed by k-s zeros. A balanced list with even total is the same list followed by an additional 1 and k-s-1 zeros. The layer classification proves that both distributions are the same uniform array of multiplicity 2^(k-s), while their total parities differ. Distinct addresses may have equal or zero residue labels, as distinct primes may have equal residues.

The actual prime-source pair N=42=2*3*7 and N=66=2*3*11 at q=5 realizes this ambiguity with k=3: both signed kernels are (4,4), but their phases differ. This pair is inherited from the companion packet, not newly discovered here.

## 3. The root class still matters after phase recovery

Recovering the quadratic phase does not recover the number of eligible split sources. The exact rule is

    Z_q(N)=D([sqrt(-N)])

when -N is a square, and zero otherwise. The signed kernel can lose the choice of root class even while retaining the phase.

A particularly small actual hostile is q=5:

    N=11: D=(2,0), (-N/5)=+1, Z=0;
    N=19: D=(2,0), (-N/5)=+1, Z=2.

The single split sums are 12 and 20. Thus equal nonuniform kernels, equal source size, and equal quadratic phase do not determine candidate eligibility without a full-support premise. These are arbitrary squarefree sources; q is not their smallest missing prime. This does not contradict the published smallest-missing-prime saturation mechanism.

A two-factor example is q=17, N=86=2*43 versus N=118=2*59. Both kernels, in signed representatives 1,...,8, are

    (2,1,0,0,0,0,0,1).

Both have positive phase for -N, but Z is respectively 0 and 2. Here the split sums for 118 include 2+59=61 and 1+118=119=7*17; those for 86 are 45 and 87.

## 4. Odd-order quotients lose the phase at a different step

If h is odd, the period subgroup must be trivial, since its order divides both h and 2^k. Thus C(D) has exactly one element: the complete binary subset-sum kernel determines the total T uniquely in G. The same proof works in any finite abelian group of odd order.

For q congruent to 3 modulo 4, the signed group H_q has odd order. Therefore its full divisor kernel determines the signed product class [N] exactly. Nevertheless the quadratic phase is not determined: the Legendre character does not descend to H_q because -1 is a nonresidue. Passing from N to -N preserves [N] and reverses the phase. The earlier actual sources N=330 and N=510 at q=7 illustrate this: both have the same signed kernel and signed product class, but opposite phases.

This locates two distinct failures precisely. At Fermat primes, the character does descend to the quotient, but a translation stabilizer can hide the center. At primes congruent to 3 modulo 4, the center is uniquely recoverable, while the signed quotient itself has discarded the character. In the latter case an odd-order reconstruction theorem cannot recover data that were removed before the kernel was formed.

## 5. What this changes in the proof program

The full divisor kernel does not always require an independently retained quadratic marker. At Fermat moduli a nonuniform kernel already determines that marker, and uniform multiplicity one determines the blocked phase as well. The exact blind regime is uniform even multiplicity at the abstract-list level. In contrast, the root-class location can remain invisible even outside that regime.

This is a concrete instance of recovering lost source information from an involution plus stabilizer, rather than assuming every compressed representation loses the same data. The remaining greedy-prime obstacle is still recurrence of the selected next prime outside an index-two subgroup; none of these recovery statements establishes that recurrence.
