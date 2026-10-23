# Both stripes in one space: exact partition functions and frequency regions

Status: PROVED consequences of the audited occurrence-dependent branch classification in fibonacci/BRANCH_ADDENDUM.md. A single shared expansion dictionary or shared parse tree is a stronger condition and is not assumed here.

Let a(H),b(H),c(H) count A,B,C in the canonical prefix of length H. The deleted-zero grammar permits independent B→A changes. The retained-zero grammar additionally permits independent C→B changes. A positions are forced in both.

## Exact multivariate enumerators

Sum x_A^(number of A) x_B^(number of B) x_C^(number of C) over all length-H prefixes that extend to a coherent infinite stripe. The exact polynomials are

Z_D,H=x_A^a x_C^c (x_A+x_B)^b,
Z_R,H=x_A^a (x_A+x_B)^b (x_B+x_C)^c.

The proof is the Cartesian-product branch theorem: each free position supplies its own two-term factor. This is not multiplication of unrelated marginal observations; the theorem has already proved that these particular choices are mutually compatible. A specified finite initial stripe pins its corresponding factors to monomials. Every remaining choice still extends.

Consequently the supplied old35B, when imposed on a canonical-C position, halves the retained-zero prefix count for H≥35 and leaves its asymptotic exponential growth rate unchanged. In the deleted-zero grammar it is impossible. In the fixed-block retained-zero grammar it survives only one further step; that is a different model, not a contradiction.

Setting all variables to1 gives2^b and2^(b+c). The logarithms divided by H converge to phi^-3 log2 and phi^-1 log2. The partition functions also give all color-count distributions by coefficient extraction, not merely their totals.

## Every attainable limiting frequency

Write alpha=phi^-2 and beta=phi^-3, so2alpha+beta=1. If a retained-zero stripe has limiting color frequencies, they have the form

f_C=alpha-u,
f_A=alpha+v,
f_B=beta+u-v,

where0≤u≤alpha and0≤v≤beta. Here u is the absolute density of changed canonical-C positions, and v that of changed canonical-B positions. Necessity follows by counting the changes. Conversely, every pair(u,v) in this rectangle is attainable: enumerate the eligible C positions and choose a subset with relative density u/alpha, and independently choose B positions with relative density v/beta. For example choose position j when floor(j theta)−floor((j−1)theta)=1 for the desired relative density theta. The canonical eligible-position densities then give the required absolute densities.

The deleted-zero family is the edge u=0 of this region. Finite prescribed prefix choices do not affect it. Arbitrary branch choices need not have limiting frequencies at all; increasingly long alternating blocks of decisions can make a density oscillate.

In particular every branch has lower asymptotic A density at least alpha>1/3. Even this enlarged family cannot have all three colors at equal frequency. Three available labels therefore do not by themselves supply a measure-preserving cyclic color symmetry.

Under the uniform distribution on all compatible retained-zero prefixes, each free bit is independent and fair. The limits of the expected prefix frequencies are

A:1/2, B:phi^-1/2, C:phi^-2/2.

This is a precise maximum-entropy choice on finite prefix families, not a claim that the user must choose it or that it has a physical interpretation.

## Aperiodicity survives every branch choice

Set gamma=phi² (distinct from the density beta=phi^-3 above). For every positive period q, the canonical A positions floor(k gamma) and canonical C positions floor(k gamma)−1 each meet every residue class modulo q infinitely often. This follows from density of the fractional parts of k gamma/q: an interval of length1/q specifies each possible residue of floor(k gamma).

Every compatible branch is forced to have A at the former positions and not-A at the latter. An eventually periodic A-indicator of period q would therefore assign both values to every residue class, a contradiction. Thus the A-indicator, and hence the full stripe, is not eventually periodic in ANY retained-zero branch, even if its limiting frequencies are rational or do not exist.

A finite-horizon check illustrates the need to retain the clock: the first5000 canonical A occurrences cover only95 of144 residue classes modulo144; full coverage first occurs at occurrence11001. This is a checked finite control, not a uniform quantitative density theorem.

Its ordinary generating function is consequently not rational. A bounded integer sequence with a fixed constant-coefficient rational recurrence is eventually periodic: finitely many possible consecutive state blocks evolve deterministically. This argument applies specifically to the forced A-indicator. Other color indicators can be identically zero in some branches and must not be called nonrational universally.

## Additive-number-theory guard

The canonical choice is useful because it retains a deterministic Beatty/Zeckendorf arithmetic register. The broader family is useful because it separates representation grammar from selection. Neither should erase the other.

Any subset of canonical C positions can be selected to remain C. Choosing that subset by primality or by sum-of-two-squares membership would encode the tested property by definition, not discover it. Meaningful arithmetic statements must either fix the selector independently or quantify over the whole family. The additive packet instead retains actual labels, shared Fibonacci digits, carries and quadratic norms, rather than inferring primality from an unconstrained recoloring.

Run check_branch_partition_functions.py normally and with python3 -O. It compares direct prefix enumeration with independently multiplied color-count polynomials, including the old35 pin, and checks the finite expectation identities.
