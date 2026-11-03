# Fixed-root iteration and finite-ray obstruction

## Exact profile dynamics

Use the notation and construction of `../gadgets/RESULT.md`: h=H(T), s is the number of Hamiltonian paths starting at v, q=H(T-v), and F counts admissible ordered path pairs (A,B). For k≥1 transitive new vertices X with v→X→T-v, the complete rooted profile is

    h' = h-s+q+(2^k-1)(F-q),
    s' = q,   q' = q,   F' = q+2^k(F-q).                 (1)

The first two identities are the prior gadget theorem; q'=q because deleting v leaves the ordinal sum X→(T-v). To derive F', split its ordered path pair by the old vertices belonging to A and B. If A contains no old vertex, admissibility forces A empty, since every x∈X loses to v; this case contributes q. Otherwise A's old subpath must end at an old vertex beating v. Both paths must place all their X vertices before their old vertices because X dominates T-v. Every partition of X between A and B is possible and determines both transitive blocks uniquely, so each original pair counted by F−q yields 2^k new pairs. This proves the last identity.

Consequently, another insertion of l vertices at the *same root* gives the same profile as a single insertion of k+l vertices. Explicitly, s′=q′=q and F′−q=2^k(F−q); the new count is h−s+q+(2^(k+l)−1)(F−q). Indeed, the newest block dominates all older inserted vertices, so the union is itself transitive and the constructions are isomorphic fixing T. Thus repeated fixed-root insertions cannot yield a new digit-by-digit Euclidean induction: they trace precisely one exponential ray. This does not obstruct changing roots, changing the gadget, or using an infinite adaptive set of seeds.

## Finite-ray simultaneous prime obstruction

Consider any finite family R of integer rays a_i+b_i 2^k, k≥1 (constant rays allowed). There are infinitely many odd primes p for which neither p nor 7p lies on any ray in R. More strongly, they contain a full reduced residue class modulo some explicitly constructible odd modulus.

Let r be the number of rays. Choose m=12r+1 and put q=2^m−1. This m is not divisible by 3 and satisfies 2^m−1>(2rm)^2. To verify the inequality elementarily, at r=1 it reads 8191>676. Increasing r by 1 multiplies 2^(12r+1)−1 by more than 4096, while the right side is multiplied by ((r+1)/r)^2((12r+13)/(12r+1))^2 <16. Induction proves the inequality for every r≥1. Then gcd(7,q)=1. As q is odd, φ(q)≥sqrt(q): indeed, for each odd prime power ℓ^e, φ(ℓ^e)/sqrt(ℓ^e)=ℓ^(e/2−1)(ℓ−1)≥1. Modulo q, each ray takes at most m values, since 2^m≡1. The combined forbidden residues for p, namely a_i+b_i2^k and 7^(-1)(a_i+b_i2^k), therefore number at most 2rm. Since φ(q)>2rm, there exists a reduced residue a modulo q outside all forbidden residues. Dirichlet's theorem gives infinitely many primes p≡a (mod q). For every such p, neither p nor 7p belongs to any ray. The assertion is unconditional but uses Dirichlet's theorem; it gives no claim about primes on an individual ray.

The same conclusion applies if finite exceptional seed counts are adjoined: enlarge m until φ(q) exceeds 2rm plus twice their number. If ordinal sums of these gadgets are also allowed, their counts are products of the allowed factors. A prime p can then occur only if p itself is an allowed factor. If 7 is not an allowed factor, 7p can occur only if 7p itself is an allowed factor (for p≠7), since its only positive divisors are 1, 7, p and 7p. Thus finite fixed-root seed rays plus ordinal products cannot settle either infinite lane. This remains only a restriction on this particular strategy, not on tournament path spectra.

Example: the strong seed ray 5+9·2^k from the gadget packet has residues {10,14,15,23,25} mod 31. Their multiples by 7^(-1)=9 mod 31 are {2,8,11,21,28}. Hence every prime p≡1 mod 31 avoids both the ray and its seven-prime counterpart; p=311 is an example.

## Checks

`python3 check.py` and `python3 -O check.py` use explicit exceptions and produce matching output; each also checks the explicit modular example, the finite inequality for r=1,...,100 and exhaustively checks (1), using the independent path-subset DP from the gadget packet, over every labelled tournament of orders 2–5, every root and k=1,2,3. This finite test is a control only; the proofs above are all-order. The proof's use of Dirichlet is standard and is not a claim of novelty.
