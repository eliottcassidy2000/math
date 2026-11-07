# Session 218: Smith numbers, repunits and totient scaling

**COMPLETED: proved bridges, exact classifications and independent audits. Carmichael's totient conjecture and Chan–Pak completeness remain OPEN.** No claim of historical priority or fixed-base repunit-prime infinitude is made.

## 1. Which sequence is 14,26,34,38,...?

The even nontotients begin **14,26,34,38,50,62,68,74,76**. They are numbers with no preimage under Euler's phi. The repo also uses the first four in a different sequence: exponents k for which composite 2^k−1 satisfies the Camion strong-tournament congruence. Its corrected complete list through 64 is **4,6,10,12,14,26,34,38,62**.

These sequences are not equal: **50 is an even nontotient but fails the Camion test**, at prime 251. We independently reproduced the complete small totient prefix and this exact obstruction. The older Camion census omitted nonsquarefree values63 and4095; MISTAKE-651 remains its correction authority.

[Prior-work reconstruction](provenance/PROVENANCE.md) · [Independent provenance review](audit/PROVENANCE_REVIEW.md).

## 2. A genuine bridge, rather than a prefix coincidence

If M=2^q−1 and W=(2^q+1)/3 are both prime with q>5, then **2q is an even nontotient**, while 2^(2q)−1 belongs to the inherited Camion family. If ell=2q+1 were prime, Fermat would make it divide3MW; it would equal M or W, both larger than ell. The full inverse-totient fiber of 2q, for prime q>3, is {2q+1,2(2q+1)} when 2q+1 is prime, and empty otherwise. The exceptions q=3,5 are sharp.

[Proof, support-aware scaling and boundaries](totient/PROOF.md) · [Independent core audit](audit-core/CORE_REVIEW.md).

This also explains two distinct uses of phi: phi(q)=phi(2q)=q−1 gives two INPUTS with the same output, while 2q itself can be a missing OUTPUT. Empty fibers are not counterexamples to Carmichael's no-singleton-fiber conjecture. Carmichael numbers, Carmichael lambda and Camion pseudoprimes are kept separate throughout.

## 3. Complete all-base Smith classification

**A composite n is Smith in some integer base if and only if n is not 6,8,9,12.** For every eligible n other than 4, the explicit base is

    b = n − sum_p v_p(n) p + 1.

For4, any base at least 5 works. Prime factors are counted with multiplicity. For n=2p, p>=5 prime, the formula gives b=p−1 and digits22. For all other eligible n>12, b>n/2 and the two digits are1 and F(n)−1, with every prime factor a single digit.

Thus **every even nontotient has a constructive Smith base**. In particular14,26,34,38 use bases6,12,16,18. This is a meaningful explicit map but a weak discriminator: it embeds nearly every composite sequence. It proves no fixed-base coincidence or infinitude result.

[Complete proof](smith/BASE_COMPLETENESS.md) · [Independent review](audit/BASE_COMPLETENESS_REVIEW.md). Producer checks cover every eligible composite through 100000; the independent construction check covers every eligible composite through 200000.

## 4. Repunit primes in all bases: positive seeds and sharp obstructions

For every prime p>=5, both 2p and p² are Smith in base p−1, where p is the repunit11. Allowing the base to vary therefore trivially supplies infinitely many repunit primes; one fixed base is a different question.

If R_b(k) is prime and2<=c<b, then cR_b(k) is Smith **exactly when c=k=2**. Longer repunit primes need a multi-digit multiplier. One uniform decimal construction is **1540 R_10(k)** whenever R_10(k) is prime and k>=4. By contrast, for a binary repunit prime, **no multiplier2<=c<=2^k** works; the boundary is real, since 34*7 is binary Smith.

Exact carry and block-grouping identities describe what a base change loses. In particular, **b²−1 is never base-b Smith**, so grouping2^(2q)−1 into base 2^q does not transport the Smith property. Repunit length-factorization and prime-order constraints are proved separately.

[All six Smith/repunit results](smith/SMITH_REPUNIT_BRIDGES.md) · [Independent review](audit/REVIEW.md).

## 5. When doubling repairs a missing totient

For an odd prime q and a>=1, **2^a q is a totient iff** some 2^j q+1 is prime with1<=j<=a, or q−1=2^j for such a j. The latter corresponds to a q² factor in the preimage. Complete preimages are parametrized by this unique q-supplier, distinct additional Fermat primes and a power of 2.

Examples:14 becomes28=phi(29);26 becomes52=phi(53);34 remains missing at 68, then136=phi(137). Uniform doubling therefore preserves neither missingness nor a complete inverse fiber. This criterion does not promise an eventual prime 2^j q+1 for every q.

[Dyadic classification](totient/DYADIC.md) · [Exact scaling table](totient/INTEGRATION_TABLE.md) · [Independent audit](audit/REVIEW.md).

For the full even-nontotient prefix through 76, the table retains factorizations, constructed Smith bases, digit weights and exact first dyadic repairs. In particular38 needs multiplication by 32 and62 by 128; merely doubling each term does not align the sequences.

## 6. Repair to prior tournament multiplicity claims

THM-4884 and its source had treated finite absence of some multiplicities below target 625 as refuting a global Ford-type statement. That implication is withdrawn in the canonical theorem, source and HYP-9477. Even an asymptotic growth claim would not exclude finitely many later small fibers. Exact singleton examples still refute the no-singletons tournament analogue. MISTAKE-723's numerical repair is preserved.

[Correction record](provenance/MISTAKES_ENTRY.md) · [Exact patch](provenance/SCOPE_CORRECTION.patch). The multiplicative/fiber analogy supplies no map preserving totient gaps and Hamiltonian-path gaps: halving14 and26 gives7 and13, but7 is forbidden while 13 is attained. No Chan–Pak coverage or eventual threshold is claimed.

## Reproduce

Source-path mentions inside frozen reports record original provenance. From this directory use:

    python3 totient/verify.py
    python3 totient/verify_dyadic.py
    python3 totient/integration_table.py
    python3 smith/check_smith_repunit.py
    python3 smith/check_base_completeness.py
    python3 provenance/check_distinct_prefixes.py
    python3 audit-core/verify_core.py
    python3 audit/verify_independent.py
    python3 audit/check_base_existence.py
    python3 audit/check_prefix_independent.py
    python3 audit/check_integration_independent.py

The totient producer scripts use SymPy. All independent verification scripts and the Smith/provenance producer scripts use the standard library. Receipts are beside the scripts; SHA256SUMS freezes the full packet. Finite checks support the written proofs. The historic full Camion census and full tournament multiplicity inventory were read with their correction lineage, not re-enumerated wholesale.

## Research method and limits

The live board tracked five objects: totient fibers, Camion residues, repunit factorization, digit carries and tournament-class fibers. Exact source/target maps and hostile examples replaced numerical alignment. The principal surviving sidecars are prime support, factor multiplicity, chosen base, full CRT residues and the distinction between a target count and a source object. Current concurrent sheaf work was preserved but supplies no dependency for these arithmetic results. The useful next questions concern fixed bases and controlled scaling, where the broad moving-base classification no longer trivializes compatibility.
