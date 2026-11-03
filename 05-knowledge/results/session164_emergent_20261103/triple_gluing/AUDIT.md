# Independent audit of the S164 triple-composition packet

**Verdict:** Core mathematical claims pass. The published note incorporates two domain clarifications: C,D>1 in the exact-fiber paragraph, and N>=2 in the divisor-balance classification. These are wording corrections, not failures of the substantive composition results.

## Results

- The literal terminal amalgamation preserves both old tournaments and yields exactly E(r+s). The path address inverse recovers each input, so the claimed multiplication is a bijection rather than just a count match.
- Gaussian primitivity, the missing relative-conjugation bit, the k-input 2^(k-1) fiber size, the coordinate-gcd carry, the telescoping cocycle, and the multiset two-valued associativity are correct. Multiplicities at complete cancellation are correct whenever the combined input support is nonempty.
- The explicit shared-prime p^2qr versus qr fibers and all displayed Pythagorean triples are correct. The shared-prime alignment predicate is precisely the stated balance predicate in this family.
- Cassini's product is correct for every n>=1. Its norm-coprimality proof can also be shortened to D-C=F_(2n+2), hence gcd(C,D)=gcd(F_(2n+1),F_(2n+2))=1. For n=3j+2 the parity and nonempty split-support hypotheses are satisfied. The positivity polynomial and three square-sum identities are correct.
- The guard counterexample is sound: (2+i)^1=2+i, (2+i)^3=2+11i, and (2+i)^4=-7+24i yield canonical charts (3,1), (13,9), and (31,17). The norm-625 product has carry one and loses the power-of-three denominator property.

## Domain clarifications

1. Section 3's zero-output formula reads 2^(z-1). If C=D=1 were allowed, z=0 and the actual quotient has one orbit, not one half. The introductory v=0 convention makes this worth stating explicitly. Requiring C,D>1 suffices; a more general alternative is to exclude only C=D=1 and state that exceptional fiber separately.
2. Section 4 inherits a classification valid for N>=2. Since the preceding section retains norm one, explicitly preserve this restriction. Under the given proper-divisor definitions, N=1 also satisfies F=S+U trivially.

## Reproduction

Both `python3 check.py` and `python3 -O check.py` pass 88,873 explicit checks and produce byte-identical JSON. Audited script SHA256: `239337d9a2dbbb0695bd44ef9da958163d6ca81fe1dfac467a03377659d5cc6a`.

An independently written supplemental verifier, `audit_extra.py`, also passes identically in normal and optimized modes. It checks 255 exponent-profile pairs over primes 5 and 13 with exponents 0..3 (including one unit input), 1,065 exact output multiplicities, direct retained-arc amalgamation for 36 pairs r,s=1..6, and Cassini/norm-coprimality for every n=1..1000.

## Prior-source overlap

Direct comparison with the inherited endpoint/Pythagorean note, the S163 triple-codec and strands notes at commit d9d9f3e483c54d83ebd569b91a85b7d015923cfe, the divisor-balance note, and the ternary-bridge note supports the packet's attribution. The endpoint decoder, Gaussian quotient, phase-sensitive AP continuation, divisor classification, and rational ternary guard are inherited. The new packaging is the literal composition and its relative-bit/cancellation accounting, together with the Cassini branch selection.

The external primary sources were opened and support the cited prior mechanisms: Yekutieli's rational-circle group and fixed-hypotenuse organization (https://arxiv.org/abs/2101.12166), and Buchstaber–Dragovic §2.1 Example 1's [x+y,|x-y|] two-valued law and duplicate-unit convention (https://amj.math.stonybrook.edu/pdf-Springer-final/018-0085.pdf). No historical-priority claim is needed or established by this audit.
