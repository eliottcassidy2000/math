# Independent audit: square-gap prime support

Status: PASS. The all-center mathematical statements audited below are PROVED by the elementary arguments; the census is FINITE-EXACT only for 1 <= e <= 3000. No general hourglass nonexistence or novelty claim is established.

## Findings

No mathematical correction to the supplied RESULT.md is required. Necessary conditions remain explicitly necessary; projective incidence is not mistaken for arithmetic addition. The result is suitable for publication with the source-material exclusion below.

1. **Gap/Pythagorean bijection and count.** From a²+b²=2e², a and b have the same parity, so u=(b+a)/2 and v=(b-a)/2 are integers. Positivity and distinctness are automatic for a positive gap. The inverse gives d=2uv. At stripped center n, Gaussian factorization has conjugate exponents a_i+j_i and a_i-j_i. A zero coordinate product occurs only for the balanced zero vector. Equal absolute products at fixed norm determine the unordered absolute coordinate pair; the full ambiguity is exactly units and conjugation. Thus the nonzero vectors modulo simultaneous sign, not independent coordinate signs, index gaps, yielding (product(2a_i+1)-1)/2. There is no diagonal orbit because n²=2u² has no positive integer solution.

2. **Stripping.** An inert prime q dividing the center must divide both coordinates; otherwise after removing their common q power one gets a nonzero square root of -1 modulo q. Iteration gives the complete q power. The analogous 2-adic iteration follows from squares modulo 4. Therefore all inert/2 factors divide all seven roots, not merely all gaps, and the squared stripping identity is exact.

3. **Valuations.** After extracting p^(a-|j|), a nonzero allocation j leaves only one of the two Gaussian primes over p. A coordinate divisible by rational p would force the other conjugate prime to divide the residual too, a contradiction. Since p is odd, the factor 2 in the gap contributes no p valuation. Hence the exact unbalanced formula holds. The balanced case guarantees only v_p(d)>=2a, and e=65,d=3000,p=5 disproves equality. The unique-minimum p-adic obstruction works for all sign choices. A positive unique maximum allocation gives a strictly unique minimum valuation, including in the presence of balanced coordinates. If all allocations are zero, maximum multiplicity is automatic.

4. **All-divisor heredity and exact root gcd.** For any t|e, t²|d implies t² divides both endpoint squares, hence t divides both roots. This proves the exact intersection identity without any prime-support restriction. At stripped center n, the proposed t divides all six endpoint roots. If M_i>0, a gap with exact valuation 2(a_i-M_i) rules out further common root division; if M_i=0, all center powers are already exhausted. Thus the common gcd with n is exactly t. Restoring inert/2 factors multiplies that gcd by m. After primitive reduction, maximum-twice implies at least two p-unit gaps at each surviving center prime, giving gcd(e,gcd(d_r,d_s))=1. This is not a claim that the gaps are pairwise coprime away from center primes.

5. **Projective indexing.** With squarefree split support, {-1,0,1}^k is in bijection with F_3^k, and the nonzero simultaneous-sign quotient is exactly PG(k-1,3). Only coordinate permutations/signs arise from the named ordering/orientation choices; no arbitrary projective symmetry is claimed to preserve gaps. Coordinate-zero subspaces agree with square-divisibility and center deletion, including the empty all-deleted case. The e=65 projectively dependent triple has gaps 4056,3000,3696 and no Schur relation. Integer phase addition does not reduce modulo 3 without losing a nontrivial Gaussian phase. The arithmetic sidecar is essential.

## Independent computation

Reproduce with:

    python audit.py > audit-normal.json
    python -O audit.py > audit-optimized.json

The audit uses a smallest-prime-factor sieve independently of the candidate's trial division; it independently enumerates endpoint squares and Pythagorean coordinate squares. Both enumerate exactly the candidate Gaussian gaps at every center 1..3000.

- 3,172 positive gaps.
- 24,496 exact hereditary checks, for **every divisor** of every center, including 1 and e.
- 4,501 triple checks for the exact common divisor of the center and all six endpoint roots; the primitive pair-gcd consequence is checked for every triple passing the valuation filter, without assuming a Schur relation.
- At e=65, 4 of 4 triples pass the valuation filter; at e=1105, 150 of 286 do. None forms a Schur relation.
- Independent dual projective construction via kernels of covectors finds 13 points and 13 four-point lines, with a unique line through every pair. All eight coordinate-deletion subsets at e=1105 agree with hereditary arithmetic scaling.
- Explicit negative controls include the balanced-valuation strengthening, the count-three near miss at 125, and both filter insufficiency examples. Nonpositive valuation arguments and an injected failed require are rejected with RuntimeError in both normal and optimized modes.
- The candidate is replayed from an audit-local copy, under both normal and -O Python. Both outputs match its supplied certificate exactly. Neither candidate nor audit uses assertions for validity guards.
- The full audit itself passes in both normal and -O mode, with identical JSON outputs. Candidate SHA256: ed1ab261f28aec007037c5b54a7356da25d1018bb9a82b9c32d718a272e865d1.

The all-center proofs do not rely on the absence of Schur triples in this finite range. The arithmetic primitive-hourglass consequences remain necessary restrictions on a hypothetical object, not existence demonstrations.

## Provenance and publication boundary

An independent fresh request to https://zenodo.org/api/records/21910218 confirms DOI 10.5281/zenodo.21910218, publication date 2026-08-13, and resource subtype **preprint**. The source PDF lists Valentyn Voronsbekher. Metadata's creator-name punctuation is awkward, but the paper author name agrees with the candidate. The supplied PDF MD5 matches the newly fetched primary record. The metadata title begins “On magic hourglasses...” while the PDF title begins “Magic hourglasses...”; this is a title variant, not a source mismatch.

The primary PDF's Lemma 7 supplies inert/2 stripping, Lemma 10 supplies Gaussian allocation formulas for its prime-power and semiprime cases, and Proposition 11 gives the one-prime-power valuation obstruction. The general multi-prime Gaussian argument in the candidate is self-contained. The candidate correctly warns that positive products in the source formulas need absolute values. It properly credits the known obstruction and does not elevate the source to refereed literature.

The preprint's semiprime argument is **not independently audited here and not used**. Its separate authorship/verification note describes assistance for its later elliptic-curve argument; this does not license importing that argument as verified.

Publication must use an explicit allowlist. The input directory contains MagicHourGlassPrimes.pdf and preprint.txt, which are full third-party source material, not authored session results. **Exclude both from the public packet**. A suitable candidate allowlist is RESULT.md, check.py, check.out, check-optimized.out, and certificate.json; include an authored citation/status note if desired. This audit directory likewise needs only its authored report, script, and output certificates, not fetched primary metadata. No full source article was copied into the audit directory. The candidate's temporary-source-path paragraph should be adapted or removed when publishing, so public reproduction instructions do not imply that excluded local copies ship with the packet.

Repository AGENTS.md, current guardrails, research protocol, and core-paper status conventions were inspected read-only. No repository files were changed; publication and Git coordination remain with the parent.
