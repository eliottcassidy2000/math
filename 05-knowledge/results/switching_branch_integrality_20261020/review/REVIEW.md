# Independent review: branch integrality and pairing boundary

Verdict: PASS. No substantive mathematical correction is required in either audited document.

The arithmetic theorem proves both branch counts integral for every positive n other than 2. At n = 2 both counts are 1/2. It does not construct the combinatorial pairing requested in HYP-9172(b). The two negative pairing statements have the stated narrower scopes: the matching-flip operation fails already at n = 4, and the inherited n = 8 census excludes a fixed-point-free involution on all unadorned even-Aut isomorphism types. Neither statement excludes a suitably decorated or differently weighted carrier.

## Versions reviewed

- `session12-switching/BRANCH_INTEGRALITY.md`, SHA-256 `9fbb830de4d47a16d4c961a03637e0fac0fc9be00236e9d0d9459415064f4e4a`.
- `session12-switching/PAIRING_BOUNDARY.md`, SHA-256 `a0b76d94b5f36cfdd420a0f29bc6af6eeb83fc3a970827495f725cde0a40636c`.

Only files in this review directory were written. The producer directory, repository files, and Git state were not changed by the reviewer.

## Proof audit

1. **Exact external dependency.** Independently opened the publisher's PDF of Brusyanskaya, Klyachko, and Vasil'ev, *What do Frobenius's, Solomon's, and Iwasaki's theorems on divisibility in groups have in common?*, Pacific J. Math. 302 (2019), 437–452, DOI 10.2140/pjm.2019.302.437. Printed p. 437, the second PDF page, explicitly gives the required theorem: the number of solutions of x^m = 1 in a finite group H is a multiple of gcd(m, |H|). The paper proves a generalization that contains this theorem. The citation and page locator in the producer document are accurate. This is not a normal-complement theorem and does not use the stronger exact-root-count conjecture.

   Primary source: https://msp.org/pjm/2019/302-2/pjm-v302-n2-p03-s.pdf

2. **Dyadic stabilizer weights.** If |H| = 2^s m with m odd, an element has odd order exactly when its order divides m. Thus the odd-order element count equals the x^m = 1 solution count. Frobenius implies a(H) = mt for an integer t, and a(H)/|H| = t/2^s. No solvability or normal-complement assumption is hidden here.

3. **Restricted Burnside identity.** For each orbit O with stabilizer H, counting pairs (g, x) with x in O and g an odd-order stabilizer element gives |O|a(H). Division by |G| gives a(H)/|H|. Summing over orbits proves the identity used in Step 1. For G = S_n, odd-order permutations are precisely those whose cycle lengths are all odd. THM-479 supplies the relevant fixed-class count 2^(o_2-k+1). Its conjugacy-class denominator is z_mu = |C_{S_n}(g)|. These are the same objects and action throughout; the stabilizer is the class automorphism group, not the automorphism group of an arbitrarily chosen member tournament.

4. **Termwise 2-integrality.** With k odd parts, sum floor(l_i/2) = (n-k)/2 and each pairwise gcd is at least 1, giving the displayed lower bound on e. All l_i are odd, so v_2(z_mu) = sum v_2(m_j!) <= v_2(k!). The k = 1, 2, 3 cases are correct. For k >= 4, binomial(k-1,2) >= k-1 >= v_2(k!), the latter also following from v_2(k!) = k-s_2(k). The only failure is k = n = 2 with mu = (1,1). This proves each other odd-partition summand has odd reduced denominator.

5. **Intersection and complementary branch.** A rational with both a power-of-two denominator and an odd denominator is integral. Subtracting N_odd from the actual orbit count gives integral N_lev. At n = 1 the values are (1,0); at n = 2 they are (1/2,1/2). No exceptional n divisible by four has been overlooked.

6. **General-prime criterion.** Replace the odd part by the p-prime part of |H|. The same argument proves B_p in Z[1/p]. The stated classwise fixed-count divisibility gives B_p in Z_(p), hence in Z. Zero fixed-point counts are correctly excluded from the valuation condition. The criterion really is sufficient rather than necessary: take G = C6 acting on a fixed singleton plus G/C2, and p = 2. The three p-regular elements have fixed counts 4,1,1, all centralizers have order 6, and their terms are 2/3,1/6,1/6. Their sum is 1 although two terms are not 2-integral. This optional illustrative example is verified in the independent control script.

## Source and pairing audit

Read THM-479 in full and §5 of `05-knowledge/results/procgen_tbij_20261001_natural_bijections.md`, including Theorem B and its proof. The source states the stabilizer-ratio formula and records exactly the n = 8 census quoted by the producer. Its groups and proportions are consistent:

- 19 C2-Sylow types, rho = 1/2;
- 2 C4-Sylow types, rho = 1/4;
- 2 V4-Sylow types, rho = 1/4 and 3/4;
- 2 C8-Sylow types, rho = 1/8;
- 2 D8-Sylow types, rho = 1/8 and 5/8.

Thus there are 27 even-Aut types. Their 1-rho contributions sum to 15, while half their number is 27/2. This is a parity obstruction to any fixed-point-free involution on the full unadorned carrier, irrespective of whether it mixes Sylow types. The full n = 8 census remains inherited FINITE-EXACT evidence; it was not rerun in this review.

Independently reimplemented n = 4 and n = 6 using ordered-triple switching invariants and adjacent-transposition orbit traversal, without importing producer or repository code. The normalized representatives, stabilizers, odd-order proportions, all involution-flip targets, and commuting witnesses agree:

- n = 4: representative 0 has Aut order 4 and rho = 1/4; representative 2 has Aut order 12 and rho = 3/4. Every automorphism involution's matching flip stays in the same type. In particular, (02)(13) takes mask 0 to flip mask 7 with commuting transporter (13), and (01)(23) takes mask 2 to flip mask 5 with commuting transporter (23).
- n = 6: the four even-Aut types all have rho = 1/2 and pair 0 with 40, and 4 with 41. The two other types have odd Aut order. Hence (N_odd,N_lev) = (4,2).

The smallest orders n = 1,2,3 were also enumerated. For all five tested orders, the direct odd-order restricted Burnside sum over every labelled class and permutation agrees with the independent orbit-stabilizer sum and odd-partition formula. Independent concrete group controls recover rho(C2) = 1/2, rho(A4) = 3/4, and rho(GL(3,2)) = rho(PSL(2,7)) = 5/8. For those groups, the root-count/p-regular correspondence and divisibility also pass at p = 2,3,5,7.

## Reproduction

Run from this directory:

    python3 -B independent_controls.py > controls.normal.json
    python3 -B -O independent_controls.py > controls.optimized.json
    cmp controls.normal.json controls.optimized.json

Both modes passed and produced byte-identical JSON. All mathematical checks use explicit exceptions, so optimization does not disable them. These finite controls are corroboration; the all-order proof is the audited argument above plus its accurately identified external dependency.

## Promotion recommendation

Record the numerical portion of HYP-9172(b) as PROVED, with Frobenius root-count divisibility explicitly cited. Keep the structural combinatorial-pairing portion OPEN. Preserve FINITE-EXACT/inherited provenance on the n = 8 carrier obstruction. No claim of novelty relative to the entire literature is established by this review.
