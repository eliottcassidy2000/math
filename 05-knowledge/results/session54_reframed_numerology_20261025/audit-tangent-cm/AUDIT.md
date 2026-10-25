# Text-only publication edition

This is a publication-safe derivative of the completed independent audit. Full third-party PDFs and extracted article texts are omitted; precise citations, inspected hashes and short notes remain in sources/REFERENCE_NOTES.md. Historical private-packet verification below is preserved as provenance. Publication-edition verification is recorded at the end.

# Independent audit: tangent-Legendre CM classification

## Mathematical verdict

PASS for the all-integer theorem and final closed packet. The integer parameters m>4 with singular-modulus j(tan²(pi/m)) are exactly {8}. This does not identify the MRB constant or defect, and does not rely on completeness of an idoneal/abelian-CM catalogue.

The reviewer independently checked the following logical joints.

1. Original Legendre substitution: with c=cos(2pi/m), lambda=(1-c)/(1+c), t=4c²−2, direct symbolic cancellation gives the stated rational j map. The primitive-root conductor is n=m/gcd(m,2), including odd m. For n>=3, denominators at t=±2 never occur. Clearing the nonconstant cubic gives field-index at most three.
2. CM/nonfundamental orders: the ring-class theorem applies to every imaginary quadratic order, not just maximal orders. All conjugates are real because Q(j) is contained in a totally real field. The ideal-class parametrization is injective, and conjugation inverts ideal classes. Thus every class has order at most two, and Q(j)/Q is elementary abelian. Its quotient inside the real cyclotomic Galois group kills all squares.
3. Finite conductor completeness: the producer's repaired condition is |U(n)^2| in {1,2,3,4,6}, not merely <=6. CRT makes every odd-prime square-group factor divide one of these six-bound possibilities. The prime 11 contributes factor5 and is genuinely excluded. The allowed primes are precisely 2,3,5,7,13, with the stated exponent bounds. The 2-adic a<=5 box is a valid loose bound; exact screening discards a=5.
4. Alternative independent conductor proof: projection to each prime-power quotient U(p^a)/±1 is surjective on square groups. For odd p its square subgroup has size p^(a−1)(p−1)/gcd(4,p−1), hence <=3. For 2^a,a>=3 its size is 2^(a−3), since −1 is not a square. This independently yields the stronger a<=4 and all the same remaining prime/exponent bounds. In particular 11 fails with quotient square-group size5.
5. Degree/collision: instead of repeating the producer's resultants, the reviewer recovered each real cyclotomic polynomial by triangular Laurent-coefficient subtraction from Phi_n. In the quotient field the rational j element was constructed by polynomial inversion. Its multiplication matrix characteristic polynomial matches every one of the 30 producer polynomials exactly and is squarefree. Thus no field degree is lost for these conductors. Brute screening n<=5000 returns exactly the same30; this is only a hostile sanity check, not the finiteness proof.
6. Integrality: the monic irreducible rational polynomial criterion is the correct one. The eight elementary-abelian conductors reduce to {4,8,12,24}; neither rational near-integers nor nonmonic integer polynomials are used as substitutes.
7. Discriminant bounds: the identity j(t)−1728=256(9t+14)^2/((2−t)^2(2+t)) is strictly positive on all conjugates. Equality would require rational nonintegral t=−14/9. The principal reduced CM form has a=1 for every order discriminant, and its absolute j value is below the positive trace in degrees2 or4. The cited2079 estimate therefore bounds all possible discriminants. Independent rational Taylor lower bounds with40 terms already certify the producer's thresholds15,19,32.
8. Reduced forms: an independent rectangular enumeration of a,c<=|D| and b²=D+4ac, rather than the producer's a<=sqrt(|D|/3) loop, matches every primitive reduced form for every discriminant −31<=D<=−3. Boundary signs and primitive gcd conditions agree. Nonmaximal discriminants −12,−16,−27,−28 are explicitly present. The remaining H−15 is verified directly in the primary paper and is coprime to the n12 polynomial.
9. Positive control: n4 corresponds only to m8 in the allowed parameter range. Its j=8000, CM discriminant−8 identification is inherited from Session52's exact AGM/lattice argument, which was reviewed. The inherited degree-two isogeny equation and its scaling isomorphism were also independently verified by symbolic substitution. This is positive CM identification, not integrality alone.

## Primary reference audit

- Cremona–Sutherland, Computing the endomorphism ring of an elliptic curve over a number field, sections2.2–2.3, explicitly discusses nonfundamental order discriminants, irreducible class polynomials, ring class fields, conjugation and class inversion. The adjacent catalogue assertion is not used. https://lucant.org/papers/2023/230123-Cremona.pdf
- Bilu–Luca–Pizarro-Madariaga, Rational Products of Singular Moduli, opening degree/integrality statement, Proposition2.1, Proposition2.5, and Table2 supply precisely the inputs claimed. Proposition2.6 also explicitly verifies existence of exactly one a=1 reduced form. https://arxiv.org/pdf/1410.1806

Both PDFs were independently fetched and locally converted to text. SHA256:
- cremona-sutherland.pdf: 2037856cfc69c7db0703ef3a18e34631767eea9171d6a05f41ba090cca470d7a
- rational-products.pdf: 976e516695aaadcfa0b7f9be4f92a23f035fb3f1b6c58a503021d4680b951363

## Executed checks

- Producer optimized certificate: normal and Python -O outputs identical, matching producer results.normal.json.
- independent.py: independent multiplication matrices and symbolic original Legendre identity; normal and -O outputs identical.
- audit_small_discriminants.py: independent forms and forty-term exponential bounds; normal and -O identical.
- hostiles.py: altered j8000→8001, omitted conductor168, missing nonmaximal D−12 form, and a false trace bound all rejected in normal and -O modes, eight rejections total.

The first producer replay used an older unoptimized full rational-composition expansion and was interrupted after replacement with the equivalent modular-Horner version. It produced no result and is not counted as verification. No repository files were edited by this auditor.

## Final closed-packet verdict

PASS. Producer manifest SHA256 fd64b457a39ca6c9746ecd1a86293c250d6bccae460babcaeb6234d4273dd585 pins14 payload files, including both primary PDFs and both inherited snapshots. A manifest/file-set verification, exact certificate replay and four producer mutation controls passed in normal and optimized modes after copying the complete packet to an unrelated temporary directory. Corrupted, missing, extra and nested-extra files were rejected in both modes. The earlier closure gaps (unlisted extras accepted; inherited sources recorded but not carried/verified) have been repaired. The audited mathematical files remained byte-identical.

The independent scripts read a frozen producer-results snapshot inside this audit packet and are independently relocatable. No repository edits were made.

## Text-only publication edition: final verdict

PASS. Publication producer packet: session54-tangent-cm-publication. Manifest SHA256: 66465399d9a82a4ec9afc661fc4cb93502b2f632157d08ebf4a1366f21b16eba. Its12 payload files have a verified exact file set. The proof, mathematical checker, rejection checker, manifest verifier, result/rejection JSON files and inherited sources are byte-identical to the independently audited originals; see source-equivalence.json.

The new edition was copied to an unrelated temporary directory and its verifier, certificate and rejection suite replayed successfully in normal and optimized modes. Corrupt, missing, extra and nested-extra controls were all rejected in both modes. All independent matrix, reduced-form and hostile scripts in this text-only audit edition also replayed and matched the prior normal/optimized receipts exactly.

Third-party PDF/article-text copies are absent from both publication editions. Their titles, URLs, precise locators, inspected hashes and short notes are bibliographic provenance, not redistributed payload or runtime dependencies. Original research/audit packets were not changed. This preserves the distinction between a closed executable certificate and explicitly cited theoretical inputs. No repository edits were made.
