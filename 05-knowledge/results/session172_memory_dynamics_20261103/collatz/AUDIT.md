# Independent audit of the S172 Collatz memory packet

Verdict: substantive PASS. No theorem or numerical correction required in the inspected NOTE.md/check.py. Minor research-target clarification recommended below.

Reproduction: normal and optimized runs of check.py agree byte-for-byte with each other and the submitted certificate.json. The independent checker collatz_audit_independent.py uses explicit finite sums for offsets, inverse residues rather than deriving words from representatives, and actual orbit iteration around the predicted threshold boundary. It checks 15,377 boundary representatives across all words through length 12, including 11 nondescending small exceptions; explicit surviving words through depth 18 (7,495 survivors); and all-one ascent witnesses to depths 65,128,257,1024. Its normal/optimized results agree.

Mathematical checks:

- Literal seed profiles, 17/19 readouts, ordinal multiplication and canonical ordered strong components justify the tape construction. Independent full-profile DP in this audit worker's spectrum packet already verifies these same seed identities by a different algorithm.
- Raw parity words exhaust all 2^k classes, so the warning against restricting them to a proper code is necessary and correct. A separate injective encoder preserves all words.
- Fixed-weight offset injectivity follows from the inverse residue identity; equality of full affine maps forces equal slopes and offsets. The conditional entropy calculation is exact for uniform residues.
- Integer thresholds use strict descent correctly: Dn>a is equivalent to n>floor(a/D). The minimum over contracting prefixes and the finite exception count are correct, including residue zero's least positive representative 2^k.
- The finite autonomous observer theorem correctly applies even on positive integers: equal states imply every parity prefix equal and therefore equality modulo every 2^k.
- Survivor DP counts every prefix barrier, rather than checking only terminal weight. Finite-horizon actual and coefficient failures differ by finitely many exceptional representatives, so their natural densities coincide.
- The exponential binomial upper bound and passage to density-one finite strict stopping are valid. Their quantifiers do not imply convergence for every positive integer.
- All-one cylinders prove that no bounded descent horizon works, and the claimed exact affine ascent formula is correct.
- The repeated-block rich probe really has 2^k outputs, but costs 2^k−1 storage blocks and as many extra read vertices. It does not violate the local binary-query lower bound.

Editorial refinement: Section 6 should explicitly restrict the certificate-tree domain to n>1 and treat n=1 as the induction base. Otherwise its literal instruction that all finite exceptions must have verified strict descent cannot close the node containing 1. This concerns the proposed future workflow, not the proved results. Sent to author.

Source verification: directly opened the cited IMPAN Terras record, which confirms author, title, Acta Arithmetica 30 (1976), pages 241–252 and DOI 10.4064/aa-30-3-241-252. Directly opened the cited Lagarias stopping-time exposition, which attributes the density-one theorem to Terras and independently Everett and discusses the coefficient/actual distinction. References support the historical claims; the packet's proof is self-contained. Replace temporary source paths by repository links during integration.

Publication integration: the author applied the n>1 / n=1 base-case clarification before publication; source paths now use repository links.
