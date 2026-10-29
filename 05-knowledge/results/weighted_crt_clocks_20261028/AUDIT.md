# Same-agent proof and computation audit

No independent referee or Lean verification is claimed. All statements have
their declared source, quantifiers and failure boundary in the topic notes.
The initial repository source is5ffa9ee1ecce; external results are explicitly
credited and the computational dependency hashes are in the receipt.

Run from the repository root:

```text
python -B 05-knowledge/results/weighted_crt_clocks_20261028/audit.py
python -B agents/check_docs.py
```

The first command reruns all five checks normally and under `-O`, compares
the resulting JSON byte-for-byte, checks every relative Markdown link, and
adds independent ring, necklace and weighted-variance checks. It writes
[audit_receipt.json](audit_receipt.json). Assertions are not used as gates.

| Lane | Exact universe | Independent route / hostile |
|---|---|---|
| CRT | All cyclic no-11 words at n=9,11,13,15,21,23,25 | Ring coordinate uniqueness, full reflection/translation partitions, necklace formula |
| Ring algebra | Odd n=3..301; period55 prime-power chart | Norm determinant, Smith gcd, root equation, exact order55 mod121; no all-word census at55 |
| Unit clocks | u=(5+sqrt29)/2 modulo11^h,13^h, h=1..4 | Exact proper-divisor order tests; all-h theorem uses binomial lifting |
| Cycle index | k=0..3; literal groups and all subsets through k=2 | Wreath recurrence versus directly listed centralizer permutations |
| Canonical carrier | Every binary array and every ordered tuple of up to3 ports through k=2 | Literal group orbits; actual strong-atom extensions; global converse; shared-phase hostile |
| Hostile pair | Literal twelve-vertex AAA/BAA and every4096 subset | Full exact automorphism search, subset orbits versus cycle index |
| A000568 | n=0..20 | New operator organization versus inherited odd-partition sum; no supplied count used as generator input |
| Signed affine | 1585 words, r<=5,A<=12 | 70612 literal odd residue tests when A<=8; exact signed cycle guards |
| Eta | Five prime levels, q^40 | Logarithmic recurrence versus independent truncated products and inverse convolution |

Type audit: a weighted counting polynomial is not a faithful tuple classifier;
rotation period differs from translation period; a prime-power residue ring
is not automatically a field; word reflection is not a ring automorphism;
Haar probability is not a positive-integer basin measure; the accelerated,
Terras and standard Collatz clocks are distinct; a one-source periodic affine
map is not an identity on two sources. All-size conclusions are supported
by the supplied proofs, never by a finite table alone.

The proofs preserve the inherited 81/27 control and the exact loop bound,
without interpreting that bound as an optimal attachment carrier. No
historical near-miss has been promoted to canon, and no global theorem ID
has been reserved. The result index routes to this packet rather than
inserting session prose above any maintained warning.
