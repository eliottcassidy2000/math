# Independent audit: prime-circulant power-sum boundary

## Verdict

**PASS. No mathematical defect found in Theorems 1 or 2, the arbitrary-exponent constructions, or the numerical examples.** The strongest survivor is the full classification and exact two-prime census, not merely the bounded computational evidence.

For an odd prime with `p−1=2^e m`, tracelessness detects nontrivial multiplier symmetry for every tournament connection set exactly when `m=1`, when `m` is an odd prime power, or when `e=1` and `m` is a product of two distinct odd primes.

For `p−1=2^e q^a r^b`, with distinct odd primes `q,r`, set `B=2^(e−1)q^(a−1)r^(b−1)` and `C=2^q+2^r−2`. The exact labelled counts are

`N_tr=C^B`,

`N_rigid=C^B−2^(rB)−2^(qB)+2^B`.

## Critical proof checks

- **Field and exponents.** The Fourier matrix is invertible over `F_p`; every averaging denominator is nonzero. On odd functions the product `∏(I−A_q)` retains exactly the unit Fourier indices. Periodicity modulo `p−1` covers arbitrary integer exponents, not only those explicitly computed.
- **Independent blocks.** Each prime-order subgroup preserves the lower primary digits and moves transitively through its top digit. Antipodality pairs only the two-primary slices. Thus exactly `B` independent rectangles remain, without hidden carry or compatibility conditions.
- **Boolean lemma.** The mixed difference of `(U+V)^2=1` gives `2ΔUΔV=0`. This proves that one coordinate function is constant in every odd characteristic. Each rectangle therefore has exactly `C` choices.
- **Rigidity.** A tournament multiplier stabilizer has odd order. Every nontrivial subgroup of the cyclic multiplier group contains a unique prime-order subgroup. Breaking all those subgroups excludes every stabilizer, including higher prime powers. The nonrigid two-prime population is exactly the union of the order-`q` and order-`r` invariant populations, whose sizes are `2^(rB)` and `2^(qB)` with intersection `2^B`.
- **Arbitrary primary exponents.** All selectors in 5A–5C preserve the subgroup claimed for each summand, while the displayed increments `q_j^(a_j−1)` break every prime-order subgroup. No unused-coordinate or quotient-lift stabilizer remains. The constructions exhaust the complement of the positive regimes.

## Exact computational scope

An independent rerun of `verify_boundary.py --bound 1000` produced output byte-for-byte identical to `verification.json`. Both files have SHA-256

`c95995037f145a9ae286186e314bf7c407074049a438dee0e4c3c46ee324a146`.

The script's coverage and reconstruction were inspected. It exhaustively enumerates every tournament connection set for odd primes at most 31, completely counts all `2^30` sets at 61 by an exact `15+15` meet-in-the-middle calculation, and tests one explicit construction at each of the 62 bad-regime primes at most 1,000. The last operation is not an exhaustive census of sets at those primes.

Fresh calculations from the stored residue sets confirm every unit-index sum and the trivial multiplier stabilizer at 127, 151, 271, and 331. Formula evaluations are:

| Prime | B | C | Traceless | Rigid traceless |
|---:|---:|---:|---:|---:|
| 61 | 2 | 38 | 1,444 | 360 |
| 127 | 3 | 134 | 2,406,104 | 308,448 |
| 151 | 5 | 38 | 79,235,168 | 45,648,000 |
| 271 | 9 | 38 | 165,216,101,262,848 | 130,031,594,956,800 |

The counts at 127, 151, and 271 are proved formula evaluations, not exhaustive enumerations. Supplemental construction checks at 1801, 2311, and 2731 also passed; they are summarized in [audit_extra_examples.json](audit_extra_examples.json) and are not needed for the proof.

## Scope caveats

The correction at 331 is right: `330=2·3·5·11` is squarefree. The complete criterion confirms minimal failure at 61 overall and at 127 among primes congruent to 3 modulo 4; 211 remains a valid later example.

THM-4724 was checked directly: section 3 leaves general uniform-dilation-rule non-existence OPEN, whereas section 4's last sentence overextends Theorem W to all local rules. The new note correctly confines that source claim to its derivative-family scope and does not certify Theorem W's separate arithmetic proof. No general local-rule theorem or tournament/graph bijection follows here.

The graph normalization warning is valid: at `p=7`, `E={±1,±2}` has stabilizer `{±1}` and `σ_2(E)=3`, a nonsquare. A same-square-class target repairs normalization only when a suitable coefficient exists. Its stated even-graph status also holds: its complement is a 7-cycle, and its dihedral automorphisms reverse an even number of its edges.

Optional clarification only: explicitly choose the “least unit index” from `1≤k<p−1`. This is the convention already implicit in the Fourier argument and is not a defect in either theorem.

## Publication revision

**Audit status: PASS retained.** Publication-layout changes replace file references with verified relative links; the mathematical statements, proof, construction formulas, and computational certificates are unchanged. The approved proof now has SHA-256

`8950f93635fe4d583a7f021f6e8d2ae6d1de1f991b88538425ab79a510612140`.

The least-index interval recommended above is explicit in the approved proof. Its repository reference resolves from the publication directory to the tracked THM-4724 file; its local links resolve to approved packet payloads.
