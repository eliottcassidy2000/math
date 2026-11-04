# Exact divisor-network interfaces and fourth-power clauses

**Status:** PROVED elementary interface lemmas and finite, replayable conditional clauses. Odd perfect-number existence remains OPEN. No originality is claimed for sigma chains, strong connectivity, or multiplicative cut identities.

## Inheritance and source boundary

- The relevant repository carrier is `07-reflections/lrc14-perfect-number-packet-merge-codex-s174.md`: factorization and prime/composite status are necessary sidecars; a product-only quotient loses the target. Here the extra indispensable sidecar is the **exact exponent of every supplier**.
- The repo's LRC route-purity machinery motivates exact interfaces but no LRC theorem is used or implied.
- Ochem–Rao, *On the number of prime factors of an odd perfect number*, https://www.lirmm.fr/~ochem/opn/opnf.pdf, gives the classical Euler form and factor-chain methodology. The present proofs below are self-contained apart from Euler's classical form.
- Zhu–Li, *Valuation Networks, Order-Filtered Primitive Excess, and Computer-Assisted Obstructions in the Odd Perfect Number Problem*, https://zenodo.org/records/22005604, posted August 2026, already advertises strong connectivity and exact cuts. Its landing-page status explicitly says unreviewed and not a proof of nonexistence. We independently prove the elementary facts needed here; no computational bound from that preprint is imported.

Concept board: exact component valuations; cut interfaces; strong connectivity; exponent-sensitive fourth-power clauses; hostile exponent changes. Anchor = odd-perfect constraints; niche = exact partial-interface certificates; wildcard = transfer of prime/composite and supplier-label sidecars from the repo's LRC arithmetic carriers.

## 1. Strong connectivity, with its sharp limitation

Let N>1 be perfect, and let P be its prime support. Give p in P its actual exponent a_p=v_p(N). Draw p→r when r is in P and r divides sigma(p^a_p). Omit any divisor-sum prime outside P (for an odd perfect number the only possible such prime is 2).

**Lemma 1 (PROVED).** This graph is strongly connected.

Proof. Otherwise it has a nonempty proper sink set S. Put A=product over S of p^a_p and B=N/A. Multiplicativity and sigma(N)=2N give sigma(A)|2AB. The sink condition says every odd prime factor of sigma(A) belongs to S; for even N all its prime factors belong to S (including 2 if it occurs). Because A and B are coprime, in every case sigma(A)|2A. Since A>1, sigma(A)>A. Consequently the positive integer 2A/sigma(A) is less than 2, and equals 1. Thus sigma(A)=2A. Multiplicativity now gives sigma(B)=B, impossible for B>1. Therefore no proper sink exists. QED.

The graph has no self-loops because sigma(p^a)=1 modulo p. Thus every prime lies on a directed cycle. However graph structure alone does not recognize perfection: N=120 has a strongly connected prime graph while sigma(120)=360=3N. The integer valuation budgets are essential.

## 2. Exact cut and partial-completion interfaces

Let N be odd perfect. Euler's form is N=q^alpha m^2, with q,alpha congruent to 1 modulo 4 and gcd(q,m)=1. For a nonempty proper subset S of its prime support, let A be the corresponding full prime-power block. Let

- O be the product of r^{v_r(sigma(p^a_p))} over p in S and odd r outside S;
- I be the analogous product over p outside S and r in S;
- epsilon=1 if q is in S, and 0 otherwise.

**Lemma 2 (PROVED).** O/I = sigma(A)/(2^epsilon A). Moreover O and I are coprime, since their prime supports are disjoint. Consequently, with g=gcd(sigma(A),2^epsilon A),

O=sigma(A)/g, and I=2^epsilon A/g.

Proof. At each vertex r the total incoming r-adic valuation is a_r, by sigma(N)=2N. Cancel the internal edge factors in sigma(A)/A. The only even factor is the unique 2 from the Euler component. This proves the ratio. Reduced fractions then prove the formulas. QED.

If epsilon=0, O>I; if epsilon=1, O<I because A is a proper unitary divisor and sigma(A)/A<2. These directions do not imply acyclicity: all non-Euler vertices have positive logarithmic divergence while the Euler vertex absorbs it, but every cut still has edges both ways.

**Lemma 3 (PROVED; useful for partial candidates).** Suppose A is a product of components whose exponents have been fixed exactly in a proposed odd perfect N=AB, with gcd(A,B)=1. Set

R = sigma(A)/gcd(sigma(A),2A).

Then R divides B. In particular gcd(R,A)>1 rejects the proposed exact assignments. If A contains the Euler component, B is a square, so

C = product over r|R of r^{2 ceil(v_r(R)/2)}

divides B. Therefore sigma(A)*sigma(C) <= 2*A*C. A strict reverse inequality rejects the branch; equality forces B=C, and that complete candidate can be tested directly.

Proof. sigma(A)|2AB and cancellation of gcd(sigma(A),2A) imply R|B. The square ceiling follows from the parity of exponents of B. Abundancy is multiplicative on coprime integers and strictly increasing along proper divisibility: each local factor 1+1/p+...+1/p^e strictly increases with e, and a new prime contributes a factor greater than 1. Apply this to C|B. QED.

**Critical scope:** exact assignments cannot silently be replaced by lower bounds. If only p^a|N is known, sigma(p^a) need not divide sigma(N). For example sigma(3^2)=13, whereas sigma(3^4)=121; the edge 3→13 disappears when the exponent is raised. Likewise 3 divides sigma(13^2) but not sigma(13^4).

## 3. Compact fourth-power clauses

These are conditional results about genuine odd perfect numbers, not a finite search for all of them.

**Proposition 4 (PROVED).** If 17||N, then 3^4|N or 13^4|N.

Proof. Exponent 1 makes 17 the Euler prime. sigma(17)=18 forces v_3(N)>=2 and that valuation is even. If it is not at least 4, it equals 2 and sigma(3^2)=13 forces 13|N. If v_13(N) is not at least 4, it equals 2. But sigma(13^2)=183=3*61 adds one more 3 to the two already supplied by sigma(17), giving v_3(sigma(N))>=3>v_3(N)=2. Contradiction. QED.

**Proposition 5 (PROVED).** If 5||N, at least one of 3^4,13^4,61^4 divides N.

Proof. If all three fail to divide N, sigma(5)=6 successively forces 3^2||N, sigma(3^2)=13 forces 13^2||N, and sigma(13^2)=3*61 forces 61^2||N. But sigma(61^2)=3783=3*13*97. The components 5,13^2,61^2 supply three factors of 3, exceeding the fixed exponent 2. QED.

**Proposition 6 (PROVED).** If 53||N, then 3^4|N, because sigma(53)=54 contains 3^3 and the non-Euler exponent of 3 must be even.

The script additionally generates exact clauses for every prime q<=97 with q=1 mod4 under the hypothesis q||N. Their disjunctions are:

| q | At least one p^4 divides N, for p in |
|---|---|
| 5 | 3,13,61 |
| 13 | 3,7,19,127 |
| 17 | 3,13 |
| 29 | 3,5,13,31 |
| 37 | 3,13,19,127,5419 |
| 41 | 3,7,13,19 |
| 53 | 3 |
| 61 | 3,7,13,31,331,5233 |
| 73 | 3,7,37,67 |
| 89 | 3,5,13,31 |
| 97 | 3,7,13,19,127 |

**Finite certificate logic.** Suppose every listed fourth power is absent. Start with q^1. Every newly forced ordinary prime then has exponent exactly 2, unless its demanded valuation exceeds 2, which already contradicts the hypothesis. Recompute all divisor sums and residual budgets after each insertion. The recorded termination is an excess incoming valuation, abundancy greater than 2, or a newly forced valuation at least 3 (hence at least 4). Each record has exact integers and component factorizations. This proves the displayed disjunctions; it does not allow freezing exponent 2 without the explicit absence hypothesis.

## 4. Reproduction and next useful operation

Run `python3 verify_networks.py`. It writes `network_certificates.json` and checks:

- 500 sigma values against independent direct divisor enumeration;
- all eight nontrivial unitary cuts of the positive controls 6,28,496,8128;
- the strongly-connected nonperfect hostile control 120;
- exact edge disappearance when exponents change;
- all eleven clauses above, with each transition exposed as a JSON record.

No installs, floating point, probabilistic primality, or external computation is needed. Trial division certifies all factorizations used. Odd-perfect positive controls are unavailable because no odd perfect number is known. The even-perfect controls test only the general algebra and graph statements.

The next operation is to combine these clauses with independently proved residue/order constraints while branching on exponent labels. Merely increasing the finite prime cutoff or repeatedly forcing exponent 2 would reproduce a known limitation: unresolved higher exponents can remove old edges. A universal proof would need a terminating argument covering those exponent branches; none is supplied here.
