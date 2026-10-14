# Complete positional-charge coverage for the 16 cyclic tournament classes at 17

Date: 2026-10-14. Status: exact finite theorem, certified by explicit witnesses in the actual integral quotient rings. This extends the single p=17 control in the [Golay obstruction and norm-prime note](path_charge_golay_obstruction_20261014.md). It does not assert a theorem for larger Fermat primes.

## 1. Finite theorem

Let S select one element from each pair {s,-s} in F_17*. Define the cyclic tournament by u -> v iff v-u belongs to S, and put

\[
 \eta_S=1+\sum_{s\in S}\zeta_{17}^s,\qquad
 R_S=\mathbb Z[\zeta_{17}]/(\eta_S).
\]

**Theorem (finite, certified).** For every such S, the map

\[
 (0,v_1,\ldots,v_{16})\longmapsto
 \sum_{i=0}^{16}i\zeta_{17}^{v_i}\pmod {\eta_S}
\]

from normalized Hamiltonian paths of T_S onto R_S is surjective.

**Certificate proof.** The 256 labelled connection sets partition into 16 multiplier orbits, each of size 16. For one representative of each orbit, the certificate supplies a Hamiltonian path starting at 0 for every residue of R_S. There are 48,160 witnesses in total. The verifier independently checks each path's permutation and directed-edge conditions and evaluates its charge in the exact integral lattice quotient. It also verifies that the 16 multiplier orbits are disjoint and cover all 256 sets.

Multiplication of vertex labels by a nonzero a transports paths and charges via the ring isomorphism induced by zeta -> zeta^a, from R_S to R_(aS). Thus the representative certificates prove the assertion for every S. The earlier [Fermat classification note](fermat_platonic_cyclotomic_orientation_20261014.md) identifies these 16 multiplier orbits with the 16 prime-cyclic isomorphism classes at 17.

This supplies a finite-set carrier for H(T_S)/17 >= |R_S| = disc(T_S), stronger than merely recomputing the already verified counting inequality. It is a proof by finite certificate, not a symbolic construction of a path from an arbitrary residue.

## 2. The actual quotient rings

For every representative, the certificate records multiplication by eta on the integral basis 1,zeta,...,zeta^15 as a matrix M. It also supplies an integer unimodular U and upper-triangular H with M U=H. The verifier checks both the lattice equality and det U=+/-1, independently computes det M, and checks det(I+K)=2^16 det M.

In these 16 cases, the resulting additive presentation has a single nonunit diagonal entry d (or none if d=1), with the class of 1 generating the quotient. Consequently the quotient is the ring Z/dZ, with the following image of zeta. The exact S and every witness are in the JSON certificate.

| Representative index | d=disc | Image of zeta in Z/dZ |
|---:|---:|---:|
| 0 | 1 | 0 |
| 1 | 103 | 8 |
| 2 | 307 | 273 |
| 3 | 103 | 61 |
| 4 | 24617 | 1535 |
| 5 | 647 | 221 |
| 6 | 409 | 5 |
| 7 | 137 | 72 |
| 8 | 14621 | 12634 |
| 9 | 2381 | 1689 |
| 10 | 1361 | 755 |
| 11 | 443 | 370 |
| 12 | 239 | 36 |
| 13 | 613 | 37 |
| 14 | 1871 | 973 |
| 15 | 307 | 9 |

All nontrivial orders in this table are prime except 24617=103*239. That row has

\[
 R_S\cong\mathbb Z/24617\mathbb Z
 \cong\mathbb F_{103}\times\mathbb F_{239},\qquad
 \zeta\longmapsto(93,101).
\]

The 24,617 witnesses in that row cover the **joint product residues**, rather than merely covering each projection separately. This distinction matters: surjectivity onto each factor need not imply surjectivity onto the product. No quotient of composite cardinality was replaced by a finite field.

Every observed quotient here is reduced. That is a checked feature of these 16 cases, not a general claim about Fermat-prime tournament quotients. The preceding norm-prime theorem allows higher primary valuations in general. Every prime factor in this table is 1 modulo 34, as that theorem requires.

## 3. Discovery and independent verification

A randomized insertion routine was used only to discover witnesses. It performed 478,240 trials in 6.39 seconds across the 16 classes, stopping a class when every residue had a witness. The largest quotient needed 260,440 trials. The resulting witness file is approximately 993 KB.

No probability estimate is used in the proof. The verifier does not replay the random search: it directly checks all 48,160 frozen witnesses, the quotient presentations, determinants, multiplier coverage, and prime-factor controls. It reports 160,225 active checks. Fresh normal and optimized Python runs produce byte-identical expected output, taking approximately 0.30 seconds per run in the staging environment.

The verifier reuses the standard-library integral-quotient checking functions from `path_charge_code_obstruction_check.py`, published with the predecessor packet. The discovery builder is not required for replay.

From the repository root:

```sh
python3 04-computation/experiments/fermat_p17_path_charge_check.py > /tmp/fermat_p17_charge_check.json
cmp 05-knowledge/results/fermat_p17_path_charge_check.json /tmp/fermat_p17_charge_check.json
python3 -O 04-computation/experiments/fermat_p17_path_charge_check.py > /tmp/fermat_p17_charge_check_O.json
cmp 05-knowledge/results/fermat_p17_path_charge_check.json /tmp/fermat_p17_charge_check_O.json
```

The four-file packet contains this note, the witness/quotient certificate JSON, the verifier, and its expected-output JSON. The manifest records the existing helper dependency separately.

## 4. Scope of the remaining Fermat question

Combining this finite theorem with the predecessor's p=3 and p=5 controls proves positional-charge surjectivity for every prime-cyclic tournament at the first three Fermat primes 3,5,17. The p=11 Golay counterexample still rules out an unrestricted odd-prime statement.

The structural distinction remains the exact norm-prime restriction: at a Fermat prime p, all residue characteristics divide the norm through primes at least 2p+1 and split completely in the cyclotomic field. This rules out the particular small-characteristic composition collapse at 11. Neither splitting nor this complete p=17 certificate proves surjectivity at 257, 65537, or any further Fermat prime.

The next mathematical target is therefore still a composition-and-order argument in the cyclotomic quotient, rather than a larger brute-force census. One must first identify which coefficient compositions are attainable in an ideal coset, then establish Hamiltonian admissibility for the positional composition. The finite witnesses establish that this combined target succeeds throughout the first nontrivial Fermat sector; they do not yet provide a uniform construction.
