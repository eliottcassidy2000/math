# Independent audit of the twin-prime Berggren packet

**PASS, 26 October 2026. No proof or implementation defect found.** The
producer files were read and replayed without editing them. This is an
independent review by the arithmetic/radii lane, not producer self-certification.

Reviewed frozen files:

| File | SHA-256 |
|---|---|
| `../tree/RESULT.md` | `748c64f37d09b24cb425cde4a1c193f4f293395837f789460fbc9f5dd6a8be10` |
| `../tree/check_tree.py` | `63adaee9b46c898e26fd1229552498c386cedfc2a06f290790d249157bf022c4` |
| Both producer replay receipts | `ea1ebee807b5fb0e09400bfa66634178b146ee09d94252c207bad39c8c7c233c` |

The producer runs in ordinary and optimized Python each pass **198,535
explicit gates** and match the frozen normal receipt byte for byte. The
independent [audit_tree_independent.py](audit_tree_independent.py) passes
**2,892 gates**, also byte-identical under `-O`; both independent receipts
have SHA `b415bdd91d896366fb7131998d3bb506a42c5190292c674ce859908365b82bbd`.

## Proof review

1. **Two-root normalization and total multiplicity.** For a primitive triple
   in odd/even-leg orientation, `c-b,c+b` are odd and coprime: a common
   divisor divides both c andb. Their product is `a²`, so both are squares.
   The inverse roots recover `m-n,m+n`. The height identity is relative to
   the odd leg, while the U ray uses height from the even leg. They are
   correctly distinguished. At gap two, total factor count `Omega(a)=2`
   forces both nonunit factors prime; the root excludes the unit factor.
   The hostile63 correctly separates total and distinct factor counts.

2. **All prime-leg edges, not only a selected sample.** Primality of
   `(m-n)(m+n)` forces `m-n=1`. I independently expanded all three branch
   odd legs: U gives `p+2`, A gives `p(2p+1)`, D gives `p(2p-1)`.
   For p>=3 the latter factors are both nonunits. Thus every prime-to-prime
   tree edge is on U and exactly a twin edge. The argument is exhaustive.

3. **Line graph and root offset.** The source edge indexed by
   `j=(p-3)/2` has words `U^j,U^(j+1)`; its image has word `D^(j+1)`.
   Edge-sharing iff consecutive indices proves both directions of adjacency
   transport. The image is the D ray with its root removed before filtering.
   It is a graph isomorphism of the stated induced subgraphs after taking
   the line graph, not a bijection of the original source vertices. Source
   isolated prime vertices disappear; the producer states this loss.

4. **Global component boundary.** A three-prime run p,p+2,p+4 must include
   the only prime divisible by3, hence p3. The next candidate9 fails.
   Consequently exactly one pair of twin edges meets, and exactly one edge
   exists in the target graph. This finite-looking conclusion rests on a
   complete modular proof; it does not assume there are infinitely many
   twin pairs or isolated target vertices.

5. **Semiprime and perfect-square fibres.** For distinct odd primes p,q,
   the only coprime factor splits of pq are `(1,pq)` and `(p,q)`. Both
   produce valid primitive triples, and the nontrivial gap is `(q-p)²/2`.
   The three aliquot formulas use the correct prime/square shape. In the
   lower-square map, perfect image six occurs at p5 and triangle(35,12,37).
   Allowing either square adds exactly the predecessor triangle(15,8,17),
   whose upper square is25. Thus the common square25, its opposite flags,
   and the distinct centres4 and6 are all correctly retained. This is
   distinct from the radii lane's fixed incoming collision at p3.

## Independent computational path and sources

The independent checker enumerates primitive triples directly from sides
with c<=2000, without generating Euclid pairs or tree words. All319 triples
obey the root-chart equations; their prime-leg branches obey the claimed
classification. It recovers194 distinct-semiprime odd-leg fibres a<=999 by
factoring `a²=(c-b)(c+b)`, independently of the producer's Euclid solve.
For61 twin pairs with lower prime below2000, repeated D matrices verify the
root shift, and every pair of edges is checked for adjacency equivalence.
Direct divisor sums verify the perfect-square anchor. No producer import
is used and every gate survives optimization.

The [Tripathi primary paper](https://www.fq.math.ca/Papers1/46_47-4/Tripathi.pdf),
introduction and Theorem2, was independently read. The claimed classical
fixed-leg specialization is accurate. I also read
[Janičková–Csókási, §2 and Lemma2.4](https://arxiv.org/pdf/2304.05230): their
matrices A,B,C agree with this packet's U,A,D. The source's extracted Euclid
statement omits an explicit opposite-parity condition; the producer records
this and retains the necessary condition. No reliance on that omission is
present. The unique-six theorem was independently proved and audited in
the previous packet and its scope is imported correctly here.

Reproduce the additional controls from this directory:

    python3 -I -B audit_tree_independent.py
    python3 -I -B -O audit_tree_independent.py

Later status/link-only promotion of the producer is outside these byte hashes
but does not alter the mathematical audit; any proof change needs its own
recorded delta or review.
