# One transitive order of viewpoints sees every odd residue at 2, 3 and 11

For distinct vertices set \(u\to v\) in \(D(T)\) when an odd number of the Hamilton words of \(T\) put \(u\) before \(v\). Oddness of \(H(T)\) makes this a tournament. The [pair-order theorem](ORDER-PARITY.md) computes it from a single tournament by toggling every pair in every directed cycle's vertex set. The two endpoints can be nonconsecutive in both the word and the cycle. If \(D(T)\) is transitive, its unique order canonically labels \(T\), and **\(T\) is rigid**: every automorphism preserves all the pair counts. That canonical order need not itself be a Hamilton word of \(T\).

**PROVED.** For arbitrary nonnegative \(a,b,c\) and any residue \(r\bmod 2^a3^b11^c\) admitting an odd representative (with \(r\) odd if \(a>0\)), there is a tournament of height \(H\equiv r\) with transitive \(D(T)\). It can be an ordinary member of an all-black doubled node of mass \(2H\). At fixed modulus **all** these residues can occur at one common order and at every sufficiently large order. The nonsingleton parts come from six explicit tournaments on at most nine vertices.

## A fixed six-block integer construction

For the ordinal join \(T\Rightarrow U\), every Hamilton word concatenates words of the factors. Inside either factor every precedence count is multiplied by the odd height of the other; every cross pair is seen forward in an odd number of words. Hence \(D(T\Rightarrow U)=D(T)\Rightarrow D(U)\). In particular transitivity is preserved. Here are the seeds in their own canonical precedence orders, with zero-based outgoing masks. The [replay](precedence_bank.py) checks every pair directly on all full words.

| \(H\) | Vertex count | Outgoing masks |
|---:|---:|---|
| 5 | 4 | 6, 8, 10, 1 |
| 17 | 6 | 30, 60, 56, 48, 32, 1 |
| 37 | 6 | 30, 44, 48, 52, 34, 1 |
| 77 | 8 | 222, 252, 248, 240, 96, 129, 160, 16 |
| 129 | 9 | 254, 508, 504, 496, 480, 448, 384, 256, 1 |
| 751 | 9 | 150, 504, 506, 257, 40, 73, 153, 312, 113 |

First take \(a\ge3,b,c\ge1\); smaller moduli follow by projection. The four unit heights \(5,17,37,751\) generate the entire unit group \(U=(\mathbb Z/2^a3^b11^c)^\times\). To see the 2-primary part, read four bits of a unit \(u\): its sign modulo four; the quadratic character of 2 modulo \(u\), negative on \(3,5\bmod8\); and the quadratic characters modulo 3 and 11. Encode a negative character by one. By the unit group decompositions at prime powers, these are the four independent coordinates of \(U/U^2\).

| Unit height | Mod-four sign | Character of 2 | Character mod 3 | Character mod 11 |
|---:|---:|---:|---:|---:|
| 5 | 0 | 1 | 1 | 0 |
| 17 | 0 | 0 | 1 | 1 |
| 37 | 0 | 1 | 0 | 0 |
| 751 | 1 | 0 | 0 | 0 |

These rows span all four bits over \(\mathbf F_2\), hence generate the 2-primary factor: a finite abelian 2-group quotient equal to its own squares is trivial. For the odd primary factors, \(v_3(5^2-1)=1\), \(5\bmod11\) has order five, and \(v_{11}(5^5-1)=1\). Ordinary prime-power lifting then shows that the order of 5 contains \(3^{b-1}\), 5, and \(11^{c-1}\), supplying every odd primary part of \(U\). Each is cyclic, and a subgroup surjecting on each Sylow subgroup of a finite abelian group equals the whole group. Positive powers suffice since all unit orders are finite.

To hit a nonunit, start with \(129^i77^j\), where \(i=\min(b,v_3(r))\), \(j=\min(c,v_{11}(r))\), with \(v_p(0)=\infty\). The identities \(129=3\cdot43\) and \(77=7\cdot11\) give exactly the required truncated valuations. In each prime-power component the units act transitively on elements of any fixed valuation. Select the correcting unit by CRT and express it by the four unit heights. Ordinal joining produces the resulting **integer** product. For a finite computational selector, breadth-first multiplication on residues returns the actual ordered factors, including at a nonunit. One may choose each word nonempty by inserting a positive power that represents the identity. Remove singleton SCCs from the chosen seed words if present and put singleton padding at their beginning. Their final component is then nonsingleton, and their first component is a singleton, which guarantees distinct converse classes. The same padding realizes all the words at a common sufficiently large order. Rigidity gives ordinary board mass \(H\); all \(2H\) boards in the merged node are black.

At the exponent-61 perfect number \(E=2^{60}(2^{61}-1)\), this proves simultaneous selection to arbitrarily high precision at all three primes for either neighbor. Yet \(E+1=2658455991569831744654692615953842177\) is coprime to **every** height in this bank. Its exact value is therefore outside the multiplicative construction. This exhibits the missing integer condition within one fixed family: every local test at these three primes passes even though no word realizes the target. It says nothing adverse about realizing it by other tournaments. The [line inverse](TARGET.md) gives an independent obstruction for a bounded source across exponents.

## Exponentially many classes with the same parity order

For each \(s\ge4\) let \(J_s\) be the ascending tournament with just its outermost arc reversed. From the [subset description of its Hamilton words](FIBERS.md#a-rigid-example-separating-the-three-counts-at-all-sizes), an interior pair occurs in reverse order in \(2^{s-4}\) words. An endpoint/interior pair occurs in reverse order in \(2^{s-3}\), and the outer pair in \(2^{s-2}\). Thus \(D(J_s)\) recovers the ascending transitive order for \(s\ge5\); for \(s=4\) it transposes the two interior ranks and is again transitive. These are strong blocks of heights \(1+2^{s-2}\).

Take arbitrary ordered words of singletons and the \(J_s\). The strong-component orders distinguish every word. All their precedence tournaments at size \(n\) are isomorphic to the transitive tournament. The exact number **in this subfamily** has generating function
\[
 \frac1{1-z-z^4/(1-z)}
 =\frac12\left(\frac1{1-z-z^2}+\frac1{1-z+z^2}\right).
\]
The coefficient is \((F_{n+1}+c_n)/2\), where \(F_0=0,F_1=1\) and \(c_n\) repeats \(1,1,0,-1,-1,0\). Already at size 36 this gives 12,078,909 different ordinary classes with the same transitive precedence tournament. In each such class that transitive order provides a canonical labeling. Taking \(s=2^j+2, j\ge1\), realizes every finite product of the Fermat numbers \(2^{2^j}+1\) in the same family. Taking the consecutive factors \(j=1,\ldots,k\) gives the integer \((2^{2^{k+1}}-1)/3\) by difference of squares. These identities require no primality assumption on the factors.

## Eleven is the first order with every odd residue modulo 24 in this fiber

**FINITE-EXACT.** Over ordinary isomorphism classes with transitive \(D(T)\), the counts by residue are as follows. Earlier orders embed by adding a single source.

| Order | 1 | 3 | 5 | 7 | 9 | 11 | 13 | 15 | 17 | 19 | 21 | 23 |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 8 | 15 | 0 | 14 | 0 | 13 | 0 | 10 | 0 | 11 | 0 | 10 | 0 |
| 9 | 48 | 0 | 38 | 4 | 48 | 0 | 35 | 2 | 31 | 0 | 43 | 4 |
| 10 | 186 | 2 | 123 | 8 | 170 | 0 | 131 | 4 | 134 | 2 | 141 | 10 |

Thus every count in this fiber through order eight is split modulo four. At nine, 10 of 253 classes are inert. At ten the only missing odd residue is 11. At eleven the height 5363 supplies it, and all eleven previous residues persist. The following actual tournaments have transitive precedence order; the [direct enumeration](precedence_bank.py) recovers that order and each exact height from the outgoing masks.

| Order | Residue | Height | Outgoing masks |
|---:|---:|---:|---|
|10|3|1563|782, 508, 504, 432, 353, 193, 393, 17, 160, 510|
|10|19|1243|838, 508, 504, 369, 161, 449, 400, 265, 16, 510|
|11|11|5363|1806, 1532, 504, 1456, 353, 193, 393, 17, 160, 1534, 500|

If the canonical precedence order must itself be a Hamilton word of the original tournament, the order-ten census leaves just five residues, with counts 74, 76, 16, 59 and 10 in residues 1, 9, 13, 17 and 21 respectively. This stricter physical condition is outside the eleven-vertex threshold theorem.

Neither order-ten score multiset is symmetric under \(d\mapsto9-d\), hence the two ordinary classes in each of its two newly occupied residues consist of the displayed class and its distinct converse. Rigidity makes the unique merged mass respectively 3126 and 2486. This count is within the transitive-precedence fiber. A switching class imposes a different physical condition and can meet all odd residues already at order nine. [THM-5313, nine-vertex switching selector](../../../01-canon/theorems/THM-5313-switching-derivatives-of-hamiltonian-paths-and-cycles.md)

For a direct word enumeration through nine and a bit-parallel cycle enumeration through ten, compile the two accompanying C++ programs and stream one representative per isomorphism class. The independent nine-vertex calculations agree on all bins. The order-ten exclusion uses all 9,733,056 representatives. The integer witness at eleven suffices for the upper bound. [McKay's tournament enumerations and format](https://users.cecs.anu.edu.au/~bdm/data/digraphs.html)

```sh
c++ -O3 -std=c++17 precedence_census.cpp -o /tmp/chanpak-precedence-direct
gentourng -q 9 | /tmp/chanpak-precedence-direct 9
c++ -O3 -std=c++17 precedence_batch.cpp -o /tmp/chanpak-precedence-batch
gentourng -q 10 | /tmp/chanpak-precedence-batch 10
python3 -B precedence_bank.py
```
