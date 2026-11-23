# What exactly two gaps would force

**PROVED:** the two gaps 7 and 21, the multiplicative classification below, and the stated finite anchor certificates. **OPEN:** that there are no further gaps (the bof–Royle conjecture, [Chan–Pak, Conjecture 6.2](https://www.math.ucla.edu/~pak/papers/DominoPermanent16.pdf)). Write \(\mathcal H\) for the actual Hamiltonian-path count spectrum and \(S=\mathbb Z_{>0}^{\mathrm{odd}}\setminus\{7,21\}\) for the conjectured spectrum.

## The holes occupy one triangle ray but both signs

Ordinal join multiplies counts: a Hamiltonian path must first traverse the first block and then the second. In particular a directed triangle realizes \(h\mapsto3h\) for each attained \(h\). This flips the [split/inert sign](../chanpak_signed_square36_20261123/README.md) because \(3\equiv-1\bmod4\). Thus \(7\) is inert, but \(21=3\cdot7\) is split. No classification by this sign alone explains the exclusions.

The obstruction lives at the bottom of **one** \(3\)-multiplication ray:
\[
7\ \text{absent},\qquad 21\ \text{absent},\qquad
63\ \text{attained},\qquad 189,567,\ldots\ \text{attained}.
\]
This ray statement is unconditional. The two holes certify failure of downward division even by an attained factor: \(63\in\mathcal H\) and \(3\in\mathcal H\), but \(63/3\notin\mathcal H\). Conditionally on the full spectrum, no other ray based at a \(3\)-free odd integer has missing levels. In particular \(3n\) would be attained for every odd \(n\ne7\), and \(9n\) for every positive odd \(n\). Neither conclusion can be reversed to prove the original conjecture by division.

## The multiplicative proof obligation is much smaller

A count \(h>1\) is an **atom of \(S\)** if it cannot be written \(h=ab\) with \(a,b>1\) in \(S\). The set \(S\) is closed under multiplication: the only candidate bad composite product is \(21=3\cdot7\), and one factor is absent. Its complete atom list is
\[
\boxed{\{p:\ p\text{ odd prime},\ p\ne7\}\ \cup\
       \{7p:\ p\text{ prime},\ p\ge5,\ p\ne7\}\ \cup\
       \{49,63,343\}.}
\]

Proof. Write \(h=7^a m\) with \(7\nmid m\).
- If \(a=0\), all odd divisors of \(h\) lie in \(S\), so precisely the primes are atoms.
- If \(a=1\), the allowed values require \(m\ne1,3\). For prime \(m\) the only nontrivial split uses 7, hence fails. For \(m=9\), the only possible splits use 7 or 21, hence fail. Every other composite \(m\) either has a prime divisor \(p\ne3\), allowing \((7p)(m/p)\), or is \(3^k\) with \(k\ge3\), allowing \(63\cdot3^{k-2}\).
- For \(a=2,3\), if \(m>1\) split off \(7^a\). The two pure powers 49 and 343 cannot be split without a factor of 7. For \(a\ge4\), split off 49, leaving an allowed factor divisible by \(49\).

All elements of \(S\) factor into its atoms by strict descent. If an atom of \(S\) is attained, every tournament attaining it has exactly one nontrivial strongly connected component, since the count is the product of the component counts and any nontrivial component has height \(>1\). So **the conjecture is equivalent to constructing strong tournaments for just the primes \(p\ne7\) and the twisted primes \(7p\) for \(p\ge5,\ p\ne7\)**: the three exceptional composite atoms already have the explicit witnesses below. Equivalently the least additional missing count, if any, belongs to the prime or twisted-prime family. This does not assert that the same tournament or gadget realizes both members of a pair.

The two families have opposite split/inert tags, \(\chi_4(7p)=-\chi_4(p)\). A hypothetical mechanism that *only* takes ordinal products of smaller counts cannot hit \(7p\): it would have to use the absent 7. But that is not a divisibility obstruction to a strong core. For instance the [signed shell of height \(546973=7\cdot78139\)](../chanpak_kaprekar_decimal_20261123/COMMON-BASE.md) realizes such a product without an ordinal height-7 component. This is a precise test for a constructive approach: it has to make the \(7p\) channels intrinsically, rather than transfer a residue or divide an already existing multiple.

For 49, 63 and 343, exact outgoing masks (zero-based bit positions) are:

| Height | Outgoing rows |
|---:|---|
| 49 | \((114,120,3,117,100,68,4)\) |
| 63 | \((238,236,248,240,163,192,16,64)\) |
| 343 | \((178,168,227,165,46,192,155,16)\) |

Run [`node atom-anchors.js`](atom-anchors.js) for subset DP and checks of strength. No minimal order is claimed for these particular witnesses.

## Two useful consequences, and the boundary

* The deletion of 7 and \(3\cdot7\) changes the arithmetic of factorization, not its residue support. Even restricted to the power-of-seven ray, \(49^3=343^2\) gives two different lengths of atomic factorization. Under the conjecture, no congruence test can detect an isolated hole: every progression containing one of the holes contains arbitrarily large attained odd integers in the same progression (also after refining by any finite additional list of moduli).
* In the inherited [ternary height tree](../session304_collatz_peeling_20261113/ARITHMETIC.md), the child rules \(3h-2,3h,3h+2\) isolate the nontrivial seeds \(19,23,61,63,65\) when the two forbidden nodes are deleted. Seen multiplicatively, 19, 23 and 61 are prime atoms, 63 is the new composite atom on the damaged triangle ray, and \(65=5\cdot13\) is not an atom. Thus there are two complementary notions of a seed: one for backwards ternary navigation, the other for ordinal joining. An arithmetic seed list does not prove that all three child moves are physically available on a common ready class.
