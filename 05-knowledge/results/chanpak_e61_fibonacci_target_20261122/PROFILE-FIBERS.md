# Substituting strong components in the exponent-61 lower classes

**PROVED, with FINITE-EXACT small component counts.** At least **414,387,118,080** nonisomorphic strong tournaments of order 91, and at least **30,453,618,843,648,000** of order 93, have \(E_{61}-1\) Hamiltonian paths. Within each of these two families the entire polynomial counting ordered weak path covers is identical. Every one-arc neighbor of every member has Hamiltonian-path count divisible by three. Consequently none of these lower classes can flip to \(E_{61}+1\).

## Which strong components have the same profile?

Write \(Z_C(x)\) for the number of ways to partition all vertices of \(C\) into \(x\) labeled directed paths, allowing empty paths. Write \(H=Z_C(1)\) and \(Q=Z_C(2)\). The [91- and 93-vertex lower constructions](../chanpak_e61_compact_volchkov_20261122/STRONG-93.md) use the cyclic wrap \(U\to a\to b\to U\), with \(a,b\) singletons and \(U\) an ordinal sum of strong components. Their count \(H(U)+Q(U)\) permits independent substitutions by components of the same order and the same \(H,Q\).

The table gives the exact number of strong isomorphism classes in each indicated \((n,H,Q)\) fiber through order nine. At order ten the table counts only the supplied component and its converse. All the indicated choices have the same entire \(Z_C(x)\) within their row.

| Order \(n\) | \(H\) | \(Q\) | Choices | Family order |
|---:|---:|---:|---:|---:|
| 3 | 3 | 12 | 1 | 93, twice |
| 4 | 5 | 32 | 1 | 93 |
| 5 | 13 | 104 | 1 | 93 |
| 6 | 31 | 328 | 2 | 93 |
| 8 | 217 | 4096 | 2 | 91 |
| 8 | 331 | 5368 | 4 | 91, 93 |
| 9 | 305 | 9268 | 2 | 93 |
| 9 | 1321 | 24800 | 2 | 93 |
| 9 | 1525 | 26420 | 2 | 91 |
| 9 | 1599 | 27804 | 2 | 91 |
| 10 | 451 | 21136 | \(\ge2\) | 93 |
| 10 | 1057 | 38656 | \(\ge2\) | 93 |
| 10 | 1661 | 50840 | \(\ge2\) | 91 |
| 10 | 3963 | 94224 | \(\ge2\) | 91 |

For the census generate all nonisomorphic strong tournaments at orders 3, 4, 5, 6, 8, 9 using *gentourng* with strong connectivity required. Their counts are 1, 1, 6, 35, 6,008 and 178,133. Dynamic programming on the last vertex gives \(H(C[S])\) for every subset. Summing \(H(C[S])H(C-S)\) gives \(Q\). For the surviving classes calculate \(a_j\), the number of partitions into \(j\) unordered nonempty directed paths. The path containing the least remaining vertex specifies the first block uniquely and supplies a factor \(H(C[S])\). Then \(Z_C(x)=\sum_{j=1}^n a_j(x)_j\), where \((x)_j\) is a falling factorial. This compares every slot count exactly. Each supplied order-ten component has a different out-degree multiset from its converse, and reversing every path preserves each \(a_j\). The replay additionally uses a degree-restricted isomorphism test for the converse of each seed.

## A local collision beyond reversal

For the four strong order-eight classes at \((H,Q)=(331,5368)\), the common coefficient vector \((a_1,\ldots,a_8)\) is
\[
 (331,2353,4662,3817,1500,296,28,1).
\]
There are two converse pairs. One pair has the sorted vertex-deletion heights
\((23,51,81,81,95,95,95,137)\); the other has \((37,69,73,75,81,93,109,129)\).
Within each pair the two degree sequences are \((1,3,3,3,4,4,5,5)\) and \((2,2,3,3,4,4,4,6)\), exchanged by reversal. The deletion heights are preserved by reversal and distinguish the pairs. Thus even on eight vertices, every unmarked slot evaluation can agree while the counts after a single vertex deletion distinguish the classes. The checker prints all four upper triangles.

## Why the multiplicities multiply

The unique maximal proper nontrivial module of the wrap is \(U\). Its strong components have an intrinsic linear order. For fixed multiplicities \(m_t\) of profile types, \(r_t\) eligible classes per type, and \(L=\sum_t m_t\), the exact number of permitted wrapped classes is
\[
 \frac{L!}{\prod_t m_t!}\prod_t r_t^{m_t}.
\]
Indeed a profile word determines the positions of all types, after which each position independently receives a class. Distinct resulting words of isomorphism classes cannot be identified. Both constructions have five two-choice component types and one four-choice type. Each of these six types occurs once. The multiply occurring types, singletons and at order 93 the triangle, each have just one choice. Hence the multiplier for each construction is \(2^5\cdot4=128\), giving the bounds \(128\cdot41!/35!\) and \(128\cdot34!/(24!\,2!)\).

Within each order, the common \(Z_U(x)\) follows by multiplication over ordinal cuts. For a fixed wrap the number of fillings of \(k\) nonempty runs through \(U\) is the \(k\)-th finite difference of \(Z_U\) at zero. Thus the wrapped polynomials also agree. Finally, all substitutions preserve the two disjoint three-divisible components and the vanishing cyclic-root weight used by the [mod-three module-path theorem](../chanpak_decimal_kernel_20261122/MODULE-PATH.md). This gives the asserted restriction on every one-arc neighbor.

Compile the [census worker](profile_fibers.cpp) and run the [family calculation](profile_fibers.py) with *gentourng* installed and available on the executable search path.

~~~sh
c++ -O3 -std=c++17 profile_fibers.cpp -o /tmp/e61_profile_fibers
python3 -B profile_fibers.py /tmp/e61_profile_fibers
~~~
