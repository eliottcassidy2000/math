# No transitive inflations of the exponent-three marked classes can cross exponent 61

**PROVED using a finite residue calculation.** Take **any** tournament pair at counts 27 and 29 differing at one specified arc. Replace each other vertex, independently, by a transitive tournament of **any positive order**, preserving the between-block orientations. Leave the two arc endpoints singletons. The resulting marked pair cannot have counts \(E_p-1,E_p+1\), in either order, at **any exponent \(p\equiv1\pmod{12}\)**. In particular this rules out all such lifts from exponent 3 to either 13 or 61. There is no upper bound on the new tournament size or any block size.

The inherited marked-pair grammar strips exterior singleton runs. Their replacements remain whole exterior transitive tournaments and have path count one. The remaining four primitive marked isotypes are two on six vertices and two on seven vertices, with both completions strong. From the complete lower catalogue, take every strong class of height 27 and each arc giving 29, modulo lower automorphisms. The strong cap at height 27 is seven. The following rows and arc labels are representatives; reversing the comparison of final heights is also allowed in the sieve. [Marked series and strong cores](../chanpak_arc_fiber_grammar_20261122/README.md)

| Lower rows | Edited arc |
|---|---|
| 2, 32, 11, 19, 7, 29 | (1,2) |
| 32, 5, 9, 3, 15, 30 | (1,5) |
| 2, 4, 17, 7, 75, 31, 47 | (2,3) |
| 2, 4, 17, 7, 75, 31, 47 | (3,4) |

## Why a finite table covers unbounded widths

In any fixed tournament quotient, compress each maximal contiguous run in a transitive block to its block label. An ordered skeleton visiting a block of width \(w\) exactly \(r\) times has
\[
      J(w,r)=r!\,{w\brace r}=\sum_{j=1}^{r}(-1)^{r-j}\binom rj j^w
\]
fillings in this block when \(w,r\ge1\). Assign each labelled vertex to a run, use every run, and read each fiber in increasing order. Choices of fillings in different blocks are independent. For a prime \(q\), \(J(w,r)\equiv0\pmod q\) when \(r\ge q\). Otherwise Fermat's little theorem makes \(J(w,r)\bmod q\) depend only on \(w\bmod(q-1)\). Consequently each completed Hamilton count modulo \(q\) depends only on the residues of the individual positive widths modulo \(q-1\). This applies **even to cyclic ports**, for which the unbounded skeleton language may have arbitrarily many visits. Widths 1 and \(1+k(q-1)\) are indistinguishable by this modulus.

For \(q=3,5,7,13\), all four periods divide 12. The pair of desired counts modulo those primes likewise depends only on \(p\bmod12\). It therefore suffices to enumerate the \(12^d\) width vectors in \(\{1,\ldots,12\}^d\) for each seed, where \(d\) is its number of nonendpoint vertices. For each residue vector, require both endpoint residues in **one common orientation** at all four primes.

The modular replay first counts skeletons modulo \(q\) by their last symbol, the multiplicities \(0,\ldots,q-1\) of the port symbols, and a two-bit occupancy mask for the edited endpoints. At each addition sum over permitted previous symbols. Retain the states visiting both endpoints and every port. Apply the one-coordinate Stirling transform \(J(w,r)\) along each port dimension to obtain the modular table for that completion. For the small inflations of total excess width at most three, the program also constructs the literal graphs and compares both endpoint counts by direct subset recursion.

| Exponent modulo 12 | Residue boxes passing both heights for **each** six-vertex seed | Boxes for **each** seven-vertex seed |
|---:|---:|---:|
| 1 | 0 | 0 |
| 2 | 0 | 0 |
| 3 | 1 | 1 |
| 4 | 0 | 0 |
| 5 | 1 | 0 |
| 6 | 0 | 0 |
| 7 | 0 | 0 |
| 8 | 0 | 0 |
| 9 | 1 | 0 |
| 10 | 0 | 1 |
| 11 | 0 | 3 |
| 0 | 1 | 0 |

The table entries count width residue vectors prior to any actual height constraint. Zero certifies impossibility. In phase 3 the unique passing vector consists entirely of ones modulo 12, since the original widths already give an exact exponent-three crossing. At exponent 61, the six-vertex types are ruled out by 3, 7 and 13 alone. For either seven-vertex type, **every proper subset** of these four moduli admits at least one residue vector; the joint test rejects them all. This is a simultaneous constraint on the two counts across the same physical arc. It does not presume that independently matching scalar witnesses give a flip.

## Exact arithmetic on the affine lifts

As a separate absolute calculation, prescribe widths \(1+29a_i\) at exponent 61, for arbitrary nonnegative integers \(a_i\). For a strong candidate the order bound 157 gives \(\sum_i a_i\le5\). Strongness is preserved by the inflations, since the original strong quotient gives routes between all blocks and routes leaving and returning to any given block. There are 756 normalized seed/width presentations, counting both zero-growth base cases and nonzero growth across the four seeds. In 42 the growing ports contain a cyclic triangle. All three widths are at least 30, and three rounds of visits in cyclic order give at least \(J(30,3)^3>E_{61}+1\) paths on these blocks alone. Insertion into a directed path preserves the lower bound on passing to the full tournament. For each of the remaining 714 cases, between consecutive ordinary vertices the transitive ports occur in order and no port repeats. The finite complete skeleton recurrence therefore terminates, and summing products of exact \(J\)'s gives both heights. Exactly 28 of these 714 retain the original signed difference of two, none has difference minus two, and neither requested endpoint occurs. The common marked decomposition is retained throughout; these totals count width presentations.

Run the following in this directory. The first command recovers the primitive marked classes and the absolute affine table. The remaining commands check all four infinite-width exclusions.

~~~sh
python3 -B seed_lifts.py
c++ -O3 -std=c++17 seed_width_sieve.cpp -o /tmp/chanpak_seed_width
/tmp/chanpak_seed_width 1 2 2 32 11 19 7 29
/tmp/chanpak_seed_width 1 5 32 5 9 3 15 30
/tmp/chanpak_seed_width 2 3 2 4 17 7 75 31 47
/tmp/chanpak_seed_width 3 4 2 4 17 7 75 31 47
~~~

[Exact affine calculation](seed_lifts.py) · [Modular width sieve](seed_width_sieve.cpp)

## Every exponent-61 straddle needs more than four transitive blocks around its arc

**PROVED using finite-exact input.** In any hypothetical strong central pair of counts \(E_{61}-1,E_{61}+1\), keep its marked endpoints as singleton blocks. The other vertices cannot be partitioned into four or fewer nonempty transitive modules. Here each module is homogeneous towards every vertex outside it in the shared missing-arc graph. This concerns the internal structure of the unknown lower witness. It applies regardless of the path counts or signed gap of the quotient obtained by contracting the blocks.

Contract such blocks while keeping the two endpoints. If both expanded completions are strong, both quotient completions are strong. Enumerate unoriented missing-arc quotients by labelling the endpoints 0,1, assigning all the other arcs, and identifying codes under permutations of the endpoints and independent permutations of the other vertices. Discard a code if either completion fails strong connectivity. Every permissible quotient occurs once. The archive's exterior-singleton theorem and Busch bound ensure the expanded central order is at most 157. Write \(d\) for the number of other blocks, so that the widths sum to at most 155.

For each quotient use the four-prime period-12 test above. For each passing residue box write the widths as \(w_i=b_i+12t_i\), where \(1\le b_i\le12\) and \(t_i\ge0\), and enumerate those within the central order bound. Test the two counts jointly at the extra primes 11, 17, 19 and 23, stopping once every candidate is gone. There are no both-strong quotients on at most three vertices capable of such a crossing. The complete larger table at the prescribed target is:

| Other blocks | Both-strong partial quotient types | Types with a passing width box | Passing boxes modulo 12 | Width vectors within the size bound | Survive all modular tests |
|---:|---:|---:|---:|---:|---:|
| 2 | 2 | 0 | 0 | 0 | 0 |
| 3 | 19 | 0 | 0 | 0 | 0 |
| 4 | 217 | 17 | 26 | 24,882 | 0 |

The two- and three-block obstructions already hold for every exponent congruent to 1 modulo 12, independent of order. The four-block conclusion here is specific to exponent 61. Removing the requirement of homogeneous transitive blocks enlarges the search space.

After compiling the modular replay above, use `python3 -B transitive_ear_bank.py 6 | /tmp/chanpak_seed_width --bank` for the four-block row. Replace 6 by 4 or 5 for the other rows. [Complete partial-quotient generator](transitive_ear_bank.py)
