# Why the same module cannot fill the six-path deficit

A proposed repair of the b=2 column would produce `7(3·2^a−1)` for every a≥1 by substituting a transitive block of size a into one marked vertex of a fixed tournament. Testing the first two widths would require a single-root homogeneous duplication from **35 to 77**.

**Finite-exact check with an all-order scope.** No tournament with H=35 has a vertex whose homogeneous duplication yields H=77, regardless of the original order. Here duplication means two vertices joined by one arc with identical contacts to every external vertex; this does not constrain arbitrary cloning with different external contacts.

The SCC counts multiply, and no SCC has count 7. Since 35=5·7 and every nonsingleton SCC has count ≥3, there is exactly one nonsingleton SCC, with H=35. [Busch's lower bound](https://www.combinatorics.org/ojs/index.php/eljc/article/view/v13i1n3) puts its order at most seven. Duplicating a singleton SCC has no effect. For a marked vertex in the nonsingleton SCC, singleton padding likewise has no effect before or after the split. So it suffices to test all labeled order-seven tournaments with a single fixed marked vertex: pad smaller cores by ordinal singletons and relabel the mark to vertex 1.

The independent forward subset DP in clone_probe.cpp computes all order-seven counts, then literally constructs and counts the order-eight duplication of the marked vertex whenever the count is 35. It finds zero hits. In particular, the same low-rail sequence `3·2^a−1` cannot be lifted exactly by **any** fixed quotient plus one uniform transitive block for all a≥1, even if the new quotient has more return routes. This does not exclude changing the quotient with a, a family starting at a≥2, or a multi-block encoding.

Run:
```sh
g++ -std=c++17 -Wall -Wextra -Werror -O3 clone_probe.cpp -o /tmp/rail-clone
/tmp/rail-clone | diff -u clone_probe.out -
```
