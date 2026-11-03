# A short exact certificate for the two tournament gaps

## Statement and dependencies

For every finite nonempty tournament T, the number H(T) of directed Hamiltonian paths is neither 7 nor 21. This proof uses Busch's all-order lower bound for a strong tournament on n≥3 vertices, H(T)≥5^((n−1)/3), and a tiny exact finite computation below. Source: Arthur H. Busch, *A Note on the Number of Hamiltonian Paths in Strong Tournaments*, Electronic Journal of Combinatorics 13 (2006), N3, https://doi.org/10.37236/1141. Busch's Corollary 1 in fact gives the exact minimum

```
              3·5^((n−3)/3)  if n≡0 (mod 3),
min H(T) =       5^((n−1)/3)  if n≡1 (mod 3),
              9·5^((n−5)/3)  if n≡2 (mod 3),
```

for n≥3, with first values 3, 5, 9, 15, 25, 45, 75, … for n=3,4,5,6,7,8,9, …. The weaker cited bound, not its incorrect reading as equality for every n, is enough.

The tournament condensation is a transitive tournament, hence its strong components have a unique total order and every Hamiltonian path traverses them as consecutive blocks in that order. Thus H(T) is the product of the component counts. Every component is a singleton with H=1 or a strong tournament on at least three vertices with H≥3. If H(T)=7, a strong component has H=7; Busch's bound forces its order at most 4. If H(T)=21, because 7 is impossible, component multiplicativity forces one strong component of count 21, whose order is at most 6 by Busch's bound. The finite verification below excludes both.

## Exact finite check

Any tournament has a Hamiltonian path: successively insert each new vertex at the first position where it beats the existing path vertex, or append it if none. Relabel its vertices so that 0→1→⋯→n−1. Thus it is enough to inspect orientations of the remaining (n−1)(n−2)/2 arcs. The number of codes for n=1,…,6 is 1, 1, 2, 8, 64, 1024. This is an exhaustive labeled cover, not a claim about the number of isomorphism classes.

`python3 certificate.py` computes H by directly inspecting all n! vertex permutations, with integer arithmetic and no pruning. Independently, it evaluates the odd-cycle collection formula using all directed 3- and 5-cycles and all disjoint pairs of triangles. These are all possible collections of vertex-disjoint odd cycles at n≤6. The program fails on any disagreement. The second method is a cross-check, not a necessary assumption for the finite exclusion: the direct permutation computation itself certifies it. Every code's strongness is determined separately by transitive closure. The exact sorted strong spectra are:

```
n=3: 3
n=4: 5
n=5: 9, 11, 13, 15
n=6: 15, 17, 19, 23, 25, 27, 29, 31, 33, 37, 41, 43, 45
```

There are no strong tournaments of order 2. Singleton strong components have count 1. The full code histogram is recorded in `certificate.out`. `python3 -O certificate.py` has the same output: checks use explicit exceptions rather than assertions. An independent canonical-label calculation `classes.py` yields all 35 strong isomorphism classes of order 6, with canonical code and count in `classes.out`. The canonical code uses the ordered pairs (0,1),(0,2),…,(4,5), bit k=1 if the first vertex dominates the second, and the minimum over all 720 relabelings. This supplemental calculation establishes a compact second route to reproducing the spectrum, although verifying completeness of the representatives still requires checking the finite enumeration.

These results establish the two exclusions for tournaments of all orders; they do not prove that every other positive odd integer occurs.
