# Two return arcs: square counts, a square obstruction, and Pythagorean data

Status: **PROVED** all-order counting and symmetry statements; **FINITE-EXACT**
independent small controls. This is the family repaired by the three-black
construction in [detours.md, section 3](detours.md). The square in its path
count and the failure of the arc-disjoint square mechanism have different
proofs; neither is inferred from matching numerical values.

## 1. The product mechanism

For a,b>=2, let T_(a,b) be the transitive tournament on 0,...,a+b, with
exactly the pairs (0,a) and (a,a+b) reversed. Write A={0,...,a-1}, c=a,
B={a+1,...,a+b}. Removing c leaves a transitive tournament.

Every directed Hamiltonian path has the unique form

```
A_I increasing, B_J increasing, c,
A_(A\I) increasing, B_(B\J) increasing.                 (1)
```

The admissible choices are exactly

```
I subset A\{0}, or I=A;
J empty, or J={a+b} union a subset of B\{a+b}.            (2)
```

Indeed, each of the two portions separated by c must increase. If any A
vertices follow c, the first must be 0, its only out-neighbour in A. If any
B vertices precede c, the last must be a+b, its only in-neighbour in B.
These are necessary conditions. They are also sufficient: when J is empty,
a nonempty I ends at a vertex other than 0; when A\I is empty, a nonempty
B\J begins at a vertex other than a+b. Here a,b>=2 handles the full/empty
boundaries. All other joins follow the transitive ordering.

Thus the **raw number of directed Hamiltonian paths** is

```
H(T_(a,b)) = (2^(a-1)+1)(2^(b-1)+1).                    (3)
```

The map is a bijection from paths to two independent subset choices; it
retains the entire path if I,J are kept. Multiplying just their cardinalities
loses endpoint and arc-incidence information needed for metagraph moves.

## 2. Balanced blocks give literal squares and a blue diagonal

Put T_k=T_(k,k), k>=2, and q=2^(k-1)+1. Reflection sigma(i)=2k-i is an
anti-automorphism. The only repeated scores are those of {0,1}, equal to
2k-1, and {2k-1,2k}, equal to 1. All intervening vertices have distinct
scores. An automorphism must preserve each repeated pair, and cannot swap
its ends because they have a directed edge. Consequently Aut(T_k) is trivial.
Any two anti-automorphisms differ by an automorphism, so sigma is unique.

A path in (1) is blue exactly when its reversal is sigma, equivalently

```
J = sigma(A\I).                                         (4)
```

This gives one admissible J for each admissible I. Therefore

```
H(T_k)=q^2,  b(T_k)=q,  Delta(T_k)=q(q-1).               (5)
```

Here b is the number of blue Hamiltonian paths and Delta is the black
path-orbit deficit h-b. Trivial Aut makes h=H and raw counts equal orbit
counts. Equation (4) identifies a genuine diagonal in a q-by-q choice
array. The off-diagonal entries are precisely the black paths.

The first values are:

| k | order | q=b | H | Delta |
|---:|---:|---:|---:|---:|
|2|5|3|9|6|
|3|7|5|25|20|
|4|9|9|81|72|
|5|11|17|289|272|
|6|13|33|1089|1056|

The inherited deficit-two series D_2(z)=z^4 is unchanged: these are
deficit-six and higher examples, not additional deficit-two centres.
The ordinary path-count series for this family alone is

```
sum_(k>=2) H(T_k) z^(2k+1)
 = 4z^5/(1-4z^2) + 4z^5/(1-2z^2) + z^5/(1-z^2).
```

In particular H(T_(k+1))/H(T_k) tends to 4 and b/H=1/q tends to zero.
This doubling mechanism is different from the Fibonacci transfer matrix in
[the earlier golden-ratio packet](../spectrum_basins_20261020e/arithmetic.md).
The appearance of 17 is explained here by a binary subset count, not by a
claimed identification of those two mechanisms.

## 3. An all-order obstruction to the square *method*

Let P=(0,...,2k). There is exactly one directed Hamiltonian path Q in T_k
arc-disjoint from P, and it is blue with the same positional reversal sigma.

To prove uniqueness, delete the arcs of P. Vertex 1 has indegree zero and
vertex 2k-1 has outdegree zero, forcing the endpoints. Vertex 2k can only
be followed by k, while vertex 0 can only be preceded by k. Thus every path
contains the consecutive block (2k,k,0). Outside k all arcs still increase,
so the remaining vertices form one increasing run before this block and
one after it. Adjacent labels in each interval A and B cannot lie in the
same run: their arc was deleted. Hence the two intervals must alternate
between the two runs. The locations of 0 and 2k fix the alternation uniquely:

```
Q = (odd labels in A, even labels in B, k,
     even labels in A, odd labels in B),                (6)
```

each sublist increasing. This is a legal path, and reflection reverses it.
The exact square colour law in [detours.md, section 4](detours.md) therefore
makes the P,Q square all blue. There is no black arc-disjoint replacement
at any order in this family. Nevertheless every marked P edge has a black
detour of length at most three by the two-block repair. This is an obstruction
to a proposed universal method, **not** a bridge or disconnected black core.

## 4. A precise, limited Pythagorean map

Equation (5) supplies the integer q=sqrt(H)=b. Since q is odd, define

```
(L_odd,L_even,C) = (q, (H-1)/2, (H+1)/2).                (7)
```

These are primitive Pythagorean triples: substituting H=q^2 proves the
identity, gcd(q,(q^2-1)/2)=1, and the even leg and hypotenuse differ by one.
This is Euclid's formula with consecutive parameters (q+1)/2,(q-1)/2.
The first triples are (3,4,5), (5,12,13), (9,40,41), (17,144,145).

The source is this marked tournament family, the target is primitive
triples, and the preserved data are the exact arithmetic H=q^2 and blue
count q. The map forgets the tournament's arcs and all detour certificates;
recovering those requires the construction of T_k. Its triangle ratios tend
to 0:1:1, not to the isosceles-right limit. Thus this explains an actual
difference-one family without supplying a Kuratowski, Tutte, or snark-minor
reduction. The marked torus-grid map in [torus.md](torus.md) preserves
different data, namely winding and Euler length.

## 5. Reproduction and scope

Run `python3 05-knowledge/results/strong_blue_games_20261020f/double_return_counts.py`
and the same command with `python3 -O`. The checker compares a direct DFS
path generator with the independent subset construction for every a,b in
2..5; it also checks balanced k=2..7, all blue predicates, the unique
arc-disjoint hostile, and literal automorphism counts through order nine.
Larger k=2..100 controls check only the arithmetic identities, not graph
enumeration. All infinite claims above have structural proofs.
