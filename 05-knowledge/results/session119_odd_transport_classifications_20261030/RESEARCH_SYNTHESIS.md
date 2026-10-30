# Synthesis: what transport unifies, and what a major advance still needs

**Session 119, 30 October 2026.** Anchor: exact bridges between odd cycles, odd
functions and odd arithmetic. Niche: complete character-covariant tournament
families. Wildcard: square-class dimension and the four odd eighths. The worktree
began at db3daffd7c4999cc931595731b65e9d983d784d5, preserving earlier shared-gate
results and all concurrent main-branch work. Classical facts and inherited
homogeneous classifications are credited in the individual proofs.

## 1. A precise common theorem, not a common adjective

The [transport theorem](synthesis/TRANSPORT_AND_EIGHTHS.md) says an affine cycle
closes exactly when (id-M)t=B. The signed, arithmetic and character cases specify
what M is and which values t are permitted. The map to a determinant keeps M but
can lose B. The map to an orbit quotient can lose the return sign. The map to a
rank can lose the nonzero-sign requirement. Each loss has an explicit hostile in
this packet; none is repaired by calling all the objects odd.

The completed dictionary retains these different data:

- Signed graph: edge signs, actual right-hand sides, spanning-tree paths and pins.
- Odd-function lattice: reflection pairs, additive-order strata, signed orbit
  closure and coefficient ring, including its fixed-point two-torsion.
- Arithmetic cycle: ordered exponents, ordered carry, denominator and positivity.
- Character family: conductor, stabilizer loops, support strata and the distinction
  between labelled symbols, multiplier orbits and arbitrary isomorphism classes.
- Eighth phases: the four characters of U(8), with the real/even and imaginary/odd
  parts retained. Unit differences of an even cyclic group give a bipartite graph,
  not the full tournament required by the odd-cycle Hamiltonian formula.

## 2. What is established and what was already known

[INHERITANCE](INHERITANCE.md) recovers the old homogeneous translation and incidence
counts, prime-power Legendre frames, the Moree/THM-4792 order criterion and the
Collatz carry formula. They are not new results of this session. The exact
inhomogeneous/integral interface, full operator classification, composite support
classification and family comparisons are the proved synthesis pursued here.

The fixed-character prime-power family has only one tournament isomorphism class
per (p,a), despite its 2^a labelled valuation frames. Separate characters on strata
escape the mixed-prime obstruction. At N=21 the resulting sixteen symbols form
exactly four isomorphism classes within that family, certified by multiplier
orbits and distinct fourth spectral moments. No arbitrary 21-vertex census is claimed.

The final square-class lift removes the character assumption altogether while
retaining invariance under square units. Its exact count exponent is
C=(product_p(1+2a_p)-1)/2. At three distinct prime factors, odd functions on the
square-class quotient need no longer be characters: the complete three-coordinate
list has eight signed characters and eight signed even-character multiples of
majority. At231 this doubles the labelled family from4096 to8192. No claim that
all added symbols have new abstract isomorphism types is made.

Equivalently, in Boolean coordinates on an r-dimensional square-class quotient,
choose a linear coordinate t that changes under the antipodal vector and r-1
coordinates z that do not. Every antipodal-odd sign is uniquely
(-1)^(t+g(z)) for an arbitrary Boolean function g. This is a complete uniform
formula, not a gate lower bound. The quotient dimension r, graph order N and
original numerical input-bit length are distinct resource parameters.

### Concurrent odd-descent work retained

While the final checkpoint was being published, commit
89d77408612b6db64a5b76bffc080e75bc28d6ae added the separate
[odd-descent and local-relaxation packet](../odd_descent_local_relaxations_20261030/README.md).
Its signed-affine spanning-tree reduction overlaps the mechanism here and is
credited as concurrent work, not an exclusive discovery of this session. Its
README and initial proof sections were read; the full median/affine-relaxation
certificates were not independently replayed here. The incoming commit and its
index changes are preserved in our final parent chain.

That packet also states precise constraint-preserving median and XOR descent
results. The distinction matters for our majority lift: an odd majority function
can define a new square-class symbol, but combining arbitrary valid tournament
objects by majority would require a separate constraint-preservation theorem.
No median-descent claim is imported merely from functional oddness.

## 3. A decisive Ramsey boundary from the family theorem

Let t(T) be the maximum transitive subtournament order. For any tournaments X,Y,

    t(X[Y])=t(X)t(Y)

for the uniform lexicographic product. Indeed a transitive subset meets copies of
Y whose indices form a transitive subset of X; within each copy it uses at most
t(Y) vertices. Conversely choose maximizing chains in a maximizing chain of copies.
Therefore every prime-power valuation frame satisfies

    t(T)=t(P_p)^a.

For p=3 this is 2^a on 3^a vertices, the inherited ternary-tower case. Thus merely
changing frame signs cannot manufacture a logarithmic-transitivity Ramsey family:
it changes only a labelling. This is a concrete obstruction, not a numerical failure
of a search. A constructive Ramsey advance must leave that one-isomorphism-class
family or use a genuinely different objective.

## 4. Bounded research targets, not claimed breakthroughs

**Ramsey target.** Use the mixed-prime, stratumwise and square-class lifts to search
for smaller transitive subtournaments while retaining literal arithmetic labels.
A useful uniform theorem must account for prime exponents as well as the number
of distinct primes; the ternary-power boundary refutes a blanket logarithmic bound.
Counting symbols does not count nonisomorphic tournaments or establish Ramsey bounds.

**Arithmetic target.** Use the explicit cyclic cokernel coordinate as a certified
filter for candidate periodic words. Any positive-cycle exclusion still has to
control the ordered carry over an unbounded family. Equal determinant, exponent
sum and exponent multiset do not suffice, as the negative-cycle hostile proves.
Periodic-word classification alone also does not prove every orbit is periodic.

**Algorithmic target.** Use spanning-tree transport to solve signed inhomogeneous
systems in linear arithmetic-operation count, and use the signed multiplier
cycles to solve the corresponding odd-function equations. These are exact linear
subproblems. They do not solve unrestricted coloring or Boolean circuit complexity.

**Circuit boundary.** Odd-weight Boolean functions can have no constant nonzero
translation derivative and still have n-1 full-B2 gates. The Session 117 examples
supply that counterbenchmark. Antipodal-odd square-class functions are different:
their prescribed derivative is constant and they have even truth-table weight.
Neither parity property alone supplies a gate lower bound.

No major open conjecture, new unrestricted Ramsey-number bound, Collatz convergence
claim, LRC closure or explicit-NP gate lower-bound record follows from this packet.
The results identify exact solvable interfaces and concrete next obstructions.
