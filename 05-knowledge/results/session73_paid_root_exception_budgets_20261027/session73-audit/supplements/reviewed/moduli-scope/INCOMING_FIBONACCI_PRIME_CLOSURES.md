# Incoming pattern-clause closures: Fibonacci and primes

Session73 additive supplement to the bounded incoming HYP-9419(4) correction.
**PROVED elementary:** Fibonacci obstruction for every integer modulus m>=2.
**CITED corollary of Dirichlet:** prime-pattern obstruction for every m>=2.
These are not conjectures requiring an unbounded numerical census. The broad
polynomial-degree exclusion is separately **REFUTED**, as established in the
frozen core by x^4 modulo7 (and modulo3). No other HYP-9419 clause is settled.

The target is pinned to commit `64a7985e7054821c5958dc85f6dd489c844111f2`:
`05-knowledge/hypotheses/HYP-9419-unifying-candidates-three-layer-composition-rank-laws-dissection-obstructions-pattern-clocks.md`,
clause4 and its title/status. Its nonzero pattern S_m defines arcs
x->x+s for s∈S_m. The inherited skew criterion requires S_m∩(-S_m)=empty
for even an oriented graph. Having both1 and−1 violates this already.

## F. Fibonacci has both signs of1 at positive indices

Let F_0=0,F_1=1 and F_(n+2)=F_(n+1)+F_n. For every m>=2, the pair update

    A(a,b)=(b,a+b) on (Z/mZ)^2

is a permutation, with inverse A^(-1)(a,b)=(b-a,a). The orbit of(0,1)
is therefore periodic from time0, not merely eventually periodic. Write
its period as P. Since(0,1),(1,1),(1,2) are distinct when m>=2, P>=3.
Backward twice gives

    A^(-1)(0,1)=(1,0),   A^(-2)(0,1)=(-1,1).

Periodicity now implies `F_(P-2)=-1 modm`, with P-2>=1. Together with
F_1=1, the nonzero Fibonacci pattern contains ±1 at positive indices.
Thus its connection graph is not oriented and hence not a tournament for
**every m>=2**. At m=2, the two residues coincide but are self-opposite,
which still violates disjointness; this case is not an exception. QED.

This concerns the full residue pattern, not successor determinism or a
finite prefix. The earlier clock-even Fibonacci scalar-successor hostile is
a different statement. No general second-order sequence claim is made.

## P. The prime pattern is already obstructed by Dirichlet

Dirichlet's theorem ensures infinitely many primes in every residue class
coprime to a fixed modulus. Apply it to1 and−1 modulo m>=2. Both are units,
so the nonzero prime pattern contains both, proving that its connection graph
is not oriented and not a tournament. At m=2, any odd prime supplies the
self-opposite residue1. This is a **CITED corollary**, not a new proof of
Dirichlet's theorem or a numerical inference. Primary source:
[Dirichlet's 1837 paper, R. Stephan translation, v2](https://arxiv.org/abs/0808.1408v2).
The formulation was also checked in
[Andrew V. Sutherland, MIT18.785 Lecture18, Theorem18.1](https://ocw.mit.edu/courses/18-785-number-theory-i-fall-2021/mit18_785f21_lec18.pdf).

More precisely, Dirichlet gives the full prime pattern

    S_m = (Z/mZ)^* union
          {p modm : p prime, p divides m, p modm != 0}.     (P1)

Any prime not dividing m is a unit and every unit class occurs; the only
other prime values are the finitely many prime divisors of m. This retains
the nonunit residues rather than silently replacing the full pattern by its
unit subset. The ±1 obstruction needs only the first term of(P1).

## Scope of the repaired incoming clause

- Polynomial degree>=3: **REFUTED** as a universal exclusion.
- Fibonacci: **PROVED for all nontrivial moduli m>=2**, by the pair permutation.
- Primes: **CITED corollary for all m>=2**, using Dirichlet, not PNT or GRH.
- m=1: the nonzero connection set is empty and the graph has one vertex.
  Under the standard vacuous convention this is a tournament, so “at any
  modulus” must exclude m=1 in the surviving nonexistence clauses.
- Nothing here treats a finite observed prefix as the whole infinite pattern
  or proves a uniform numerical cutoff in the modulus.

`check_patterns.py` enumerates complete Fibonacci pair orbits for2<=m<=512,
checks both inverse steps and the positive index P-2, and independently tests
pair-update invertibility on every pair for2<=m<=32. For primes it finds and
trial-divides explicit ±1 witnesses for the same511 moduli, within a declared
sieve bound10^6. These are finite controls, not the proof of the all-m claim.
Normal and -O receipts are byte-identical. External source metadata is recorded
in `EXTERNAL_SOURCES.json`; repository snapshot pins retain their exact bytes.
