# Joint incidence, transport and guards: an operational composition theory

2026-10-26. **PROVED elementary interfaces:** guarded Fourier formula, response
transport, precision-changing clock maps, and the two-stack simulation below.
**FINITE-EXACT:** the displayed censuses and 81/27 second response.
**OPEN:** the scalar cyclicity conjecture and economical full attachment carriers.
The general organizing framework is a proposal, not a universal classification.

## Main advance

The tournament/even-graph correspondence becomes useful for composition as
a **change of observable basis**, with operators and realization data retained.
For a named Hamilton path P, its legality guard multiplies the complement
operator. In signed graph coordinates this is convolution, not a diagonal
edge toggle. This supplies a concrete extension of the source-independent
operation envelope in the concurrently received
[Session72 packet](../session72_even_graph_unifying_carriers_20261026/README.md).
That packet arrived at `00f1c9aa39` while this work was running; its Fourier
envelope and positive-cone results are inherited, not claimed again here.

The second response already separates the 81/27 hostile. Let

    (L f)(T) = sum_{P a directed Hamilton path of T} f(T^op xor E(P)).

Then L1=H and L^2 1 counts two successive legal path-complement choices:

| source | Aut | H | L^2 1 |
|---|---:|---:|---:|
| C3[A+,A+,A+] | 81 | 462915 | 202898835819 |
| C3[A+,A+,A-] | 27 | 462915 | 202869738135 |

Here A+=C3->1, A-=its converse. The difference is 29,097,684. Under one
uniform re-marking/complement step, their expected partner H values differ
by 1077692/17145. This is a new finite response test on the inherited control,
not a new construction of the pair or a replacement for its full response law.
The earlier [total-variation calculation](../unnamed_carrier_dissection_20261026/CONNECTIONS.md)
remains stronger for comparing entire target laws.

## Three exact bridges

1. **Tournaments -> signed graph observables.** Relative to the group that
   fixes the required ports, admissible signed graph-orbit sums form a
   complete basis. A gluing operation becomes a tensor in these coordinates;
   sequential operations multiply their matrices. Positivity is carried by
   the transported cone, not by requiring the new coefficients to be positive.
   [Proof, formula and four-vertex example](FOURIER.md).
2. **Integer sequences -> clock arrows.** Reduction modulo m is a legitimate
   deterministic state only if it determines the next arrow and its guard.
   Affine integer maps descend; Fibonacci values need a second coordinate;
   division consumes precision. Exact arrows from modulus qm to modulus m
   replace false autonomous arrows on one clock. [Clock and machine proofs](CLOCKS-MACHINES.md).
3. **Turing configurations -> two stacks -> marked tournament tapes.** Binary
   stack push/pop are multiplication and division by two. The inherited
   TT3/C3 tape encoding then gives controlled tournament rewrites, retaining
   head, program, origin and read guards. This is an exact simulation with
   those data, not universality of native metagraph complement.

## A precise new classification target

The elementary converse symmetry gives

    span_Q{1,L1,L^2 1,...} subset functions invariant under tournament converse.

Exact rational arithmetic through six, and a prime-modulus rank certificate
at seven together with the converse upper bound, give:

| n | tournament types | converse orbits | response-space rank |
|---:|---:|---:|---:|
|2|1|1|1|
|3|2|2|2|
|4|4|3|3|
|5|12|10|10|
|6|56|34|34|
|7|456|272|272|

**OPEN CONJECTURE C:** equality holds for every n. Equivalently, the constant
vector is cyclic for the counted path-complement operator on converse-invariant
functions. It would classify types up to converse by all successive total
response counts. It is stronger than pairwise separation by those sequences.
The 81/27 result passes its first nontrivial hostile test at n12, not an n12
class census. Next decisive test: exact or certified modular rank at n8.

This conjecture concerns raw path counts. The presentation histogram divides
each source row by its automorphism order; the probability kernel divides by H.
They are different operators. In particular a stochastic kernel always fixes1,
so replacing L by the normalized kernel makes this conjecture false from n3.
No all-size strong-atom response theorem is inferred from the table.

## The larger proposal

Classify objects by the responses of their **permitted future contexts**.
Joint incidence specifies what can be matched, transport specifies how frames
change, and guards specify whether the proposed match is realized. The exact
quotient is relative to this operation language. [The framework, history and
falsifiable obligations](THEORY.md) distinguish:

- state separation from linear dimension of distributions;
- symmetry compression from erasure by a genuine reset;
- finite observability at fixed size from a uniform finite carrier;
- existence of a finite quotient from an algorithm that can find it.

The last distinction is especially sharp for machines: every single deterministic
run has a finite quotient for the sole observer "has halted", but no uniform
algorithm can construct such a quotient for every program/input. Full tape-read
observations instead have growing, and across unbounded tapes infinite, index.
These statements are proved in CLOCKS-MACHINES; they do not decide a cryptid.

## Validation and reproduction

`check.py` independently enumerates all labelled tournaments and graphs through
six, uses every vertex permutation, tests signed-orbit orthogonality, computes
exact rational operator matrices and ranks, checks all four-vertex named-path
character truth tables, and checks ordered joins. It also checks arithmetic
clock compositions and native tape edits against the independent stack formulas.
There is no inherited atlas filter or external package dependency.

`rank.cpp` independently generates the complete tournament atlases through
seven and certifies rank over the prime 1000000007 (checked by trial division).
At n7 it covers all 2,097,152 labelled tournaments. Full rank272 over that
field gives a lower bound272 over Q; converse invariance gives the matching
upper bound. The earlier ranks agree with the Python rational calculation.

`hostile.cpp` independently enumerates all Hamilton paths and their normalized
presentations for both 12-vertex sources, then uses subset dynamic programming
to compute every partner H. `audit_hostile.py` uses recursive memoization in
Python for both sources and twelve recorded partner controls. The full second
moment has one implementation; the sampled checks are not advertised as a
second exhaustive sum. All Python checks agree normally and under -O. New
proofs have a local mathematical audit; no separate agent review is claimed.

Run from the repository root:

```text
python -B 05-knowledge/results/clock_even_guarded_composition_20261026/verify.py
python -B 05-knowledge/results/clock_even_guarded_composition_20261026/verify.py --replay
```

Default verifies seals, both receipt pairs, and cross-checks. `--replay` compiles
both C++ programs with g++, repeats them in a temporary directory, and runs
both Python programs normally and optimized. `SOURCES.json` pins the recovered
mechanisms and literature scope. `SHA256SUMS` uses LF-normalized text bytes.
No theorem ID is reserved and no pre-audit incoming claim is silently promoted.
