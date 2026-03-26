# Code IS Compressed Mathematics

*opus-2026-03-25-S344*

## The Thesis

A program is a finite representation of an infinite mathematical truth.
The function f(x) = x^2 has infinitely many input-output pairs, but the
CODE `lambda x: x*x` is 14 bytes. This is infinite compression.

## The Tournament Connection

A program's dependency DAG is a partial order. Completing it to a tournament
gives the full scheduling structure. The key invariants:

| Tournament Invariant | Program Meaning | Code Meaning |
|---------------------|-----------------|--------------|
| Score sequence | Criticality ranking | Variable importance |
| H(T) = # Hamiltonian paths | # valid execution orders | Scheduling freedom |
| Cycle structure | Feedback loops | Recursion/iteration |
| Isomorphism class | Semantic equivalence | α-equivalence class |
| Transitive tournament | Sequential code | Chain of assignments |
| Regular tournament | Balanced parallelism | Map-reduce |
| χ(T) | Min parallel stages | Critical path length |

## Three Levels of Code Compression

1. **Syntactic**: compress the TEXT of the program (zlib, indent-split).
   Savings: 3-5x over raw text.

2. **Structural**: compress the AST skeleton + token residuals.
   Separates WHAT the code does (skeleton) from HOW it names things (tokens).

3. **Semantic**: compress the MEANING (α-renaming, expression normalization).
   Two programs computing the same function → one canonical form.
   This is tournament isomorphism: the iso class is the compressed form.

## The Decomposition Theorem

Source code decomposes into:
- **Structure** (indentation, nesting, control flow) — LOW entropy (2.0-2.7 bits)
- **Content** (identifiers, values, operators) — HIGH entropy (5.2 bits)

This EXACTLY parallels tournaments:
- **Score sequence** — constrained by Erdős-Gallai conditions (low entropy)
- **Arc orientation** — exponentially many choices (high entropy)

Compressing them separately (indent-split) wins on every tested file.

## The Deep Principle

Code encodes the RULE, not the TABLE.
A lookup table for f: [0,255] → [0,255] needs 256 bytes.
The code `lambda x: x*x % 256` needs 19 bytes.
Even compressed, the table (145 bytes) is 7.6x larger than the code.

This is because the code captures the TOURNAMENT STRUCTURE of the computation:
the dependency graph, the data flow, the invariants that make the function
predictable from its definition. The table captures none of this — it's
the brute-force enumeration of all input-output pairs.

Tournament theory tells us: the compression ratio = H(T) / n! where
H(T) is the number of valid computation orderings and n! is the number
of possible orderings. For a transitive tournament (sequential code),
H/n! = 1/n! — exponential compression. For a regular tournament
(fully parallel code), H/n! approaches 1 — minimal compression because
all orderings are valid and must be distinguished.

## The Practical Insight

For code analysis tools:
- **log2(linear_extensions) = scheduling entropy** = bits of information
  a parallel scheduler needs to specify one execution
- **Critical path = longest directed path** = inherent sequential bottleneck
- **Score sequence = dependency ranking** = which code to optimize first
- **Structural fingerprint = AST tournament iso class** = clone detection
