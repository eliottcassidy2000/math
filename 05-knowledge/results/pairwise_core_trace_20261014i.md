# Guarded three-matrix compression of signed graph sums

2026-10-14. **PROVED:** the exact trace identity and its arithmetic-operation
bound. **FINITE-EXACT:** the explicitly listed graph/prime controls.
**CITED:** the attached preprint's counting mechanisms. These are finite
evaluation tools, not an all-depth Collatz or Mahler certificate.

## Inheritance and the source of the transfer

The closest proved mechanism is the injective low-degree reduction in
[the wheel/run interface](wheel_run_interface_20261014h.md). Its hostile
examples are collision terms and a signed sum of zero with nonempty support.
The corrected near miss is replacing extension-field Frobenius traces by
Fourier coefficients at prime-power indices. The least-used sidecar is the
collection of inequalities between vertices in different blocks, including
pairs that are not edges of the graph.

The supplied paper is Sanjay Jain, Junqiang Peng, Frank Stephan and Mingyu
Xiao, *Solving Restricted Versions of #SAT Faster than 2^m: Algorithms and
Lower Bounds*, arXiv:2610.15426v1 (11 October 2026), supplied as
`C:/Users/Eliott/Downloads/2610.15426v1.pdf`. We use its techniques, without
claiming independent verification of every complexity result. Pages11 and16
were rendered and inspected after text extraction.

Three mechanisms matter here:

1. Its clause branching retains a small pivot while branching on a neighbour.
   A measure `m+eta/k` pays both ordinary progress and pivot replacement.
   This suggests a precise question for adaptive arithmetic controllers:
   **which retained witness supplies a bounded debt, and which authenticated
   transition pays its replacement?** Merely changing patterns does not pay.
2. Lemmas5--6 and Theorem5 express signed sums of unary and pairwise factors
   as the trace of three matrices. The occurrence-two restriction is what
   makes this factorization applicable; unrestricted ternary factors do not
   acquire pairwise structure just by grouping variables.
3. Sections6.2.1--2 preserve cancellation-sensitive data by a mixed-radix lift.
   If nonnegative bounded scores are p,q and the sign is (-1)^h, retain
   `E=p+Bq+B^2 h`, with B larger than either score. Coefficients sharing E
   have the same sign. This retains the signed histogram, not the identities
   of all assignments inside a histogram cell.

Our concept board is: carry / quotient / collision / sign / retained witness /
terminal core. Each is a coordinate that can be removed only after proving
that the desired predicate survives its removal.

## Exact graph identity

Let G be a directed loop-free labelled multigraph on n>=2 vertices, and let
q be an odd prime. Define

\[
 S_G(q)=\sum_{x:V(G)\hookrightarrow\mathbb F_q}
              \prod_{(u,v)\in E(G)}\chi(x_v-x_u).
\]

Repeated edges are permitted. If n>q there are no injective maps. If the
number of edges is odd, multiplication by a nonsquare reverses every
summand, so S_G(q)=0. Suppose henceforth that it is even.

Every injective assignment has a unique affine normalization x_0=0,x_1=1.
The edge product is unchanged because its scaling factor is
`chi(x_1-x_0)^|E|=1`. Thus the full sum is q(q-1) times the normalized sum.
The two fixed vertices need not be joined by an edge.

Split the remaining vertices into three blocks I,J,K. A state in each block
is an injective assignment from that block into F_q minus {0,1}. Put all
within-block and fixed-to-block edge factors in a unary U_i. Put cross-block
edge factors, **and every cross-block distinctness indicator**, in V_ij.
The fixed-to-fixed edges contribute a constant c. Set

\[
 A_{ab}=U_I(a)V_{IJ}(a,b),\quad
 B_{bc}=U_J(b)V_{JK}(b,c),\quad
 C_{ca}=U_K(c)V_{KI}(c,a).
\]

Then, exactly,

\[
 \boxed{S_G(q)=q(q-1)c\operatorname{Tr}(ABC).}
\]

Expanding the trace proves the statement: each unary appears once, each
cross factor once, and the three internal plus three cross inequality sets
are exactly global injectivity. No approximation or limiting argument occurs.
The same identity holds over an odd finite field if field arithmetic is
supplied; the present implementation deliberately accepts only prime fields.

For block sizes a>=b>=c and Q=q-2, padding to powers gives matrices of sizes
at most Q^a by Q^b, Q^b by Q^c, Q^c by Q^a. Blocking into square Q^c matrices
gives an arithmetic-operation upper bound

\[
 O^*(Q^{a+b+c-3c+\omega c}).
\]

Here omega may be any admissible square-matrix multiplication exponent;
factor construction and final contraction fit this bound since omega>=2.
This is an operation count with integer bit costs additional. For fixed G,
intermediate integers have O(n log q) bits. The supplied prototype uses
classical Python-integer matrix multiplication, not a fast-omega algorithm.

For the new six-vertex core, splitting4 remaining vertices as2+1+1 gives
O^*(q^(omega+1)), versus direct O(q^4) normalized enumeration. With the
classical backend both have exponent4; no measured asymptotic speedup is
claimed. For W4, the inherited specialized second-moment formula is already
O(q^2), better than this generic trace construction. Always retain that
specialized evaluator when it applies.

## The new core and independent controls

The [degree-ten reduction](paley_next_cores_20261014i.md) leaves, besides the
wheel, a six-vertex core V: K6 with the path edges01,12,23,34,45 removed,
every retained edge oriented from smaller to larger label. Its full sums are:

| q | V(q) | Affine-normalized sum |
|---:|---:|---:|
|7|336|8|
|11|-12320|-112|
|19|-76608|-224|
|31|52080|56|

All were computed twice, by the trace and an independent direct permutation
sum. At31, block state counts are812,29,29, requiring47,937 stored matrix
entries. The experiment also checks W4 and K4 at3,5,7,11,19,31, including
empty injective universes. Reversing one edge multiplies the sum by chi(-1).
The total unsigned number of maps is(q)_n, so the positive and negative
counts are exactly((q)_n+S)/2 and((q)_n-S)/2.

Two hostile controls explain why the retained coordinates matter:

- Atq=7, dropping cross-block inequalities changes the wheel sum from-1176
  to-1344. Absence of an edge does not license a vertex collision.
- K4 atq=7 has signed sum0, but420 positive and420 negative embeddings.
  Cancellation cannot certify absence. Recording the sign separately is
  already a two-cell lossless histogram for this question.

The checker rejects composite q=9, Boolean integer aliases, loops, malformed
guard flags and exceeded matrix-entry caps. It uses exact object arrays;
there are no floating-point arithmetic or modular-reconstruction assumptions.

```text
python -B 04-computation/experiments/pairwise_core_trace_20261014i.py
python -B -O 04-computation/experiments/pairwise_core_trace_20261014i.py
```

The saved [JSON](pairwise_core_trace_20261014i.json) contains288 checks.
Of these,222 independently enumerate every full injection for every simple
labelled graph on2..4 vertices overq=3,5,7, without affine normalization or
matrix factorization. Two further controls check parallel and reversed
parallel edges. The exhaustive audit was proposed and first run by a
separate reviewer, then incorporated into the reproducible checker.
Normal and optimized execution produce identical saved JSON; LF-normalized
SHA256 is `d91ffd69a0ec73f97598373975ffc887081e271b7ef985d521c9967500508f19`.
An independent reviewer audited the normalization, three factor assignments,
bit/injectivity guards, resource cap and arithmetic-operation bound.

## What this changes in the larger recursion picture

The source/target map here is from injective assignments to three block
states. It preserves the complete signed sum because its sidecars restore
all removed interactions. The [triangular Mahler quotient](recursion_carry_bridges_20261014i.md)
has the analogous one-bit sheet, but a different predicate: strict real
safety. Its two complementary lifts have opposite safety despite the same
entire quotient orbit. This is a concrete common mechanism, not an
identification of graph counting with integer dynamics.

The immediate next targets are to evaluate the six-vertex core V by a
specialized moment identity, and to find a source-dependent Mahler witness
whose replacement has a paid finite bound. The matrix trace solves a finite
aggregation problem. It does not establish uniform control of the infinite
Hamiltonian expansion, ordinary positive stabilization of Mahler codes, or
universal Collatz descent. Those remain separate obligations.
