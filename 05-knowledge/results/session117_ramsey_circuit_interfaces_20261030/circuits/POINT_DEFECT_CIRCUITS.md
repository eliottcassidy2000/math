# Extreme-Walsh symmetry is a small-circuit family, not a hard-circuit amplifier

**Session 117 — PROVED with a finite-exact base certificate.** No novelty claim
relative to the circuit-complexity literature. This gives a precise obstruction
for importing the repository's session115 odd-eighth phase classification into
CircuitBreaker. It does not improve the unrestricted asymptotic record.

## Model and inherited object

`C(f)` counts gates in a deterministic acyclic, one-output circuit with all 16
binary Boolean truth tables, unrestricted depth, fan-out and sharing. Each
retained gate costs one; input wires and output selection are free. A retained
NOT/projection gate is not free. The constructions below absorb literal signs
into adjacent binary truth tables, rather than introducing a free-NOT basis.

For `n>=3`, `v in F2^n`, `r in F2^n`, and `c in F2`, let

    f(x) = c + r.x + 1[x=v]                 (all additions in F2),
    k = wt(r).

The inherited session115 Walsh classification says these are exactly the
functions with a Walsh coefficient of magnitude `2^n-2`, with unique affine
phase and defect point for `n>=3`. At n=3 they exhaust odd-weight truth tables;
at higher n they do not. The present statements concern this classified family,
not all odd-weight functions, all cubic phases, or all NP languages.

Input translations and permutations, and output complementation, preserve the
minimum gate count for functions depending on at least two inputs: relabel the
input pins' truth tables and complement the final gate's truth table. This
argument remains valid with shared fan-out. Thus it suffices to take `v=1^n`,
`r=(1^k,0^(n-k))`, `c=0`.

## Theorem 1: complete all-dimension gate classification

For every member above, all n>=3:

| affine support k | gate complexity |
|---|---|
| 0 | exactly n-1 |
| 1 | exactly n-1 |
| 2 | exactly n |
| 3 | exactly n+1 |
| k>=4 | exactly n+k-2 |

In particular the complete extreme-Walsh family has a uniform `2n-2` upper bound.
For an arbitrary sequence of such functions, `limsup C(f_n)/n <= 2`. No member
sequence can establish a B2 linear lower-bound coefficient greater than 2.
This rules out this particular synthesis route even before comparing it with
the 3.1 or 3.5 literature coefficients. It does not constrain different families
constructed using these phases as components, whose additional work must be
counted separately.

### Explicitness and input-length accounting

Here n is the number of Boolean input bits. The earlier Walsh cube has2^n labels; that truth-table length is not the circuit input length. A concrete fixed total P language is obtained by taking v=1^n,r=1^n,c=0 for every n>=3 and defining the two shorter positive lengths and the empty input separately (for example, output0). It is evaluated in linear time and has exact gate size2n-2 at every n>=3. An empty witness makes it an NP family as well. More generally any polynomial-time choice of the parameters gives a P family; arbitrary nonuniform parameter choices are not automatically an NP language.

### Upper bound proof

For k=0 the function is the conjunction of n literals and uses n-1 gates. For
k>=1, set `P=x2...xn`. Then

    f = (x2 + ... + xk) + (x1 AND NOT P).

Compute P in n-2 gates, combine it with x1 using the single B2 gate
`a AND NOT b`, and, if k>=2, compute the parity of k-1 inputs in k-2 gates
and take one final XOR. This costs n-1 when k=1 and n+k-2 otherwise.
Output constants and input translations are absorbed as above. The formula is
verified on all inputs, all k, for n=3,...,10 by the accompanying checker.

### Lower bound proof: dependency and gate elimination

Every variable is essential. For an inactive affine coordinate, the derivative
is the indicator that all other coordinates equal their defect values. For an
active coordinate it is one plus that indicator. Both are nonzero for n>=2.
A fan-in-two one-output DAG depending on n distinct inputs has at least n-1
gates: its connected ancestral graph has at most 2s incoming gate edges and
at least n+s-1 edges when it contains s gates. This proves the k=0,1 bounds.

We use a one-gate restriction lemma: if f depends on an input y and the
restriction f|y=b depends on at least two other inputs, then

    C(f) >= C(f|y=b)+1.

Take an optimal circuit and fix y. At least one gate receiving y becomes
constant or unary. Remove it and propagate its constant, projection or
negated projection through the descendants' binary truth tables, continuing
as needed. At least one original gate disappears. A negated projection at the
terminal output can be absorbed into the remaining output gate; the assumed
at-least-two essential variables ensure the remaining output is neither a
constant nor a lone input literal. Thus no terminal NOT gate needs to be
added. This is ordinary B2 normalization, not a free-NOT assumption.

### Pure-parity defect: a forced branching input gives the sharp induction

Write `F_m(x)=x1+...+xm+x1...xm`. For m>=3, some input of every
normalized circuit for F_m must feed at least two distinct gates. To prove this,
suppose every input feeds only one gate. Choose the earliest gate; after
removing unary/constant gates it combines two distinct inputs xi,xj. Both inputs
can affect the remaining circuit only through this one output bit, even if that
bit has arbitrarily many descendants.

Fix all other inputs to one. Up to an additive constant the resulting pair
function is `p+q+pq`. Fix the other inputs to an assignment not all one. Up to
an additive constant the pair function is `p+q`. (Such an assignment exists
because m>=3.) Both functions would have to factor through the same one-bit
gate `g(p,q)`. This is impossible: the joint value pairs of `(p+q,p+q+pq)`
on `00,01,11` are `(0,0),(1,1),(0,1)`, three distinct values, whereas a
single bit has only two values. Fixed additive constants cannot alter this
count. Therefore a branching input exists.

Set that input to one. Each of its at least two distinct receiving gates becomes
constant or unary and can be removed by the normalization argument above.
Repeated pins cannot fake the two receivers: a gate with two copies of the same
input was already unary and removed. The restricted function is
`1+F_(m-1)`. For m>=4 it still depends on at least three inputs, so terminal
negation creates no extra gate. Consequently

    C(F_m) >= C(F_(m-1)) + 2       (m>=4).

The exhaustive base certificate gives `C(F_3)=4`; hence
`C(F_m)>=2m-2` for every m>=3. The upper bound matches.

For general k>=3, successively fix all n-k inactive variables at their defect
values. Each deletes at least one gate, giving
`C(f)>=n-k+C(F_k)=n+k-2`. For k=2, fix n-3 inactive variables and use the
three-variable k=2 base cost3, obtaining `C(f)>=n`. The k=0,1 dependency
bounds above finish the complete classification.

The branching-input argument and all-k induction were supplied by the session117
root research lane and audited here. No claim of literature priority is made.

## Finite-exact base certificate, independently duplicated

`check_point_defect.py` exhausts **all** B2 circuits on three inputs with at most
three gates, including shared subgraphs. Two independent representations agree.

1. Complement-normalized available-function sets. Initial signals are the
   three inputs. Each useful gate on distinct available functions has one of
   five forms, up to complement: AND, `a AND NOT b`, `NOT a AND b`, OR, XOR.
   Constant/unary/redundant gates may be discarded when seeking a minimum
   for an at-least-two-essential-variable function. Input-sign choices are
   accounted for by these five forms. Sets are legitimate DAG states because
   fan-out is free and future gates may use any prior wire; they are not
   formula states. The successive state counts are 15, 180 and 2168.
2. Raw ordered signal sequences. This path keeps both signs, every degenerate
   gate, both ordered input pins, repeated pins, and all 16 truth tables.
   It makes no normalized-state or nonredundancy assumption. There are 38
   one-gate sequences and 1744 two-gate sequences; all possible third gates
   are then evaluated.

Both yield 116 truth tables up to output complement, or 232 actual functions.
For the four k strata the costs are 2,2,3,4; the fourth is absent through three
and has the explicit four-gate construction above. Every translated defect
point and output complement is separately checked. Input-coordinate
permutations transport other r of the same weight.

This is an exhaustively reproducible finite proof component, not a SAT/DRAT or
proof-assistant certificate. The checker has no optimization-removable
assertion statements. Normal and `python -O` outputs agree byte-for-byte.

## Additional four-variable audit

`check_four_inputs.cpp` independently tests the induction at its first new
dimension. It enumerates normalized available-function sets through four gates,
then tests every possible fifth gate for the target `F_4`. State counts are
1,30,735,18037,473767. The target is absent through five gates and has the
six-gate construction already verified by Python. This is corroboration, not
a replacement for the all-dimension induction. Compile with a C++17 compiler:

    g++ -O3 -std=c++17 check_four_inputs.cpp -o /tmp/check_four_inputs
    /tmp/check_four_inputs

## Boundaries, uniformity, and what remains open

- The n>=3 hypothesis is essential: for n=2, k=2 the function is OR (after
  translation/sign changes), with one gate, not n or n+k-2 as an equality.
- An arbitrary parameter sequence (v_n,r_n,c_n) need not define an NP language.
  If the parameters can be computed in time polynomial in n, the family is in
  P by direct evaluation and hence NP with empty witness. Define n=0,1,2
  separately, e.g. constant zero, for a total all-length language.
- A complex symmetry group, extreme Walsh spectrum, or large orbit does not
  by itself yield high general circuit complexity. Here the explicit
  affine-plus-single-defect representation is a short circuit description.
- The Ramsey bridge should seek a property stable under the richer
  independent-parameter substitutions used in current lower-bound work,
  not infer hardness merely from homogeneous-set avoidance or orbit size.
  That is a research criterion, not a new lower-bound theorem.

## Reproduction

    python3 check_point_defect.py
    python3 -O check_point_defect.py

The recorded output is `check_point_defect.json`. Scope: all three-variable B2
DAGs with <=3 gates; all four three-variable translated/complemented strata;
20,464 truth-table evaluations of the general construction through dimension10.
