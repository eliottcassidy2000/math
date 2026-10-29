# Marked and unmarked mass are different enumeration targets

**Status:** PROVED from the inherited aligned-pair/self-converse bijection and
orbit-stabilizer. Small values are independently FINITE-EXACT in
`check_transport.py`. The groupoid distinction is classical; no novelty claim
is made for orbit-stabilizer or inverse-automorphism weighting.

For a self-converse tournament W of odd order2n+1, write

    A(W)=Aut(W),
    i(W)=number of anti-involutions,
    c(W)=order of the centralizer of any anti-involution.

All anti-involutions are A(W)-conjugate, so c(W) is well-defined and

    i(W)c(W)=|A(W)|.

The aligned-pair class corresponding to W has automorphism group of order
c(W), not |A(W)|. Because there are2^(n²) labelled aligned-pair inputs on[n],

    sum over self-converse classes W of 1/c(W)
      = sum_W i(W)/|A(W)|
      = 2^(n²)/n!.

This is the **marked mass**, not the ordinary unmarked mass
sum_W1/|A(W)|. The class bijection preserves the number of classes, but not
ordinary groupoid cardinality.

For n=0 the singleton has one signed anti-involution, whose underlying
permutation is the identity. All three numbers are1. This uses the signed
symmetry convention, not the false claim that the identity permutation has
order2. For n=1 the transitive triple has(A,i,c)=(1,1,1) and the directed
triangle has(3,3,1). The marked mass is2 and the ordinary mass is4/3.

Complete small values from the checker:

| n | self-converse order | class count | marked mass | ordinary mass |
|---:|---:|---:|---:|---:|
|0|1|1|1|1|
|1|3|2|2|4/3|
|2|5|8|8|88/15|
|3|7|88|256/3|23776/315|

A direct labelled count gives the same formula. There are
(2n+1)!/(2^n n!) involution permutations with one fixed vertex; each admits
2^(n²+n) anti-invariant tournaments. Thus the number of labelled pairs
(W,chosen anti-involution) is

    (2n+1)! 2^(n²)/n!.

This is not the number of labelled self-converse tournaments, since W is
counted i(W) times. It is also not an algorithm for removing that multiplicity
without recovering the relevant symmetry information.

## Exact sampling consequence

A uniformly chosen n×n bit matrix, interpreted under the **natural aligned
pair action**, produces a lift class W with probability proportional to

    1/c(W)=i(W)/|A(W)|.

Uniform labelled self-converse tournaments induce class weights proportional
to1/|A(W)|. Therefore the raw lift sampler is biased by the number of
anti-involutions. At order3 it chooses transitive/cyclic classes with
probabilities1/2,1/2, whereas uniform labelled tournaments give3/4,1/4.

If i(W) is supplied exactly, rejection with probability1/i(W), followed by
uniform vertex relabeling, gives an exact uniform labelled self-converse
sample. This is a conditional procedure: computing i(W), and its rejection
rate, must be charged. It is not a proven faster sampler or a free symmetry
count. With n=1 the mean acceptance probability is2/3.

## What the sidecar buys

Knowing c(W) and i(W) recovers the full group order, but not the full group or
its action. The companion transport note supplies an exact coset transversal:
retain an actual base anti-involution t0, every alternate mark s, and the
actual centralizer permutations. Then

    Aut(W)=disjoint union_s sqrt(s t0) C_A(t0).

This is the constructive upgrade from a scalar mass correction to a complete
symmetry carrier. It applies when these source data are already known; it
does not hide the cost of discovering them.
