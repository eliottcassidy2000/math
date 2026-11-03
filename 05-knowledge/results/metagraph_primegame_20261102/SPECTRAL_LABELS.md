# Algebraic labels: specify the operator before using its spectrum

**OPEN identification.** The user's proposed labels for the three pure-blue
five-vertex merged-metagraph nodes are grouped as

    {0}, {+i sqrt(2),-i sqrt(2)}, {(-1+sqrt(3))/2,(-1-sqrt(3))/2};

the other five self-converse nodes have proposed labels

    2, i, -i, (-1+sqrt(5))/2, (-1-sqrt(5))/2.

The inherited `groups_of_five_primegame_merged_metagraph_20261031.md`
identifies the three pure-blue classes as TT5, 1->C3->1, and regular T5;
the eight self-converse classes and two converse pairs give ten merged
nodes. That counting statement is not a definition of an algebraic operator.
No assignment rule for the displayed labels was located in the bounded
recovery search. Clarification of their origin has been requested.

There are precise conditional constraints, useful for choosing that rule.

1. An eigenvalue of an integer matrix is an algebraic integer. The numbers
   beta=(-1+-sqrt(3))/2 have irreducible monic polynomial
   X^2+X-1/2 and conjugate product -1/2. They are not algebraic integers, so
   cannot be unscaled eigenvalues of an integer adjacency matrix.
2. For an ordinary tournament adjacency matrix A,
   A+A*=J-I. If Av=lambda v, then

       Re(lambda)=(|1* v|^2/||v||^2-1)/2 >= -1/2.

   Thus both (-1-sqrt(3))/2 and (-1-sqrt(5))/2 are excluded from that
   particular spectrum. A different signed or rationally normalized operator
   need not satisfy the same constraint.
3. The first five numbers have sum -1 and the second five sum +1. Neither
   list by itself can be the full spectrum of a five-by-five loopless
   adjacency matrix, whose trace is zero. This does not rule out assigning
   blocks or individual roots to metagraph nodes, which may be what is meant.

The repaired research question is therefore to name a port-sensitive operator
and a normalization, then prove its behavior under native complement and
re-marking. Choosing companion matrices for the desired polynomials would
realize the labels tautologically, but would not explain the tournament map.
The older [THM-505 spectral/packing correction](../../../01-canon/theorems/THM-505-the-ocf-non-spectral-defect-H-equals-spectral-skeleton-plus-witt-defect.md) is relevant: spectral data
alone can lose the overlap information needed for Hamilton-path counts.

In this session the PRIMEGAME-to-Fano map is already defined and tested;
this proposed spectral map remains a separate open lane. Cardinality or
root-list similarities are not used as dependencies of the birth theorem.

Concurrent integration supplies a defined positive alternative, credited to
[Session143's exact tournament spectra](../session143_magic_square_cyclic_bridges_20261103/tournaments/RESEARCH.md)
and its [integral reconstruction theorem](../session143_magic_square_cyclic_bridges_20261103/cyclic_integral_order/PROOF.md).
For the regular T5, A=P+P^2 and chi_A=X^5-5X^2-5X-2. Ordinary integer
polynomials in A have index11 in Z[P], but P=A-A circ A^2 reconstructs the
chosen shift integrally, including modulo11. This locates a genuine
operation-dependent loss and its repair. It does not assign the proposed
radicals to the n5 metagraph nodes; that separate interpretation remains OPEN.
