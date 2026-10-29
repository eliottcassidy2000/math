# Stabilizer-sensitive attachment recurrences

**PROVED** S1-S3. The binary A_k specialization is inherited from
[Session68 CLASSIFICATION, section5](../session68_tournament_composition_classification_20261026/session68-classification/CLASSIFICATION.md),
including its exact polynomial and the two rigid palettes. Existence of
regular subsets for arbitrary tournament automorphism groups is classical,
as documented in the [Gluck recovery](../signed_tournament_palette_rigidity_20261021/session23b-stabilizers/REPORT.md).
We make no novelty claim for distinguishing-colouring methods; see also
[Chan, wreath-product distinguishing actions](https://arxiv.org/abs/math/0601414).

## S1. Weights and stabilizer order close together

Let W_k be the tournament on k ternary digits, oriented by cyclic order at
the first differing digit. A_k is its full automorphism group; C_k centralizes
digit negation. A vertex palette has nonnegative integer cost, encoded by
Q(x)=sum_a x^cost(a), where colours of equal cost remain distinct. Define

    P_k(Q;x,z)=sum_(A_k-colouring orbits [f]) x^cost(f) z^log_3|Stab(f)|,
    B_k(Q;x,z)=sum_(C_k-colouring orbits [f]) x^cost(f) z^log_3|Stab(f)|.

With P_0=B_0=Q(x), the all-order recurrences are

    P_(k+1)(Q;x,z)
       =[P_k(Q;x,z)^3-P_k(Q;x^3,z^3)]/3
         +z P_k(Q;x^3,z^3),                            (S1)
    B_k(Q;x,z)=Q(x) product_(j=0..k-1) P_j(Q^2;x,z).   (S2)

Q^2 denotes the palette of **ordered pairs** of original colours, with
additive cost, not independent orbits. In(S1), substitution x->x^3 scales
the costs of all colours; it does not cube their multiplicities.

Proof: a nonconstant triple of child types has a free cyclic root action,
so divide its count by3 and multiply its child stabilizers. A constant
triple has stabilizer H wr C3 if its child's stabilizer is H, adding1 to
three times the exponent. For C_k, the fixed child carries C_(k-1), while
the folded opposite pair carries one A_(k-1) acting on the paired palette.
Independent group factors multiply their stabilizer orders. This proves
(S1),(S2) and explains why marginal side stabilizers cannot be sampled
independently.

The labelled-mass checks, with N=3^k, are

    3^((N-1)/2) P_k(Q;1,1/3)=Q(1)^N,
    3^((N-2k-1)/4) B_k(Q;1,1/3)=Q(1)^N.                (S3)

They follow directly from orbit-stabilizer, independently of the recurrence.
For binary Q=1+x, the C_2 stabilizer distribution after x=1 is160+32z:
160 free orbits of size3 plus32 fixed colourings account for all512 inputs.

## S2. Primitive necklaces occur locally, at every node

Extract the coefficient of z^0. For a weighted palette put

    R_0(Q;x)=Q(x),
    R_(k+1)(Q;x)=[R_k(Q;x)^3-R_k(Q;x^3)]/3.             (S4)

This is the length-three primitive-necklace transform applied recursively
to **free child types**. For a q-colour unweighted palette,

    r_0(q)=q,   r_(k+1)(q)=(r_k(q)^3-r_k(q))/3.

Thus r_k(2)=2 for every k, the inherited Session68 result. These are not
the ordinary primitive necklaces of length3^k: at length9 there are56
primitive binary necklaces, but only2 hierarchically free colour orbits.
Local child rotations are additional symmetries. The classical
[formal-series Witt transform](https://doi.org/10.1016/j.disc.2005.03.004)
provides context for(S4); the distinct recursive action is stated explicitly.

For the centralizer, the weighted regular-subset orbit polynomial is

    R_C,k(x)=(1+x) product_(j=0..k-1) R_j((1+x)^2;x).   (S5)

Here regular means **trivial setwise stabilizer**, not merely a pointwise
base. The first nonconstant paired level is

    R_1((1+x)^2;x)=2x+5x^2+6x^3+5x^4+2x^5.

At every later level its lowest exponent triples, and the leading
coefficient remains(2^3-2)/3=2. Therefore for k>=1 the least selected size is

    b_k=sum_(j=0..k-2)3^j=(3^(k-1)-1)/2,
    number of minimum regular-subset orbits=2^(k-1).     (S6)

The empty sum at k=1 gives b_1=0 and one orbit, since C_1 is trivial. At
k=0 the minimum is also0, with one orbit. The first sizes are0,1,4,13.
This coincides with the inherited [minimum pointwise base size](../session93_reversal_lucas_clocks_20261029/ports/PORT_ORBITS_AND_LCA.md), but(S5)
independently proves attainability by sets with no nontrivial set permutation.

| k | all binary C_k orbits | regular binary C_k orbits | minimum size | minimum-size orbits |
|---|---:|---:|---:|---:|
|1|8|8|0|1|
|2|192|160|1|2|
|3|887808|425600|4|4|
|4|29258450074632192|2670085775168000|13|8|

The k=4 row comes from an exact proved recurrence, not enumeration of2^81
inputs. Literal binary and three-colour stabilizers through k=2 check it
at the smaller scales in [check_stabilizers.py](check_stabilizers.py).

## S3. Future attachments must use the actual residual stabilizer

For a finite G-set V and named palettes P,Q, partition joint colourings
by their first coordinate. There is an exact orbit decomposition

    (P x Q)^V/G  = disjoint_union_([f] in P^V/G) Q^V/G_f. (S7)

Choose one representative per first orbit; its stabilizer G_f is exactly
the residual equivalence on the second colouring. Changing representative
conjugates the action, so the count is independent of that choice. With
palette sizes p,q this gives

    a_G(pq)=sum_([f]) a_(G_f)(q).                       (S8)

For C_2, two binary attachments have88064 joint orbits, whereas192^2=36864
independent marginal pairs. The complete conditional action supplies the
missing information. The script verifies(S8) by listing each stabilizer
and independently applying its cycle index.

Stabilizer **order** alone is still insufficient. The regular actions C9
and C3 x C3 have the same order and fixed-point moments, but10 versus12
three-subset orbits, as independently checked in
[Session94's audit](../session94_recursive_ports_crt_clocks_20261029/audit/AUDIT.md).
For conditional counting retain the stabilizer cycle action; for actual
joint composition retain its named transporter/relative frame. The recursive
key can supply a stabilizer portrait: product for a nonconstant child triple,
wreath product for a constant triple, and folded diagonal action for C_k.

This is the group-action version of the recovered
[conditional boundary-weight multiplication](../context_leaves_prime_weights_and_union_20261011d.md):
independence holds after the shared boundary/source is fixed, not before.
