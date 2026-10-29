# The smallest response carrier for rational attachment queries

Status: **PROVED elementary**, relative to the exact query language below.
This is not the smallest carrier for arbitrary tournament attachments,
integral cokernels, or path-native complement.

## The allowed future operations

An object is a nonsingular even-dimensional integral alternating matrix M
with named integer port columns P. Allowed operations are:

1. attach an integral alternating coupling E by M -> M+P E P^T;
2. re-mark the retained span by P -> PR with a specified integer R;
3. take global sign reversal M -> -M;
4. ask for the resulting determinant and whether the next matrix is
   singular. A singular result stops a branch of this query language.

No new port outside the retained span is introduced. For this language the
exact response carrier is

    (delta,G) = (det M, P^T M^-1 P).

Its update laws are P3-P4 of [TWO_PORTS](TWO_PORTS.md), re-marking
G -> R^TGR, and sign reversal G -> -G with delta unchanged. These laws
prove that equal carriers have equal observations along every allowed
finite operation sequence, including singularity decisions.

## Separation proves minimality

The empty query reads delta. For ports i,j, use the two couplings

    E_+ = e_i e_j^T-e_j e_i^T,     E_-=-E_+.

If g=G_ij, their determinant ratios are `(1-g)^2` and `(1+g)^2`.
Consequently

    g = (det(M+P E_- P^T)-det(M+P E_+ P^T))/(4 det M).

Thus the empty query and two one-step tests per port pair recover the
whole carrier. Two objects are response-equivalent **iff** their
`(delta,G)` agree. Any sufficient representation must distinguish these
equivalence classes. This is information minimality relative to these
queries, not a claim that rational entries admit a bounded bit encoding.
If a tested matrix is singular its determinant zero is a valid terminal
observation, so this case does not invalidate separation.

Using oriented Pfaffians would recover g in one signed query, but requires
the vertex-order sign sidecar. The Pfaffian sign correction in this
session is a literal example of why that sidecar matters.

## What remains invisible

Take

    M_a = J + 9J + J + J,
    M_b = J + 3J + 3J + J,

where plus means block-diagonal sum, and retain just the two coordinates
of the first J as ports. Both have delta=81 and G=-J. Every allowed
future determinant query therefore agrees. Their paired cokernels are
H_9 and H_3+H_3, respectively. This is an explicit hostile to using the
minimal response carrier as a full integral or Selmer carrier.

Likewise the tournament 81/27 hostile concerns automorphism stabilizers
and cyclic triples, which these determinant queries do not observe. Its
full transporter data must remain whenever that is a future operation.
The equality of some numbers in the two controls creates no equivalence
between their mathematical predicates.

## Keeping a clock changes the observation language

If T is an isometry of M and the ports are independent with
`T^T P=PC`, the coupling preserves T exactly when

    C E C^T=E.

Indeed the update term transforms to P C E C^T P^T. This is a clock on
port **functionals**; using TP=PC instead would put the transpose on the
wrong side. Under an invertible integer re-marking R, transport the clock
as C -> R^-1 C R. Under global sign reversal it remains the same.

For dependent ports one must retain ker P and test
`P(CEC^T-E)P^T=0`; G alone need not reveal these dependencies. Also the
two separating probes above may fail the clock constraint. Therefore
the full G is sufficient but is not automatically minimal for a
clock-restricted query language. The exact H9/H27 fixed-clock census is
a concrete example of fewer legal observations. Determining its smallest
carrier under the full intended re-marking operations remains a separate
problem.

The bridge to Selmer computations is precise but limited: local
character matrices also retain only what their future tests see. The
August rank30 certificate is sufficient for independence, but it loses
point signs and cannot determine a full Mordell-Weil group or Sha.
