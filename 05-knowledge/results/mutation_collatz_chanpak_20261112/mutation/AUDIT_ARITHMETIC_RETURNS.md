# Independent review of the clock guard and rooted triangle response

**VERIFIED proof audit.** Reviewed [arithmetic A1--A5](../arithmetic/PROOF.md)
and [returns R0--R3](../returns/PROOF.md) at their stated scopes. No
mathematical correction is required. This review does not rerun the older
mutation census or prove the inherited all-order7/21 exclusion anew.

## Arithmetic clock and factor guard

For m=4x+epsilon, one odd branch followed by j actual legal halvings ends
at4(3x+epsilon)/2^j+epsilon. Substitution gives exactly
[3m+(2^j+1)epsilon]/2^j. Modulo3 its numerator vanishes precisely when
2^j=-1, hence precisely for odd j. This calculation is independent of the
initial m modulo3, and uses a recorded ordinary-clock halving depth. The
two single-step branch formulas both take a3-divisible modulus out of that
residue class. No accelerated-clock or convergence conclusion follows.

The factor criterion N=3h with h attainable is both necessary and
sufficient for existence of an ordinal-triangle presentation: deletion
of that ordinal block yields an actual residual graph, and ordinal joining
a C3 gives the converse construction. It is not sufficient for a fixed
source to have the factor or to reach such a presentation by mutation.

The H63 hostile is therefore valid at every graph order: the required
quotient is21, which the inherited exclusion forbids. The cited H63 source
is an actual positive witness;43 and127 lie in the attained small interval
of the inherited exclusion/coverage theorem. The minus trajectory
11->32->16, hence43->127->63, obeys the stated ordinary signed clock.
It illustrates a failure of the physical guard, not a purported graph
operation between those three counts.

For N=7*3^r with r>=2, adjoining r-2 triangles to an H63 source gives the
lower bound. A presentation with r-1 or r triangles would have quotient21
or7; more than r is arithmetically impossible. This proves the exact
maximum r-2 and does not revive the refuted claim that this seven-tower
consists of gaps. Any admissible7p count cannot have an ordinal triangle:
divisibility would force p=3, whose target21 is excluded.

The nonzero-digit obstruction is also exact: both3H+2 and3H+4 are units
modulo3. Every finite composition of short return blocks preserves H, so
it cannot rearm an ordinal-triangle guard there. This is not an assertion
that all ordinary mutation words preserve H, or that other reset guards
are impossible.

## Full one-root polynomial at the front triangle

The triangle polynomial P_C(m)=m^3+2m follows from one-, two- and
three-bin occupancy. Requiring one named endpoint in the first bin gives
1+(m-1)+(m-1)^2=m^2-m+1. For C->K, binwise concatenation multiplies cover
counts; a first-bin terminal point in C forces that bin to contain no K
vertices, hence P_K(m-1), whereas a first-bin starting point in C allows
all m bins of K. For K->C these entrance and exit factors exchange.
The first-bin convention and the empty-K boundary are explicit and valid.

At m=1,2 the displayed root state is(H,F,s)=(3h,3h,h), and the deletion
count is h because deleting one triangle vertex leaves an ordinal TT2.
The inherited birth update gives(5h,6h,2h) and then(9h,11h,3h). Thus
the newly available three front roots have zero defect, and the next two
newborn defects are h and2h. Degree-one availability is correctly kept
separate from positive defect1/2.

The suffix factor argument is literal: every birth rooted inside A of
A->K keeps the whole K suffix untouched and dominated by the enlarged A.
The induced old triangle forces H(A)>=3 by insertion monotonicity. If
H(K)>1, neither prime nor7p factorization can result without an impossible
factor7. For a first mark inside K, cyclic rotation of the untouched C3
still acts freely on Hamilton words, so the output remains3-divisible.
These arguments do not cover operations crossing the ordinal interface.

## Two-root readiness classification

In the exposed-triangle normal form, a,b are the two successive newborns
from(W,c). Root a has start and deletion count d and defect zero. At b,
deleting b leaves the first birth, of count d+F_W-s, while paths starting
at b have count F_W. The difference is d-s. This agrees with the
indegree-one defect identity. Inserting c after the last in-neighbour on
each Hamilton path of R injects H(R) paths into W-paths not starting at c.
The proof explicitly assumes I nonempty, as required for exactly two
roots; otherwise c is a source and the defect vanishes.

If d-s is1 or2, H(R) must be1 because tournament path counts are positive
odd. Therefore R is transitive. Its unique order makes s either0 or1,
depending on whether its first vertex beats c. Defect1 forces(d,s)=(1,0),
so W is transitive and I is a nonempty initial segment of length j. The
cover count U_I is2^j-1, giving H(T)=1+2^(j+1).

Defect2 forces(d,s)=(3,1). The H3 SCC classification gives a C3 with
transitive padding, and c being a possible first vertex removes all prefix
padding. Since R itself is transitive, c joins exactly its first two
vertices to form the initial triangle. Thus I={r2}, with all later vertices
in the suffix; U_I=2 and H(T)=9. Both converse constructions satisfy the
stated conditions. The first family's values never have factor7 because
the powers of2 modulo7 are1,2,4.

Finally, READY section4 was checked directly. Its no-re-marking theorem
assumes a strong input of order at least4 and concerns the output of two
successive newborn births. That is exactly the scope used in R3. Its C3
exception produces the small H9 positive boundary and is not suppressed.

## Scope verdict

The arithmetic clock, existential factor guard, positive-cut normal form
and attachment defects remain distinct. The comparison supplies precise
obstructions and small real constructors without equating Collatz labels
with tournament counts or claiming recursive exact-target coverage.
