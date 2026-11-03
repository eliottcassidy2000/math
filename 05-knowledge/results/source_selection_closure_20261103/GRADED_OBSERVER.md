# Enlarging the operation: eighteen polynomial attachment coordinates

**PROVED: an exact positive differential carrier for two-polarity births;
no bounded-dimensional linear observer exists even on strong sources.**
This is a wildcard operation test, not a claim of prime-lane coverage.

## G1. The cheap enlargement and the expensive enlargement differ

The inherited [common observer](../adaptive_state_induction_20261103/COMMON_OBSERVER.md)
already handles all twelve operations in which a newborn dominates every
forgotten vertex, has an arbitrary incoming mask on the two retained ports,
and replaces either retained port or is itself forgotten. In particular a
newborn losing to both ports costs no extra observation coordinate. The
branch repair in this packet uses precisely that inexpensive enlargement.

Allowing the opposite polarity as well is different. A negative-polarity
newborn loses to every forgotten vertex and still has an arbitrary incoming
mask on retained ports. There are now 24 commands, with physical names and
slot replacement retained. Forgotten vertices remain in the tournament.
Every old vertex has the declared arc to the newborn, so each command is an
actual tournament construction on any source with two distinct ports.

The old seventeen queries no longer determine every response. TT2 with its
two vertices marked and TT3 with its first two vertices marked have identical
old query vectors: the additional last singleton is an ordinal sink block
under every positive-polarity continuation. A negative-polarity birth that
beats only port1 has counts 1 and 3 respectively. This pair concerns actual
sources, one of which is weak; the stronger linear obstruction below uses
only strong sources and strong children.

## G2. Exact finite list of polynomial coordinates

A spanning directed path cover is an unordered set of nonempty paths covering
the source. Record each component's ordered pair of endpoint labels, using
0 and 1 for retained ports and * for forgotten endpoints. Remove all
(*,*) components and call the remaining multiset its skeleton s. Let m be
the number removed. A port can occur in at most one component; if it occurs
at both endpoints, that component is the singleton port.

There are exactly eighteen skeletons. Assign independently to each port one
of absent/internal, start-only, end-only, or singleton: sixteen possibilities.
The two additional skeletons are the single components (0,1) and (1,0).
For every skeleton retain the ordinary generating polynomial

    P_s(z) = sum_m c_(s,m) z^m,

where c counts actual path covers with that signature. No upper bound on m
is silently imposed. At order n there are at most 18(n+1) coefficients;
many are zero. The one-port version has four skeleton polynomials.

For each of the 24 commands there are explicit nonnegative integral
polynomial matrices A(z), B(z) such that

    P_new(z) = P(z) A(z) + P'(z) B(z).                  (G2)

The frozen receipt lists all coefficient entries, with source/target skeleton,
derivative order, power of z, and integer weight. Their polynomial degrees
are at most two. This claim is an exact carrier, not a minimality assertion
about eighteen polynomials or a bounded number of scalar coordinates.

**Proof.** In a new path cover, the newborn can have at most one incoming
and one outgoing old path component. Delete it to recover a unique old
cover and those choices. Joining the same old component on both sides would
make a directed cycle and is excluded. After joining, forget the selected
port and retain the specified names.

A (*,*) component can be selected only on the outgoing side in positive
polarity, and only on the incoming side in negative polarity. Thus at most
one such component is selected. Selecting none contributes a monomial times
P_s. Selecting one gives exactly m choices and contributes a monomial times
P_s', including the correct multiplicity of indistinguishable signatures.
Forgetting may create new (*,*) components but never deletes an unselected
one. At most two are created, since there is one new component and at most
one old port is forgotten. All coefficients count actual choices and are
nonnegative. The newborn-deletion bijection proves (G2) for all actual
sources and all finite command words, not just tested histories.

Hamilton count is recovered by summing c_(s,m) with |s|+m=1. The carrier
keeps paired endpoints and the complete forgotten-component multiplicity.
Marginal start/end counts or merely the total number of components would
lose the permitted joins. Actual source existence, strongness, and a legal
construction word remain separate obligations.

## G3. Why the unbounded degree is necessary for linear observation

For m>=1 let T_m be the strong tournament with two vertices a,b and a
transitive homogeneous module M of size m, with

    a -> M -> b -> a.

Mark only a. Let w_r consist of r repetitions of the following two commands:
add y of outdegree one, beating only the current root, and keep y;
then add z of indegree one, losing only to y, and keep z.
These are actual one-port commands. Every intermediate graph is strong:
adjoining a vertex with nonempty in- and out-neighborhoods to a strong
tournament preserves strongness.

Let f_r(m)=H(T_m after w_r). Then

    f_r(m) = sum_(j=1)^(r+2) c_(r,j) j! S(m,j),
    c_(r,j)>=0,                 c_(r,r+2)=(r+1)!.      (G3)

Here S(m,j) is a Stirling number of the second kind. To see the expansion,
contract each maximal run of module vertices in a Hamilton path to a block.
For a fixed exterior word with j block positions, assign the m transitive
vertices surjectively to those ordered blocks; each block is forced into
its transitive order. There are j! S(m,j) such assignments.

Every interior gap between blocks starts with one of b,y_1,...,y_r and ends
with one of a,z_1,...,z_r. Hence there are at most r+1 gaps and r+2 runs.
At the maximum, all exterior vertices lie in two-vertex gaps. In the possible
matching of gap starts to ends, b must pair with a; successively y_i must
pair with z_i (its only other possible z partner is z_(i-1), already used).
The r+1 pairs can appear in any order, separated by module runs. This proves
both existence and the exact top coefficient (r+1)!.

Since j!S(m,j) is a linear combination of 1^m,...,j^m with coefficient one
at j^m, the functions f_0,f_1,... are linearly independent: each has a new,
strictly larger exponential base. Therefore no finite-dimensional linear
future-H representation can handle both polarities uniformly over all
actual strong one-port sources. This excludes a fixed matrix enlargement
of the old seventeen-state observer, not a nonlinear encoding, a restricted
source family, or a finite-dimensional observer for one fixed source.

The first coefficient rows are

    r=0: (3,1)
    r=1: (13,12,2)
    r=2: (123,163,70,6)
    r=3: (1789,2946,1996,492,24).

They are exact finite controls of (G3), not the basis of its all-r proof.

## G4. Reproduction and connection boundary

Run `python3 -I -B graded_carrier.py`; its stdout is `graded_carrier.json`.
The independent implementation enumerates all eighteen skeletons and all
24 commands, tests 2,592 signature/multiplicity identities at m=0..5,
and compares all 288 one-birth transitions from twelve literal sources
against complete path-cover enumeration and a separate Hamilton-path DP.
The sources are all labelled tournaments of orders two and three, plus TT4
and the strong H5 order-four source. Strong clone controls use r=0..3,
m=1..6; the seventeen-query hostile uses TT2/TT3. No large support scan is
involved. The only inherited receipt read supplies the seventeen query words.

The connection contract is explicit: actual named sources map to their
eighteen endpoint-cover polynomials; all 24 future Hamilton responses are
preserved; source isomorphism and ancestry are lost; actual words and names
must accompany the coefficients. The cheapest hostile is TT2 versus TT3.
The practical conclusion is to try the already closed mask3 operation first.
If two-sided births become necessary, the forgotten-component polynomial is
a constructive attachment coordinate, with an unavoidable unbounded degree
for universal linear observation.
