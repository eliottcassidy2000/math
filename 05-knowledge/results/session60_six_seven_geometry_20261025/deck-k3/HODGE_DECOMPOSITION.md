# Seven K3 transcendental summands, not merely seven two-forms

**Status: PROVED over C.** Standard inputs are rational Hodge theory,
Lefschetz (1,1), finite-quotient transfer, and the surface facts and explicit
projection lifts proved in `DECK_CHARACTERS_AND_K3_QUOTIENTS.md`.
No integral decomposition, motivic isomorphism, or arithmetic descent is asserted.

Write T(V)_Q=NS(V)_Q^perp inside H^2(V,Q), with the cup-product form.
For a smooth projective surface its rational (1,1) classes are exactly
NS(V)_Q. This latter space is nondegenerate by the Hodge index theorem.

## 1. A quotient-transfer lemma with the resolution ambiguity removed

Suppose a finite group A acts on a smooth projective surface V, and W is
any smooth projective model of C(V)^A. The quotient correspondence induces

    T(W)_Q -> T(V)_Q^A

as an isomorphism of rational Hodge structures. Its pullback multiplies
the intersection pairing by |A|.

**Proof.** Resolve the graph to obtain V' -> V (a composition of point
blow-ups) and a generically finite morphism h:V' -> W of degree |A|.
Blowing up adds only rational divisor classes, so it leaves T unchanged.
Pullback preserves T: for t in T(W) and a divisor D on V',
(h^*t,D)=(t,h_*D)=0. The projection formula gives h_*h^*=|A|.
For the other composition, restrict to a nonempty open locus on which the
quotient is a finite etale A-cover. There the trace identity is
h^*h_*=sum_(a in A) a^*. The difference between the corresponding global
cycle identity and this open identity is supported over the complementary
curves and points. Such terms act on H^2 through divisor classes, or pair
against divisor classes, hence vanish on T. Thus the trace identity holds
on T. Its restriction to T(V)^A is |A| times the identity, proving
surjectivity and injectivity over Q. Finally
(h^*x,h^*y)=|A|(x,y) by projection formula. Changing the smooth model W
adds or removes only divisor classes and leaves this conclusion unchanged.

This lemma is deliberately on T, not on all H^2: exceptional divisors
would make the unqualified full-cohomology statement false.

## 2. Decomposition theorem

Let Y=Y_7, let sigma_j be its seven coordinate involutions, and let
f_j:Y -> S_j be the generically double-cover morphisms to the seven
six-root K3 resolutions. Then

    T(Y)_Q = orthogonal_direct_sum_(j=0)^6 f_j^*T(S_j)_Q.  (H1)

The j-th summand has character c_j=all indices except j, and its pairing
is twice the pairing on T(S_j)_Q.

**Proof.** The deck group G_7 has only rational-valued characters, so its
idempotents decompose H^2(Y,Q) into rational Hodge substructures. The
canonical-character theorem shows that only c_0,...,c_6 have any (2,0)
part. Because the characters are real, every other character also has
zero (0,2) part and is entirely of type (1,1), hence algebraic. Therefore
T(Y) is supported on the seven characters c_j.

The quotient-transfer lemma for <sigma_j> identifies T(S_j) with
T(Y)^{sigma_j}. Among the seven supported characters, c_j is the only
one trivial on sigma_j. It follows that this invariant subspace is exactly
the c_j summand. Distinct character summands are orthogonal: an element
on which their product character is -1 negates their pairing while
preserving the cup product. This proves (H1) and the scale factor two.

In particular, this is an actual correspondence of transcendental Hodge
structures, stronger than the equality of p_g values or the decomposition
of the seven-dimensional space H^{2,0} alone.

## 3. Picard-rank consequence for the hourglass

For Y_7, chi(O)=8, K^2=16. Noether's formula gives c_2=80; q=0 gives
b_2=78. Since each S_j is K3, b_2(S_j)=22. Thus (H1) gives the exact identity

    rho(Y_7) = sum_(j=0)^6 rho(S_j) - 76.                 (H2)

For the hourglass the central deletion S_0 has no branch triple points;
each of the six outer deletions has two. Formula (9) of the deck packet
therefore gives rho(S_0)>=16, rho(S_j)>=18 for j!=0 and already implies
rho(Y_7)>=48.

With the independently established central Kummer/CM result
rho(S_0)=20 and the explicit D12 transitivity on the six outer deletions,
write rho_o for their common Picard number. Then (H2) sharpens to

    rho(Y_7) = 6 rho_o - 56,
    rho_o in {18,19,20},
    rho(Y_7) in {52,58,64}.                              (H3)

The packet does NOT determine which value of rho_o occurs. An exact
calculation for one outer-deletion K3 would determine the whole hourglass
Picard rank. The central rank20 and D12 statements are supplied by the
separate hourglass-symmetry/CM packet, not assumed from a numerical match.

## 4. Six-root symplectic transfer

The same argument applies to the group H=(C2)^4 on a six-root K3 S.
Only its total-parity character can support T(S), and H lies in its kernel;
therefore H acts trivially on T(S). If Z resolves the double-plane quotient,

    T(Z)_Q ~= T(S)_Q,

under a correspondence whose pullback pairing is multiplied by sixteen.
This transfers Picard rank, but by itself does NOT transfer the property
of being a Kummer surface or give an integral lattice isometry. The
separate explicit Gale transformation in `KUMMER_GALE_BRIDGE.md` is what
places the central six-root model itself on the classical Kummer locus.

## 5. Ledger / audit boundary

Source -> target: seven-root surface -> seven six-root K3 surfaces.
Map: the sum of the actual deletion correspondences f_j^*.
Preserved: rational transcendental Hodge structures; cup form up to factor2.
Lost: algebraic divisors, integral discriminants, and arithmetic twists.
Sidecar: character projectors and the resolution graphs.
Hostile test: replacing T by H^2 adds exceptional and other divisor classes
and breaks the dimension identity; replacing <sigma_i,sigma_j> by the
single pair flip likewise changes the invariant space.
