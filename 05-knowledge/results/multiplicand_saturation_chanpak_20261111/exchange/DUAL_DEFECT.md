# Two positive decompositions of the same quantitative slack

2026-11-11. **PROVED consequence of two independent proofs.** The incoming
Session261 commits40eba89fd5 and d42d89d3f9 were read during close-out.
Their [even-component proof](../../session261_creative_chanpak_20261111/bridges/EXTENSION4.md)
proves quantitative domination without an exchange graph. Our
[star proof](STARS.md) separately supplies the literal simultaneous
injection. Neither proof is represented as a verification of the other.

Keep U,A,B,z,G and the path–cut incidence graph of PROOF.md. Define

    sigma=H(U+z)-H(U)-d(G)+1.

For an even nonempty vertex set W, let p(W) be the number of Hamilton
paths of U[W] starting in A and ending in B. Let epsilon(W)=1 exactly
when every connected component of G[W] has even order, and let kappa(W)
be its number of components. Session261 proves, by grouped
inclusion-exclusion and the marked odd-cycle formula,

    sigma = sum_W [2*p(W)*H(U-W)-epsilon(W)*2^kappa(W)]. (1)

Each bracket is nonnegative. For admissible W, a canonical endpoint path
in each component can be concatenated in all kappa(W)! orders, because
there are no cross-component A-to-B edges. Thus p(W)>=kappa(W)! and
2*kappa(W)!>=2^kappa(W). Nonadmissible W has no subtraction.

Our incidence identity and star theorem give simultaneously

    sigma = u+(1-tau)+E,                              (2)

where u counts unused newborn Hamilton paths, tau<=1 counts closed tree
components, and E is the nonnegative surplus of all other components.
These are exact formulas for the same integer. The coordinate change
preserves total slack; it does not identify individual W fibres with
individual incidence components.

## Equality and obstructions can now be transferred

Equality in (1) holds exactly when every admissible W has kappa<=2,
p(W)=kappa!, and H(U-W)=1, while every even nonadmissible W has p(W)=0.
Equality in (2) holds exactly when there is one closed star, u=0, and
every other component is either closed unicyclic or a tree with one stub.
Consequently these two descriptions characterize the same rooted
tournaments. No claim that either description is a small computed list
is made.

For example, an induced3K2 in G supplies slack at least4 in (1), so
the incidence resources in (2) must also total at least4. More generally,
an induced kK2 and all subsets of its edges give the inherited bound

    u+(1-tau)+E >= sum_(j=1)^k binom(k,j)*(2*j!-2^j).   (3)

The first values for k=1,...,5 are0,0,4,48,408. Distinct edge subsets give
distinct W, so no fibre is counted twice. This locates forced surplus
without pretending that the inequality selects an exact path increment.

The incoming
[finite equality-basis theorem](../../session261_creative_chanpak_20261111/bridges/FINITE_EQUALITY_BASIS.md)
therefore applies to the incidence equality architecture as well. It
reduces positive slack either to a six-vertex triangle-complement witness
or to a word over eight letters on a fixed triangle plus a transitive
spine. Subsequence monotonicity gives a finite basis; exact weighted
automata make synthesis terminate. The entire basis and its maximum
order remain uncomputed. This is a second way to make an infinite class
manageable, complementary to the fixed-grade arithmetic compiler.

## The incoming two-port kernel uses the same retained module data

The independent
[bilinear kernel](../../session261_creative_chanpak_20261111/constructors/BILINEAR_KERNEL.md)
constructs an actual tournament from modules A,B with

    H_out=149hg+72(cg+hd)+37cd,

where h=H(A), g=H(B), c=C2(A), d=C2(B). In the coordinates used by our
offset carriers, q=Q2(A)=c+2h and r=Q2(B)=d+2g, this is exactly

    H_out=9hg-2(qg+hr)+37qr.                          (4)

The first expression has positive run coefficients; the second shows
that both constructions retain the same two input statistics. Swapping
modules preserves the count but need not preserve selected-port state.
Composing (4) with our one-port carriers still requires the output's Q2;
H_out alone does not close the recursive state. A useful next constructor
must control that missing output coordinate or preserve a marked return.
This is the precise remaining transfer gate, not a scalar saturation
claim or a completed proof of Chan–Pak.
