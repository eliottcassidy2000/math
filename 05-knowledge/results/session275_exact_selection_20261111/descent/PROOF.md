# Projective descent distinguishes the two holes: seven is rigid, twenty-one is not

**Status: PROVED elementary classification from the inherited proper-graft law;
FINITE-EXACT literal witness and small-boundary controls. Chan–Pak remains OPEN.**
No claim of priority for the scalar count 609 or the seven-vertex tournament.
The new target is the previously unclosed *normalized cofactor 21* question,
not another prefix of the ordinary count spectrum.

## Inheritance and the change of direction

The five-core inverse of Session251 retains two unbounded binary addresses.
Its coefficient pullback grows the slope, so it is not a finite-state selector.
The older `inverse65_zeckendorf_20261103/ZECKENDORF.md` makes the analogous
point concretely: a numeral transducer needs a legal matrix cycle and its
physical boundary, not just the right scalar recurrence. The distant
`decoder_mahler_catalan_20260925.md` supplies the useful positive lesson:
retain the *ordered carry* and distinguish an exact route certificate from
an existence theorem. Here the analogue of that carry is gadget content.

The two-prime-lane reduction and its atoms 49,63,343 are already
THM-4094 (`01-canon/theorems/THM-4094-hamiltonian-matching-deficit-and-two-prime-lane-completeness.md`),
not new here. A targeted provenance search also recovered THM-2500's
finite-hole prime-ray theorem. This prompted a sharper question: can a
physical construction have a forbidden *normalized* cofactor even though
its actual count is legal?

THM-5107 (`01-canon/theorems/THM-5107-physical-graft-content-excludes-cofactor-seven.md`)
excludes normalized cofactor seven for every proper graft and explicitly
leaves twenty-one unexcluded. This note completely classifies that latter
case and supplies an actual witness. It does not undermine either
ordinary exclusion H=7 or H=21.

## 1. Definitions and inherited law

For a nonempty gadget B put H=H(B), and

    Q=sum_(X subset V(B)) H(B[X]) H(B[V(B)\X]), H(empty)=1.

For a source (A,a), require a to be a nonsink. Set D=H(A), d=H(A-a),
E=number of Hamilton paths ending at a, U=number of ordered spanning
2-path covers whose first path is nonempty and starts at a (the second path
may be empty), and f=U-D-d+E.
The proper graft has B→a and A-a→B. The inherited Session229/242 identity is

    N=Hd+Q(D-E)+(Q-H)f,
    D,d positive odd, 1<=d<=D, 0<=E<=min(d,(D-1)/2),
    f>=0 even.

Let g=gcd(H,Q), h=H/g, q=Q/g, and r=N/g. Thus

    r=hd+q(D-E)+(q-h)f >= h+q.                    (1)

If |B|>=2 then q is divisible by four, h is odd, and gcd(h,q)=1.
Moreover cutting a Hamilton path at each of its n+1 positions injects
(n+1)H objects into the Q ordered covers. Consequently

    q >= (n+1)h,       n=|B|.                    (2)

The law, inequalities and cut injection are inherited, not new.

## 2. Exact cofactor-21 classification (PROVED)

**Theorem.** A proper graft has N/g=21 if and only if

1. Q(B)=20H(B); and
2. A is transitive and a is its penultimate vertex in the transitive order.

Arbitrary vertices before a are allowed. Every such output has N=21H(B),
and every gadget occurring in this classification has at most 19 vertices.

**Proof.** A singleton B has g=1, so r=21 would contradict the established
ordinary exclusion. Assume n>=2. If n>=5, (1)-(2) give 7h<=21.
The possibility h=3 would require 18<=q<=18, incompatible with 4|q.
Thus h=1. For n<=4, the complete pairs (H,Q) are

    n=2: (1,4);
    n=3: (1,8), (3,12);
    n=4: (1,16), (3,24), (5,32).

These follow from the four elementary order-four isomorphism types, or
independently from the finite complete enumeration in check.py. The last
pair has h+q=37>21; all others have h=1. Hence in every case h=1 and
q belongs to {4,8,12,16,20}.

For q=4,8,16 the normalized pair (1,q) is the actual packet of a transitive
gadget of order 2,3,4. Replacing B by that gadget in the unchanged source
would give an actual tournament of H=21, impossible. This uses physical
normalization, not division of a scalar.

For q=12 or20, D>=3 would force D-E>=2 and r>=1+2q>21.
Therefore D=1, d=1, E=0, and

    21=1+q+(q-1)f.

With nonnegative even f this forces q=20 and f=0. A tournament with H=1 is
transitive. Let t>=1 be the number of its vertices after the nonsink a.
Its first path in the definition of U can contain any subset of those t
vertices, and the remaining vertices form the second path uniquely. Thus
U=2^t and f=2^t-2. The equation f=0 is exactly t=1: a is penultimate.
Conversely this source has (D,d,E,f)=(1,1,0,0), so its graft is H+Q=21H
when Q=20H. Finally (2) gives n+1<=20, hence n<=19. QED.

The ordinary H=1 iff transitive fact also follows from the inherited
positive odd-cycle expansion: a directed triangle already adds paths.
No all-order census is used.

## 3. An explicit non-dyadic integral packet (FINITE-EXACT witness)

Use outgoing row masks

    B=(108,41,2,4,47,12,62).

Direct subset DP gives (H,Q)=(29,580). Independent permutation counting
checks H(B)=29. Grafting onto the source of TT2 gives

    T=(2,508,433,165,9,17,189,49,249),
    H(T)=609=21*29,  gcd(H(B),Q(B))=29.

The nine-vertex permutation count independently gives 609, without using
the graft law. This is a literal counterexample to extending THM-5107's
normalized-seven exclusion to normalized twenty-one. Its normalized
packet (1,20) is nonphysical: H=1 forces a transitive tournament, for which
Q is a power of two.

This also clarifies a dangerous apparent descent: 609/29=21 is arithmetically
smaller, but not realizable. A correct descent must keep the integral lift
29 of the nonphysical packet (1,20). The target itself is not new; the
certified obstruction concerns the proposed normalization operation.

## 4. An infinite controlled lift family (PROVED)

Ordered joins multiply both H and Q: each of the two paths in an ordered
cover independently splits across the two consecutive components. Therefore
an ordinal join of m>=1 copies of B and a transitive k>=0 block has

    H=29^m,     Q=29^m * 20^m * 2^k.

The same TT2 graft gives strong outputs

    N_(m,k)=29^m(1+20^m*2^k),
    gcd(H,Q)=29^m,   N_(m,k)/g=1+20^m*2^k.       (3)

Strongness is literal: the quotient is the directed triangle a→v→B→a,
so any module vertices can reach each other through a and v.
Only (m,k)=(1,0) has normalized cofactor21. Formula (3) is not a claim
that all normalized cofactors are physical counts, nor that they cover
all allowed odd values. It exposes an indefinitely reusable *lifted*
interface in which forgetting the lift is demonstrably unsound.

## 5. What changed for exact selection

For any proper graft and normalized r, (1)-(2) imply

    |B| <= r-2.

Thus projective target descent still has an ordinary-height obligation,
but the gadget side is bounded by the *normalized target*, not by N.
At r=21 the equality classification is much stronger: all counterexamples
are exactly the finite gadget packet fibre Q=20H with the penultimate-root
transitive source. There are finitely many such gadget isomorphism classes;
we do not enumerate that fibre through order19 here.

Seven and twenty-one now play different roles in this construction:

- 7 is excluded even after dividing by gadget content (inherited theorem).
- 21 is excluded as an actual count but occurs after division, in exactly
  the packet/source family above (this note).

The useful new principle is **lift-sensitive descent**: an inverse may
reduce the target only if it transports an actual normalized interface,
not merely a rational ray or its gcd. This is the same kind of loss as
ordered carries in the Mahler decoder, but the map and counterexample are
now explicit tournament constructions.

The constructive next obligation remains positive: find a return that
changes the lift while keeping the desired absolute count, or prove that
an appropriate prime/7p target always has a content-compatible predecessor.
Neither the finite fibre bound nor (3) discharges that obligation.

## Reproduce

    python3 05-knowledge/results/session275_exact_selection_20261111/descent/check.py

The script validates B and its graft by independent permutation counting,
recomputes all labelled packets through order4, and checks four instances
of the all-order ordinal-product family by ordinary subset DP. Assertions
are not used for required checks. There is no larger coverage experiment.
