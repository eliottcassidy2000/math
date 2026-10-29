# A ternary root recurrence and a strict completion obstruction

Status: **PROVED, with FINITE-EXACT certificates for the native source.**
The proof and its alternative verification paths were developed in this
session; there is no claim of an independent-agent audit.

The full native ternary minimum is strictly greater than the sharp
radical bound. For the saved Paley-seven three-run source

    V = C3^2 + C9^6 + C81,        |V| = 3^18,
    Gamma = Aut(Paley7) x S3,

every equivariant embedding into a finite abelian group with a perfect
Gamma-invariant alternating Q/Z-valued pairing has index **at least
3^10**. Hyperbolic doubling still gives index3^18. Thus the only remaining
possible minimum indices are3^10,3^12,3^14,3^16,3^18. Attainment of3^10
is **OPEN**. This is an injective completion problem, not a claim of
ordinary-quotient compatibility or an arithmetic Selmer realization.

## Inheritance and the concept board

- Anchor: the native symplectic completion with all original210 port maps.
  The closest mechanism is [N4](FULL_NATIVE_COMPLETION.md): the common
  invariant radical is R=C3^3+C9+C27, of order3^8, and one form has exactly
  this radical. The pairing quotient is C3^6+C9^2.
- Niche: the source-fixed root recurrence [N7/N8](FULL_NATIVE_COMPLETION.md)
  at2. We extend its carry-and-character description to two order-three
  clocks and an involution. Roots and characters are retained together.
- Wildcard: the older carry obstruction
  [THM-2657, odometer-carry-root-lift-nonsplit-extension-and-cech-cocycle](../../../01-canon/theorems/THM-2657-odometer-carry-root-lift-nonsplit-extension-and-cech-cocycle.md)
  and [THM-2779, bockstein-symplectic-decoder-frame-torsor-and-heisenberg-root-degree-gate](../../../01-canon/theorems/THM-2779-bockstein-symplectic-decoder-frame-torsor-and-heisenberg-root-degree-gate.md).
  Their preserved coordinate is a carry or central phase. Their sources
  and actions are different from the present native module.

The six live concepts are: native ports; the radical and its dual;
root carries; joint clock relations; ordinary quotient maps; arithmetic
realization. The corrected near miss is that an optimal restriction form
need not have an equally small equivariant perfect completion. The
Heisenberg example in N5 is the small hostile;81/27 remains the native
tournament hostile. The least-used useful sidecar is the actual lift of
a quotient generator, including its carry under EACH clock.

## T1. Exact source-fixed recurrence at every finite size

Let A be ANY finite abelian3-group with an alternating pairing b and
pairing-preserving clocks F,R,S satisfying

    F^3=R^3=S^2=1,  FR=RF,  SF=FS,  SRS=R^-1.

Every Gamma-stable index-three overgroup, with chosen quotient generator,
has a presentation

    A' = (A + Zx)/(3x-v),
    Fx=x+w_F,    Rx=x+w_R,    Sx=epsilon*x,
    epsilon in {1,-1}.                                    (T1)

The involution can always be put in this form because2 is invertible on
A. The required and sufficient equations on v,w_F,w_R are

    Sv=epsilon*v,
    3w_F=(F-1)v,          3w_R=(R-1)v,
    (1+F+F^2)w_F=0,      (1+R+R^2)w_R=0,
    (S-epsilon)w_F=0,    (1+epsilon*RS)w_R=0,
    (F-1)w_R=(R-1)w_F.                              (T2)

These follow by applying each defining relation to x; conversely they
make the displayed lifts well-defined automorphisms with every required
relation. In the native setting the seven-cycle acts trivially on A and
on the quotient. Its lift must be trivial too: if T=1+N then N^2=0 and
T^7=1+7N, so N=0 on a3-group.

Choose v from the simultaneous F,R-fixed, epsilon-S eigenspace in A/3A,
then choose an exact S-eigenvector lift. These choices cover all source
extensions. Replacing x by x+a changes

    v -> v+3a,   w_F -> w_F+(F-1)a,   w_R -> w_R+(R-1)a.

For the normalized S-action one uses Sa=epsilon*a. Distinct marked
presentations may describe the same extension; no unlabelled counting
claim is made.

Write z(a)=b'(a,x) for the extra character of an extended alternating
pairing. Its exact necessary and sufficient conditions are

    3z(a)=b(a,v),   z(v)=0,
    epsilon*z(Sa)=z(a),
    z(Fa)+b(Fa,w_F)=z(a),
    z(Ra)+b(Ra,w_R)=z(a).                          (T3)

The first line makes the form descend through3x=v. The remaining lines
are precisely invariance on all pairs(a,x); pairs from A and(x,x) impose
no further constraints. Thus T1--T3 are a full recursive construction,
not just necessary tests. On elements(a,t),0<=t<3, addition carries v
whenever t+t' reaches3. The saved integer Smith chart transports the
new source injection, root, form, and clocks together.

For fixed b,epsilon,v these equations form an affine linear system over
a sufficiently large3-power modulus. A full local Smith elimination
retains nonunit pivots and supplies both a particular solution and the
entire homogeneous kernel. The checker audits this solver against every
literal fibre for all2-by2 matrices over Z/3 and Z/4. A positive control
constructs C3^2 -> C3+C9 -> C9^2, with a perfect pairing at the end.

The recurrence has the immediate odd-prime version with px=v and the
norm sums1+F+...+F^(p-1),1+R+...+R^(p-1), when F,R have order p and S
has the same conjugation relations. The implemented native census uses
p=3. No computational census for every odd prime is claimed.

## T2. Why a minimum-radical completion needs a decreasing first step

If A embeds into a perfect paired W with index |rad(b)|, then

    A-perp = rad(b) is contained in A.

Thus A is coisotropic. Every intermediate subgroup A<=B<=W is also
coisotropic, since B-perp is contained in A-perp and hence in B. For an
index-three intermediate B, its radical therefore has order |rad(b)|/3.
Equivalently the new character z must be nonzero on rad(b).

A Gamma-stable first intermediate always exists here. The normal
seven-cycle acts trivially on the native source. In any paired ambient,
its averaging projector(1+T+...+T^6)/7 is self-adjoint. Its fixed summand
is therefore perfect, contains the source, and is Gamma-stable because
the seven-cycle subgroup is normal. Discarding the other summand cannot
increase the index. Thus a minimum completion has T=1, and its action
factors through C3 x S3. The normal3-subgroup generated by F,R has a
nonzero fixed vector on the3-torsion
of W/A; its fixed space is stable under S. Since2 is invertible, that
space has a one-dimensional S-eigenline. Its inverse image is the needed
index-three extension. This also proves that repeated T1--T3 steps are
complete for this coisotropic completion question. More generally, the
same fixed-vector argument can be repeated on each intermediate quotient
without assuming coisotropy. Including every invariant starting form
and every legal carry/character choice therefore gives a finite exhaustive
search scheme for the minimum, stopping no later than hyperbolic doubling
at index |V|. The present checker implements the exact step equations
and first-step exclusion, not that full potentially enormous search.

The condition that z read the old radical is exact: if z is nonzero
there, its image has order3 and the extended radical has one-third the
old order. If z vanishes there, the extended radical has three times
the old order. In the latter case z lies in the image of the old form,
so after subtracting an element of A from x it is orthogonal to A; the
new radical maps onto the new quotient. This gives a useful search
budget: an index3^(8+2s) completion starting with radical3^8 has exactly
s radical-increasing steps and8+s decreasing steps.

## T3. The native first-step obstruction, allowing EVERY source form

Use the nine saved heterogeneous coordinates, of orders
(3,3,9,9,9,9,9,9,81). We do not treat their clock matrices as arbitrary
9-by9 matrices over Z/81. Each row has its own target order.

1. For each sign epsilon the root-class space in T1 has dimension2
   over F3. There are18 sign-marked root classes in total. All16 nonzero
   classes fail T2 before any pairing is imposed. The certificate for
   each failure is a literal left annihilator y of the group-equation
   matrix, with y*rhs nonzero modulo81. Consequently every allowed
   first overgroup splits as an ABELIAN group, A'=A+C3. Its clock action
   may still have nonzero carries.
2. With v=0, w_F,w_R belong to A[3], and z has order at most3. The carry
   cocycle spaces have dimensions6 and7 for epsilon=+1 and-1. Generator
   changes in A[3] contribute coboundary dimensions2 and3. There are
   therefore3^4=81 clock-carry classes for each sign.
3. The complete invariant-form module was retained in
   `primary_native_clocks.json`. Its restriction to A x A[3] has
   dimension3 over F3. All its36 saved module generators are checked to
   lie in the span of the three retained torsion-pairing images. No
   nonunit coefficient or higher-precision form is silently discarded:
   T3 with v=0 only reads this restricted image.
4. The common radical R maps to a SINGLE line in A/3A, spanned by the
   sixth source generator e_6. Every other saved generator of R is
   divisible by3 in A. Thus z can read R iff z(e_6) is nonzero.
5. For all162 carry classes, solve T3 with BOTH the nine character
   coefficients AND all three invariant-form coefficients variable.
   In every case the equations force z(e_6)=0. The receipt stores an
   explicit row combination equalling this coordinate functional.
   A second path checks row ranks on all729+2187=2916 literal carry
   pairs, including every gauge representative.

The gauge change alters z by b(-,a), which vanishes on the COMMON radical
for every invariant form. Thus passing to the81 carry classes does not
hide a radical-reading character. The exhaustive literal-carry audit
also verifies the conclusion without relying on that reduction.

If index3^8 were attainable, the restriction form would have radical
exactly R and T2 would require a radical-decreasing first step. T3 rules
out that step for every invariant restriction form. This contradiction
excludes index3^8. Since both |V| and the order of a perfect alternating
group are squares, the next possible index is3^10, proving the stated
bound. The obstruction is native; we do not identify it with the simpler
single commutator formula in the Heisenberg hostile.

## Preserved maps, consequences, and remaining obligations

The source is the native three-run cokernel, the target is a perfect
finite alternating module, and the map is an injective homomorphism
intertwining every named clock. The preserved predicate is the original
pairing after restriction. Forgetting the root or clock carries destroys
the extension equations. The cheapest decisive test was the first
coisotropic step, not enumeration of ambient paired groups of order3^26.

This gives a strict refinement of the earlier doubling/squaring picture:
square orders and even elementary-divisor multiplicities constrain the
object, but compatible source embeddings and clock lifts constrain its
composition. For the complete native primary product, the established
minimum has the form

    2^40 * 3^a * 11^6 * 79^3 * 5861^6 * 392809^3,
    a in {10,12,14,16,18}.

The factors at5 and17 require no enlargement. The value of a remains
unknown. Primary splitting is characteristic, so these bounds combine
without coupling different primes.

The source-fixed carries are structurally analogous to the retained
integer carry in the cheaper ternary return controller and to the older
13-adic odometer in THM-2657. There is no map here from those dynamics
into this finite module, and no new exclusion of Collatz exceptions.
Likewise an arithmetic Selmer transition would have to supply the actual
source diagram before applying this recurrence. This result proves
neither finite prime support nor bounded primary exponent for Sha.

Next decisive test: permit ONE radical-increasing step and test whether
nine later decreasing steps can complete the native source with its
saved optimal form. A successful chain proves index3^10 attainable.
Failure of a greedy chain is not an exclusion: retain all carry,
character, and source-form choices needed for a complete search or find
a new invariant obstructing every chain. A complete index3^10 exclusion
must ALSO consider invariant source forms with radical order3^10 and
ten decreasing steps. **Subsequent closure:** the complete
[729-form census](TERNARY_FORM_STRATA.md) now proves that no such forms
exist, and reduces the324 optimal forms to20 sufficient representatives.
The pure full-clock minimum and native attachment realization remain
separate problems. Sha finiteness, Collatz, and LRC(14) remain open.

Reproduce with `python -B ternary_root_extensions.py`, and repeat with
`python -B -O ternary_root_extensions.py`. The file writes
`ternary_root_extensions.json`; all119,542 checks remain active under -O.
Then run `python -B ternary_obstruction_verify.py` and its -O variant.
This separate verifier imports no extension solver, checks979 gates
including two hostile certificate mutations, and replays81/27. The
normal/optimized receipts are byte-identical; see
`ternary_extension_mode_receipts.json`. It is a read-only verification
of the mathematical certificate and source data, with a separate audit
receipt as its only output.
