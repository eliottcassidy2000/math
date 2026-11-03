# Synthesis: repair the object before transporting its spectrum

**Session143 / checkpoint1.** Mathematical statuses below are separate from
publication status. The main square-of-squares conjecture remains **OPEN**.
The session began2026-11-02 23:55UTC; checks continued2026-11-03UTC.

## Correct count, incorrect spectral identification

The exact universe is all `2^10=1024` labeled tournaments on vertices0..4,
quotiented by all120 permutations. The results are12 classes,8 self-converse,
10 converse classes,9 adjacency characteristic polynomials, and1 regular
class. A Burnside argument proves the first three counts independently of
the enumeration. Determinants, power traces and cycle counts independently
agree on each characteristic polynomial. The primary McKay catalogues match
all three relevant lists, not just their counts.

A loopless tournament matrix has `tr(A)=tr(A^2)=0`. Any proposed adjacency
spectrum summing to1 or-1 fails before deeper group theory is considered.
The regular class instead has eigenvalues

    2, -1/2 +/- (i/2)sqrt(5+2sqrt(5)),
       -1/2 +/- (i/2)sqrt(5-2sqrt(5)).

Every tournament automorphism group has odd order, so a faithful order20
`5:4` cannot be an individual tournament automorphism group. The repair is
not to discard the group: retain the four translation-invariant orientations
on a chosen `C5` vertex action. The affine group acts on that family, with
translation kernel `C5`; its `C4` quotient realizes cyclotomic Galois action.

## A positive bridge with a specified operator

For a complex group representation and `g` of order5, Fourier decomposition
of its restriction to `<g>` proves the spectrum of `rho(g)+rho(g)^2` exactly.
When the four nonidentity character values equal `c`, its polynomial is

    (x-2)^((d+4c)/5) q(x)^((d-c)/5),
    q(x)=x^4+2x^3+4x^2+3x+1.

The cited Baby Monster character of degree4371 has values21 and-4 on `5A`
and `5B`. The respective multiplicities are891+4*870 and871+4*875. These
operators have the same eigenvalue support, but different multiplicities.
The spectrum alone forgets exactly the information that distinguishes those
classes. The order-five operator construction is generic, not uniquely
sporadic, and needs no unconstructed five-point Baby Monster action.

No explicit Baby Monster matrices or generator/intertwiner have been
constructed here. No character value is identified with a nonreal tournament
eigenvalue. The cited minimum faithful ordinary degree4371 excludes a
nontrivial24-dimensional characteristic-zero Baby Monster representation;
modular and nonlinear actions are outside that statement.

## Inherited mechanisms, remaining obligations

The theorem-level proofs and exact source identities are in the lane files.
This session uses the inherited odd-group inverse-pair orientation principle
(THM-3447), cyclic Cayley-algebra mechanism(THM-2532), and trace obstruction
(THM-1858), with full slugs in the tournament source-pin file. Historical
sources are provenance, not an automatic status promotion.

A literal eigenvalue-field inclusion into `Q(zeta17)` is impossible;
`Q(zeta5)` and `Q(zeta17)` intersect in `Q`. An `X0(13)` relation would need
an explicitly chosen geometric object and map, not a number-field/function-
field conflation. The 13-point `M13` puzzle retains a hole endpoint and legal
concatenation, so its move system is not silently a tournament.

Still open in the user's broad proposal: a specified lambda-cube referent;
a magic-square/octonion map that preserves all square and line-sum constraints;
any actual modular/Hecke correspondence; and, above all, a magic square of
nine distinct integer squares. Checkpoint1 does not announce that later
lanes are audited or that any of these obligations has been paid.
