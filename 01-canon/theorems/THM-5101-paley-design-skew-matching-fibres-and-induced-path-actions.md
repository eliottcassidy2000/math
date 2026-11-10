---
id: THM-5101
title: "Paley design skew-matching fibres and the exact induced Hamilton-path action"
status: "PROVED general restoration/transport statements; independently audited FINITE-EXACT fibres at7,11,19 and all32 eleven-circulants. No extremality or Chan-Pak completeness claim."
session: Session240, 2026-11-10
source: 05-knowledge/results/session240_seven_eleven_symmetry_20261110/paley/REPORT.md
audit: 05-knowledge/results/session240_seven_eleven_symmetry_20261110/audit_design/AUDIT.md
---

# THM-5101: eight and twelve are the missing orientation coordinate

[Proof and exact census](../../05-knowledge/results/session240_seven_eleven_symmetry_20261110/paley/REPORT.md)
· [Independent audit](../../05-knowledge/results/session240_seven_eleven_symmetry_20261110/audit_design/AUDIT.md).

For an incidence structure with n points and n distinct blocks, a compatible
skew matching f is a point/block bijection with x not in f(x) and exactly one
of y in f(x), x in f(y) for distinct x,y. These matchings correspond exactly
to tournaments with the specified out-neighborhood block set. The side-preserving design
group G acts by f -> beta f alpha^(-1), and G_f=Aut(T_f).

For Paley designs on7 and11 points, complete independent matching searches
find exactly8 and12 matchings, respectively. Each is one faithful
2-transitive orbit under its design group, of orders168 and660, and has
stabilizer21 or55. The PSL identifications and exceptional original point
actions are inherited from corrected HYP-9357, not inferred from orders.
The new finite content is the complete matching fibres and their actions.

If the matchings form one orbit, the disjoint union of all pairs(f,h), with
h a Hamiltonian path of T_f, is the induced set G x_(G_f) HP(T_f). Its G-action
is free, and its quotient is canonically HP(T_f)/G_f. At7 and11 this gives
1512/168=189/21=9 and1141140/660=95095/55=1729. The enlarged group therefore
produces no improved divisibility denominator for a fixed raw H.

Endpoint counts in P7 are1 versus8, and in P11 are628 versus1101, according
to the square/nonsquare difference from the initial vertex. These directly
witness the information lost by the bare2-transitive design action.

Hostiles: Paley19 has only one compatible matching and group171. Among all32
eleven-circulants,20 have one matching,10 have two, and only the two Paley
orientations have12. The p+1 pattern cannot be extrapolated to all primes.
No global Paley maximality or Chan-Pak coverage claim follows.
