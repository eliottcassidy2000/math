# Fixed-baseline charge rigidity with a moving apex

At fixed ordinary Hamilton count prior work identifies all Q₂ values
as a finite union of dyadic rays. For an apex cut, stronger anchoring
is natural: keep the full assembled count N fixed, and try to move
the offset without changing that baseline. Even allowing the ordinary
spine count to change does not create an additive fine-offset actuator.

Let v be the apex, T the assembled tournament, B=T−v, h=H(B), and
N=H(T). Mark the vertices of B beaten by v. As usual S and d are
the total and odd-minus-even weighted starts in marks, including
Hamilton counts of complements. Factor the strong components into
T=L→C→R, where C is the component containing v, and let D=C−v.
Let S_D,d_D be the weighted starts on D for the inherited marks,
and Q_R the two-cover count of R. Empty parts have ordinary count
and two-cover count one; empty D has S_D=d_D=0. Then

    h = H(L) H(D) H(R)
    N = H(L) H(C) H(R)
    S+h = H(L) Q_R (S_D+H(D))
    d = h                         if R is nonempty
    d = H(L) d_D                    otherwise.

Indeed B=L→D→R. All marks in L are zero and all marks in R are one.
Paths starting in D contribute H(L)Q_R S_D, allowing an optional
R path. Paths starting in R contribute H(L)H(D)(Q_R−H(R)).
When R is nonempty, signed optional paths in R cancel, and its
nonempty signed starts contribute H(R). Otherwise only D contributes.
These arguments also cover C singleton and either empty side.

For fixed N, every nonsingleton strong component has count at most N,
and its order is bounded by Busch's strong-tournament lower bound.
There are only finitely many possible ordered reduced lists at this
N: every nonsingleton component has count at least 3, so the list
has at most floor(log₃N) such components. At C use a finite list of
vertex positions (with a separate option C singleton). Unbounded
singleton multiplicity in L has no effect on the formulas; singleton
vertices in R contribute a factor 2 to Q_R. Their relative positions
within R are irrelevant to Q_R and h. Thus the full set of h,d,S+h
is a finite union of exact dyadic rays in S+h with constant h,d on
each ray, plus the finite boundary for empty R. Conversely their
inflation is literal: append sink singletons after the apex
component, retaining H and h and doubling S+h each time. If R was
empty, split off the zero-length member to account for the d change.

A wedge of finitely many dyadic rays cannot contain arbitrarily
long consecutive integer intervals in the charge S+h (in fact it
has only O_N(log X) points through X). This is not a statement about
profiles with only h fixed, disconnected external observers, or
multi-vertex selective recontact. To escape it, let the physical
baseline drift with the target, or retain a richer interface whose
choices are not localized by the SCC of one embedded apex.
