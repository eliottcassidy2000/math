# The count 63 occurs at rooted gap one, not at rooted gap two

Write (d=R_+-H) for the absolute gap of the apex. The [contact-pair
budget](CONTACT_PAIR_DEFECT.md) gives (d≥ r(r-1)) for (r) apex
outcontacts. At gap two only one or two contacts are possible:
zero contacts give gap zero. The two-contact branch is already
classified by five rays, none containing 63. We close the remaining
branch here. This is not a new scalar exclusion.

## One contact and a power-of-two gap: an exact terminal reduction

Let (v) be the outer apex and let (u) be its unique outcontact.
By [the finite-support identity](ORDINAL_DEFECT.md), (d) counts
full puncture paths ending at (u). If (d>0), (u) belongs to
the terminal strong component of the puncture. Any earlier component
multiplies that count by its odd path count. If (d) is a power of
two, all these factors are one, so the earlier components are
singletons. They form a prefix that dominates both the remaining
graph and (v), and can be deleted without changing either (H)
or (R_+). The remaining puncture is strong; the singleton case contributes
only the special pair stated below.

In the nonsingleton puncture delete (u) and let (B) be the last
strong component. Paths ending at (u) correspond to paths of the
deletion whose last vertex beats (u). Again, any earlier component
contributes an odd factor dividing (d), hence contributes one. The
earlier components therefore form a transitive wing before (B).
Let (x) be its binary outcontact word from (u) (1 means (u)
beats that vertex), interpreted as an integer. Leading zeros would
dominate all of the puncture, so the minimal strong presentation has
no leading zeros. Address (x=0) means an empty wing.

Let (U) be the vertices of (B) beaten by (u). The endpoint
counts at the vertices of (B-U) must sum to (d). Every vertex
of a strong tournament can end a Hamiltonian path, including
a singleton, so (|B-U|≤ d). Conversely, choose (B,U,x) with
that endpoint sum, and require (U ≠ ∅) if (x=0);
this gives a strong inner puncture: (u) reaches the first wing
vertex (or some vertex of (B)), and a vertex of (B) beats (u).
The outer augmentation is also strong, since (v) beats (u)
and every other vertex beats (v). The only singleton-puncture
exception gives ((H,R_+)=(1,2)) and (d=1).

## Read the prescribed integer, not a residue

Use the inner profile of (B) with apex (u): (h) is the original
path count, (K) the augmented one, (f) the ordered two-color
count, (P,Q) the exposed start and dual end counts, and (T) the
number of full paths starting outside (U). Let (g(x)) be the
[transitive-wing augmented count](compile_absolute.py), with
(g(0)=1). Then the outer readout is exactly

```
H = h g(x) + (f+Q)x + K+P−d + T(x mod 2),     R_+ = H+d.
```

To prove it, write (D=wing→ B). Every subset path counted
by the outer start response starts at (u); stripping (u)
leaves either an empty path (with h_D choices for the complement)
or a path starting at an inner contact (with P_D total choices). Thus (H=K_D+h_D+P_D-d).
The ordinal equations give (P_D=fx+P) and
(K_D=h(g(x)-1)+K+Qx+T(x mod 2)). Both equations remain valid
for an empty wing.

Write (b=|B|). All nonempty subsets of (U) contribute at
least one to (P), so (H+d≥ P+1≥ 2^max(b−d,0)).
If (x>0), then (H+d≥ fx≥2^b x), so the number of bits
of (x) plus (b) is at most ⌊log₂(H+d)⌋+1.
In particular any one-contact pair with a positive power-of-two
gap has a witness of order at most
⌊log₂(H+d)⌋+d+2, if it has a witness
at all. The possibly much larger transitive prefix was neutral.
This bound can be used before any search; it is not claimed for
arbitrary numbers of contacts.

For (d=2), (B) is not a singleton, therefore (b≥3).
It has a nonempty (U) and a nonempty complement. Hence both
(B) and (B+u) are strong. The rotations of a Hamiltonian
cycle give (h≥ b), (K≥ b+1). The search bound used
here is then (H≥ 2^{b-2}+2b-2). For (H≤63), this reduces
the search to (b≤7), and each wing is tested only up to
((f+Q)x+K+P+h-2≤63).

## The exact comparison at 63

The [standalone replay](replay_gap_two_63.js) generates all
isomorphism classes of terminal blocks through order seven by
one-vertex extension and canonicalization within degree classes.
It finds 162 eligible marked terminal frames for gap two and
checks all surviving wings by literal subset-path dynamic
programming on the augmented tournament. The one-contact values
through 63 are exactly **9, 23, 29, 37, 51, 53, 59**. With the five two-contact rays, the complete gap-two list through 63 is
**5, 9, 11, 13, 19, 23, 25, 29, 31, 37, 39, 47, 49, 51, 53, 55, 59**.
In particular **(63,65) has no tournament witness**, at any order.

Yet **(63,64) does**: an eight-vertex strong tournament with
outgoing masks (190,252,248,176,160,192,25,64), the last vertex
being the outer apex. It is the singleton (B) with (U) empty
and wing address (x=19); literal counts give
((h_puncture,P_puncture,S,H,R_+)=(25,39,1,63,64)).
Thus even a difference of one in the requested response can separate
a realized scalar count from a forbidden joint pair.

The terminal reduction and comparison at 63 came from a background
pass; the explicit logarithmic order bound and standalone independent
replay are included here.
