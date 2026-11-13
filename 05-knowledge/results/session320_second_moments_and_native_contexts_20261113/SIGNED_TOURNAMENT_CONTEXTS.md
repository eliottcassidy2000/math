# Opposite signs from common strong tournament contexts

The [Rule 30 context theorem](RULE30_CONTEXT.md) gives both signs at each
sufficient depth. The [Session 319 neutral switch](../session319_arc_usage_and_signed_cycles_20261113/RESERVOIR_POLYNOMIAL.md)
cannot behave this way under common homogeneous contexts, because its
nonempty-cover profiles are componentwise ordered. Here is a literal
incomparable seed that escapes that particular obstruction.

The outgoing row masks of two strong nine-vertex tournaments are
```
A = (452,349,136,385,141,351,412,34,148)
B = (16,413,209,277,160,15,315,105,181).
```
Their (H,Q(2),Q(3)) are (931,19556,189849) and (931,19388,191097).
So they agree at one slot, and the respective differences are
ΔQ(2)=168, ΔQ(3)=−1248. In particular
ΔF_2=168 and ΔF_3=−1752: they are genuinely incomparable in
nonempty-cover coordinates.

Let W(T)=C3[T,1,1] be the inherited cyclic wrap. Its slot identity is
Q_W(z)=z²Q(z)+z[Q(z+1)−Q(z−1)]. Thus at equal H a single wrap
reads ΔQ(2), whereas two wraps read 5ΔQ(2)+2ΔQ(3).
The two literal signs are +168 and −1656. This is not a comparison
of two unrelated contexts on the two members: within each comparison
the context is identical. Both contexts give strong tournaments.

More generally first ordinally pad each source with TT_m, m>=0.
This leaves H fixed, and multiplies ΔQ(2), ΔQ(3) by 2^m,3^m
respectively. Consequently the one-wrap differences are +168·2^m,
and the two-wrap differences are 840·2^m−2496·3^m, negative for
every m>=0. Thus the sign switch persists for arbitrarily large
literal strong tournaments. If desired compare at equal final order
by padding the single-wrap inputs by two extra transitive vertices
before their wrap; the positive sign persists. No exact common
strong equalizer is implied by the two signs.

For a quotient with marked indegree two, its three-run coefficient
and two-run coefficient c_3,c_2 are enough to state the precise
next test: substitution in the two sources gives equal H exactly if
168 c_2−1752 c_3=0, or c_2:c_3=73:7 if both are positive.
The scalar H of the starting pair is not an interface: a different
run response changes which continuation can cancel it.

## A literal exact return, and all three signs at the same added size

For a second pair the required ratio is realizable. Use
```
X = (388,449,202,195,335,351,129,48,204)
Y = (284,33,506,418,138,17,59,355,114).
```
Their profiles (H,Q(2),Q(3)) are (1047,21492,204453) and
(1047,21252,205173); thus ΔQ(2)=240, ΔQ(3)=−720 and
ΔF_3=−1440. A single wrap gives difference +240, two give −240.
Now use the common strong five-vertex quotient with outgoing masks
`(24,9,3,4,14)` and substitute in its vertex 0. It has marked
indegree two, and literal layout counts (c_1,c_2,c_3)=(13,12,2).
The marked block can have at most three runs since each entry uses
a distinct incoming exterior vertex. The output count is exactly
```
13H+12F_2+2F_3 = -5H+6Q(2)+2Q(3).
```
Both strong outputs therefore have **532623** paths, and direct
subset DP checks the count independently of the layout formula.
They are not isomorphic or converses: their Q(2) values are 14875056
and 15009456. The same ending was achieved, but the next two-slot
response is still different.

The two-wrap and exact-return contexts each add four vertices.
A positive context of the same added size is one wrap after
padding the source ordinally by TT_2, giving +960. Thus there
are three *physical strong* common contexts of identical output
order 13 realizing positive, zero and negative discrepancies
on one pair of starting tournaments. This is a tournament analogue
of one Rule 30 three-way comparison, not Rule 30's universal
all-depth statement. In particular the current proof supplies
neither a repeatable equalizer with restored readiness nor a
construction at arbitrary prescribed targets.
