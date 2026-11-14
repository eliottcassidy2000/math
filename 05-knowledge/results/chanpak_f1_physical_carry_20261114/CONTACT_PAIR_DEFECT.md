# Contact-pair contraction and exact two-contact gap-two rays

Use the [finite-support defect](ORDINAL_DEFECT.md) of a puncture
B with apex v. The apex beats precisely the contact set U.
Here H is the augmented count, h the puncture count, P the
unsigned start response, and R_+=h+P. The identity
R_+−H=E+S counts rejected complementary path pairs E and
full puncture paths ending in U, counted by S.

## Two contractions per contact pair

For each unordered pair {u,w} ⊆ U, E_uw counts the summands of E
whose two endpoints are u,w in either order. Delete the two
endpoints and substitute a single marked vertex. In the AND
contraction, an old vertex beats the mark iff it beats both
contacts; in the OR contraction, iff it beats either contact.
Then

```
E_uw = H(AND) + H(OR) + J.
```

Indeed concatenate the two complementary paths after identifying
their two endpoints into a marked position. An ordering of the
remaining vertices around the mark lifts by a distinct choice of
contacts (a,b) from A×B', with A the contacts beaten by the
predecessor and B' the contacts beating the successor. With no
predecessor or successor the corresponding set is full. Count one
lift for |A|=2 and |B'|≥1 (AND), one for |B'|=2 and |A|≥1 (OR),
and one if A and B' are distinct singleton sets (J). These count
all the lifts including two when both sets are full. In particular
the contact-to-contact arc is invisible to E_uw, though not to H.

Each contracted tournament has at least one Hamiltonian path,
hence E_uw≥2. The second difference of augmented H over the four
apex masks empty, {u}, {w}, {u,w} equals −E_uw: h+P−S is affine
in the mask and E has one term per selected pair. Rédei parity
then makes E_uw even. Since E is the sum over unordered pairs,
the **whole tournament** satisfies

```
R_+ − H ≥ |U|(|U|−1).
```

No ordinal assumption is needed. This is a second independent
budget from the source-weighted ordinal reader bound.

## The minimum contact-pair frame

Equality E_uw=2 requires both contractions to be transitive
(the only tournaments with one Hamiltonian path) and J=0.
Both are transitive iff the remainder is transitive and its
order has three successive groups: vertices beating both
contacts, vertices beating exactly one, vertices beating neither.
In that order J=0 iff the middle has at most one vertex of each
mixed type. Indeed a placement with distinct singleton sets
needs two mixed endpoints of the same type. Conversely two
such vertices give a placement: end the left piece at the
later one, start the right piece at the earlier one, and put
intermediate vertices in either piece in their inherited order.
There is no old cross-edge condition at the marked cut.

For exactly two contacts, gap two further requires S=0.
In the minimum frame this is equivalent to a nonempty last group:
if nonempty, it is an ordinal suffix with no contact. Otherwise
every vertex of the remainder beats a contact, so the terminal
strong component contains a contact. A Hamiltonian path can end
at any vertex of the terminal strong component (break a Hamiltonian
cycle there; the singleton case is immediate). Thus S>0 in that case.

## The complete five rays, not just a finite census

A prefix beating both contacts would beat every later vertex and
the apex, so it is absent in a strong witness. Even if permitted
in a nonstrong witness it leaves both H and R_+ unchanged: it
has h=1 and P=0. The common-outneighbor suffix is transitive of
length s≥1. Let A consist of the two contacts and the zero, one,
or two mixed vertices. From the no-contact reader law,
R_+=h_A+2^s P_A and H=R_+−2.

The entire finite middle calculation is below. Outgoing masks
are in bit-index order with the two contacts at indices 0 and 1.
Interchanging the contact names handles the other mixed order;
each listed core has E=2. The optional arc between contacts
is represented by the distinct rows.

| Mixed vertices | Core outgoing masks | h_A | P_A | H |
|---:|---|---:|---:|---|
| 0 | (2,0) | 1 | 3 | 3·2^s−1 |
| 1 | (0,5,1) | 1 | 5 | 5·2^s−1 |
| 1 | (2,4,1) | 3 | 6 | 6·2^s+1 |
| 2 | (8,5,9,2) | 5 | 13 | 13·2^s+3 |
| 2 | (10,4,9,2) | 5 | 14 | 14·2^s+3 |

Each ray is physically attained for every s≥1; the apex reaches
both contacts, every mixed vertex is beaten by a contact, and
all core vertices reach the suffix, which beats the apex.
Given H, subtract the ray intercept and test that the quotient
by its coefficient is an integral power of two at least two.
The union has O(log X) counts through X and contains none of
7, 21, 63. Thus the known count 63 needs a different interface;
this two-contact classification is not a scalar exclusion for 63.
Gap two with **one** contact is not covered by these rays: it
requires exactly two full puncture paths ending at that contact.

For comparison, a strong diagonal witness (R_+=H) must have one
contact u, and u must lie in the first strong component A of the
puncture, with a nonempty ordinal remainder D. Necessity follows
from the contact bound, reachability from the apex, and S=0.
Conversely these conditions suffice, by the ordinal orientation
and because no full puncture path ends at u. Its exact diagonal
formula is h_A h_D + P_A(u) f_D.

The two contractions and minimum-frame criterion came from a
background pass. The exact finite core table and all-order ray
decoder are a further consequence. Run
**python3 replay_contact_pair_defect.py** for 10,650 literal
marked pairs (1,922 minima) and 25 literal ray witnesses.
