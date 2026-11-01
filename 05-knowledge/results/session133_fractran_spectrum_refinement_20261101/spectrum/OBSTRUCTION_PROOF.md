# Sharp local witnesses for bounded tournament Hamiltonian counts

**Status: PROVED structural reduction from cited classical results;
FINITE-EXACT complete forbidden bases in the stated bounded range.
No literature-priority claim. Ordinary spectrum completeness remains OPEN.**

## 1. Object, inherited tools, and the precise new direction

For a finite tournament T, let H(T) count ordered directed Hamiltonian vertex
words, with H(empty)=1. For an integer B>=1, define the hereditary class

    C_B = { T : H(T)<=B }.

An induced-minimal obstruction is a tournament F outside C_B all of whose
proper induced subtournaments belong to C_B. These are isomorphism classes,
not labelled tournaments, cycle-conflict graphs, or automorphism-normalized
path counts.

The inherited machinery is:

1. H is induced-subtournament monotone. This follows from the Rédei insertion
   injection, or from OCF: deleting vertices only removes available literal
   odd-cycle packings. Consequently minimality is equivalent to checking
   **every one-vertex deletion**.
2. Strong components have a unique linear dominance order and their H values
   multiply. Singleton components contribute one.
3. Every strong tournament of order at least four has an induced strong
   subtournament on one fewer vertex. This is the case k=n-1 of Moon's
   [*On Subtournaments of a Tournament*, Theorem 2 (1966)](https://doi.org/10.4153/CMB-1966-038-7).
4. Busch's strong bound is H(C)^3>=5^(|C|-1). For |C|>=4,
   `5^(|C|-1)>3^|C|`: start with 125>81 and multiply the two sides by 5
   and 3. For a strong triangle, H^3=27=3^3. Multiplying over nonsingleton
   SCCs proves that a tournament without singleton SCCs satisfies

       H(T)^3 >= 3^|T|.                                    (1)

   [Busch, EJC 13 (2006), N3](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v13i1n3/pdf).

Items 1--4 and the deletion-only kernel bound implicit in (1) already occur
in the inherited `rotation_owners_spectrum_20261020e/spectrum.md`. What is
proved below is a **sharp induced-obstruction-order cutoff**, followed by an
explicit finite forbidden-family characterization. It is not a new proof of
spectrum completeness or a claim that the classical ingredients are new.

## 2. Sharp obstruction-size theorem at every threshold

**Theorem.** Let B>=1 be an integer, and let

    r = min { j>=1 : 3^j > B }.

Every induced-minimal obstruction to C_B has at most **3r vertices**. This
bound is attained for every B: the ordinal sum of r directed triangles is
an induced-minimal obstruction on exactly 3r vertices.

**Proof.** Let F be minimal, of order n. F has no singleton SCC: deleting
such a vertex would leave H unchanged and contradict minimality.

If every SCC is a triangle, let their number be s. Then H(F)=3^s>B, while
deleting any vertex turns its triangle into a transitive two-vertex block
and leaves H=3^(s-1)<=B. Hence s=r and n=3r.

Otherwise some SCC C has order at least four. By Moon choose v in C with
C-v strong. Every SCC of F-v is still nonsingleton: C-v has order at least
three and the other SCCs are unchanged. Minimality and (1) give

    3^(n-1) <= H(F-v)^3 <= B^3 < 3^(3r).

Thus n-1<3r, and integrality gives n<=3r. For sharpness, the ordinal sum of
r triangles has H=3^r>B and every one-vertex deletion has H=3^(r-1)<=B.
Monotonicity handles every smaller induced subtournament. QED.

**Consequences.** C_B has a finite forbidden induced-subtournament basis,
computable by enumerating all tournaments through order 3r and retaining
exactly the minimal bad ones. Every tournament with H>B contains an induced
witness on at most 3r vertices, and no smaller universal bound is possible.
This is a local recognition/certificate theorem for a fixed threshold; it
does not construct a tournament with a prescribed H.

For B=21, r=3 and the sharp bound is nine. Its extremal witness is
C3 ordinal-sum C3 ordinal-sum C3, with H=27 and deletion counts all nine.
For B=7, the sharp bound is six, witnessed by two ordered triangles.

### Exact-core sharpening and uniqueness at the top order

The parallel [exact core reduction](../core_reduction/THEORY.md), Section 3,
proves that the largest possible number of nonsingleton-SCC vertices among
tournaments with H<=B is

    K(B) = max_(b=0,1,2; 5^b<=B)
             [4b + 3 floor(log_3(B/5^b))].                 (3)

The logarithm means the largest integer a with 3^a*5^b<=B; the checker uses
only integer arithmetic. Replacing (1) by this exact bound in the same
strong-deletion argument sharpens the theorem:

    every minimal B-obstruction is either the ordered r-triangle chain,
    or has at most K(B)+1 vertices.                         (4)

Indeed F-v is singleton-free and H(F-v)<=B, so n-1<=K(B). The pure-triangle
case was already classified in the proof. This corollary deliberately keeps
the exceptional chain: deleting one of its vertices leaves two singleton
vertices, so its order is not its cyclic-core order.

Here K(21)=7 and K(7)=4. Thus the only possible order-nine minimal
21-obstruction is the three-triangle chain, and the only possible order-six
minimal 7-obstruction is the two-triangle chain. All other obstructions have
orders at most eight and five, respectively. The producer checks (4) for
every obstruction at all eleven thresholds. It also independently recovers
the exact-core floor controls at orders six and seven: minimum H9 with one
class/80 labelled tournaments, and minimum H15 with two classes/3360 labelled
tournaments, respectively.

## 3. Complete B=21 basis: 56 classes

**Theorem (CITED reduction + FINITE-EXACT enumeration).** A tournament has
H<=21 iff it contains none of the 56 tournaments listed in `basis21.json`
as an induced subtournament. The basis consists of:

| Kind | Order | Number of classes |
|---|---:|---:|
| Strong | 6 | 29 |
| Strong | 7 | 13 |
| Nonstrong | 8 | 13 |
| Nonstrong | 9 | 1 |

The nonstrong basis is explicitly:

* `C3 ordinal-sum S` and `S ordinal-sum C3` for each of the six strong
  five-vertex classes S (12 classes);
* `S4 ordinal-sum S4`, where S4 is the unique strong four-vertex class;
* the ordinal sum of three C3s.

### Strong completeness

If a strong minimal obstruction has n>=4, Moon gives v with F-v strong.
Since H(F-v)<=21 and `5^6>21^3`, Busch forces n-1<=6, hence n<=7.
The order-three tournament is not bad. The producer exhaustively generates
all classes through seven by one-vertex augmentation, giving
`1,1,2,4,12,56,456`. Degree-cell canonicalization is exact: an isomorphism
preserves outdegrees; all possible orderings within equal-degree cells are
checked. Orbit masses independently recover every labelled tournament at
orders one through seven, ending at 2^21.

Among these classes, exactly 29 strong order-six and 13 strong order-seven
classes satisfy H>21 and every one-vertex deletion H<=21. Actual adjacency
and all deletion values are frozen. This is finite exhaustive data with a
proved all-order cutoff, not inference from a sample.

### Nonstrong completeness, without a finite-order leap

A minimal bad tournament has no singleton SCC and at least two nonsingleton
components. Each component is a proper induced subtournament, hence H<=21.
Each has at least three vertices and H>=3.

There can be at most three components: if there were at least four,
deleting a whole component leaves a proper induced subtournament with at
least three nontrivial components and H>=27.

* With three components, if any component has order at least four, delete
  a vertex preserving that component's strong connectivity. The resulting
  three nontrivial components still give H>=27, violating minimality. So
  all three are C3, giving exactly one class.
* With two components, no component can have order at least six. A strong
  deletion would retain order at least five, whose H is at least nine;
  multiplied by the other factor, at least three, the deletion has H>=27.
  The order-five lower bound follows from Busch plus oddness, or the exact
  small atlas.
* If one component has order five and the other has order at least four,
  delete a vertex in the latter preserving strong connectivity. The product
  is again at least 9*3=27. Thus a five-vertex component must be paired with
  C3. All six such components have H in {9,11,13,15}; deleting within them
  gives an order-four tournament with H<=5, while deleting within C3 removes
  its factor three. Thus all twelve ordered joins are indeed minimal bad.
* Otherwise both component orders are three or four. The combinations
  (3,3), (3,4), (4,3) have H at most 15 and are not bad. The unique (4,4)
  join has H=25 and every deletion H<=15, so it is minimal bad.

These descriptions distinguish ordinary isomorphism classes because the SCC
word has a canonical order. Converse pairs are not silently merged.

### Exact basis plateau criterion, not a cardinality test

Write F_B for the set of isomorphism classes in the minimal obstruction basis.
For any integer thresholds 1<=a<b,

    F_a = F_b  iff  no tournament has H in (a,b].             (2)

Indeed a hereditary class is exactly the tournaments avoiding every member of
its minimal obstruction basis. Thus equality of bases implies equality of the
classes C_a=C_b, which is equivalent to emptiness of that H interval.
Conversely equal hereditary classes have the same intrinsically defined
minimal obstructions. This also proves the converse implication in (2).

The producer and harness compare the **entire adjacency-class lists**, not
just their lengths, and establish

    F_5 = F_7,             F_19 = F_21.

Consequently no H occurs in (5,7] or (19,21]. This independently recovers the
already-known exact omissions 7 and 21 (even counts are already excluded by
Rédei parity). Equal cardinalities b(a)=b(b) alone would not establish (2);
minimal forbidden patterns can be replaced without changing their number.
The structural deliverable here is the explicit local finite certificate,
not priority for the omissions and not coverage of any unbounded range.

## 4. A further precisely defined tournament sequence

Let b(B) be the number of ordinary isomorphism classes of induced-minimal
obstructions to H<=B. `basis_sequence.json` gives complete values:

| B | 1 | 3 | 5 | 7 | 9 | 11 | 13 | 15 | 17 | 19 | 21 |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| b(B) | 1 | 2 | 7 | 7 | 11 | 16 | 26 | 48 | 53 | 56 | 56 |

These are proved finite-range terms, not a conjectured recurrence. Even
thresholds have the same class as the preceding odd threshold because H is
odd. The all-threshold size theorem makes b computable, but no efficient
uniform algorithm or asymptotic law is claimed. The repeated entries at
5/7 and 19/21 record the known missing H values, not a numerological identity.

## 5. Independent verification and hostile controls

`check_basis_independent.py` does not import either producer. It checks each
listed obstruction and every one-vertex deletion by literal permutation
counting, including the order-eight and order-nine nonstrong witnesses.
It then independently covers all strong classes at orders six and seven by
fixing the directed Hamilton cycle 0->1->...->n-1->0 and enumerating every
remaining arc. Every strong tournament occurs in this cover by a Hamilton
cycle relabelling (Moon's vertex-pancyclic theorem, or Camion). This gives
512 order-six and 16,384 order-seven labelled presentations, not the
producer's augmentation universe. A backwards recursive subset count is
used there; the resulting minimal-bad class lists must match exactly.

The checker explicitly rejects an altered H value and source/sink padding
of a bad tournament. The latter is essential: padding preserves H but
**destroys induced-minimality**. Normal and `python -O` receipts must agree.
All code uses explicit exception checks rather than optimization-sensitive
assertions. No external census, nauty installation, or repository checkout
is needed to reproduce the finite certificates.
