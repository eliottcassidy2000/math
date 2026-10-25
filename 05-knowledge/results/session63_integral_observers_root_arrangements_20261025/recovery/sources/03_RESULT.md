# Palette orbits versus weighted response rank

**Status: PROVED** elementary statements below; **FINITE-EXACT** bounded controls.
This is a continuation/application of existing repository interfaces, not a
historical-priority claim. It proves no new Hadwiger–Nelson bound, bounded-width
exterior theorem, or geometric equality gadget.

## 0. Inheritance and the new contribution

Read against repository HEAD `fde0dc637d81` on 2026-10-25:

* `ports_regular_actions_hn_20261025/hn/README.md`, §§2–5: the minimal
  `q+2` Boolean path automaton, anchored C6 closure, and the full named
  boundary-count vector.
* The same directory's `RELATIONAL_LIFT.md`, §§1–3: one common exterior
  colouring, exact joint separator tables, and the `q^w(q+2)` row bound.
* `boundary_palette_compression_20261011.md`, §1: complete equality-partition
  relations and their Bell-number bound; **this quotient is inherited**.
* `spindle_tournament_interfaces_20261011g.md`, §§3–4: the induced unit-distance
  three-port higher-correlation hostile and complete cross-interface
  composition; **the partition-amalgamation rule is inherited**.
* `boundary_context_weights_20261011d.md`, §2: conditional weights and the
  multiplicity of a forgotten new colour.
* `session61_hadwiger_six_seven_20261025/palette-deck/RESULT.md`, §3: full
  ordered `F2^3` pair frames versus their parity character.

All paths above are relative to `05-knowledge/results/`. The new results are
(1) precise orbit row budgets for the relational path/C6 machine, including
partially named palettes; (2) an exact three-way distinction between processed
path rank `q`, whole-word rank `q+1`, and two-port response rank `q^2`.

Anchor: compress the actual correlated exterior interface. Niche: audit its
weighted response algebra. Wildcard: compare the lost gauge with full-frame
holonomy. The concept board is **joint relation / equality partition / named
pin / Boolean future / linear response / frame holonomy**. The inherited
meta-pattern used is “Pay algebraic rank, then audit the physical action.”
The closest proved mechanism is relational join-and-project; the canonical
hostiles are correlated port lists and the induced three-port graph below;
the corrected near miss is compression of independently renamed pieces;
the least-used sidecars are named pins and the empty-prefix flag.

## 1. Exact orbit budgets for the relational local machine

Fix a palette C of q>=2 colours and k distinct named colours, 0<=k<=q.
All constraints and observations must be invariant under simultaneous
permutations of the other q-k colours. If arbitrary named lists are allowed,
every colour mentioned by name must be included among these k constants.
Lists formed as complements of the colours of retained neighbouring vertices
are equivariant without naming any colours.

Let B_{<=r}(n)=sum_{j=0}^r S(n,j), including B_{<=r}(0)=1 for r>=0.
Define

    N_(q,k)(w) = sum_(j=0)^w binom(w,j) k^(w-j) B_{<=q-k}(j),       (1)

with 0^0=1. This is the number of orbits of w ordered colour variables
under the group fixing the k named colours pointwise. An orbit specifies
which variables take each named colour, and an equality partition of the
remaining variables into at most q-k unnamed blocks. To count them, choose
those j remaining variables, assign the others to the k names, and partition
those j variables. This proves (1), including k=0 and k=q.

**Theorem 1 (relational automaton orbit budgets).** For a separator containing
w exterior colour variables, the complete joint Boolean relation with one
active path state has at most

    2 N_(q,k)(w) + N_(q,k)(w+1)                              (2)

orbit rows. With one retained cycle anchor it has at most

    2 N_(q,k)(w+1) + N_(q,k)(w+2).                          (3)

**Proof.** The dead and free states are fixed by palette permutations, giving
two copies of the w-variable orbit space. A forced(c) state contributes one
additional colour coordinate, giving the (w+1)-variable orbit space. A cycle
anchor adds a colour coordinate to both cases. These are ambient row budgets;
reachability may reduce them. No assumption that all orbit rows are reachable
is needed. QED.

For k=0 these are `2 B_{<=q}(w)+B_{<=q}(w+1)` and
`2 B_{<=q}(w+1)+B_{<=q}(w+2)`. At q=6,w=6 they are 1,282 and 5,863,
respectively, instead of the named-table budgets 373,248 and 2,239,488.
For k=q, N_(q,q)(w)=q^w, so (2) and (3) become exactly the inherited named
bounds. Naming all colours eliminates this symmetry saving.

**Exact operations, not independent representatives.** Store the *complete
set of joint orbits*. For a join of invariant relations on variable sets U,V,
retain every admissible partition on U union V whose restrictions occur in
both input tables. Marked constants stay distinct and fixed. The number of
unnamed blocks is capped at q-k. This is the inherited complete
cross-interface composition rule, now applied also to state colours and cycle
anchors. Projection restricts the partition and recanonicalizes it, retaining
all projected rows. The local DFA and its closing-edge filter are equivariant,
so they descend to these joint orbits. Induction over joins and projections
therefore gives exactly the named-table empty/nonempty answer.

For sufficiency at a join, choose one actual assignment of the union partition.
Each input witness can be renamed, fixing the k marked colours, to agree with
that assignment on its boundary. These witnesses agree on their common
variables. Necessity is restriction. This is why one representative per
*joint orbit* is enough, while independently fixing a representative for each
piece is not. For disjoint single variables x,y at q>=2, each unary universal
relation has one orbit. The join has two orbits, x=y and x!=y. Setting each
piece's canonical colour to zero and keeping only that alignment loses the
second orbit and falsely rejects after adding x!=y.

These budgets do not bound the cost of constructing the exterior relation,
its separator width, or the number of simultaneous active components. They
count rows, not the number of distinct tables (which can be exponential in
the row budget). They also do not yield a three-state automaton on arbitrary *named* list letters:
relative alignment with the next named letter must still be retained.

### A geometric hostile to all pairwise marginal tables

The inherited exact Eisenstein points are

    (z0,z1,z2,p0,p1,p2) = (0,1,alpha,-alpha,2,2alpha),
    alpha=(1+i sqrt(3))/2.

Their induced unit graph has the triangle z0z1z2 and only the three spokes
zi--pi. At q=3 its boundary relation on p0,p1,p2 is exactly “not all equal”.
It has 24 labelled boundary assignments, each with two extensions. Every
pairwise projection is the full nine-row relation, just as for three isolated
vertices. The complete partition table has four of the five possible
partitions: it excludes 000. Thus full joint orbit compression is exact while
all pairwise marginal relations are not. The checker independently verifies
the six exact unit edges with the norm a^2+ab+b^2 and enumerates all 729
colour assignments. No new geometry is claimed here.

## 2. The path transition algebra is the full matrix algebra

Work over Q, with q>=2. Put A=J-I, d_L the indicator column of L, and

    M_L = diag(d_L) A.

For a nonempty processed prefix, v_c counts its colourings ending in c.
Appending L gives v'=M_L v, the inherited recurrence. For singleton {i},
write M_i=e_i r_i, where r_i is row i of A.

**Theorem 2 (full response algebra).** The Q-linear span of products of
singleton-list matrices of lengths one and two is Mat_q(Q).

**Proof.** For i!=j,

    M_i M_j = e_i (r_i e_j) r_j = e_i r_j.

For i=j the same outer product e_i r_i is M_i itself (not M_i^2, which is
zero). A has eigenvalues q-1 on the all-one line and -1 on the sum-zero
space, so its rows r_j form a basis. The q^2 outer products e_i r_j therefore
form a basis of Mat_q(Q), and every one is a legal length-one or length-two
product. QED.

This is a linear-span statement, not a claim that every rational matrix is
realized by a list word. Linear combinations can have negative coefficients;
actual transition products have nonnegative integer entries. The coefficient
field matters: at q=3 over F2, A has the all-one vector in its kernel. Every
nonempty product then annihilates that vector, so their span is not Mat_3(F2).
The Q hypothesis is not cosmetic.

### Three different minimal linear dimensions

A “dimension” below means the dimension of a vector space over Q supporting
linear transitions and linear output queries. It is not a lower bound against
arbitrary nonlinear encodings, and not the number of Boolean DFA states.

**(a) Nonempty processed-prefix continuation: exactly q.** The vector v is
sufficient. Prefixes consisting of the singleton lists {i} give v=e_i.
Appending singleton {j} observes r_j v. The q-by-q response matrix is A,
which is invertible, so no smaller linear representation supports every
such continuation. The empty continuation has output 1^T v and needs no
extra coordinate once at least one vertex has been processed.

**(b) Whole path-counting series, including the empty path: exactly q+1.**
Let F(w) be the number of proper list colourings of the word w, with
F(empty)=1. A dimension-(q+1) representation has state (s,v), starts at
(1,0), outputs s+1^T v, and updates by

    (s,v) -> (0, s d_L + M_L v).                             (4)

For the lower bound, use prefixes and suffixes indexed by the empty word
and by the q singleton lists. Their concatenation-count matrix is

    H = [ 1    1^T ]
        [ 1     A  ].                                      (5)

A linear representation of dimension d factors every such matrix through
Q^d, so its rank is at most d. Since A^(-1)1=1/(q-1)1,

    det H = det A [1-q/(q-1)] = (-1)^q != 0.

Thus d>=q+1, proving sharpness. This catches the tempting false claim that
the full counting series has rank q. Indeed starting the q-vector recurrence
at v0=1/(q-1)1 gives all nonempty counts correctly but outputs q/(q-1), not 1,
on the empty word. The empty-prefix flag is essential.

**(c) Two-port response: exactly q^2.** For a list word w on newly added
vertices define T_w=M_(L_m)...M_(L_1), with T_empty=I. Its entry
(T_w)_(b,a) counts extensions with old endpoint a and new endpoint b.
Store T_w, update by left multiplication, and query its entries by named
endpoint pins. By Theorem 2, attainable matrices span Mat_q. The entry
queries span its dual. The response matrix between the q^2 basis words
and the q^2 endpoint pins is therefore invertible, so any linear interface
supporting these queries has dimension at least q^2.

This lower bound includes arbitrary singleton list letters and named endpoint
queries. It does not apply unchanged to a restricted symmetric context bank.
Closing a genuine cycle of length m>=3 with lists gives tr(T_w); carrying
T_w handles every anchor at once. The trace convention at the empty word is
q, not the empty graph's one colouring; it must not be silently substituted
for F in (b).

The three ranks and the Boolean q+2 states answer different questions.
There is also a precise failure of weighting the inherited deterministic
Boolean machine: a **single homogeneous linear scalar register carried along
its deterministic state** is insufficient. At q=3 the one-letter prefixes
L={0,1} and L={1,2} both end in free. Their responses to the two singleton
suffixes {0},{2} are respectively (1,2) and (2,1). In a one-register linear
model, every fixed suffix acts from the same free state by a fixed scalar
multiplier, so the response vectors from these two prefixes would be
proportional. Their determinant is -3, a contradiction over Q.

This excludes only that deterministic-state, one-linear-register model.
It does **not** exclude an arbitrary q+2-dimensional weighted automaton:
the q+1-dimensional construction (4) already fits in one. Nor does it exclude
nonlinear scalar encodings. The inherited proportional-vector multiplicity
example distinguishes tag-only counting, but alone would not refute a scalar
sidecar; the nonproportional example above supplies the missing test.

## 3. Relative frames: what a scalar character forgets

This is a scoped connection, not an identification of palette colours with
surface sheets. If pieces carry full bijective palette frames f_v, and an
oriented seam prescribes permutation g_uv by f_v=g_uv f_u, then full frames
exist iff the ordered product of seam permutations around every cycle is the
identity. A spanning tree proves sufficiency and necessity is telescoping.
Independent local renaming changes g_uv by h_v g_uv h_u^(-1); cycle gains
are conjugated, so being the identity is invariant.

Applying only the sign character of S_q is insufficient for q>=3. A triangle
with one seam gain (012) and two identity gains has trivial sign around every
cycle but nonidentity full holonomy. This obstructs a *full frame*. It does
not by itself obstruct a single colour section, which only needs a common
fixed colour for the holonomy subgroup. For example at q=4 the same 3-cycle
fixes colour 3, so a one-colour section exists despite the frame obstruction.

Session61's ordered three-pair model is the abelian specialization A=F2^3.
A triangle with gain (1,1,0) on each edge has zero parity shadow but nonzero
full cycle gain. The lost space there has dimension 2*b1 because pair labels
are fixed. Allowing pair permutations changes the coefficient system; no
unchanged 2*b1 claim is made here.

The precise bridge is **relative alignment**: partition amalgamation keeps
all alignments compatible with the shared boundary; explicit full-frame
coordinates keep their permutation gains. A scalar sign or parity quotient
cannot replace either full joint boundary relations or full seam gains.
There is no inferred Euclidean realization or chromatic obstruction.

## 4. Verification and boundaries

Run `python3 check.py` and `python3 -O check.py` here. Always-active checks
cover orbit counts, marked-palette endpoints, orbit joins versus literal
assignments, equivariance of list transitions/anchor closure, the geometric
higher-order hostile, exact integer path counts, determinants in (5), full
matrix span, endpoint-query rank, the characteristic-two scope hostile, and the two frame hostiles. Universes and
actual counters are printed in the saved JSON outputs. Targeted mutants
remove an amalgamation or erase the empty-prefix flag; the checker must fail.

Only the bounded controls are FINITE-EXACT; the displayed formulas have the
proofs above. Complete partition composition and its known geometric hostile
are explicitly inherited. Novelty beyond this repository is not asserted.
The next useful search target remains a realizable exterior with an empty
*joint* acceptance relation; the new orbit budgets can reduce an exact search
without relaxing that target or exchanging its quantifiers.

`check_independent.py` does not import the producer: it compares permutation
Burnside counts with onto-block inclusion-exclusion, recomputes the published
orbit budgets, evaluates the Hankel determinants by the Leibniz formula, and
checks the explicit rational inverse of J-I. Normal and optimized outputs
are retained separately. `NEGATIVE_CONTROLS.json` records nonzero exits for
both targeted mutants in both interpreter modes.
