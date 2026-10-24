# Guarded packet atoms: a conditional all-size decomposition rule

24 October 2026. **PROVED:** Theorems 1–5 below, conditional on their explicit
native-input hypotheses. **FINITE-EXACT:** the accompanying exhaustive small
models, full twenty-option native reservoir checks, and inherited six native
primitive circuits. No unconditional native supply theorem, all-order constant,
or rewiring of the complete canonical atlas is asserted.

## 1. Keep three objects separate

1. The actual native record contains a source matrix, an indexed Hamiltonian
   path, its actual owner map, all declared shared pins, the chosen full kernel
   lift, the literal complement partner and reflected converse.
2. The capacity carrier is the normalized black resource, not a raw source
   name, a charge label, or an alternative lift of the same resource. A selected
   resource carries its entire complement/converse closure. The immutable
   remainder has already been subtracted from the available inventory.
3. The charge label belongs to the authenticated, owner-marked paired-charge
   group. Equality of group orders or invariant factors does not identify owner
   labels, construct native records, or transport resources between graphs.

Throughout Sections 2–5 the supply consists of N binary slots, at distinct
source owners s_i outside a common target pair U,V. Each slot has two actual
records, targeting U and V respectively. All required local predicates hold
for each record. The complete pin conditions leave exactly those two choices
if completeness of the entire pinned family is claimed. More weakly, the
connectivity statements hold inside a declared authenticated two-choice family
without claiming it exhausts the pin fiber.

All choices are jointly compatible in the actual shared coordinates. Their
2N normalized resources are pairwise distinct, their full closed port sets
are disjoint, and none meets the immutable remainder. Every retained shared
owner/kernel condition is either already fixed identically across the choices
or explicitly included in the remaining acceptance predicate. No common-host
tournament amalgamation is inferred from agreement of declared observer pins.

Fix the original finite connected weighted graph L, its critical group K, its
2-primary part K_(2), the actual paired multiplicity 2^b, and a marked map

    Q = 2^b K_(2),       ell(v) = 2^b pi_2([e_v-e_root]).

Assume ell(U)=ell(V). Let g_i=ell(s_i). The bit x_i=1 selects the V record,
and x_i=0 selects the U record. The unchanging ordinary per-slot guards must
be separately authenticated; equality in Q alone need not preserve them.

The exact remaining constraints are

    sum_i x_i = r,             sum_i x_i g_i = gamma in Q.       (1)

These are target margins and the full original nodewise dyadic current
signature, not the assertion that the selected packet is harmonic in isolation.
They are sufficient only when there are no omitted varying constraints.

**Why (1) is the exact current test here.** For a balanced difference d=y-x,
put u=sum_i d_i e_(s_i) and v=e_U-e_V. Up to the harmless global orientation
sign, the paired perturbation is 2^b(uv^T+vu^T). Every source row has charge
d_i ell(v)=0. The two target rows have opposite charges ell(u), so all
original dyadic harmonic banks are preserved exactly when ell(u)=0. This
uses the authenticated row-image theorem, including its converse. Actual
sources and targets are disjoint here. Section 6 handles overlapping roles.

## 2. Conformal decomposition inside the endpoint union

For accepted x,y, let P={i:x_i=0,y_i=1} and M={i:x_i=1,y_i=0}.
Then |P|=|M| and sum_P g_i=sum_M g_i. A **balanced subtrade** is a nonempty
pair (P',M') with P' contained in P, M' contained in M, equal cardinalities,
and equal charge sums. It is **primitive** if it has no proper nonempty
balanced subtrade. Indices remain actual slot identities even when labels
coincide; a label is not a resource.

**Theorem 1 (conformal native decomposition).** Under Section 1's hypotheses,
every two accepted assignments admit a partition of their changed slots into
primitive balanced subtrades. Switching these subtrades in any order gives a
legal route from x to y. Every changed slot switches exactly once. Every
intermediate assignment uses only endpoint resources, retains the selected
actual owner/kernel records and complement/converse closure, and has no
simultaneous normalized-resource or closed-port reuse.

**Proof.** If P=M=empty, use the empty route. Otherwise choose an inclusion-
minimal nonempty balanced subtrade; the full pair (P,M) ensures one exists.
Remove it. Equal cardinality and equal charge survive on the remainder, so
induction terminates after at most |P| subtrades. Their supports are disjoint
and their sum is y-x. Any partial sum preserves both equations (1), and its
bits lie coordinatewise between x and y. The fixed actual records therefore
give exactly the required native intermediate. The global resource-disjointness
hypothesis, not charge cancellation, proves injection. Each record includes
its original lift and equivariant closure, so no hidden lift is re-counted
as capacity. The immutable remainder stays unchanged. QED.

This theorem permits atomic simultaneous subtrades; it does not assert that
their individual slot changes are legal. The number of changed slots along
the route is exactly the endpoint Hamming distance, although the largest
simultaneous event can be large. It is not a theorem about an unknown occupied
remainder or the already complete canonical atlas.

**Constructive extraction.** For each side construct a finite subset-sum
table indexed by (cardinality, group sum), retaining an actual index witness.
Find the smallest positive cardinality appearing with the same sum in both
tables. The resulting subtrade is primitive: any proper balanced subtrade
would have smaller cardinality. Remove it and repeat. Witness-retaining
dynamic programming is polynomial in the explicit number of group elements,
not necessarily in the bit length of a compressed group description. Native
record authentication and any table-generation cost are charged separately.

## 3. A group-dependent upper bound on atomic support

Let D(G) be the least positive integer such that every sequence of D(G)
elements of a finite abelian group G has a nonempty zero-sum subsequence.
Equivalently it is the maximum length of a minimal zero-sum sequence. The
equivalence follows by removing one term from a minimal zero-sum sequence,
and, conversely, appending the negative sum to a longest zero-sum-free one.

**Theorem 2 (Davenport bound).** A primitive balanced trade with k slots on
each side satisfies k<=D(G). Consequently the route in Theorem 1 can use
events of at most 2D(G) changed slots, with no external resources.

**Proof.** Pair the positive and negative indices by any bijection, and set
d_j=g_(p_j)-g_(m_j), j=1,...,k. Their total is zero. A proper nonempty zero-sum
subsequence would give a proper balanced subtrade on the paired indices,
contradicting primitivity. Thus this is a minimal zero-sum sequence, of length
at most D(G). The result applies to each atom in Theorem 1. QED.

The argument allows repeated charge labels but never repeated resources. A
single pairing being minimal zero-sum is **not sufficient** to establish
balanced primitivity: positive labels (0,1), negative labels (1,0) in C3 give
paired differences (2,1), a minimal zero-sum sequence, yet cross-pairing gives
two one-against-one subtrades. Actual primitivity checks must allow all
cross-subsets, as the checker does.

For every finite group, D(G)<=|G|: the |G|+1 partial sums of any |G|-term
sequence repeat, producing a nonempty consecutive zero sum. For a finite
abelian p-group the sharper exact answer is the classical Olson formula.

**Lemma 2.1 (self-contained p-group computation).** If

    G = C_(q_1) + ... + C_(q_t),    q_i=p^(a_i),

then D(G)=1+sum_i(q_i-1).

**Proof.** Write d=1+sum_i(q_i-1). The group algebra over F_p is

    F_p[G] = F_p[X_1,...,X_t]/(X_i^(q_i)-1)
           = F_p[T_1,...,T_t]/(T_i^(q_i)),  T_i=X_i-1.

The second equality uses the characteristic-p identity
(1+T_i)^(q_i)=1+T_i^(q_i). Its augmentation ideal I=(T_1,...,T_t) has I^d=0:
every monomial of degree d exceeds at least one exponent bound q_i-1.
For any g_1,...,g_d in G, every 1-X^(g_j) belongs to I, so their product is
zero. If the sequence had no nonempty zero-sum subsequence, the coefficient
of X^0 in this product would be 1, supplied only by the empty subset, a
contradiction. Hence D(G)<=d. Conversely take q_i-1 copies of basis element
e_i for each i. No nonempty subsequence sums to zero, so D(G)>=d. QED.

Appending sum_i e_i to that lower-bound sequence gives a minimal zero-sum
sequence of length d. This proves sharpness for abstract zero-sum sequences,
not existence of a native primitive trade of support 2d.

More generally the bound 2D(G) is sharp in the **unrestricted abstract labelled-
slot model**: take a minimal zero-sum sequence of length D(G) on the positive
side and D(G) zero labels on distinct negative slots. Any proper balanced
subtrade would be a proper nonempty zero-sum subsequence on the positive side.
Assign formally distinct resources to all slots and choices. This construction
does not authenticate native source owners, paths or resources; repeated group
labels need not have the required native multiplicity. It cannot establish
native sharpness.

**Original order-eight corollary.** The separately authenticated marked
group is Q=(C2)^4+C128, so

    D(Q)=1+4(2-1)+(128-1)=132,       2D(Q)=264.              (2)

Thus 264 is a valid conditional upper bound for every primitive balanced
common-target binary packet of this precise type in the original graph.
It is NOT shown sharp for native packets, does not improve the known exact
maximum 16 in the fixed twenty-option reservoir, and does not bound arbitrary
multi-target or shared-coordinate systems. If the original graph/order changes,
Q changes: for Q_b=direct-sum C_(2^max(a_i-b,0)), the corresponding bound is
2[1+sum_i(2^max(a_i-b,0)-1)], omitting trivial factors. No order-independent
constant follows. A matching numerical value in another construction is
not an owner-marked map between the constructions.

## 4. What the sixteen-source obstruction actually forbids

The six inherited complete pinned circuit families have exactly two accepted
assignments at full dyadic precision. In particular the sixteen-slot family
is one primitive atom, not a concatenation of legal smaller atoms on those
same slots. Theorem 1 returns that whole atom. It cannot contradict its
minimality: a decomposition theorem is not a theorem that all atoms are small.

Let W be an enlarged legal state space and pi its map to original literal
slot choices. If every legal W-state projects to an original accepted state,
and each event changes at most b projected slots, every legal W-route projects
to a walk in the b-Hamming graph of the original accepted relation. For the
two-state sixteen-slot family this requires b>=16. This elementary projection
argument is independent of the number of fibers, neutral resource layers,
or retained kernel decorations. It does not apply without an authenticated
projection retaining the original native predicates.

## 5. Charged buffers can escape that projection, conditionally

**Theorem 3 (partial-sum buffer).** Let a primitive trade have positive slots
p_1,...,p_k and negative slots m_1,...,m_k. Put

    d_j=g_(p_j)-g_(m_j),  s_0=0,  s_j=d_1+...+d_j.

For k>=2 suppose there are k additional binary slots b_0,...,b_(k-1), outside
the original slots, with actual labels h-s_j, j=0,...,k-1. Assume the ENTIRE
enlarged supply satisfies Section 1's native, pin, owner/lift, guard, common-
target and normalized-resource-disjointness hypotheses. Initially and finally
select b_0 and no other buffer in the bit-1 target class.

Then the original trade admits a k-step enlarged legal route with at most
four changed slots per event. All original slots change once; the buffer
assignment and its complete native resource/current ledger return to their
original values at the end. Every state is injective in the enlarged resource
carrier. For k=1 no buffer is needed and the event has two changed slots.

**Proof.** Theorem 2's minimal-zero-sum argument shows s_0,...,s_(k-1) are
distinct: a repetition yields a proper consecutive zero-sum subsequence.
At step j switch m_j from 1 to 0 and p_j from 0 to 1, and move the selected
buffer from b_(j-1) to b_j, with b_k=b_0. The cardinality changes cancel
separately in the original pair and the buffer pair. The charge change is
d_j+(h-s_j)-(h-s_(j-1))=0. Equation (1), all fixed native predicates, and
resource injection therefore hold at every step. Since s_k=0, the buffer
returns to b_0. QED.

This is a supply criterion, NOT a claim that the necessary native buffer
labels, paths or free resources exist. A buffer label missing from the actual
owner image is an immediate obstruction; an existing label may still have
no eligible source/path with the required pins or available normalized tile.
Even distinct labels do not prove native supply. The joint shared-coordinate
conditions must hold for the enlarged family, not only each buffer separately.

Deleting the buffer from an intermediate enlarged state leaves debt s_j,
which is nonzero for 0<j<k. Thus projection no longer lands in the original
accepted relation: this is exactly why Theorem 3 does not contradict Section 4.
No new port is manufactured from a lift. Restoring a borrowed resource later
is not simultaneous reuse and does not create a second capacity unit.

## 6. The all-size incidence/lattice version

For a general native transition column t, retain its ACTUAL integer symmetric
zero-row-sum Laplacian perturbation Delta L_t at the actual physical owner
indices. If a common multiplicity 2^b divides every entry of every column,
write Delta L_t=2^b F_t, with F_t the DIVIDED integer perturbation. Using the
already defined ell_b(w)=2^b pi_2([e_w-e_root]) in Q=2^b K_(2), its guard
column is the vector indexed by all physical rows v,

    c_t(v) = sum_w (F_t)_(v,w) ell_b(w) in Q.            (3)

Equivalently, use the actual, undivided Delta L_t against the undoubled marked
map ell_0 in K_(2). These are the same charge because Delta L_t=2^b F_t.
DO NOT use the actual Delta L_t against ell_b: that would apply the paired
multiplicity twice and can silently erase a genuine guard. If there is no
common positive b, use b=0 and the undivided convention throughout.

Margins, shared-coordinate constraints and resource capacities remain
separate obligations. When source and target role sets overlap, add all
contributions in each same physical row before applying (3). Neither a scalar
total nor separate source/target copies is an exact replacement. The fixed
marked Q and the single, explicit multiplicity normalization are part of the
input. This clarification was made during final review before first publication;
the common-target formulas and all direct native image certificates already
use the correct single factor.

Let an authenticated finite encoding have integer variables z, integral
equalities Az=a, modular equations Cz=c mod m (rowwise positive moduli), and
coordinate bounds l<=z<=u. Suppose decoding its integer solutions gives
exactly the accepted native assignments with their chosen actual records,
and every such assignment has a retained encoding. No omitted predicate is
allowed. Define the lattice

    Lambda={d in Z^n: Ad=0, Cd=0 mod m}.

Write e conformal to d if every coordinate has the same sign and
|e_i|<=|d_i|. Its **Graver atoms** are the nonzero conformally minimal members.

**Theorem 4 (finite exact-fiber decomposition).** Lambda has finitely many
Graver atoms. Any two encoded solutions z,w can be joined by a conformal sum
of these atoms, with every intermediate coordinate between the corresponding
endpoint coordinates. Every intermediate is an exact accepted encoded state.
The atom set depends on the full incidence encoding, not just |Q| or D(Q).

**Proof.** In each of finitely many orthants, absolute values identify the
conformal order with coordinatewise order on N^n. Dickson's lemma says every
subset has finitely many minimal elements. For completeness, any infinite
sequence in N has an infinite nondecreasing subsequence: either some value
repeats infinitely or unbounded tails allow successive larger choices.
Applying this extraction to each coordinate gives comparable pairs in any
infinite sequence in N^n; an infinite antichain of minimal elements is
impossible. Thus the atom set is finite. For any nonzero d in Lambda choose
an atom e conformal to d. Then d-e is in Lambda, remains in its orthant,
and has smaller l1 norm. Induction decomposes d. Any partial sum remains
coordinatewise between 0 and d, proving the box and congruence claims. QED.

**An explicit sufficient lossless encoding.** Use one-hot variables for each
slot's authenticated actual records, with exactly one chosen record per slot.
For every normalized resource, add its occupancy inequality at the appropriate
unit capacity, after subtracting the immutable remainder; turn it into an
equality with a nonnegative bounded slack. Include all actual target margins
and the complete row charges (3). For varying shared owners, use one-hot
variables for each actual owner value, not just residues, and equality of
each occurrence's chosen value indicators. General finite residual factors
can be encoded by one-hot variables for every exact allowed factor tuple,
with occurrence consistency against slot and owner indicators. Integral
one-hot consistency is exact even on cyclic factor graphs; this is NOT
an assertion that the fractional local polytope is exact. The tables must
retain the full kernel lifts and may be exponentially large.

Applying Theorem 4 to this full encoding preserves injection and actual
shared coordinates. Slot one-hot coordinates stay within the endpoint union,
so each physical slot is unchanged or switches once directly from its old
record to its new record. Auxiliary factor/slack coordinates are retained
as well. The result can be as coarse as a single event involving every slot.
It neither promises a compact solver nor a uniform small local generating set.

**Lattice generation is not bounded-fiber connectivity.** The full Graver
atom set in Theorem 4 cannot be replaced merely by a list generating Lambda
as an abelian group. Integer expressions in those generators can cancel and
leave the permitted box. For example Lambda=Z(1,-1) is generated by
2(1,-1) and 3(1,-1), since their difference is (1,-1). In the binary fiber
z_1+z_2=1 neither listed move nor its negative is feasible from either
endpoint, although the fiber has two states and the lattice is fully generated.
The missing Graver move (1,-1) connects them. A native lattice-span certificate
and a capacity-bounded Markov/connectivity certificate are therefore distinct
deliverables. This distinction also applies when small native circuits span
the complete integer current kernel.

**Theorem 5 (why a charge-only bound fails beyond common targets).** There
are explicit finite binary assignment interfaces with trivial Q, pairwise
disjoint resources and constant per-slot guards whose sole two accepted
states differ in arbitrarily many slots.

**Proof.** For n>=3, slot i has target i in choice 0 and target i+1 mod n in
choice 1. Require exactly one selected record at every target. Target j then
has occupancy (1-x_j)+x_(j-1), so its margin equation gives x_j=x_(j-1).
The only states are all 0 and all 1. Assign a separate resource to every
slot/choice and use trivial charge. Thus the necessary atom support is n,
although D(Q)=1. QED.

This is an abstract exact-interface hostile, not an asserted realization by
actual tournament ports. It rigorously disproves a general incidence theorem
based only on Q; a native universal claim would still require native examples.

## 7. Obstruction taxonomy and what a repair must supply

- **Supply:** a group vector or divisibility solution does not create paths,
  pin completions or free normalized resources. Authenticate the actual menu.
- **Capacity:** raw lifts, relabelings and repeated source names can project to
  one tile. Subtract the remainder and check normalized injection, not counts.
- **Shared-coordinate:** two kernel lifts can agree modulo a quotient and
  disagree as actual shared owners. Keep the full values and joint factors.
- **Electrical:** margins or one scalar phase can hold while an actual row
  has nonzero full Q debt. Keep (3), including every dyadic coordinate.
- **Atomicity:** an accepted endpoint difference can be one primitive trade.
  Neutral vacancies do not split its current barrier; charged buffers must
  carry the exact intermediate debt with their own native certificates.
- **Bounded connectivity:** a small menu may generate the whole integer
  lattice without connecting a fixed capacity-one fiber. Conformal atoms or
  a separate legal-state connectivity certificate are needed, not only rank,
  index or Smith-normal-form agreement.
- **Incidence:** a common-target simplification is invalid after changing
  target pairs or identifying source and target owners. Rebuild physical rows.
- **Projection:** equal cardinalities, group types or numerical invariants
  are not maps. A usable comparison preserves actual owner coordinates,
  current predicates and resources and explains its fibers and capacity.
- **Precision:** genuine higher dyadic harmonics scale beyond saturation;
  an unchanged integer potential tested at a new modulus is not that lift.
- **Scope:** partial packet replacement is not changing the deterministic
  complement involution of the fully occupied canonical atlas.

The all-size wiring rule established here is therefore conditional and exact:
authenticate a lossless native finite interface, retain all incidence and
capacity constraints, and decompose in its exact fiber. Common equal-charge
targets reduce that fiber to balanced finite-group trades and provide (2).
Anything stronger needs a new supply theorem or additional structure in the
full incidence system, rather than a numerical analogy.

## 8. Verification and references

The accompanying checker has no optimization-sensitive assertions. It tests
all endpoint pairs in every cardinality/charge fiber for the stated small
groups and label lengths, checks every order of the extracted primitive atoms,
and replays the inherited actual native checker. It verifies all forty
normalized resources of the twenty-slot supply are distinct, rechecks the
six exact primitive families, constructs a native decomposable endpoint
pair, and exercises loss-of-margin, loss-of-charge, lift/quotient, capacity,
pairing-primitivity, neutral-buffer and multi-target-cycle hostiles. The
partial-sum buffer controls are explicitly abstract finite models, not new
native buffer certificates. Finite tests support the implementations; the
all-size claims rely on the proofs above.

See REFERENCES.md for the original Olson reference and an author-hosted
primary Graver-basis paper. SOURCE_PINS.json identifies the exact inherited
inputs. Normal and optimized receipts must agree byte-for-byte.
