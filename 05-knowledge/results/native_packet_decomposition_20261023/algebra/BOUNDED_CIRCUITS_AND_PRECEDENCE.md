# Balanced charge blocks and the cost of retained precedence

Status: PROVED, conditional physical decomposition theorem and exact abstract
obstruction. FINITE-EXACT checks are in `check.py` and `receipt.json`. No claim of
universal native supply or connectivity is made.

## 1. Inheritance and the actual interface

The closest proved mechanism is
`session45_injective_multisource_moves_20261023/balanced-circuits/THEOREM.md`
§§2,5: when the two targets have the same paired charge, a balanced source
relation preserves every retained dyadic node-current guard. That file's
middle-binomial argument produces SOME compatible exchange. The result below
instead decomposes ANY prescribed feasible endpoint difference under a stated
independence hypothesis.

The canonical hostile is the same session's
`injective-completion/THEOREM.md` §3: distinct raw paths can consume the same
normalized resource. Resource identities are the full F/R orbit keys; hidden
kernel lifts do not create extra capacity. The corrected near miss is replacing
the full higher-dyadic relation by its binary projections. The sixteen-slot
mod-8 menu has 196 accepted states although every pair of choice bits is free.
The least-used sidecar here is the physical source/target row: before evaluating
a charge, identify repeated owner rows, as required by
`dyadic-constraints/THEOREM.md`. No abstract source-copy index may hide that merge.

Our small concept board is: paired source charge; exact degree margin;
normalized resource; retained owner/clock record; endpoint-mixture precedence.
The maps below preserve the first two additively. The last three require their
own hypotheses and are not consequences of an additive code.

Fix a finite abelian group G and a finite supply of actual binary U/V slots.
Each slot retains its source owner, both labelled Hamiltonian paths, literal
indexed pins, kernel lift, converse partner, chosen original per-slot phase,
and normalized resource identities. Assume:

1. the two targets have equal paired charge in the selected harmonic quotient;
2. the complete source-row criterion is the equality of the signed source
   charge sum in G, after all repeated physical rows have been combined;
3. both alternatives pass every required per-slot predicate, including the
   SAME requested original phase in that slot; every mixture is compatible
   with the declared shared records and local path/pin predicates;
4. every mixture is resource-injective and avoids the fixed remainder's occupied
   resources.

A sufficient, stronger form of 4 is that the union of all option-resource
reservoirs is disjoint between slots and from the fixed remainder. We do NOT
require the two choices of one slot to consume disjoint resources, since they
are never occupied simultaneously. All-mixture compatibility is the exact
hypothesis used in the proof; checking only the two endpoints is insufficient.

## 2. Every endpoint difference admits bounded monotone blocks

Let D(G) be the least integer such that every sequence of D(G) elements of G has
a nonempty zero-sum subsequence. Repetitions are allowed. A minimal zero-sum
sequence has length at most D(G): apply the definition to any D(G) entries of a
longer one.

**Theorem.** Under the preceding supply hypotheses, any two feasible assignments
with the same target margins and node-current signature are connected by a
sequence of actual moves changing at most 2D(G) slots each. Each differing slot
changes exactly once; every intermediate assignment uses only the two supplied
endpoint records. The number of blocks is at most half the number of differing
slots. Every block preserves the degree margins, retained shared records,
resource injectivity, chosen per-slot phases and selected harmonic bank.

**Proof.** Let A be the slots changing U to V and B those changing V to U.
The equal-margin condition gives |A|=|B|=m. Equality of currents gives

    sum_(a in A) label(a) = sum_(b in B) label(b) in G.

Pair A and B arbitrarily, giving differences
g_i=label(a_i)-label(b_i), i=1,...,m. Their sum is zero. Extract an
inclusion-minimal nonempty zero-sum subsequence, delete it and repeat; the
remaining total is always zero. Each piece has at most D(G) pairs. Change both
slots in each pair of that piece. The changed sources balance in cardinality
and charge, so the exact current criterion applies. Conditions 3 and 4 preserve
the full physical records at every intermediate mixture. No slot changes twice.
There are at most m nonempty pieces. This proves the theorem.

If desired, refine each block into inclusion-minimal balanced source relations;
this never increases support. Minimality of the initially extracted paired
sequence alone does not imply minimality among ALL differently paired balanced
subrelations. The latter refinement makes that distinction explicit.

With nonnegative fixed per-slot costs c_s, the same construction has block cost
at most the sum of the 2D(G) largest costs among changed slots (or all costs if
there are fewer), hence at most 2D(G) max c_s. Costs count actual changed records;
owner-clock distance must be supplied by the individual path comparison.

**Constructive extraction.** For the first min(D(G),number remaining) pairs,
subset-sum dynamic programming finds a minimum-cardinality nonempty zero-sum
witness. A shorter window still has a witness because it then contains the
whole remaining zero-sum sequence. Minimum cardinality implies inclusion
minimality inside the window.
Repeated extraction costs O(m D(G)|G|) group transitions with straightforward
witness storage. This bound is deliberately conservative; authenticating the
actual native supply is a separate cost.

## 3. An elementary prime-power bound and its sharp scope

For a finite abelian p-group

    G = Z/(p^a_1) + ... + Z/(p^a_r),
    D(G) = 1 + sum_i (p^a_i - 1).

This is the classical p-group theorem of John E. Olson,
[A combinatorial problem on finite Abelian groups, I](https://www.sciencedirect.com/science/article/pii/0022314X69900213),
Journal of Number Theory 1 (1969), 8–10, DOI 10.1016/0022-314X(69)90021-3.
The new application here is the scoped physical endpoint decomposition, not
the value of the Davenport constant. Here is a self-contained proof. In the
group algebra F_p[G], write t_i=[e_i]-[0]. Its
presentation is

    F_p[t_1,...,t_r] / (t_i^(p^a_i)).

The presentation follows from (1+t_i)^(p^a_i)=1 in characteristic p and the
matching basis of |G| monomials. The augmentation ideal has nilpotence index
at most d=1+sum_i(p^a_i-1), since every degree-d monomial contains a forbidden
power. Each (1-[g_j]) belongs to that ideal. Thus the product of any d such
factors is zero. If no nonempty subsequence sums to zero, however, its identity
coefficient is 1, contributed only by the empty subsequence, a contradiction.

For sharpness, take p^a_i-1 copies of each basis vector e_i and one final vector
e_1+...+e_r. Their total is zero. A zero-sum subsequence omitting the final vector
would need each basis count divisible by p^a_i and hence every count zero.
A subsequence including the final vector needs every basis count p^a_i-1,
so it is the whole sequence. Its length is d. For the trivial group D=1,
represented by the one-term zero sequence.

The physical support bound 2D(G) is sharp in an ABSTRACT independently compatible
supply. Give D(G) distinct positive source owners the labels of this minimal
zero-sum sequence, and give D(G) distinct negative source owners label zero.
Assign all slots distinct resources and the same fixed shared-record values.
Equal margins force the same number of positive and negative choices. A
nonempty accepted mixture must include the whole minimal positive sequence,
then all negative slots. There are exactly two feasible guarded states; the
only move changes 2D(G) slots. Distinct source owners may have repeated charge
labels. This does not assert that an actual tournament atlas realizes these
multiplicities, paths or resource reservoirs.

The original order-eight paired code has

    G_8 = Z/4 + (Z/2)^4,       D(G_8)=8,
    G_dyadic = Z/128 + (Z/2)^4, D(G_dyadic)=132.

Consequently independently compatible original supplies have monotone
changed-slot budgets at most 16 and 264, respectively. For a fixed difference
with 2m slots, replace these by min(2m,16) and min(2m,264).
These are not maximum primitive sizes of the actual finite atlas. In particular
the inherited complete twenty-source supply has sharp all-dyadic primitive
maximum 16 by its own parity exclusion and finite certificate. There is no
conflict between an atlas-specific ceiling and the general group bound.

Only the p-group formula is asserted. No formula for D(G) for arbitrary mixed
prime type is inferred; the theorem in §2 still uses the definition of D(G).
The elementary general estimate D(G)<=|G| follows by taking partial sums of any
|G| terms, including the empty sum, and applying the pigeonhole principle.

## 4. A vacancy path defeats every charge-only bound

The independence hypotheses cannot be replaced by feasibility of the two
endpoints. This failure does not require a resource cycle or a large implication
SCC.

For any m>=2, introduce m slots numbered 0,...,m-1 and m+1 unit-capacity resources
p_0,...,p_m. Slot i chooses old resource p_i or new resource p_(i+1).
At the old endpoint p_m is vacant; at the new endpoint p_0 is vacant.
Both endpoints are injective. A mixture with change bits b_i is injective iff

    b_i <= b_(i+1)  for i=0,...,m-2.

Indeed the only possible collision is new p_(i+1) in slot i versus old p_(i+1)
in slot i+1. Thus feasible mixtures are precisely suffixes of changed slots,
and the implication graph is an acyclic chain with singleton SCCs.

Use G=Z/2. Give resource p_0 and p_m charge0 and each internal resource charge1;
retain the total occupied resource charge as a guard. The slot-change charges
are (1,0,...,0,1). Every proper nonempty suffix has charge1. Only the all-old and
all-new mixtures have the original total charge. Hence the only guarded move
changes all m slots. This is unbounded while D(Z/2)=2. It is also an unrestricted
lower bound INSIDE these complete binary menus, since there are only two
accepted states, not just a monotone-path obstruction.

The source of failure is visible: middle slots have individually zero charge,
and the first/last pair balances, but those separate blocks violate the resource
precedence. The guard couples a vacancy's travel to its boundary label. Nothing
here claims these arbitrary charges are genuine tournament harmonic charges.
It is an exact full-record/resource compiler hostile, with its native
realization left separate.

For comparison, discard the resource collisions while retaining the same
charges. All middle slots may move singly, and the first/last pair may move
together; the required budget becomes2. Thus the missing information is exactly
the physical occupancy relation, not a hidden arithmetic carry.

There is also a balanced U/V version in the ABSTRACT record interface. Take m
even, tag p_j with U for even j and V for odd j, and keep the same capacity
chain. Each slot now genuinely switches its target, and the two endpoints
contain equal numbers of U and V records. Give the first and last source slots
label g=1 in G=Z/2, and every middle slot label0. The signed target-margin and
source-charge changes are

    q_i = ((-1)^i, g*1_{i in {0,m-1}}) in Z plus Z/2.

The full sum is zero. Every proper nonempty suffix of odd length has nonzero
integer margin; a proper suffix of even length has zero margin but group
charge g. Thus only the two endpoints pass BOTH exact margins and group charge.
The forced budget is m, now inside the formal balanced binary interface; m=6
is the first even member exceeding2D(Z/2)=4. If capacity conflicts are removed,
the first/last pair and adjacent pairs of middle slots each balance both
coordinates, so budget2 suffices. Source labels, target tags and resource names
in this variant are declared abstract record data, not claims that the actual
tournament orbit carrier realizes these combinations.

## 5. What is preserved, and the next obligation

Source: two complete supplied endpoint packet assignments. Target: an ordered
list of actual balanced changed-slot blocks. Preserved predicates: exact target
margins, selected harmonic current signature, authenticated local phases and
pins, shared records and resource capacity. Discarded information: none from the
endpoint records; the group code is used only after its exact guard criterion
has been authenticated. Sidecar: the all-mixture compatibility/resource test.
The three-slot vacancy path already exposes the lost precedence. Five slots
first exceed2D(Z/2)=4 in the generic resource/charge interface; six slots do so
in the balanced U/V variant with exact integer margin. Arbitrary length proves
that no charge-only uniform bound survives the retained precedence.

The theorem proves a route inside a chosen endpoint mixture cube. It does not
preclude a cheaper route using records outside that cube, temporary source
vacancies or additional compensation slots. Conversely, the vacancy-path
counterexample does not show any obstruction in an independently compatible
native supply.

The root packet's endpoint-implication compiler supplies the correct interface
when independence fails: resource and retained shared-record predicates first
produce a poset, and the additive guard then selects zero-charge filters. This
file gives the edgeless-poset quantitative bound and an acyclic-chain boundary.
The next native question is how much additional compatible source supply lowers
the actual bottleneck, rather than how many formal zero-sum relations exist.

## 6. Reproduction and audit scope

From the repository root:

    python 05-knowledge/results/native_packet_decomposition_20261023/algebra/check.py
    python -O 05-knowledge/results/native_packet_decomposition_20261023/algebra/check.py

The checker uses explicit exceptions, not Python assertions. It exhausts stated
small p-group sequences and balanced labelled supplies, verifies the constructive
pairing decomposition, and compares the vacancy-path resource test against
literal occupancy and exact group sums. It also verifies source-code sizes and
hashes the inherited theorem/receipt dependencies. Finite checks support the
proof; they are not the all-size argument. `--write` regenerates the receipt.
