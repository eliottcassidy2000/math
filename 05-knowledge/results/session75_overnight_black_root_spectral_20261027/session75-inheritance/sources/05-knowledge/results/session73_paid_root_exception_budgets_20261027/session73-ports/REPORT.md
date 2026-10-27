# Session73 — Shared-port grammar beyond the forest: exact debt and a Paley11 barrier

**Status:** PROVED within the explicitly stated interfaces; finite controls are
FINITE-EXACT. No claim of historical novelty for weighted feedback-arc
substitution or the graph cycle-space identity. Native all-word reachability,
a uniform relative-exception grammar, paid ROOT positivity, and the economical
coloured aggregate remain OPEN.

## Useful results first

1. **One legal complement cannot turn any forest-framed tournament into
   Paley11.** Every such target on n≥4 has minimum feedback-arc count at most
   2n−4. Paley11 has at least19, but 2·11−4=18. This contrasts with the inherited
   explicit forest→Paley7 transition. Exact finite subset DP gives20 for Paley11,
   but that stronger computation is not a premise of the proof.
2. **A genuine uniform linear total-rank exception budget fails on elementary
   substitution families.** C3 towers with n=3^k have minimum exception-edge
   count n(n−1)/6. For k≥2 their minimum exception-cycle debt is exactly
   (n²−7n+6)/6. They nevertheless have linear-size constant-arity substitution
   descriptions. This is not a proof they are reachable by the narrower native
   J/star/relocation/complement constructor.
3. **The inherited family P7[C3^[k]] links orders7 and21 structurally:** its
   exact edge debt is n(n−1)/6 and exact cycle debt (n²−7n+6)/6. At orders7,21,63
   these cycle debts are1,50,589. This is a tournament substitution statement,
   not a bridge to raw-H gaps, tc gaps, or prime untouchables.
4. **Exact weighted substitution and exact shared-port debt compose.** The
   former retains positive vertex costs; the latter retains actual component
   incidence and parallel port edges. Supplied source masses, source names,
   active paths and integer clocks remain separate sidecars.

## 1. Inheritance, scope and source pins

The local source snapshot manifest is `SOURCE_PINS.json`. The starting shared
checkout was `92e7f12f3b95`; no shared repository file is edited by this lane.

- `f4ef2ab273c1442aa88a1079a92fb6a2c2aabd80`,
  `05-knowledge/results/weighted_renewal_moving_modulus_20261026/README.md`,
  T1/T2: exact rank/exception decoder; relocation, J/star and converse;
  complement D→D△E(P); forest→Paley7; exactly r independent cycles at7r
  vertices. The forest-to-Paley7 construction and the P7 cycle obstruction are
  inherited, not new here.
- `05-knowledge/results/native_return_necklaces_20261026/README.md`, N1–N3:
  shared anchor/source archive, separate controller/primitive clocks,
  fixed-port synchronous no-reset, separable LCA grammar and its order-four
  relocation escape. The old two-current-endpoint interface does not survive
  arbitrary relocations.
- `00-navigation/S69.md`: ordinary uncoloured all-arity multi-rooted cover
  composition is CLOSED; a small economical coloured aggregate is not.
- `00-navigation/S72.md`: deterministic full cyclic-frame relative phases
  reconstruct a common-shift orbit; arbitrary native Paley markings can have
  indistinguishable named pair types but different triple guards. Neither law
  supplies the missing port incidence in the other setting.
- `THM-4743` plus full slug
  `01-canon/theorems/THM-4743-converse-extended-gluck-for-tournaments-primitive-classification-paley-tournaments-and-c3-towers.md`:
  the C3 towers and P7[C3^[k]] are inherited named substitution families. Its
  converse-palette classification is not used to prove our debt formulas.

Anchor: paid exception debt beyond P7. Niche: exact component-incidence
assembly. Wildcard: the same quadratic debt on the inherited7→21→63
substitution family. Live board: **rank cost / cycle debt / named ports /
substitution baseline / source masses / physical clock**.

## 2. Definitions: two different debts, and two different weight types

Let T be a tournament on a nonempty finite named vertex set V. For a total
order r, let D_r(T) be the undirected set of backward pairs. Equivalently,

    x→y in T iff (r(x)<r(y)) XOR ({x,y}∈D_r(T)).

Set

    f(T) = min_r |D_r(T)|,
    b(T) = min_r beta(D_r(T)),
    beta(D) = |E(D)|−|V(D)|+components(D).

All vertices, including isolated vertices, count in components. Thus

    max(0,f(T)−n+1) ≤ b(T) ≤ f(T).

The left inequality uses components≥1. A connected exception graph in an
f-optimal order attains f−n+1. Minimizing edge debt does NOT by itself minimize
cycle debt: the natural C3-tower order is a concrete hostile below.

For optional strictly positive **vertex costs** a_v define

    f_a(T) = min_r Σ_{backward {u,v}} a_u a_v.

These pair-factor costs are not the **source masses** μ_s in a supplied
source-labelled law. No theorem below replaces native μ_s by vertex sizes,
uniform probabilities, or freshly chosen masses. The probabilities used in
the proof of G1 are an auxiliary rounding argument, not a dynamical law.

## 3. G1 — exact weighted substitution law

**PROVED.** Let Q be a tournament on indices1,…,q; T_i be nonempty tournaments
on pairwise disjoint named vertex sets V_i; and T=Q[T_1,…,T_q] be their ordinary
lexicographic substitution. Write A_i=Σ_{v∈V_i} a_v. Then

    f_a(T) = Σ_i f_a(T_i) + f_A(Q).                    (G1)

In the unweighted case a_v=1, the quotient edge cost is |V_i||V_j|.
No primality, equal-size, transitivity or stabilizer hypothesis is needed.

**Proof.** Fix any global order r. Its internal cost is at least Σ_i f_a(T_i).
Independently choose a vertex X_i from each V_i with probability a_v/A_i.
The chosen vertices inherit a random order σ on the quotient indices. For a
quotient pair i,j, the expectation of its weighted backward indicator
A_iA_j·1_{i,j backward under σ} equals the actual sum of a_u a_v over backward
cross-pairs u∈V_i,v∈V_j in r. Therefore the total cross cost is the expected
weighted feedback cost of σ, hence at least f_A(Q). Conversely concatenate
individually optimizing orders of the T_i in an A-weighted optimizing order
of Q. This attains both lower bounds simultaneously. QED.

**Equal-sized specialization.** If every |V_i|=m, then

    f(Q[T_i]) = Σ_i f(T_i) + m² f(Q).                 (G1a)

**Perturbed interface.** If a named set H of cross-pairs is reversed after
substitution, then

    |f_a(T△H)−f_a(T)| ≤ Σ_{{u,v}∈H} a_u a_v.         (G1b)

Every fixed order changes cost by at most that sum; optimize in both
directions. This is a robust bound, not an exact replacement for G1.

**Boundary.** Ordinary substitution preserves module blocks; legal path
complement can shatter those modules. Hence G1 is an exact composition theorem,
not a statement that the modular decomposition survives native complement.

## 4. G2 — exact component-incidence law for shared ports

**PROVED.** Take pairwise disjoint child vertex sets V_i with arbitrary simple
exception graphs D_i. Let H be a simple set of exceptional cross-pairs whose
endpoints are the actual named ports. Repeated mentions of the SAME physical
pair must first be combined by XOR if they denote toggles. Distinct physical
pairs remain distinct. Let D be their union.

Form a multigraph B with one vertex for EACH connected component of every
D_i, and one edge for EACH member of H, using the components containing its
two actual endpoints. In particular retain parallel edges. Then

    beta(D) = Σ_i beta(D_i) + beta(B).               (G2)

**Proof.** Put c_0=Σ_i components(D_i). Contracting all child components
preserves the final component count c=components(D)=components(B).
Then |E(D)|=Σ_i|E(D_i)|+|H| and |V(D)|=Σ_i|V_i|, so

    beta(D)−Σ_i beta(D_i)=|H|−c_0+c=beta(B).

QED. Bridges add no new debt, recovering the inherited J/star law.
Multiple hinges between already connected pieces can create debt even when
no new child is cyclic. Shared endpoints themselves need not create a cycle.

**Small hostile.** Each of two children has one internal edge. Attach two
cross-edges at distinct endpoint pairs to form C4. Both children have beta0;
B has two vertices and two PARALLEL edges, beta1. Collapsing B to a simple
edge wrongly reports0. Conversely, three singleton children joined at one
shared central port form a star with beta0; treating each use of the shared
port as a new independent vertex changes the object.

**Operation-faithful extension.** More generally retain any explicitly decoded
substitution baseline B_0 plus a global named exception set D. Relocation
changes only the active path. Full converse changes B_0 to its converse and
keeps D. Legal path complement changes B_0 to its converse and D to D△E(P).
Converse of a substitution baseline is obtained by conversing every template.
J uses a forward two-block root; star uses a backward two-block root plus its
actual hinge. This is an exact decoder/compiler with a paid exception set.
No uniform linear bound on that relative exception set is proved here.

A total-rank baseline is a very special baseline. The quadratic examples below
show why “linear exception edges against a total rank” is much stronger than
“linear-size exact hierarchical grammar.” No module-faithfulness is asserted
for the perturbed decoded tournament.

## 5. G3 — improved native complement tariff and the Paley11 wall

For a legal Hamilton path P of T, tau_P reverses every arc not in P and keeps
all arcs of P. Define

    c(n)=n−1−min(2,floor(n/2))  for n≥2.

Thus c(2)=0,c(3)=1, and c(n)=n−3 for n≥4.

**PROVED.**

    |f(tau_P T)−f(T)| ≤ c(n).                       (G3)

More generally, from any source frame with k exception edges one obtains
f(tau_P T)≤k+c(n).

**Proof.** Start with a source frame (r,D), |D|=k. The inherited exact update
is r'=reverse r and D'=D△E(P). Put s=|D∩E(P)|. If s≥1, then
|D'|=k+n−1−2s≤k+n−3, which is at most k+c(n) for n≥2 where that case occurs.
If s=0, every P-edge belongs to D'. Since P is still a legal directed path
in the target, r' decreases strictly along P. Thus r' is exactly reverse P.
Swap floor(n/2) disjoint adjacent pairs of this total order. Each swapped
pair was a backward P-edge, and an adjacent swap changes no other relative
pair. The new order has k+n−1−floor(n/2) backward edges. Taking the weaker
of the two guaranteed reductions proves k+c(n). Optimize the source rank,
and then apply the same result to the involution tau_P for the reverse
inequality. QED.

The proof for n2 has s≥1 impossible in an optimal source frame; for arbitrary
frames the displayed s≥1 bound still implies the claimed k+c(n) inequality.
For n1 define c(1)=0 separately.

**Forest one-complement corollary.** If some source frame is a forest, k≤n−1;
therefore for n≥4,

    f(tau_P T) ≤ 2n−4.                              (G3a)

**Paley11 corollary (PROVED, finite incidence input verified explicitly).**
For P11 with x→y iff y−x is a nonzero square modulo11, there are55 directed
triangles and each edge belongs to3. Every total order must make at least
one edge of each directed triangle backward. Hence f(P11)≥ceil(55/3)=19.
Since19>18=2·11−4, P11 cannot be a one-complement target of ANY forest-framed
source, with ANY legal path and ANY source/target ranks. Full converse and
renaming do not remove the obstruction.

For clarity the generic doubly-regular triangle quotient q(q−1)/6 is
110/6 at q11, so its ceiling is19, NOT20. Exact subset DP independently finds
f(P11)=20. The proof uses only19.

This is a one-step barrier. It does not exclude two or more complements,
growth, a richer baseline grammar, or a different paid source construction.
Likewise it says nothing about ordinary all-arity response composition, which
was already closed in Session69.

## 6. G4 — quadratic exact edge AND cycle debt on compact families

Let C=C3 and U_0 be one vertex, U_k=C[U_{k−1},U_{k−1},U_{k−1}], n=3^k.
The family itself is inherited; this section computes its frame debts.

**PROVED.** For all k≥0,

    f(U_k)=n(n−1)/6.                                (G4a)

For k≥2,

    b(U_k)=(n²−7n+6)/6.                             (G4b)

At k=0,1 the cycle debt is0.

**Edge proof.** Write m=3^(k−1). G1a and f(C)=1 give
f(U_k)=3f(U_{k−1})+m². The displayed solution follows by induction.
Equivalently each of the m³ cross-block directed triangles needs a backward
edge, and each cross edge belongs to m of them, forcing m² cross debt; a
block order attains it. This alternative lower bound is unweighted and
independent of the rounding proof of G1.

**Exact cycle proof.** The universal bound gives b≥f−n+1. At k≥2 the first
child A is nontransitive. Choose a minimizing order of A and split it into
nonempty prefix L and suffix R across a backward edge. Keep minimizing
orders of children B,C, and order all vertices as

    L, B, C, R.

The only cross backward pairs are L×C and R×B, exactly m|L|+m|R|=m².
All internal relative orders survive. Their cross graph has two connected
pieces L∪C and R∪B; the retained backward edge of A across L,R joins them.
Thus the entire optimal exception graph is connected, giving b=f−n+1. QED.

**Hostile to naive cycle minimization.** Natural lexicographic orders have
k+1 exception components and cycle debt f−n+k+1. At n9, their beta6 is not
minimal: the split-block rank gives beta4. A scalar edge optimum alone is
not the missing shared-port connectivity coordinate.

### The genuine7→21→63 continuation

Let V_k=P7[U_k], n=7·3^k. The inherited natural P7 order attains f(P7)=7
and its exception graph is connected (T2's unicyclic graph). G1a gives

    f(V_k)=7·m(m−1)/6+7m²=n(n−1)/6.

In a block-concatenated optimal order, every P7 exception edge becomes a
complete bipartite defect graph. The connected P7 defect graph touches every
block, so its lift connects every vertex, regardless of internal defects.
Therefore

    b(V_k)=f(V_k)−n+1=(n²−7n+6)/6.                 (G4c)

Values: n7 gives(f,b)=(7,1); n21 gives(70,50); n63 gives(651,589).
The same argument applies to iterated P7 substitution, giving the same
formulas at n=7^ell, ell≥1.

**Connection contract.** Source: the inherited exceptional converse-palette
family P7[U_k]. Map: ordinary substitution followed by the exact rank decoder.
Preserved predicate: labelled tournament arcs and the exact cost optimization.
Discarded by the scalar: marks, port incidence, active paths, stabilizers,
source masses and trajectories. Needed sidecar: actual named substitution
and exception graphs. Decisive test: at21 the exact structural cost is70
edges and50 cycles, not a continuation of raw-H's number21 by analogy.
There is NO map here to the proper-divisor-sum image predicate or to prime
untouchability. Those arithmetic questions remain independently scoped.

### What fails, and what remains open

Any class containing all U_k or all V_k has no uniform linear bound on
minimum rank-exception edges or cycles. The family descriptions themselves
have O(n) nodes in a fixed ternary/seven-ary alphabet (O(n log n) bits once
arbitrary names are included, excluding supplied numerical weights/clocks).
Thus rank-debt and description size are genuinely different resources.

We have NOT proved that the native J/star/relocation/complement closure of
singletons contains these families. Accordingly the all-word linear-budget
question for that narrower reachable class is NOT refuted here. G3 supplies
a necessary paid complement tariff for any proposed native compiler; the
parent synthesis separately accounts for its construction clock.

## 7. Source weights, marks and integer clocks

All constructions are deterministic on each named source record

    (s, μ_s, named vertices, arc decoder, active Hamilton path,
     named shared/external ports, integer operation clock).

Reframing changes the decoder, not the labelled tournament or any μ_s.
A legal relocation changes the active path, not the source law. A physical
complement increments its specified native operation clock; merely optimizing
a rank or contracting a proof graph incurs no physical move and must not be
charged as a dynamical step. If two source records acquire equal displayed
states, their names and masses remain separate in the archive; only an
explicit pushforward may sum them. Refining observation modulus is not an
extra graph move or a source resampling operation.

G2 uses the full joint port incidence. Neither pairwise cyclic phase
reconstruction nor weighted quotient costs authorize identifying shared
ports. The native Paley triple hostile inherited from Session72 remains
active; no pairwise-sufficiency claim is made for this grammar.

## 8. Reproduction and explicit universes

From this directory:

    python3 verify.py > receipt.json
    python3 -O verify.py > receipt.optimized.json
    cmp receipt.json receipt.optimized.json
    sha256sum -c SHA256SUMS

The checker uses exact integers/Fractions and explicit RuntimeError guards,
not removable assertions. Fixed seed73073 specifies the sampled tests.

- Exhaustive marked fixed-path tournaments n2…7, counts
  1,2,8,64,1024,32768. This covers all Hamilton-marked tournaments modulo
  relabelling their path, not merely randomly selected paths. Maximum observed
  |Δf| is0,1,1,2,2,3, below/at G3's0,1,1,2,3,4. No sharpness claim at n≥6.
- Exact independent subset-DP/permutation comparisons: all labelled tournaments
  n2…4 plus20 weighted order-five controls;94 cases total.
-72 deterministic-seed weighted substitution tests with arbitrary quotient
  and uneven child sizes, total order≤8, positive rational vertex costs.
- Paley3,7,11,19 exact triangle and edge-incidence counts. Subset DP at3,7,11
  gives1,7,20. No DP optimality claim for Paley19.
- Literal inherited Paley7 predecessor tree and involution replay.
- C3 towers n3,9,27,81: exact structural upper/connected-rank certificates;
  independent full subset DP at n3,9. Larger optima follow the proof, not an
  exhaustive permutation census.
-200 deterministic-seed component-incidence gluing controls, with parallel
  quotient-edge and repeated shared-endpoint hostiles.
- Unequal named source masses1/7,2/7,4/7, shared named ports, and an integer
  clock preserved by a deterministic complement update. This is a carrier
  sanity check, NOT a new native ROOT estimate.

The source snapshots and code/output hashes make the local claims replayable;
no artifact is silently treated as an all-order existence or ROOT certificate.
