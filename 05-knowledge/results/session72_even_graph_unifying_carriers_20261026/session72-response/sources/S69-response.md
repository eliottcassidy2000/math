# Session69 — a trace-polynomial prime action and a revisiting-port convolution

Status: **PROVED** elementary derivations relative to the explicitly **CITED** canonical decomposition/realization inputs; **FINITE-EXACT** controls independently audited. No proof-assistant verification or unrestricted census is claimed. The inherited 12-class Paley decoration census is explicitly recovered, not new. The new contributions proposed here are its complete field-coordinate classifier and marked move law, and a recursively closed multi-rooted response convolution for literal complement lines under arbitrary canonical substitution. The latter has a large finite state space; it is not a bounded/minimal-state theorem or an unrestricted prime atlas.

Input: read-only repository `61d7cca5d8df80169f939268269239c2f67ea724`.

## 0. Inheritance, concept board, and scope

Closest mechanisms: THM-3121 `path-cover-walk-content-substitution-kernel`; canonical strong-module classification; the rooted ordered-join convolution M7 in `05-knowledge/results/tournament_composition_classification_20261026/metagraph/RESULT.md`. Canonical hostile: the order-four same-blue-loop half-edges whose singleton extensions have H=11 and H=13. Corrected near miss: THM-4803's false all-order automorphism invariance. Least-used sidecar: the ordered component endpoint marks in a spanning path cover, retaining coincidences from singleton runs.

Board: (i) actual prime permutation action; (ii) finite-field coordinates; (iii) complete C3 realization fibre; (iv) ordered run allocations; (v) literal complement targets; (vi) endpoint flux as an intentionally lossy observer. The anchor is unrestricted prime composition; the niche is a complete action stratum; the wildcard is whether full endpoint flux distinguishes its three weight-three types (it does not).

**Recovered prior.** `cyclictriple_composition_20261026/README.md`, Theorem G and its noncyclic control, already gives 12 Paley7 decoration classes, the stabilizer strata 2/6/4, weight multiplicities, and the cube-quotient spectrum. THM-211 already identifies the two Fano systems; THM-4553 gives the order-21 group. Session68's incoming-c9a88234 copy contains this material. Some files are absent from the sparse local filesystem but are present as Git objects; this is not missing repository history. No priority claim is made for these recovered facts.

The cyclic-triple reconstruction input is **CITED**, not re-proved: Boussaïri–Chergui–Ille–Zaidi, *3-uniform hypergraphs: modular decomposition and realization by tournaments*, https://arxiv.org/html/1805.04607, Theorems 10 and 24, Proposition 51 and Corollary 57; primary HTML checked 2026-10-26. The local derivations below are elementary relative to that structural input.

## 1. The precise F8 tournament and the 1,2,4 orbit

Let F=F2[α]/(α³+α+1), and label its seven nonzero elements by α^i, i in Z/7. Let Tr(z)=z+z²+z⁴. Define, for distinct x,y in F*,

    x → y  iff  Tr(y/x)=0.

The trace-zero nonzero elements are α,α²,α⁴, since α⁴=α²+α. Their exponents {1,2,4} contain exactly one element from each inverse pair in Z/7. Thus this is a tournament, and exponent coordinates give the Paley7 rule j−i in {1,2,4}. This trace test is specific to this field; it is not claimed for all F_(2^r). In F32=F2[α]/(α⁵+α²+1), z=α and z^−1=α⁴+α both have trace0. The same rule would orient this pair both ways, so the generalization already fails at degree5.

The maps x↦u x^(2^j), u∈F*, j∈{0,1,2}, preserve the relation and form ΓL1(8)=C7⋊C3, of order21. This is the **full** automorphism group: an automorphism fixing exponent0 acts on its three out-neighbours, which induce a C3, so there are at most three choices; the in-neighbours are individually distinguished by their incidence with these out-neighbours. Translation is transitive, so the group has at most7·3 elements, all already exhibited.

Paley7 is modular-prime. Any proper module M lies wholly in the in- or out-neighbourhood of a vertex outside it, so |M|≤3. Global regularity makes its internal scores equal, so a nontrivial M would have size3 and induce C3. Each member then dominates two common exterior vertices. Its out-neighbourhood contains its unique internal out-neighbour dominating both those exterior vertices, contradicting the fact that every Paley out-neighbourhood is C3.

The translates of {1,2,4} and of {3,5,6} are two disjoint seven-line Steiner triple systems. Each base triple has all six nonzero ordered differences exactly once, so each system is a Fano plane. All these triples are cyclic in Paley7. A regular seven-tournament has 35−7·3=14 cyclic triples, proving that they exhaust C3(Paley7). The entire cyclic-triple hypergraph is **two** Fano planes, not one.

Independent direct path enumeration gives H(Paley7)=189 and tc=H/21=9. Neither 7 nor21 is an attained raw H here; 21 is a group order. THM-1370's all-order raw-H exclusions use strong-component multiplication and strong-order floors. The normalized-tc exclusions have different scope and cannot be inferred from that proof.

## 2. A complete trace-polynomial classifier, not a necklace substitution

Every function w:F*→F2 has a unique expression

    w(x)=a+Tr(bx+c x³),      a∈F2, b,c∈F.                    (F1)

Indeed the right side is

    a + bx+b²x²+b⁴x⁴ + cx³+c²x⁶+c⁴x⁵.

If it vanishes on all seven nonzero field elements, a polynomial of degree≤6 has seven roots, so every coefficient is zero. There are128 parameter triples and128 functions, proving the bijection. Also a=sum_(x∈F*) w(x) in F2: it is weight parity.

Pullback by x↦ux gives (a,b,c)↦(a,bu,cu³); pullback by x↦x² gives (a,b,c)↦(a,b⁴,c⁴). For b≠0, normalize b to1 by u=b^−1. The remaining invariant is

    r=c/b³, modulo r↦r².                                   (F2)

Consequently, for each fixed a there are exactly six ordinary action orbits:

| parameters | residual parameter orbit | stabilizer order |
|---|---|---:|
| b=c=0 | — |21|
| b=0,c≠0 | normalize c=1 |3|
| b≠0,r=0 | {0} |3|
| b≠0,r=1 | {1} |3|
| b≠0,r∈{α,α²,α⁴} | degree3 orbit |1|
| b≠0,r∈{α³,α⁶,α⁵} | degree3 orbit |1|

For b≠0 the surviving group is Gal(F/F2(r)), of order3/[F2(r):F2]. This is a field-action classifier for the **full** 21-element action. The older normal-basis necklace model retains only the C7 subgroup and has20 classes; it is not the right quotient here.

The function module is F2⊕F8⊕F8, of dimensions **1+3+3**, corresponding to the doubling cyclotomic cosets {0},{1,2,4},{3,6,5}. It is not a representation-theoretic 1+2+4 decomposition. The literal 1+2+4 exponent orbit and the divisor identity s(8)=1+2+4=7 must not be confused with invariant subspace dimensions or a canonical group-invariant flag.

Let A0=1→C3 and A1=A0^op, with a fixed cyclic orientation on the same triangle labels (changing that orientation changes no ordinary block type). Put T_w=Paley7[A_(w_i)]. The canonical maximal strong blocks force

    Aut(T_w) ≅ C3^7 ⋊ Stab_ΓL1(8)(w).                       (F3)

The semidirect product uses the actual seven-slot action. The cited C3 realization theorem and the self-converse root imply that all ordinary realizations of this fixed cyclic-triple hypergraph have such a representative. Thus the recovered 12-type census is now classified by F1–F2. The complete labelled fibre has2^15=32768 elements; Aut(C3(T)) has order6^7·42, giving respectively256,1792,5376 labelled lifts for each class with stabilizer21,3,1. The class counts2,6,4 give the required total.

Global converse, returned to the fixed root by x↦x^−1, sends

    (a,b,c) ↦ (a+1,c²,b⁴).                                  (F4)

There are six converse pairs and no self-converse decorated class, since the number of A1 blocks changes from k to7−k. These pairs are not merged in the ordinary atlas.

## 3. Pointed field coordinates and the exact modular-flip port

Mark one quotient vertex z∈F*. Normalize it to1 by x=z y, and write

    B=bz,   C=cz³.

A **block-pointed** ordinary class (one marked quotient slot, not an arbitrary marked vertex inside the28-vertex tournament) is exactly (a,B,C) modulo simultaneous Frobenius (B,C)↦(B²,C²). For each a, the64 pairs have four fixed pairs in F2² and20 orbits of length3. Hence there are48 block-pointed classes, equivalently48 directed modular-flip port orbits. Their point stabilizers have order3 for8 classes and order1 for40 classes. If an actual vertex of T_w is marked instead, there are96 ordinary types: each marked block has two invariant inner vertex orbits, its singleton and its cyclic triangle.

The indicator of z is

    δ_z(x)=1+Tr(z^−1 x+z^−3 x³).

This follows from the geometric sum Σ_(j=0)^6(x/z)^j: it is1 at x=z and0 otherwise. Flipping the block at the marked point is therefore the exact involution

    (a,B,C) ↦ (a+1,B+1,C+1),                                (F5)

which commutes with Frobenius. It has24 unoriented port orbits. Forgetting the mark before applying F5 is precisely the loss of which block is changed. This is the modular-cut operation; it is **not** the literal Hamilton-path-complement line operation treated next.

## 4. Multi-rooted cover responses close the exact prime-substitution interface

Let Q be a prime tournament on q≥3 slots, and T=Q[S1,...,Sq] with nonempty factors. Let K⊆Aut(Q) preserve their ordinary isomorphism types. These substituted blocks are the canonical maximal proper strong modules, so Aut(T)=(∏Aut(Si))⋊K after choosing coherent child models.

For c≥1, let C=(P1,...,Pc) be an **ordered** spanning cover of S by nonempty directed paths. Preserve the ordered path words at the source. Put

    U_C = S^op xor E(C),
    ρ(C) = [ U_C ; (first(P1),last(P1)),...,(first(Pc),last(Pc)) ].   (P1)

The bracket denotes a multi-rooted ordinary isomorphism type. Repeated marks are retained (a singleton path has equal first and last). Define R_(S,c)(ρ) to be the number of Aut(S)-orbits of ordered covers with responseρ. Equivalently it is the raw cover count forρ divided by|Aut(S)|. This division is integral: an automorphism fixing an ordered spanning path cover fixes every ordered vertex of every component.

For a directed quotient walk w with content c=(c1,...,cq), assign the rth occurrence of slot i to the rth ordered child component. Form Glue_Q(w;ρ1,...,ρq) by:

1. taking the local tournaments U_i fromρ_i;
2. orienting all cross-block pairs by Q^op;
3. for each two consecutive letters in w, toggling the single hinge from the end mark of the first component to the start mark of the next.

There are no equal consecutive letters in a quotient walk. Every marked component is used exactly once, so hinge pairs are distinct. This construction depends only on the multi-rooted types, not the chosen representatives. Write [Glue] for its ordinary target class.

**Prime-substitution response theorem.** The number M_T(Z) of ordinary Hamilton-path presentations of T whose literal complement target has ordinary class Z is

    M_T(Z) = (1/|K|) Σ_c Σ_(w∈Walk_Q(c))
                    Σ_(ρ1,...,ρq) (∏_i R_(Si,ci)(ρi))
                    ·1{[Glue_Q(w;ρ1,...,ρq)]=Z}.             (P2)

All sums are finite:1≤ci≤|Si|. This is an exact target-incidence convolution, not merely an H or flux formula.

*Proof.* Split a Hamilton word of T into maximal constant-block runs. Each child receives an ordered spanning path cover in order of occurrence; the block word is a directed quotient walk. Conversely those data concatenate to a unique Hamilton word. Inside child i, complementing nonpath arcs produces U_i. Between children, every arc is reversed except the consecutive inter-run hinges, exactly the gluing construction. This proves the target formula on labelled words. Quotient by the independent child automorphisms: their free action produces the product of the R multiplicities. The remaining K action is free on these underlying kernel-orbits, since a stabilizer would lift to an automorphism fixing the entire global ordered Hamilton word. Each target-class category is K-invariant. Dividing its aggregate count by|K| gives P2. Equivalently K acts freely already on the retained labelled quotient word, since every slot occurs. A response marginal after discarding that word need not be a free K-carrier; multiplicities and the walk must remain in the interface. ∎

For an ordered join only one quotient word exists and every ci=1. P2 specializes to the inherited endpoint-rooted hinge convolution M7. For a cyclic prime quotient, repeated visits are allowed and the full multi-rooted cover response is needed. The original modules generally shatter in the target; the gluing construction does not pretend they survive.

Loss ledger: a child path-cover profile remembers only Σ_ρ R; it forgets the target relation and endpoint marks. A score-labelled endpoint tensor remembers one flux observer but not the multi-rooted target. Fullρ is sufficient to reconstruct target incidence. The theorem does not establish smallest possible state space, a bounded-size state family, efficient arbitrary-prime enumeration, or a new prime census. Blue/black reflection colour additionally requires retaining the full marked-presentation reflection predicate; P2 is ordinary uncoloured target incidence.

### 4.1 The response family is recursively closed (P3)

The same construction computes every parent response R_(T,d), not only the d=1 ordinary line ledger. Replace w in P2 by an **ordered d-tuple of nonempty quotient walks**, of total content c. Assign the child component occurrences in the concatenation of the tuple, maintaining a single occurrence counter for each child. Toggle hinges only between consecutive letters **within the same walk**. The output is the resulting tournament with the ordered pair of outer endpoints of each of the d walks. Replace the indicator in P2 by equality with this output multi-rooted responseρ. Then

    R_(T,d)(ρ) = (1/|K|) Σ_c Σ_(ordered walk tuple w1,...,wd)
                          Σ_(ρ1,...,ρq) (∏_i R_(Si,ci)(ρi))
                          ·1{multi-rooted Glue=ρ}.          (P3)

There is **no d!**, since the global path covers and walk tuples are ordered. Maximal block-run splitting separately on each global path gives the same bijection and exact target reconstruction; the union of the tuple visits every quotient slot, so K acts freely. The argument also applies to a canonical linear node, where the child order is fixed and K is trivial. At a singleton the only response is its single vertex marked twice, with multiplicity1. Induction up the canonical strong-module tree therefore computes the whole response family from these singleton bases and the retained quotient relations/actions. Child responses need not be supplied by a separate path enumeration.

This is a closed exact recursive **ordinary, uncoloured target-incidence carrier**, not a theorem that its state space stays small. The response does discard internal path words, retaining their multiplicity and the full local target relation with endpoint marks; the non-unit example below proves this is genuine aggregation rather than merely renaming every path. Canonicalizing the multi-rooted targets and listing arbitrary prime quotient inputs may still be expensive. No field or numerical shortcut replaces those tasks.

### 4.2 Why the extra marks and weights are necessary

The independent audit supplied two concrete controls, then rechecked directly by the producer. Use upper-triangle bits in lexicographic pair order, bit1 for i→j.

1. Let S4 have bitcode7 (a transitive tournament). The ordered covers

       C=((0,1),(3,2)),       D=((0,2,1),(3,))

   have identical local complemented target isotypes and even identical outer-first/outer-last rooted target, canonically `((0,1),4)`. But in C3[S4,1,1], the quotient walk `(0,1,2,0)` produces literal target canonical codes81 and83. The two-root ordered-join response thus **fails under revisits**. Retaining every ordered component endpoint pair repairs this explicit failure; no assertion of absolute minimality is made.
2. Let S5 have bitcode10 and Aut order1. Some local response types have multiplicity2. In C3[S5,1,1], P2 with those multiplicities gives H=93 ordinary source ports; replacing every positive local response multiplicity by1 gives88. Thus the weights are not bookkeeping that can be erased. The exact non-unit response records are in `receipt.json`.

## 5. Endpoint specialization and a new Paley28 hostile

For each child S, let F_S(c) count ordered c-path covers. Let L_S(c;s) count those with first component starting at local score s; let E_S(c;t) count those whose last component ends at local score t; and let B_S(c;s,t) impose both marks. The complete endpoint Hamilton response follows from the same quotient-walk bijection:

- for different first/last slots i,j, use L_i(c_i;s)E_j(c_j;t)∏_(h≠i,j)F_h(c_h);
- for the same slot i, use B_i(c_i;s,t)∏_(h≠i)F_h(c_h).

Multiply by the number W_Q^(i,j)(c) of quotient walks of that content and those endpoints, then sum. Add each block's exterior score to its endpoint's local score. This specializes P2 only to endpoint observables; it does not reconstruct target isomorphism classes.

For A0 and A1,

    (F(1),F(2),F(3),F(4))=(3,18,36,24).

Their endpoint profiles are obtained independently by permuting four vertices and choosing cuts, requiring a cut at each backward adjacency. The Paley root walk table needs only5^7 content cells and seven last-slot states for a fixed first slot. Full transitivity and the two ordered-pair orbitals reduce endpoints to (0,0),(0,1),(0,3).

The common raw Hamilton count of all128 decorations is

    H = 5,307,704,039,306,577,787,515.

The ordinary normalized presentation counts are:

| root stabilizer | full automorphism order | classes | tc=H/Aut |
|---:|---:|---:|---:|
|21|45927|2|115568272243050445|
|3|6561|6|808977905701353115|
|1|2187|4|2426933717104059345|

Every decoration has c3=7+64·14=903. A literal partner has flux Δc3=score(first)−score(last)−1. Since every block's exterior score is12, only the local scores enter. The **entire raw flux histogram depends only on k(7−k)**, where k is the number of A1 blocks. This is proved, not fitted: same-colour root arc counts are binom(7−k,2),binom(k,2), and each directed mixed-colour count is k(7−k)/2 by regularity; the two endpoint orbitals have fixed response tensors; the diagonal responses agree between A0 and A1; and the same-colour responses00 and11 agree on each orbital. For the last equalities, reverse every cover/path and apply global converse, returning Q to itself by inversion: exchanging endpoints and the root anti-isomorphism together preserve the ordered-pair orbital, while the score difference and hence flux are unchanged.

In particular the three nonisomorphic weight-three classes have exactly the same raw flux histogram, although their root stabilizers are3,3,1. Thus even H plus the complete raw Δc3 distribution does not classify these ordinary objects. The two cyclic-three-subset orbits (the two Fano systems) have the same H, full automorphism order and flux histogram but remain different ordinary classes. This is the stronger same-group-order hostile. Converse sends weight3 to4, so this pair is not merely an ordinary-versus-converse collision.

For k=3 (and4), the exact histogram is:

| Δc3 | raw paths |
|---:|---:|
|−4|67938585917119459830|
|−3|468402208566561437355|
|−2|1166655991737358177398|
|−1|1555120616167066763790|
|0|1347099003419914838700|
|1|597993950323183023630|
|2|104493683175374086812|

The zero-flux category does not imply equal cyclic-triple hypergraphs or membership in the decoration fibre. Normalizing by the correct source automorphism order counts ordinary source ports; it does not by itself produce a converse-merged degree.

## 6. Boundaries

The trace polynomial supplies a complete action classifier for one modular-prime quotient, and P2–P3 supply an exact recursively closed canonical-substitution response grammar **once the prime/linear quotient relations and their child-type actions are supplied**. This closes a sufficient exact run/allocation/port interface at the full multi-rooted response level, but leaves unrestricted structural prime classification, minimal closed state compression, reflection-coloured/relocation-compatible extensions, and efficient full metagraph generation open. The quoted {7,21} raw-H gaps do not follow from the number21 appearing in ΓL1(8), and no ABC or untouchable-prime theorem is inferred from these field coordinates.


## 7. Exact controls and the nine-to-seven literal prime port map

Run without network or repository access:

    python3 -I -B check.py > replay.json
    python3 -I -B -O check.py > replay.optimized.json
    cmp replay.json receipt.json
    cmp replay.json replay.optimized.json

The standard-library checker exhausts all7! permutations for the root action, every proper subset for primeness, all128 functions and896 marked flips, the two Fano systems, and the complete nine ordinary Paley source presentations. It compiles every labelled Hamilton path in36 substitutions (quotients an ordered edge or C3; child types K1,T2,C3), totaling6160 paths, and checks the literal target relation against direct complement. In26 cases of total order≤6 it additionally computes the full compressed P2 target-class convolution and compares with independently enumerated ordinary target multiplicities. Source H is independently checked by subset DP. Normal and optimized receipts agree. These are producer finite checks. The independently supplied Session69 tournament audit accepts F1–F5 and P2–P3, recomputes the entire Paley28 endpoint table by determinant/cofactor rational-walk expansion rather than this producer DP, verifies P2 against direct target-class enumeration, and verifies all-cover P3 through order6, including a nontrivial top stabilizer K=3:7367 actual ordered covers and2942 response categories. See the separately packaged audit report and exact reviewed-source manifest. No formal verification is claimed. Additional all-cover P3 controls compare every response channel in five parent substitutions through order5, and the two missing-port/multiplicity hostiles are directly reproduced.

For Paley7 the nine source ports land at only **seven** ordinary target classes. The canonical target encoding lists bits for pairs(i,j), i<j, lexicographically, least significant bit first, and minimizes over all vertex permutations.

| source Hamilton word | target canonical code | target H | target Aut |
|---|---:|---:|---:|
|0123456|1592|151|1|
|0124563|5682|151|1|
|0124635|9400|153|1|
|0126345|9528|153|1|
|0134526|9585|153|1|
|0135246|9649|151|1|
|0152346|9400|153|1|
|0152463|9585|153|1|
|0152634|9272|155|1|

All targets have13 cyclic triples. The repeated target classes9400 and9585 are prime; target codes1592,5682,9649 have pair modules. Thus target-class multiplicity cannot be omitted even for this explicit prime source. This does **not** refute HYP-9414's carousel injectivity claim or the narrower conjecture for sources whose every partner is prime: Paley7 has nonprime partners. The old census already knows that universal injectivity is false; this is a concrete field-action-stratum control, not a global novelty assertion.
