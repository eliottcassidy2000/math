# Return-join selectors: controlled overlap beyond matchings

**PROVED:** an exact return-join cycle/packing recursion, its labelled selector decoder, and a nonmatching six-address initial gate. **FINITE-EXACT:** 225 small return joins and three nested constructors. **OPEN:** arbitrary requested predicates, regenerative selection, and Chan–Pak. All claims below are local theorems, not reserved canon IDs.

## 0. Inheritance and the change of viewpoint

Closest mechanism: Session262 `defect/EXACT_SELECTORS.md`, §§7–8, realizes a diamond mutex and products of matching guards. Hostile: three primitive cycles with conflict graph P3 cannot be the complete odd-cycle list without leakage, under the inherited strong-floor argument there. Corrected near miss: a shared-edge fan with three tips has an extra five-cycle, not merely three triangles. Least-used sidecar: cycle length and the residual tournament after deleting an endpoint path.

Anchor: exact selected source addresses. Niche: endpoint-path return joins. Wildcard: Boolean idempotents translated through an actual cycle-intersection graph. Live board: protected zero source; overlap versus disjointness; newly born long cycles; endpoint paths; all-order cover profiles.

**Change:** do not discard the long cycles. Give each its own address, then compile the complete conflict graph. This is not a way to realize a preassigned P3 with no extra addresses. It is a physical grammar that outputs its exact admissible predicate and permits composition without untracked leakage.

The one-apex tournaments used below are inherited objects, already studied in `session26_random_bridges_20261022/session26-tournament/PROOF.md` and `APEX_CLASSIFICATION_AND_SPECTRAL_LAW.md`. Their scalar H spectrum or continuation automaton is not new here. The new work is their address-level decoder and the return-join composition with full residual sidecars. Similarly, 17/224/1449 is an inherited profile from Session261/262; the claim here is its explicit nonmatching address geometry, not discovery of that profile.

## 1. The literal return-join operation

Take disjoint nonempty tournaments A,B with distinguished vertices α∈A, β∈B. First direct every cross edge A→B. Reverse just α→β, so its direction is β→α. Write R(A,α;B,β) for the resulting tournament. No other arc changes.

A path below is nonempty, directed, simple, need not span, and includes its specified endpoint. Its size means number of vertices. A cycle means a directed simple cycle modulo rotation, not modulo reversal. Only odd cycles of length at least three are considered.

### Theorem RJ1: complete new-cycle list (PROVED)

Every odd cycle of R is either an old odd cycle entirely in A or B, or uniquely

    γ(p,q) = p followed by q followed by β→α,

where p starts at α in A, q ends at β in B, and |p|+|q| is odd. In particular the latter sum is at least three. Every such pair is a new odd cycle.

**Proof.** A cycle meeting both shores must cross B→A, whose only available edge is β→α. Since the cycle is simple it uses that edge once. Cutting immediately after this edge gives a path from α to β. Before its first crossing A→B it stays in A; after that crossing it cannot return to A without repeating β→α. Thus it has the unique displayed decomposition. Conversely the seam from the end of p to the start of q exists: the only missing forward cross edge would have both paths singleton, but their combined size would then be even. The return edge closes the cycle. ∎

This elementary proof includes possible internal nontransitivity of A and B. The entire new-cycle list is a list of actual labelled objects, not just counts.

## 2. Exact packing polynomial and induced-subgraph interface

Let

    Z_X(t,u) = Σ_F t^{number of cycles in F} u^{number of used vertices in F},

summing over vertex-disjoint collections F of odd cycles, including empty. Define Z_empty=1. The signs used by the odd-cycle formula are not yet included: the number of signed states is Z_X(2,1).

### Theorem RJ2: return-join recursion (PROVED)

For R=R(A,α;B,β),

    Z_R(t,u) = Z_A(t,u) Z_B(t,u)
       + Σ_(p starts α, q ends β, |p|+|q| odd)
           t u^{|p|+|q|} Z_(A−V(p))(t,u) Z_(B−V(q))(t,u).       (RJ)

**Proof.** All new cycles contain both α and β, so a packing contains at most one. If it contains none, it is an arbitrary old packing on each shore. If its new cycle is γ(p,q), its remaining cycles avoid V(p) and V(q). They lie entirely within the shores since another cross cycle would reuse α,β. This gives the displayed bijection and inverse. ∎

The identity holds multivariately with an independent variable x_C for every odd cycle: replace t u^{|γ|} by x_γ and use the inherited variables on each deletion. This records source addresses exactly, unlike a scalar H identity.

**Closed sidecar under nesting.** Keep the same data for all induced subsets. If both α and β survive, an induced subgraph is the return join of the restrictions. If either mark is deleted, the restriction is an ordinal join. Thus (RJ), ordinal multiplication, and the endpoint-path lists form a finite exact recursively closed interface. It may be exponential in the number of vertices; no efficient uniform algorithm is claimed. A scalar H or a handful of unmarked cover coordinates cannot substitute for these residual records.

**Conflict-graph grammar.** Preserve the two old conflict graphs, with no edges between them. Add one new vertex for every eligible (p,q). The new vertices form a clique. Join γ(p,q) to exactly those old cycles meeting V(p) or V(q). This is an explicit physical realization rule for new incompatibility graphs, not an assertion that all abstract graphs are realizable.

## 3. Exact address selection and preserved core

For a compiled module M with complete labelled odd-cycle list C1,…,Ck, let Γ(M) have an edge ij iff Ci and Cj intersect. Replace k labelled triangle reservoir blocks by M, ordinally after an unchanged protected core (T,z). The new signed packing states have the explicit decoder

    empty choice at address i ↔ coordinate 0,
    selected Ci with sign ε ↔ coordinate ε∈{+,−}.

The image is exactly the sign-saturated address set whose support is independent in Γ(M). The decoder is the identity on the core packing and on every untouched shell. It is a bijection, and preserves the all-zero source copy literally on the core. Ordinal direction prevents cycles mixing the core and module. The old triangle-shell vertices are not retained, and the new addresses need not retain their old geometric ports. An arbitrary fixed ± convention is part of the decoder, not a claimed automorphism.

The exact guard in the Boolean quotient x_i²=x_i is

    f_Γ(x)=∏_(ij∈EΓ)(1−x_i x_j).

It is idempotent and its range on Boolean points is exactly the allowed support set. Algebraic idempotence alone does not build a tournament: the return-join recipe supplies the physical realization certificate. For a desired graph isomorphic to an output Γ(M), the labelled isomorphism supplies its address decoder.

**Protected source boundary.** The original core and mark z remain untouched; all old core cycle states occur in the zero-address copy. Any stronger marked-source injection established for ordinal reservoir attachment in Session262 continues to apply with the same ordinal direction. This theorem does not identify the internal module marks α,β with the protected core mark. Further ordinal recharge remains possible.

Thus the grammar extends matching selectors by actual constructors, while explicitly declaring the additional cycle addresses created by overlap.

## 4. Finite one-apex clique bank

For a binary word w=w1…wm, take transitive vertices v1→…→vm and an apex a with a→vi iff wi=1. Every cycle contains a. An odd cycle is uniquely

    a, vi, (an odd or even selected increasing interior subsequence), vj, a

with i<j, wi=1, wj=0 and an even number of selected interior vertices (its length is three plus that number). Hence its exact length enumerator is

    L_w(u)=Σ_(i<j, wi=1,wj=0) Σ_(0≤h≤j−i−1, h even)
                 binom(j−i−1,h) u^{h+3}.                         (A)

Every pair of cycles intersects at a, so Z=1+t L_w(u). Put c(w)=L_w(1). Each inversion at distance1 contributes1; each at distance d≥2 contributes 2^{d−2}. This gives a literal K_c(w) mutex selector with multiplier 1+2c(w), with one address for each actual cycle.

Initial words: 10 gives c=1,H=3; 100 gives c=2,H=5; 1000 gives c=4,H=9; 1100 gives c=5,H=11. In particular K4 and K5 guards are new nonmatching initial selectors in this interface, not merely products of pair gates. The shared-edge fan on m tips gives c=2^{m−1}, recovering its formerly unwanted longer cycles as valid additional addresses. No K3 claim is made; the inherited seven exclusion forbids an exact three-address clique selector.

The formula (A) is an address enumeration proof, not a claim of a new one-apex H recurrence. All those scalar counts are compatible with the prior one-apex work.

## 5. Two one-apex modules and a non-cluster initial selector

Suppose every cycle of A contains α and every cycle of B contains β. Let p,q be their respective odd-cycle counts. Put

    a_o,a_e = counts of odd/even-sized α-start paths in A,
    b_o,b_e = counts of odd/even-sized β-end paths in B,
    s = a_o b_e + a_e b_o.

All residual graphs in the second term of (RJ) are transitive, since the removed paths include the apexes. Consequently

    Z_R(t,1)=(1+pt)(1+qt)+st,
    Γ(R)=(K_p disjoint-union K_q) joined to K_s.                  (B)

The definition of joined means every vertex of the final clique is adjacent to every old vertex. This follows because each old A-cycle contains α and each old B-cycle contains β, while each cross cycle contains both. For p,q,s>0 the resulting conflict graph contains an induced P3, but its other vertices supply the required extra states. The inherited no-three-only-cycle P3 obstruction is respected, not bypassed.

### Literal six-address gate

Take A and B each a directed triangle, with their respective apexes as marks. In vertex order 0,…,5, one literal output has outgoing masks

    (50,60,57,17,32,8).

Its complete cycle list in canonical rotation is

    C1=(0,1,2),         C2=(0,1,2,5,3),
    C3=(0,1,3),         C4=(0,1,4,5,3),
    C5=(0,5,3),         C6=(3,4,5).

Only C1 and C6 are disjoint. Thus Γ=K6 minus edge16, and the allowed signed supports are empty, any singleton, and {1,6}. This is a six-address gate with 1+6·2+4=17 states, rather than the 729 states of six independent triangle shells. The predicate has all pair constraints except16, genuinely outside the matching grammar.

Each triangle has start-path parity counts (a_o,a_e)=(2,1) and end-path counts (b_o,b_e)=(2,1), so s=4, agreeing with the displayed four new cycles. The exact length-sensitive polynomial is

    Z_R(t,u)=1+t(4u³+2u⁵)+t²u⁶.                                (C)

This is also a small concrete explanation of why promoting leakage to addresses can be productive: two disjoint old cycles survive as one two-address support, while four new shared cycles offer mutually exclusive replacement states.

## 6. All-order cover sidecar from the same polynomial

Let Q_r(X) count ordered spanning r-path covers allowing empty slots, as in Session262. For n=|X|,

    Q_r(X)=Σ_F 2^{|F|} r^{n−|V(F)|+|F|}.                       (Q)

**Proof from inherited OCF.** Assign each vertex one of r path slots. For a fixed assignment, the number of directed paths spanning each color class is the product of its induced H values (H(empty)=1). Expand each by the inherited odd-cycle formula. A resulting packing F must have each of its cycles monochromatic; the |F| cycles and the n−|V(F)| unused vertices can receive colors independently. Summing colors gives exactly (Q). No new symmetric-function result is invoked.

Equation (Q) transforms the computed Z sidecar into every Q_r. Nonempty ordered-cover coordinates follow by binomial inversion. For the six-address gate,

    Q_r=r⁶+8r⁴+8r²,
    (H,Q₂,Q₃)=(17,224,1449).

For the word constructor, Q_r=r^{m+1}+2Σ_cycles r^{m+2−length}. Ordinal products multiply Q_r and multiply the independent-set packing polynomials. Thus combinations of these gates retain complete unmarked cover data for downstream kernels, not merely H.

**What is not recovered.** Even the entire polynomial Q_r need not recover cycle labels, intersections, endpoint path lists, or deletion sidecars. Keep those explicitly to iterate return joins. The inherited Session262 full-profile trade has this same 17/224/1449 triple; no equality of marked path/deletion data follows.

## 7. Reproducible finite cases and hostile probes

Run from repository root:

    python3 05-knowledge/results/session263_spin_flip_20261111/selectors/selector_probe.py

The standard-library script constructs literal adjacency masks. It checks all 63 binary words of lengths0…5 against direct odd-cycle enumeration and independent Hamiltonian-path subset DP. It constructs every pair of one-apex words of lengths0…3: 225 return joins, orders2…8. It compares (RJ) with independently enumerated cycle packings and compares (Q) for r=1,2,3 against ordered-cover subset convolution of direct path counts. This is an exact finite constructor test, not a tournament census or a verification project for inherited canon.

It also checks the displayed signed decoder and all cycles of the six-address gate. Three nested return joins on9 vertices retain nontrivial residual cycles, guarding against the false generalization that every new cycle leaves a transitive residual. Their exact H values are149,97,215; full masks and packing data are in `selector_receipt.json`. Singleton shore/path cases test the missing α→β forward edge: its sole problematic seam has even total length and never enters (RJ).

The examples distinguish the two live mechanisms: apex pairs permit formula (B); arbitrary return joins require the deletion sidecars in (RJ). Applying (B) to arbitrary marked tournaments is not justified.

## 8. Source/novelty ledger and next gate

- **Inherited/cited internally:** odd-cycle formula (THM-002, `01-canon/theorems/THM-002-ocf.md`); Session262 marked ordinal source/diamond/matching interface; Session26 one-apex family; Session261/262 profile17/224/1449. No re-audit or priority claim.
- **New in this packet:** complete literal return-join cycle decoder; multivariate residual packing recursion; recursively closed induced-subset interface; the explicit conflict-graph grammar; exact address selectors K4,K5 and K6 minus one edge; length-sensitive sidecars attached to those decoders. Formula (Q) is an elementary derived sidecar, not claimed as new literature.
- **Snippet inspiration only:** overlapping/disjoint odd cycles and algebraic idempotents suggested the native cycle-conflict carrier. No graph-minor/Petersen theorem is imported or asserted to transfer to tournaments.
- **Destroyed information:** old shell vertex incidence; any claim that a scalar multiplier determines the requested source addresses.
- **Preserved:** exact labelled source support through the decoder; all signs; zero source; unchanged external core; exact cover sidecars.
- **Cheapest decisive test:** enumerate actual cross cycles and compare them to endpoint-path pairs before searching a target predicate.

**Actual gain toward Chan–Pak:** the physical selector library now includes correlated supports with overlapping forbidden pairs, including an induced-P3 conflict geometry with fully controlled extra addresses. Its compatible doubleton yields a quadratic packing term rather than a pure clique sum. Nested constructors give further exact gates without another guess about longer-cycle leakage. This is a supply/compiler advance, not a target-coverage theorem: selecting an arbitrary desired number of states, synchronizing an infinite recursive source, and closing Chan–Pak remain OPEN.
