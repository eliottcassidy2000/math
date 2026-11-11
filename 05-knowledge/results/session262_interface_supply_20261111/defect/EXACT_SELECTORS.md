# Exact support selectors, their finite guards, and cover-profile blindness

**PROVED:** an exact classification of the reservoir address sets selectable by induced deletion; an all-order cover-profile transport formula; and a precise indistinguishability theorem for profile-only constructors. **FINITE-EXACT:** small selector and profile controls. **PROVED in §7:** a legal two-address correlation gate retaining the reusable core source. **OPEN:** arbitrary correlated selectors and arbitrary targets. This is a new constructive wave, not a re-audit of the inherited kernel or source modules.

## 1. The operation really can select addresses, within an exact class

Keep a nonempty core `(T,z)` and its protected fibre injection from `PROOF.md`. Its k-shell reservoir has the labelled packing product

    N_W(T,z) × Ω_k,     Ω_k={0,+,−}^k.

For a word w, let supp(w) be its nonzero coordinates. For J⊆[k], physically delete r_i, the designated vertex of triangle i, for every i outside J. The other two vertices of those triangles remain. Call the resulting tournament T_J.

### Theorem D: exact principal-ideal selector (PROVED)

The natural inclusion of odd-cycle packing states of T_J into those of the full reservoir has image exactly

    N_W(T,z) × P_J,
    P_J={w:supp(w)⊆J}.

In particular it preserves the distinguished all-zero source copy and has `|P_J|=3^|J|`. The selected tournament has the original core unchanged and still admits a new ordinal recharge shell.

**Proof.** An intact triangle has its empty packing and two labelled cyclic states. A triangle with its designated vertex deleted is a directed edge; its only odd-cycle packing is empty. Ordinal block directions prohibit any cycle spanning blocks. Inclusion therefore acts identically on old states and sets exactly the deleted-shell coordinates to zero. ∎

This is an **exact source-address selector**, not a mere comparison of counts. It costs only explicitly specified vertex deletions and retains all protected source addresses. A deletion is a mathematical construction, not an instruction to alter any user's unrelated data.

### Converse, including arbitrary shell-vertex deletions

Allow deletion of any collection of shell vertices, but no core vertices and no arc changes. A shell contributes its nonzero cycle states iff all three vertices remain. If zero, one or two remain, it contributes only the empty state. Thus every such induced deletion has image `N_W×P_J`, where J is the set of completely retained triangles. Conversely each P_J is realized by the designated-vertex operation. Hence these and only these selectors are available in this operation class.

The selector semigroup is exact:

    P_J ∩ P_K = P_(J∩K),
    delete(J) followed by delete(K) = delete(J∩K)

at the state-image level. Recharging appends an unconstrained coordinate: `P_J×{0,+,−}`. Unselected coordinates can instead be appended as directed edges, giving `P_J×{0}`.

## 2. A finite, complete certification interface

Given a desired set X⊆Ω_k, a core-independent protected selector exists by induced shell deletion iff X=P_J for some J. There is a finite certificate without any tournament search:

1. **Sign saturation:** membership depends only on supp(w), not its ± choices.
2. Let F⊆2^[k] be the accepted supports. Require F nonempty and downward closed.
3. Require F closed under unions.

For finite F, these conditions imply that `J=⋃F` belongs to F and `F=2^J`; conversely a principal ideal satisfies them. The certificate outputs the explicit deleted vertices `{r_i:i∉J}`. Equivalently, a downward-closed support family must have exactly one maximal member.

This interfaces directly with a finite initial-bank compiler: a proposed finite support rule can be accepted with a physical selector, or rejected with a short witness. A rejection witness is a pair of words of equal support with unequal membership; a missing subset of an accepted support; or two accepted supports whose union is rejected. The empty requested set is also impossible: the empty packing survives every induced deletion.

In Boolean variables x_i indicating nonzero coordinates, P_J has the exact membership polynomial

    f_J(x)=∏_(i∉J)(1−x_i)

in the quotient by `x_i²−x_i`. A proposed Boolean polynomial guard f describes a lawful selector exactly when its reduced multilinear representative equals one f_J. These statements are finite identities in the Boolean quotient; they do not require an infinite vanishing-ideal decision procedure. With an explicitly supplied truth table, the finite test above is complete. With a compact guard representation, finding J can still be computationally expensive; no polynomial-time claim is made.

**First correlated hostile.** For k=2, request `supp(w)∈{∅,{1},{2}}`, meaning “use at most one shell.” This is sign-saturated and downward closed but not union-closed. It contains five states. The realizable selector sizes are 1,3,9, so it cannot be realized by any induced shell deletion. Its predicate is `1−x₁x₂`, not one of `1`, `1−x₁`, `1−x₂`, `(1−x₁)(1−x₂)`. The union-failure certificate is `{1},{2}`.

Count five is not asserted absent from the tournament spectrum. What is impossible is retaining **these five particular inherited states**, with this protected core, by the specified operation.

**Second hostile.** “Use shell1 with its + label” violates sign saturation and excludes the zero word. It is impossible even if its cardinality happens to equal that of an unrelated realizable tournament. Labelled cycle bits are counting coordinates, not arcs that can be independently removed.

Thus capacity and target selection already separate on two shells, with a finite certificate rather than an open-ended failed search.

## 3. Complete ordered-cover state of the stabilized bank

For a tournament X and integer r≥1, define Q_r(X) as the number of ordered r-tuples of vertex-disjoint directed paths whose vertex sets partition V(X), allowing empty paths. Set Q_r(empty)=1. Let C_j(X) count ordered spanning j-path covers with every path nonempty. For nonempty X,

    Q_r(X)=Σ_(j=1)^r binom(r,j) C_j(X).

In particular `Q₁=H`, `Q₂=2H+C₂`, and `Q₃=3H+3C₂+C₃`.

### Theorem E: all-order cover transport (PROVED)

For ordinally joined tournaments X→Y,

    Q_r(X→Y)=Q_r(X)Q_r(Y).

**Actual map.** Given ordered covers `(p₁,…,p_r)` of X and `(q₁,…,q_r)` of Y, concatenate p_i q_i in each slot, omitting empty pieces. Every seam points forward. Conversely restrict each output path to X and Y. No output path can move from Y back into X, so these restrictions are directed paths and the maps are inverse. This includes empty slots and empty blocks.

For a directed edge E and directed triangle C,

    Q_r(E)=r²,
    Q_r(C)=r³+2r.

Indeed C has `C₁=3,C₂=6,C₃=6`, while E has `C₁=1,C₂=2`; substitution gives the formulas. These are explicit initialized modules for *every* finite cover hierarchy, not just H and Q₂.

Consequently, for the designated-vertex selector T_J on k shells, writing j=|J|,

    Q_r(T_J)=Q_r(T) (r³+2r)^j (r²)^(k−j).             (E)

In particular,

    H(T_J)=3^j H(T),
    Q₂(T_J)=3^j 4^k Q₂(T),
    C₂(T_J)=3^j[4^k Q₂(T)−2H(T)],
    C₃(T_J)=33^j 9^(k−j)Q₃(T)
             −3·3^j4^k Q₂(T)+3·3^jH(T).

This supplies the missing higher coordinate for this reservoir family if a constructor's U₀ requires C₃. It does not assert that C₃ is redundant for general modules. Notice the r=1,2 ratios of an intact versus a deleted triangle are both3; the r=3 ratio is11/3. Thus H and Q₂ scale together under activation at fixed k, whereas C₃ need not scale with them.

## 4. Exactly what the other lanes can see

### Theorem F: profile-only address blindness (PROVED)

For fixed T and k, all designated-vertex selectors T_J with the same |J| have identical Q_r for **every r≥1**, hence identical C_r for every r by binomial inversion. Any constructor whose output observables are fixed functions of finitely many unmarked ordered-cover coordinates cannot distinguish these selectors. This applies to a polynomial kernel once all of its required inputs are among those coordinates.

The proof is formula (E), not a finite coincidence. The operation can distinguish the size of J, but these statistics discard *which* shell addresses were kept. At k=2, the actual patterns `triangle→edge` and `edge→triangle` are different selected subgraphs; their complete ordered-cover profiles agree by commutativity of the product.

For the inherited bilinear K(H,C₂), substitution of T_J in one module slot at fixed k multiplies both input coordinates by3^j relative to T_∅. Therefore it multiplies that bilinear output by3^j when the other module is fixed. For the stated zero-tail anchor `E₀=7 Q₂(A)H(B)`, selector activation in its A or B slot has the same scalar effect. These statements use only the inherited formulas' declared dependence, not any re-verification of their kernel coefficients.

A constructor depending on C₃ can see additional dependence on j, exactly supplied above. Nevertheless it still cannot see the identity of J. Higher unmarked cover orders never repair that address loss in this family. A marked incidence statistic, or a constructor that couples specified shell positions, is necessary to do so.

**Guard on applicability.** The theorem concerns the designated one-vertex deletion menu, with two vertices retained in each inactive shell. Arbitrary deletion of additional shell vertices changes Q_r through factors1 or r as well as r²; it still selects the same principal support ideals in OCF coordinates, but the unmarked cover profile can then also detect numbers of residual blocks of each size. Marked or endpoint-restricted kernels are outside the profile-only theorem unless their required sidecars are separately shown to factor. No conclusion about an arbitrary physical graft is inferred just from its H output.

## 5. Relation to paid sources and the exact remaining gate

The other lane's paid sources `(D,U,V)=(289^m,3808^m,3808^m)` with independently interchangeable `(H,Q₂)=(17,224)` modules give another genuine source bank. A kernel declared to depend only on H,Q₂ cannot select between those alternatives using its scalar count. This is an immediate conditional consequence of that declared interface, not a claim that their marks or higher cover profiles coincide. If a downstream U₀ uses C₃, those interchangeable modules need a C₃ sidecar before that downstream interchange is justified. Our formula (E) supplies such sidecars for triangle/edge reservoirs, but not automatically for the paid-source modules.

The finite guard compiler can now carry a particularly sharp certificate:

    source = protected N_W×Ω_k;
    desired address predicate = f;
    physical selector class = induced shell deletion;
    success witness = J with f=∏_(i∉J)(1−x_i);
    rejection witness = sign/downward/union violation;
    output anchor = H(T_J)=3^|J|H(T);
    retained recursive data = core labels, block labels, and Q_r(T_J)
                              for every cover order the next kernel uses.

This is an actual finite certification interface with explicit initial objects and output tournaments. It is not another capacity bound.

**Precise stopping boundary.** The two-shell predicate `1−x₁x₂` is the first useful target for a new operation: admit either triangle alone while forbidding their simultaneous use, retain the all-zero protected source, and return a reusable port. Induced deletion cannot do it; a profile-only constructor cannot identify which individual source address it is spending; and interpreting its five selected states as an existing tournament is not a construction. A richer marked module can cross this boundary: the next creative step in §7 constructs one. The impossibility proved here remains exact for induced deletion, but is not an impossibility for module replacement.

## 6. Finite-exact controls

Run `python3 05-knowledge/results/session262_interface_supply_20261111/defect/selector_probe.py`. The script directly enumerates labelled words and all per-shell residual-size deletion patterns for k≤4, checks every sign-saturated support predicate for k≤3 against the complete certification test, and counts ordered covers of the initialized edge and triangle for r=1..5 by literal permutations and assignments. It also compares all designated two-shell deletion profiles by an independent subset-convolution cover counter. The extension also directly enumerates every odd cycle of the four-vertex gate and shared-edge fans through five tips, checks its literal ordered covers for r≤5 and a two-gate product, and checks the explicit five-word decoder. Receipt: `selector_receipt.json`.

## 7. Crossing the boundary: an actual two-address mutex gate

The induced-deletion obstruction suggests making the two optional cycles share vertices, rather than deleting one of them. This works exactly for two addresses.

Let M have vertices a,b,c,d and arcs

    a→b, b→c, c→a, b→d, d→a, c→d.

Its outgoing masks in order a,b,c,d are **(2,12,9,1)**. Its only directed triangles are abc and abd: triples acd and bcd are transitive. A four-vertex tournament has no longer odd cycle, and the two triangles intersect. Thus its complete odd-cycle packing set consists of

    empty, (abc,+), (abc,−), (abd,+), (abd,−).

### Theorem G: exact correlated pair selection (PROVED)

Replace two consecutive triangle-shell blocks by M, ordinally joined to the unchanged core and the other shell blocks. There is an explicit bijection from the new packing states onto the old reservoir states satisfying

    not(x_i=x_j=1),

where i,j are the replaced address coordinates. Send the empty packing of M to `(0,0)`, the two labelled abc cycles to `(±,0)`, and the two labelled abd cycles to `(0,±)`. Take the identity on all core and other-shell packing states. Ordinal block directions prevent cross-block cycles, so this is a complete bijection, not merely an injection or count coincidence.

The selected output has H multiplier5 in place of9. All core labels, z, and the protected zero-word source copy survive. Its cut graph at z still differs from G only by forced isolated vertices. The output remains available for further ordinal recharge. This is a lawful physical graph replacement selecting a correlated address family that no induced deletion could select.

**What is sacrificed.** The six old shell vertices are not retained, and the old shell marks are not preserved as independently deletable triangle ports. Both output cycles share a,b. The operation preserves the **abstract source addresses through the displayed decoder**, not their old vertex incidences. The allowed interface therefore consists of the protected core and labelled packing addresses, not a demand to retain all original shell marks literally. The gate comes with a complete new cover-profile sidecar below.

### Theorem H: a finite matching-selector compiler (PROVED)

Take k labelled shell addresses. Specify some forced-zero coordinates and a matching of forbidden simultaneous pairs among the remaining coordinates. Let p be the number of pairs and s the number of unpaired, unforced coordinates. Replace each paired pair by M, retain a triangle for each free coordinate, and retain a directed edge for each forced-zero coordinate. Put the resulting blocks in a fixed labelled ordinal order after the unchanged core.

Its packing states are exactly the original source states whose supports avoid every forbidden coordinate and every matching edge. The multiplier is exactly

    h=5^p 3^s.

For a fixed core, the whole new interface is

    H_out=hH(T), a_out=ha, d_out=d,
    Δ_out=hΔ, σ_out=hσ+(h−1)(d−1).

The guards have the explicit Boolean form

    f(x)=∏_(i forced)(1−x_i) ∏_({i,j} matching)(1−x_ix_j).

This upgrades the finite certificate compiler from principal ideals to a genuinely correlated class. For a given sign-saturated downward-closed support family, compute its inclusion-minimal forbidden supports. It belongs to this class iff they all have size1 or2 and the size2 ones form a matching. The disjointness condition is load-bearing: the proof is an ordinal product of separate physical gates.

The initial bank is finite and explicit: directed edge E, directed triangle C, and the four-vertex gate M. The decoder proves exact address selection at every product size. Independent gauge choices for the ± labels must remain part of the certificate; they are not asserted to be tournament symmetries.

### Complete cover-profile sidecar for the new gate

For M the ordered nonempty cover counts are

    C₁=5, C₂=22, C₃=36, C₄=24.

The first follows from the five packing states and the inherited OCF identity. For C₂, a singleton/triple split contributes `2(3+3+1+1)=16`, and the six ordered two/two splits contribute6. For C₃, choose the doubleton in6 ways, its path slot in3 ways, and order the two singleton slots in2 ways, giving36. Four singleton paths give24. Hence for all r≥1,

    Q_r(M)=5r+22 binom(r,2)+36 binom(r,3)+24 binom(r,4)
          =r⁴+4r².

For p pair gates, s free triangle blocks and t forced-zero edge blocks,

    Q_r(output)=Q_r(T)(r⁴+4r²)^p(r³+2r)^s(r²)^t.

In particular the replacement's H,Q₂,Q₃ are **5,32,117**, not the **9,144,1089** of two independent triangles. A downstream K(H,C₂) or U₀(C₃) must read the new sidecar, not reuse the old one. This supplies every finite ordered-cover coordinate the stated unmarked kernels might require, by binomial inversion.

### The first stronger boundary is sharp in a different way

A hypothetical standalone three-address mutex module preserving exactly empty and one of three signed cycle states would have H=7. The inherited global H=7 exclusion (THM-1370, cited by the sampled `THM-5107-physical-graft-content-excludes-cofactor-seven.md`) rules out such a module. This is a genuine spectral obstruction to that precise standalone seven-state selector, not merely failure of the present construction. It does not rule out a core-dependent construction with total7H(T) when that product is allowed, nor a selector that changes the state semantics.

The naive shared-edge generalization explains why the two-address gate stops. Add ordered vertices c₁→…→c_m with b→c_i→a and a→b. Every odd nonempty subset of the c_i gives an odd cycle through a,b in its unique increasing order. These are all cycles, since deleting a or b leaves a transitive tournament. All cycles intersect, so H is

    1+2 Σ_(j odd) binom(m,j)=1+2^m  (m≥1).

For m=2 this is the desired5. For m=3 it is9 rather than7: the additional five-cycle supplies exactly the unwanted two states. Thus overlapping shared-edge mutexes introduce a concrete higher-cycle leakage, not a mysterious numeric failure.

The new stopping problem is now specific: compile a nonmatching conflict family with a marked decoder and full cover sidecar, while controlling these longer-cycle states. The gate above already resolves the first correlated selector and expands exact multipliers from powers of3 to the physical semigroup `3^s5^p`.

**Inherited module cross-link.** The coordinator identifies the same four-vertex diamond as the Session261 nine-grid s2 source, with singleton selected mark2 and U=9. That marked-interface identity is inherited here, not re-derived or used in the selector proof. It offers an additional route from this new cycle-overlap decoder to the existing offset constructors; the selector claims themselves retain the original core mark and use the complete unmarked cover sidecar above.

### A second sharp obstruction: the three-address path of conflicts

Assume a primitive-cycle selector has exactly three distinguished odd cycles C₁,C₂,C₃, with C₁ and C₃ vertex-disjoint, C₂ meeting each, and no additional odd cycles. Its proposed independent-set state count is

    1 + 3·2 + 2² = 11.

But the union of intersecting directed cycles is strongly connected, so the induced tournament on `V(C₁)∪V(C₂)∪V(C₃)` is strong. The two disjoint odd cycles force at least six vertices. The inherited strong floor `H≥15` for every strong tournament of order≥6, as stated in **THM-1370, `01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md`, proof steps3–4**, contradicts H=11 for this induced tournament. Therefore this exact primitive-cycle model necessarily acquires additional odd-cycle packing states, totaling at least four beyond the proposed eleven.

This is stronger than the three-way mutex's scalar-seven obstruction: eleven is an allowed scalar tournament count, but it cannot realize **this primitive cycle-intersection geometry without extra cycles**. It does not prohibit tournaments containing three cycles whose intersection graph is P₃, does not prohibit an unrelated H=11 tournament, and does not rule out a different abstract decoder for eleven states. No census or re-verification of the inherited strong floor is used.
