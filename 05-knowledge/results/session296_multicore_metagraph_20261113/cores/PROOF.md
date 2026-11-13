# Exact endpoint-ready handoffs past the mod-45 wall

**PROVED finite entrance ladder, with actual named interfaces.** The
all-clock continuation and Chan–Pak remain **OPEN**.

## 1. What barrier is removed, and at what scope

The inherited [entry resource theorem](../../entry_resources_20261113/MACHINES.md#M5)
and [uniform return](../../cycle_opening_six_return_20261112/compiler/MODULE_RETURN.md)
have source (A_k=45·2^k+7) and ready return (B_k=135·2^k+23). Both
have an indegree-one root with defect (1); the latter is **not** another
member of the former family because their residues modulo (45) differ.
The next demand is (C_k=405·2^k+71), which is 26 modulo 45.

We obtain physical strong defect-one sources at **every (C_k) for
1 ≤ k ≤ 9**. These are literal marked witnesses from several
different recursive core addresses, not arithmetic solutions standing in
for a tournament. For (k=1) there is a longer chain of ready handoffs:

    97 → 293 → 881 → 2645 → 7937 → 23813 → 71441.

Every arrow realizes (3h+2) by two indegree-one births followed by a
same-count named replacement (the first arrow uses the inherited uniform
return). The same old root and predecessor can survive the replacement
as physical vertices. No ready return at the next target 214325 is
claimed. These are not infinite rays, and nine entry
clocks do not prove the entry formula at any other clock.

## 2. Actual recursively composed cores

The construction is the inherited [one-start funnel](../../polarized_path_calculus_20261103/CALCULUS.md)
and [six-coordinate compiler](../../berggren_target_return_20261111/reset/funnel_six_state.py).
Its use here supplies targets outside the preceding fixed charts. Start
with a singleton, marked at its only vertex. For each binary tail word
(w), form (F_w(T,r)) with a transitive tail on the letters of (w),
uniformly beaten by the old core (T). A new vertex (q) beats only
(r) among old core vertices, and beats precisely the zero-letters of
the tail. Mark (q). Each word is empty or ends in 1. The singleton
uses a nonempty word in all our addresses.

Every marked funnel root has one start: a Hamilton path starting there
cannot enter the tail first (there is then no way to visit the old
core). It enters the old root, then takes its unique full old-core
path before the transitive tail. The new rooted core is strong: from
the old core one can reach the new root, either through an old
non-port vertex or through the terminal tail vertex; the new root
enters the old root and the old core reaches the entire tail. All
later empty-tail cases have non-singleton strong old core. Finally
add (z) beaten only by (q) and mark (z). The strongness persists
and (z) has defect (H(T-z)-s_z(T)=1), since a Hamilton path starting
at (z) can be followed by every old Hamilton path **except** those
starting at (q).

| Desired ready H | Funnel addresses, chronological (ε = empty) | Order of ready source |
|---:|---|---:|
|293|1, 10001|10|
|881|11, 1, 11|10|
|2645|111, 11, 1|11|
|7937|01, 01111101|14|
|23813|1011, ε, ε, ε, ε, 1|13|
|71441|000011010110001|18|
|1691|1011, 101|11|
|3311|1, ε, 01, ε, 1|11|
|6551|101, ε, 01101|13|
|13031|010101111011|15|
|25991|011, ε, 0011001|15|
|51911|1, ε, 11, ε, ε, 1, ε, ε|14|
|103751|1, 001111, 1, ε|14|
|207431|11, ε, ε, 00011, 1|15|

The mark and predecessor are always the final two names in construction
order. The (293) row is the inherited (Q(TT_1)), with exactly its
construction order. This table is checked both via the exact funnel
compiler and via independent literal subset path DP in
[`verify_ready.py`](verify_ready.py). The literal formula
`D[S,j] = sum D[S−{j},i]` over predecessors i in S−{j} with singleton
base cases counts every path once by its last vertex. Running it on
the full set, the root-deleted set, and with a prescribed initial root
checks both the absolute count and the defect.

## 3. Same-target physical replacements

Let (T) be a displayed ready source at (h), and (R) the raw
two-birth extension with newest names (x,y) (sole predecessor of (y)
is (x)). The inherited two-birth identity gives (H(R)=3h+2).
At every declared target select its displayed ready (Q), with
predecessor (q) and root (z). Keep (x=q) and (y=z) physically.
Of the remaining labels, match the first (min(|R|,|Q|)-2) in
construction order; delete or add unmatched vertices and orient all
comparisons as in (Q). Flip exactly the mismatched matched-pair arcs
recorded in the receipt. No retained arc incident to the root flips:
its only incoming arc is from its matched predecessor, and every
other old retained vertex loses to it. Every new vertex also loses to
it. Both endpoint graphs are strong, their whole counts coincide,
and the resulting named root has defect (1), proving it really is
an eligible input to the next two births. This is a macro, not a
claim of count preservation at the intermediate single edits, nor
of optimal old-arc retention.

For k = 2, …, 9, the exact inherited (B_k) graph is the address
(1,1^k0001) followed by the same terminal birth. Therefore the same
explicit relabel/replace recipe certifies all nine transfers from
the inherited family to the new finite ladder. Counting these larger
raw graphs is not needed for the proof: the all-height inherited
birth identity fixes their scalar, and the incidence recipe fixes their
root star. For (k=1), the verifier independently recomputes raw
full counts at all five handoffs of the displayed chain after (293).

## 4. Complete first-target ready fibre, not a higher-clock exclusion

The inherited [complete ready inverse](../../ready_return_pentagon_20261103/OBSTRUCTION.md)
shows that any strong defect-one source at (H=881) is an indegree-one
birth from a strong one-start parent with count at most (881+1)/2 = 441.
Busch's strong-order bound (f(13)=625) forces parent order at most
12. The inherited one-start grammar is recursive over all sizes,
including exactly the allowed words of §2; the parent order profile
through 12 is (1,3,8,21,55,144,377,987,2584,6765) (orders 3–12).
The complete deterministic inversion finds precisely the two addresses

    (11,1,11)          H(parent)=75,  E_parent(2)=807,
    (1,ε,00101)        H(parent)=219, E_parent(2)=663.

The child count is (H(parent)+E_parent(2)-1), so both equal (881).
They have distinct orders (10) and (11), hence are not isomorphic.
The source graph determines the rooted parent by deleting its root
and marking its sole predecessor, so the complete parent inversion
is also complete for the rooted strong children. The weak ready
sources are precisely these strong sources followed by a transitive
tail, by the cited SCC theorem.

In particular an empty ready fibre at the very first new target cannot
be the reason this iteration has resisted a proof. Neither readiness
nor a construction for finitely many addresses, however, supplies
the missing **uniform** entrance rule. The adjacent finite-atlas
obstruction suggests variable/non-dyadic clocks as the next place
to look; the literal ladder is a collection of test cases for such
a rule, not a substitute for one.
