# Strong-subtournament endpoint floors block induced rearming

Session 194, 2026-11-05. **PROVED elementary theorem, independently reviewed, with FINITE-EXACT controls.** The unrestricted Hamilton-path spectrum outside `{7,21}` remains **OPEN**.

## 1. Question and inheritance

The endpoint-matched fan in [the published ternary-fan proof](https://github.com/eliottcassidy2000/math/blob/801be8dc61a35076678363c2ead0d6d835b57ee9/session194/ternary-fans/PROOF.md) gives genuine induced one-vertex children with counts `3N+2` and `3N+4`. Its local audit proves that the literal parent and children lack an endpoint-one mark. The question here is whether swapping the varied join operand, adding a few vertices, or rearranging exterior interfaces can repair that missing mark.

The closest mechanism is the elementary injection of Hamilton words from an induced tournament into a larger tournament. The necessary sidecar is **reachability of the prescribed final endpoint**, which allows the insertion to preserve that endpoint. The hostile is a nonstrong source with previously impossible endpoints. The connection is from an induced strong bank to a lower bound on every positive endpoint count in every ambient tournament; no module hypothesis is needed.

The live concepts are prescribed-endpoint insertion, an exterior directed bridge, induced strong cores, endpoint-matched join coefficients, and an interface-return obstruction. No literature novelty claim is made. The theorem is proved directly below.

## 2. Prescribed-endpoint insertion

Write `e_x(T)` for the number of directed Hamiltonian vertex words of T ending at x.

**Lemma 1.** Suppose every vertex of T can reach x. Let W be any nonempty subset containing x. Every Hamilton word of T[W] ending at x extends, injectively, to a Hamilton word of T ending at x.

**Proof.** While W is proper, some exterior vertex z has an outgoing arc into W. To see this, follow a directed path from any exterior vertex to x and take the last vertex before its first entrance into W. Choose such z deterministically from T and the current set W, independently of the Hamilton word.

For a current word `p_1,...,p_k=x`, let p_j be its first vertex beaten by z. If j=1, prepend z. Otherwise p_(j-1) beats z by minimality of j, so insert z between p_(j-1) and p_j. This is legal and leaves x last, even when j=k. Replace W by W union {z} and repeat.

The chosen sequence of new vertex names depends only on the growing sets, not on the word. Deleting the inserted names recovers the original word, proving injectivity. QED.

In particular, if `e_x(T)>0`, every vertex reaches x along a Hamilton word, and hence

    e_x(T[W]) <= e_x(T)       whenever x belongs to W.

The positivity/reachability condition cannot be omitted: adjoining a universal sink makes every old endpoint count zero.

Applying the lemma starting from W={x} also proves that every vertex of a strong tournament can be a Hamilton endpoint. This uses no Hamilton-cycle theorem.

## 3. All-order induced-core floor

**Theorem 2.** Let S be an induced nonempty strong subtournament of an arbitrary tournament T, and put

    q = min_{s in V(S)} e_s(S) >= 1.

Then, for every vertex x of T,

    e_x(T)=0  or  e_x(T)>=q.                         (1)

If T is strong, the zero alternative never occurs. By reversing every arc, the analogous statement holds for counts of Hamilton words starting at a vertex, using the minimum starting count of S.

**Proof.** Fix an x for which `e_x(T)>0`. If x belongs to S, Lemma 1 applies directly to at least q Hamilton S-words ending at x.

If x is outside S, every vertex of S reaches x. Choose a simple directed path from a vertex of S to x. Let s be the **last S vertex on that path**, and retain the suffix starting at s. This is an actual directed bridge

    s -> z_1 -> ... -> z_k=x,

whose vertices after s are all outside S. We do not delete interior S vertices and pretend that the remaining vertices are adjacent: choosing the suffix is essential. Fix this bridge using only T,S,x, before considering any Hamilton word.

Append `z_1,...,z_k` to each Hamilton S-word ending at s. There are at least q such words. The appended paths are distinct, end at x, and span the same set W consisting of S and the bridge's exterior vertices. Lemma 1 then inserts every remaining vertex while preserving x. Deleting all vertices outside S recovers the original Hamilton S-word, so the entire construction is injective. This proves (1). QED.

An alternative proof handles exterior endpoints by induction using `e_x(T)=sum_{u->x} e_u(T-x)`. The explicit bridge proof above shows the actual injection and avoids any ambiguity about endpoints outside S.

The same conclusion with `q=min_all e_s(S)` is formally true for nonstrong S, but q is then zero and the statement is vacuous. Replacing this q by the minimum *positive* endpoint count of a nonstrong S is false.

**Sharpness and zeros.** Adjoining a universal source preserves every old endpoint count and adds one zero endpoint count, so q is attained in a proper extension whenever it was attained in S. Adjoining a universal sink gives zero at every old vertex and H(S) at the new sink. Thus the positive-count qualification is necessary.

## 4. Exact hostile to dropping strongness

Take

    S=[10,12,9,0],

which is a directed triangle on 0,1,2 followed by a sink 3. Its endpoint vector is

    [0,0,0,3].

Add vertex 4 with sole outgoing old arc `4->0`. The induced extension is

    T=[10,28,25,16,1].

It is strong and has endpoint vector `[1,1,1,3,3]`. Its unique Hamilton word ending at 0 is

    (1,2,3,4,0).

Thus the source minimum positive endpoint count is 3, while the extension has endpoint count 1. A previously impossible old endpoint becomes possible. This does not contradict Lemma 1, which only supplies the lower bound `e_0(S)=0` at that root. It also explains why the strong-core hypothesis in Theorem 2 cannot be replaced by a min-positive convention.

The independently reviewed relabelled hostile supplied by the construction worker is S'=[6,12,0,5], with endpoints `[0,0,3,0]`, and T'=[22,28,16,5,8], with endpoints `[1,1,3,1,3]`.

## 5. The ternary banks and arbitrary ambient surgery

The three special banks, rooted at b=1, are

    B   =[2,4,97,39,15,83,27],
    B_2 =[2,4,225,39,15,83,155,59],
    B_4 =[2,4,225,39,15,211,155,27].

They are strong. Their complete endpoint and starting vectors are:

    B ends   =[21,15,5,5,4,3,4]
    B starts =[3,3,3,11,13,11,13]

    B_2 ends   =[64,43,9,15,15,9,7,9]
    B_2 starts =[7,7,7,22,27,22,38,41]

    B_4 ends   =[65,41,9,13,16,7,8,12]
    B_4 starts =[7,7,7,21,23,33,43,30].

Consequently any tournament containing an induced copy of B **or its converse** has every positive endpoint count and every positive starting count at least 3. Containing B_2 or B_4, in either orientation, raises the bound to 7.

This applies to **arbitrary ambient tournaments**, not only a sequence of vertex births. Adding any number of vertices, changing arcs outside the retained bank, changing arcs between the bank and its exterior, deleting exterior vertices, changing roots, or moving the bank between operands cannot create an endpoint-one or start-one mark while an induced copy of the bank survives.

There is an even smaller retained obstruction: all three banks contain the induced regular five-tournament

    U=[24,9,3,20,6],

on the old labels 2,...,6. Every endpoint and starting count of U equals 3. Thus destroying the full seven/eight-vertex bank but leaving this five-vertex core unchanged still cannot restore a unit endpoint. Any successful repair must break every such induced strong obstruction, for example by changing or deleting vertices/arcs inside this retained core, or reconstructing a different witness. Those repairs are not proved to exist here.

## 6. Swapping join operands does not repair the fan

Recall that J(A,B) takes `A->B^op` and reverses only the cross arc from a to b. Its exact count is

    H(J(A,B))=h_A h_B-e_a e_b+R_a R_b.

Swapping the operands gives the converse graph, after the evident block relabelling:

    J(B,A) is isomorphic to J(A,B)^op.

The formula is symmetric. Replacing either operand's triple `(h,e,R)` by `(3h,3e-d,3R)` changes the joined count from N to

    3N+d*e_opposite,   d=2 or 4.

Therefore the same endpoint-matched family needs an endpoint-one mark on the opposite operand in either arrangement. The previous output has no such mark, and no induced enlargement retaining its bank can supply one. At a positive mark of an enlarged special child, the two increments would be at least 14 and 28 rather than 2 and 4; a zero mark gives two zero increments. Using an unconversed join presentation merely trades ending counts for starting counts; the dual floor excludes that apparent escape too.

This excludes rearming the endpoint-one interface, not producing literal 3N+2 or 3N+4 values: the existing fan already produces those values. This is an exact blocker for these interface-preserving proposals. It is not a proof against different coefficient changes, new physical banks, noninduced arc repair, deleting part of the retained obstruction, or scalar reconstruction with a different tournament. Nor does returning an endpoint-one mark alone automatically prove scalar ternary iteration; the correct marked join presentation would still have to be retained.

## 7. Strongness of the complete joins

For these strong banks of order at least two, J(A,B) is strong **if and only if** the marked vertex a reaches every vertex of A.

Every A vertex can reach B^op directly, since all cross arcs except one point toward B^op and the bank has at least two vertices. The bank is strong and can reach a through b->a. Conversely, b->a is the only cross arc into A, so every A vertex reached from the bank must be reachable from a inside A. This proves the equivalence.

The canonical one-apex T(v) used in the fan satisfies this condition for every v: when v>0 its first bit is one, so the apex reaches the first chain vertex and hence the entire chain; v=0 is a singleton. In particular, the even-v parent and both children are all strong and every endpoint/start count obeys the positive bound 3 or 7, with no zeros.

For a general rooted A with e_a=1, strongness need not hold. For example, A can be a transitive pair rooted at its sink. The join then has an initial singleton component and a strong terminal component containing the bank. Zero endpoint counts outside the terminal component are permitted and do not provide a unit interface.

## 8. Verification and independent review

`verify.py` uses terminal-removal dynamic programming for all endpoint counts and independent directed reachability for strongness. Its universe is every one of the 33,867 labelled tournaments of orders one through six, every nonempty induced subset, and every final endpoint choice. It checks both Lemma 1's inequality and Theorem 2, including proper induced strong cores with q>1. Explicit endpoint-preserving and exterior-bridge injections are additionally constructed through order four. Further bridge injections on each special bank include both a new source and a new sink, testing nontrivial floors 3 and 7 with exterior vertices left to insert. The precise join strongness criterion is tested on all 29 rooted A through order three against all three banks (87 joins). The counterexamples, sharpness/zero controls and all bank vectors are checked exactly.

Run `python3 verify.py` and `python3 -O verify.py`. Gates use explicit exceptions, not assertions. The all-order results follow from the proofs, not this finite universe. The construction worker independently approved the original induction proof and nonstrong hostile; the reusable-ternary worker independently reviewed the direct bridge and insertion proof, including the zero-endpoint caveat and arbitrary-ambient scope.

Final verification: the complete verifier passed normal and optimized Python after the report-variable repair. Both stdout and receipt files are byte-identical across the two runs. The final census contains 150 additional literal special-bank exterior-bridge injections, in addition to the controls stated above. The proof and applications have independent mathematical PASS reviews; the companion audit records the frozen file hashes.
