# A prime branch requires changing the old source relation

**PROVED relative extension reduction and explicit named construction;
FINITE-EXACT complete attachment and distance certificates.** The specified
balanced tournament C_2 has H=59 and admits no induced extension with H=181
at any order. A source-changing native excursion repairs this branch while
retaining all eight physical names, including the root and child. It preserves
H at its endpoints, not at intermediate steps. The full spectrum remains OPEN.

## T1 The ternary extension question has a finite strong target

Let T be a fixed strong tournament with H(T)=h>=3, and let N be either
3h+2 or 3h+4. Then T has an induced extension of count N if and only if it
has a **strong** induced extension of count N. The latter has order at most
L(N), where L is the inverse of the Moon–Busch strong minimum.

Proof. In any extension U, the strong old tournament T lies within one SCC K.
The inherited insertion injection gives H(K)>=h, and SCC multiplicativity
gives H(K)|N. Here N is odd and not divisible by three. Every proper positive
divisor of N is therefore at most N/5<h, since h>=3. Hence H(K)=N. Delete
all other SCCs; no old vertex is deleted. Conversely, a strong extension is
an extension. The cited strong floor bounds its order. QED.

This extends the same SCC/divisor mechanism already used for p<q<3p in
[INDUCED_OBSTRUCTIONS](../../metagraph_primegame_20261102/INDUCED_OBSTRUCTIONS.md).
It is a specialization of inherited finite H-decidability to an actual
old source, not a new general decision theorem. It applies to both ternary
branches regardless of whether N is prime or 7p.

For reference the exact floor is

    f(3k)=3*5^(k-1), f(3k+1)=5^k, f(3k+2)=9*5^(k-1), k>=1.

The primary theorem and its scope are pinned in
[Session133 core reduction, section 1](../../session133_fractran_spectrum_refinement_20261101/core_reduction/THEORY.md).
The bound does not limit a native excursion that changes old edges.

## T2 A complete obstruction at the first balanced prime source

Use the two-apex word 202032, with apices 0->1 and remaining vertices
2,...,7 in increasing order. A letter's bit i means apex i beats that
chain vertex. Its physical names are

    (c_2, x, reset_root, reset_child, y, c_1, birth_a, birth_b).

This is exactly the prior balanced C_2 after the explicit old-label map
(0,1,4,5,6,7,2,3). Its physical root 2 has profile (59,17,17,59), and the
protected root-child arc is 2->3. No equality-of-H substitution is used.

Since f(10)=125 and f(11)=225>181, T1 reduces the existence of any induced
H181 extension to zero, one, or two additional vertices. The
[interior producer](../interior/check.py) enumerates all first cuts and all
eligible second cuts, using monotonicity only to skip first children whose
count already exceeds 181. The independent checker in this directory
enumerates all **131072 ordered two-birth cuts without that pruning**.
Both find no H181 extension.

The independent method enumerates every directed simple path of each old
graph to obtain start and end counts on every vertex subset. Let W_uv count
ordered nonempty two-path covers ending at u and beginning at v. If the new
vertex is beaten by I and beats O, its count is exactly

    sum_(u in I) end_u + sum_(v in O) start_v
        + sum_(u in I,v in O) W_uv.

Deleting the newborn from a Hamilton path proves this partition into three
disjoint cases. It uses neither the producer's subset recurrence nor its
pruning. Its six H179 positive pairs reproduce the producer after converting
in-masks to out-masks; the mask convention is serialized explicitly.

Consequently arbitrary old-root changes, interior placements, fresh port
choices, and any number of births cannot give this actual C_2 an H181
descendant while retaining all old arcs. Both 59 and 181 are primes.
This refutes a universal source-preserving ternary closure on the balanced
outputs, rather than merely leaving that hypothesis unproved.

## T3 Three old arc changes repair the branch with the named ports intact

The inherited fifth exceptional ready family in
[READY section 3](../../two_entry_funnels_20261103/READY.md)
has H=14*2^k+3. At k=2 it gives a seven-vertex H59 source with a two-ready
root. Adjoining a global sink preserves both H and the root's four-coordinate
profile: the sink must be appended to the unconstrained final path bin.
This makes a source on the same eight physical names as C_2.

Choose the named alignment given in the receipt. It yields rows

    (66,212,121,243,193,211,0,69)

and the same physical root 2 and child 3, still with 2->3. The root profile
is now (59,5,5,61), so F-H=2. Precisely these three old pairs differ:

    {2,7}, {4,5}, {6,7}.

The [inherited one-auxiliary compiler, METAGRAPH section 5](../../adaptive_metagraph_reformulation_20261103/METAGRAPH.md#5-one-temporary-vertex-compiles-any-old-arc-reversal)
implements each arc reversal with one temporary birth, two legal native
K_P moves, and deletion of that temporary vertex. The certificate contains
every intermediate relation, incident cut, Hamilton word, and named deletion.
For all three gadgets the clocks are three temporary births, six native
moves, three deletions, and peak order nine. All eight old names survive.

After the successive temporary deletions the H values are

    59 -> 161 -> 125 -> 59.

Counts at the native intermediates are also saved. The root-child edge is
reversed at one intermediate native move and restored; it is not claimed
to be fixed throughout the excursion. Thus this supplies physical port
transport with updated attachment data, not an unchanged observer or a
stepwise neutral path in the H59 fibre.

Now add vertex 8 beaten only by root 2, then vertex 9 beaten only by 8.
The final rooted profile at 8 is

    (181,59,59,181).

The two original ports retain names 2 and 3; the next balanced interface
has explicitly recorded fresh names 8 and 9. Vertex 6 remains a global
sink. The construction retains that physical name; deleting it would
recover the strong H181 core, but is not part of this program.

**Sharp endpoint distance with its guard.** Three is the minimum number
of old arc differences needed to replace C_2 by a same-order H59 source
that retains the endpoint arc 2->3 and admits any H181 induced extension.
The interior certificate checks all 28 one-pair reversals and all 378
two-pair reversal sets. None of the former preserves H59. Five of the
latter do; one reverses 2->3 and is excluded by the guard. The remaining
four are strong, and their complete two-birth searches contain no H181
extension. T1 promotes those four exclusions to all orders. The three-pair
construction above attains the bound. This is not a shortest native-move
theorem, nor a distance theorem without the protected endpoint edge.

## Attachment data and scope of the repair

The receipt stores the seventeen inherited **general two-live-port static
birth queries** at the original source, after the source exchange, and at
both the old and new pairs of the H181 output. These are distinct from the
seventeen two-apex word-suffix queries. Full named graph rows remain the
carrier when native moves change old edges; a query list does not authorize
an unrecorded source replacement.

This example proves that one blocked prime branch can be repaired by a
small endpoint change compiled into actual native operations. It does not
make the ready-root method universal: the inherited H35 obstruction has
no one-ready or two-ready representative at any order. The induction must
permit other terminal gadgets or other predecessors, as discussed in
[the conditional reformulation](../conditional/NOTE.md).

Run `python3 -B check.py` and `python3 -B -O check.py` from this directory.
The output must reproduce `check.json`. The independent native replay is
in [interfaces/audit_transport.py](../interfaces/audit_transport.py).
