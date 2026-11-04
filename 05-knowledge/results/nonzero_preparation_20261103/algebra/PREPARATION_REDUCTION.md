# Exact reduction of weak-allowed native preparation to two ready spectra

**PROVED reduction — [independent audit](PREPARATION_AUDIT.md) passed.** This is an exact
necessary-and-sufficient reduction, not a proof that the required integer
spectra are closed. It uses the separately audited native single-arc compiler
and the proved complete joint guard of
[JOINT_CLASSIFICATION.md](JOINT_CLASSIFICATION.md).

## 1. The operation contract

The input is an actual finite tournament T with H(T)=h>1 and one retained
physical root r. All old names, including any other marked ports, must remain.
A preparation may make arbitrary vertex births, use legally guarded native
Hamilton-path complements, and delete only temporary helpers. Its final
tournament V must have H(V)=h, but intermediate Hamilton counts may change.
Weak intermediate and final tournaments are allowed. Other port responses
are transported on the actual named graph, not required to be constant.

For s in {1,2}, preparation must supply an actual auxiliary a satisfying

    sigma_a(V)=s, tau_(a,r)(V)=0.                       (1)

It is followed by the usual two permanent C2 births. The root response F
is always read at the same old name r. This contract deliberately distinguishes
H-neutral preparation from an operation preserving H at every intermediate.

Define the two unlabelled integer sets

    R_s = {H(Q): Q is nontrivial strong, and some a,r satisfy
                    N_Q^+(a)={r}, e_r(Q−a)=s}.         (2)

These are spectra of actual marked objects. Membership is not inferred from
an integral solution for their summary coordinates. The unique/two-start
grammars in [TWO_ENTRY](../../two_entry_funnels_20261103/TWO_ENTRY.md), read
in converse and followed by the sole-exit birth, enumerate these classes.
The five exceptional two-exit scalar-ready families are excluded by the joint
guard and do not belong to this enumeration unless they independently have
another sole-exit mark.

## 2. Exact preparation equivalence

**Theorem.** For every actual T,r with h>1 and each s in {1,2},
an H-neutral preparation satisfying the contract above exists if and only if
h belongs to R_s.

Necessity. The full joint classification forces a's sole exit to be r and
e_r(U)=s for U=V−a. A positive ending count puts r in the terminal SCC C of
U. Every Hamilton path traverses the SCCs in order, so

    e_r(U)=H(prefix)*e_r(C)=s<=2.

Each nontrivial strong SCC contributes at least three Hamilton paths.
Therefore every SCC preceding C is a singleton: the prefix is TT_l. Since
a's only exit is r, that prefix beats a and every vertex of C−{r} beats a.
Consequently V=TT_l→(C+a). If C is nontrivial, C+a is strong: a reaches all
of C through r, and C reaches a through any other vertex. If C is the
singleton r, V is transitive with a penultimate and r last, so H(V)=1,
contrary to h>1. Thus Q=C+a is the required strong ready core with the same
marks. Since H(TT_l→Q)=H(Q), it proves h in R_s. This step forgets outside
old names only to prove membership; it does not delete them in an actual
operation.

Sufficiency. Choose an actual witness Q,a_Q,r_Q for h in R_s, of order m.
Put N=max(|T|,m) and l=N−m. Form W=TT_l→Q. The inherited prefix identity gives

    H(W)=h, sigma(W)=s, tau(W)=0,
    F_r(W)=2^l F_(r_Q)(Q).                            (3)

Embed the old names injectively in W, with the old physical r placed at r_Q;
assign the remaining positions distinct fresh names. This defines a target
adjacency matrix on the named set. It is not a relabelling move on T.
If N>|T|, first birth N−|T| permanent vertices, with any explicit incidence
rule (for example, successive sinks). Then compare its actual matrix with
W's matrix. For each differing pair, use the inherited one-temporary-vertex
compiler of
[METAGRAPH §5](../../adaptive_metagraph_reformulation_20261103/METAGRAPH.md#5-one-temporary-vertex-compiles-any-old-arc-reversal).
Its direct guards are constructed in
[NATIVE_SIX §4](../../multiport_grammar_20261103/component/NATIVE_SIX.md#4-the-inherited-legal-compiler-instantiated-on-every-discrepancy).
Each compiled pair reversal deletes only its helper, so all old and permanent
names survive. After all discrepancies are reversed the graph is exactly W.
This proves the required preparation.

If d unordered pairs differ, the preparation uses N−|T| permanent births,
d helper births and deletions, and2d guarded complements, with peak order
at most N+1 (or N when d=0). No strongness guarantee for these intermediates
is supplied by the general compiler. The subsequent C2 response is exactly

    H_child=3h+2s, F_child=3*2^l F_(r_Q)(Q),
    (F−H)_child=3*(2^l F_(r_Q)(Q)−h)−2s.             (4)

Thus the changing root defect is paid. Full named adjacency remains a valid
carrier for arbitrarily many retained ports; no finite summary is claimed
to determine arbitrary future native moves.

## 3. A constructive unbounded grammar criterion

Let C be a set of integers greater than1 containing an actual seed count h0.
If

    C subset R_1 intersection R_2,
    h in C => 3h+2 and3h+4 in C,                      (5)

then the preparation equivalence followed by C2 gives a native realization
of every finite nonzero word from h0, retaining all names and the root at
each step. The number of permanent names and the root defect may grow; weak
targets are allowed. This is a genuine object-level implication because the
construction selects actual witnesses and then compiles their actual arcs.

Conversely, if every finite nonzero word from an actual seed can be followed
by both H-neutral preparations under this contract, its scalar orbit lies
in R_1 intersection R_2. That orbit is closed under the two maps in (5).
Therefore unbounded two-branch closure in this precise weak-allowed language
is **equivalent** to ready-spectrum membership of every ternary descendant.
The conclusion does not assert (5) for h0=1429 or for any other seed.

The witness selection can be made algorithmic, without an oracle for actual
objects. The cited Busch bound H(Q)>=5^((|Q|−1)/3) bounds every strong witness
of H=h. Enumerate labelled tournaments only at orders n with

    5^(n−1) <= h^3,                                  (6)

test strongness and the exact adjacency-derived guard, and choose the first
witness in a fixed naming order. This is finite and complete, though not an
efficient general integer solver. Alternatively the complete rooted one/two
start grammars greatly reduce the universe. The bound is sourced and checked
in the [induced proof §4](../induced/PROOF.md#4-the-complete-finite-universe-and-all-order-reduction).

## 4. Hostiles and what remains open

The original H1429 graph has no ready auxiliary, and its induced extension
class cannot even hit the next H4289. The explicit old-arc surgery nevertheless
prepares both branches. This is a positive control for the distinction between
an actual graph, its induced descendants, and the set R_s.

The H4289 second-round +4 preparation uses TT1→Q. It therefore illustrates
the name-padding part of the sufficiency proof, with its changed root
response. It is not evidence for an all-strong version of the theorem.
The all-strong version cannot simply use a transitive prefix; it requires
a sufficiently large strong witness. The later
[three-helper compiler](STRONG_COMPILER.md) supplies the strong intermediate
schedule for every such actual witness, so retained target order is the
remaining difference between these two preparation criteria.

The 34-core all-tail selector is only a subset of (2). Failure of that bank
at a later descendant is not nonmembership in R_s. Conversely, a finite
tree of successful memberships does not prove (5). Establishing such an
invariant integer set, or proving that no relevant orbit stays in both ready
spectra, is the new exact closure question. Global Hamilton spectrum coverage
still requires a separate coverage argument even if one such orbit closes.
