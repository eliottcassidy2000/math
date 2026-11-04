# Constant-region surplus in shared binary circuits

**Status: PROVED incidence identity and conditional deletion certificate.**
This extends the constant-antichain counting argument recovered from Session
118. It is a local structural result, not a new unrestricted circuit lower
bound, a proof of the external elimination candidate, or a literature-priority
claim. The only finite computation used below is corroboration.

## 1. Inheritance and the exact question

The closest recovered mechanism is the constant-antichain lemma in
`05-knowledge/results/session118_shared_gate_deletions_20261030/audit_ap2/AUDIT.md`:
r incomparable constant seeds in a live binary DAG force 2r distinct
constant/unary deletions when the output remains nonconstant. Its proof counts
edges after closing the seeds under two-constant-parent gates.

The canonical hostile is the comparable chain `ab -> abc -> abc XOR x`,
which has two comparable constant seeds but only three deletions. The corrected
near miss is the genuine tight GP/P2 prefix in
`session118_shared_gate_deletions_20261030/worst_cases/TIGHT_CYCLE_HOSTILE.md`:
the common child of two affine equations need not disappear. The least-used
sidecar is *distinct child incidence*, retained by that work but discarded by
the scalar sequence inventories.

This note makes that incidence sidecar exact. It identifies every surplus
term and the full equality topology. In particular, merging two constant
branches is not a license to double-count their common child.

## 2. Model and definitions

Let D be a finite acyclic single-output circuit. A gate has at most two input
pins, each reading a primary input, an allowed boundary value, or an earlier
gate. All gates are in the ancestry of the output. Count edges between
distinct gates, not occurrences of pins: repeated pins count once. The output
reference itself is not a gate child.

Let C be a nonempty set of old gates known to be constant under a specified
legal restriction. Assume:

1. the output is not in C;
2. whenever every input pin of a gate reads a member of C, that gate is in C.

For ordinary binary truth tables, condition 2 holds if C includes all gates
that are constant after the restriction, or if C is obtained from certified
seeds by constant propagation. The restriction may have a nonrectangular
allowed domain: the constant values need only be valid on that domain.

For a gate v, let `d_C^-(v)` be its number of distinct parents in C, and let
`d^+(v)` be its number of distinct gate children in the old current DAG. Set

    B = {v outside C : d_C^-(v) > 0},
    r = #{v in C : d_C^-(v) = 0},
    f = #{v in C : d_C^-(v) = 1},
    E = sum_{v in C} (d^+(v)-1).

Here r counts sources *inside the constant region*, not arbitrarily chosen
constant gates. The letter f counts one-parent constant nodes, not Boolean
functions. Such a node may be forced by an annihilating pin, by other known
semantics, or by repeated input pins. No independence of the constant values
or realizability of arbitrary seed assignments is assumed.

Because every C gate is live and is not the output, each has a gate child,
so E is nonnegative. Because a B gate cannot have both pins in C, it has
exactly one distinct C parent.

## 3. Exact surplus identity

**Theorem.** Under these hypotheses,

    |C| + |B| = 2r + f + E.                       (1)

**Proof.** Put c=|C| and b=|B|. The C-induced graph has r vertices of
indegree zero, f of indegree one, and c-r-f of indegree two. Hence

    e(C,C) = f + 2(c-r-f) = 2c-2r-f.

Each vertex of B has exactly one distinct C parent, so `e(C,B)=b`. Counting
all edges out of C gives, independently,

    e(C,C)+e(C,B) = sum_{v in C} d^+(v) = c+E.

Equating these expressions yields (1). This counts every physical gate and
every distinct child once. It allows arbitrary fanout, sharing, and repeated
pins. No Boolean function is replaced by a formula. QED.

**Why the source count is the repair.** In the comparable-chain hostile,
`C={ab,abc}`, `B={abc XOR x}`, r=1, f=1, E=0. Thus (1) gives 3, not the
incorrect 4 obtained by treating both C vertices as independent sources.

### Exact equality classification

The certificate has size exactly 2r if and only if f=E=0. Equivalently:

- each C vertex has exactly one distinct gate child;
- each nonsource C vertex has exactly two distinct C parents.

Consequently C is a disjoint union of full binary merge trees, oriented from
their source leaves toward their roots. Each tree root has one outgoing edge
to its own distinct boundary gate in B. Distinct trees cannot share that
boundary: it would then have two C parents and belong to C. Conversely every
such merge forest has f=E=0 and attains (1) with equality.

This classifies the size of the *certified set C union B*, not necessarily
the entire simplification. Further alias propagation or semantic rewriting
may delete more gates.

## 4. Deletion consequences and the output convention

For the deletion corollaries in this section, assume additionally that the
restricted output function is nonconstant. This is stronger than merely
requiring that the output gate is absent from the chosen certified subset C;
the incidence identity itself does not need semantic nonconstancy.

In the full binary Boolean basis, each C gate can be replaced by its certified
constant. Each B gate then becomes a constant or a literal of its remaining
input and can be rerouted. Negations on internal edges are absorbed into the
truth tables of their users. The replacements use predecessor values and do
not create a cycle. Thus, in a model permitting a free final constant or
literal, at least

    2r + f + E

distinct old gates can be deleted. Shared users are redirected, not copied.
If a constant region contains gates whose constancy is semantic rather than
found by the chosen algorithm, their constant values must first be separately
certified. Identity (1) is not such a certificate.

The pinned external candidate explicitly permits free final constants and
literals. That is the convention inherited by Session 118. See Section 1 of
`submission/elimination.md` at commit
`4a61f024a4993592d97b49ca7d00979029b7fece` in
`Layr-Labs/circuitbreaker-llm-verifier`.

For ordinary strict B2 gate complexity with no free terminal NOT, the always
safe conclusion is instead a net saving of at least

    2r + f + E - 1.                               (2)

There is at most one new terminal NOT, and none is needed if the reduced
output is a positive primary input or if its negation can be absorbed into
a surviving output gate. In particular, if the restricted function depends
on at least two unrestricted primary inputs, the full `2r+f+E` bound holds.
This last criterion is stated in the ordinary fixed-coordinate model, not
as an uncharged affine change of input basis.

The correction in (2) is sharp. Take `A=a AND b` and `O=NOR(A,x)` and fix
a=0. Then C={A}, B={O}, r=1, f=E=0. The signed-output model deletes both
gates and outputs NOT x for free. Strict B2 needs one gate for NOT x, so
the net strict saving is one, not two. This is a model boundary, not a
refutation of Session 118's inherited permissive-output convention.

## 5. Antichains, GP, and honest telescoping

In a DAG whose gates have exactly two distinct parents, if an ancestry-antichain
S of r seeds is closed only under gates with both parents already in C, the
induced sources are exactly S and f=0. Hence

    |C|+|B| = 2|S| + E.

This refines the recovered 2r bound: all excess of that formal closure is
exactly its aggregate distinct-child fanout. Extra one-pin forcing, when
included in C, supplies the additional f term of (1).

For a single GP constant seed G, pure two-parent closure produces no other
constant gate in a distinct-pin DAG. Its certificate consists of G and its
k distinct children; equation (1) reads `1+k=2+(k-1)`. Thus it recovers the
existing child-count surplus K, rather than manufacturing a new one. If the
actual propagation additionally contains one-parent constant gates, those
contribute f. The tight one-child case has r=1, f=E=0 and is exactly the
two-gate base case.

For a sequence of restrictions, apply the theorem to each *current* live
DAG, before its corresponding deletions. If C_j and B_j refer to objects
deleted in that step, these sets are disjoint across steps, so

    total deleted gates >= sum_j (2r_j+f_j+E_j)

in the permissive-output convention. Do not use historical outdegrees after
their children have disappeared, and do not add E again if a row's existing
child-count guarantee has already charged it. The identity exposes possible
surplus; it supplies no lower bound on its frequency in minimum, low-dispersal
circuits. The Session 118 GP/P2 bottleneck therefore remains OPEN.

## 6. Concept board and connection contract

1. Constant propagation / induced source forest: (1) retains exactly the
   incidence lost by a seed count. Cheapest test: comparable chain.
2. DAG sharing / distinct child degree: merging preserves the identity;
   branching increases E. Cheapest test: merge-then-boundary versus fork.
3. Elimination inventories / local surplus: f and E telescope only on
   current disjoint deletions. Cheapest test: one-child GP equality.
4. Boolean model / terminal sign: the incidence count is unchanged while
   strict gate savings can lose one. Cheapest test: NOR(A,x), A=0.

The connection is `restricted shared DAG -> (r,f,E,B)` by retaining its
constant-induced edges and boundary. It preserves the exact number of
certified deletable old gates. It forgets the surviving Boolean function,
input-coordinate cost, and global minimality. Those must remain separate
sidecars. The equality merge-forest and strict-terminal hostile are the
decisive tests; no transfer to an asymptotic circuit coefficient is claimed.

## 7. Reproduction and scope

Run `python3 verify.py` and `python3 -O verify.py` in this directory.
The checker is standard-library-only, uses explicit checks rather than
optimization-removable assertions, and prints its exact universe counts.

- All live distinct-pin binary DAG topologies with two external vertices
  through six gates; every nonempty constant region excluding the output
  and closed under all-C pins. The test checks (1), equality, components,
  and both repeated and unrepeated boundary sharing.
- All live repeated-pin topologies with two external vertices through five
  gates, with the same complete constant-region enumeration.
- Every all-16-table B2 circuit of at most two gates on three inputs with
  unordered pins (truth-table asymmetry is retained), including repeated
  pins; every nontrivial coordinate restriction. The truth-table test uses
  actual constants and checks signed and strict simplification counts.
- Explicit sharp merge forests, one-pin forcing, fanout, comparable seeds,
  and the strict terminal-negation hostile.

These finite universes do not establish all-size circuit coverage. The
displayed edge count proves the all-size identity, and literal substitution
proves its scoped deletion consequence. No installs, source candidate code,
startup helper, or forbidden ledger were used.
