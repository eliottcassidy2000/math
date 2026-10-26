# Recovering a tournament from the arcs used by its Hamilton paths

2026-10-26. **PROVED: P1–P4 independently audited PASS; classical separator cited.**
[Independent mathematical audit](../audit_paths/PROOF_AUDIT.md).
The all-order reconstruction below is a deduction from the classical
nontraceable-arc separator theorem and a proved-here deletion lemma; no external
priority is claimed. **CITED:** Busch's Theorem 4.12, with exact source below.
**FINITE-EXACT:** exhaustive labelled tournaments through order six, checked by
two different path enumerations. Local labels P1–P4 belong only to this note.

A labelled tournament is determined by the **union of directed arcs occurring
in its directed Hamilton paths**. Thus its full labelled Hamilton-path language
also determines it. Keeping one path or only a path count does not suffice.
The inverse is polynomial in the supplied union's size; this does not make
listing all Hamilton paths polynomial in the tournament's size.

## 1. Inheritance, representations, and the question

Closest inherited mechanism: [Session25 agreement certificate, Theorem 1](../../session25_spectrum_wiring_20261021/session25-wiring/RESULT.md),
which recovers a *particular* path from two endpoint tournaments on a common
labelled vertex set. [Strong-atom C1–C5](../../strong_atom_continuation_vitali_20261026/README.md)
keeps the relative alignment and the names of re-marking operations; its complete
legal-guard signatures were verified injective through order seven. The present
argument makes that named-guard assertion all-order. The canonical hostile to
forgetting alignment is the same packet's unbounded family of `2^k` aligned
states with equal separately rooted endpoints. The corrected near miss is
[rooted-lines R3](../../rooted_lines_frobenius_sheets_20261026/README.md): separate
entrance and exit types do not recover the joint orbit (the 289/253 sandwich).
The least-used sidecar here is a **one-vertex deletion reachability relation**.

[HYP-9416 / fresh re-marking](../../../hypotheses/HYP-9416-fresh-re-marking-two-port-quotient-joint-classes-are-minimal-and-line-responses-classify-strong-atoms.md)
keeps multisets of successor *classes* after forgetting operation names. That
is a different observer and remains open in its stated all-order forms.
[THM-4526 / shaved tournaments](../../../../01-canon/theorems/THM-4526-shaved-tournaments-path-plus-chord-and-the-largest-unavoidable-graphs.md)
asks for a single spanning digraph embeddable in every tournament; our union
depends on the particular tournament and is not a universal shaving. Its
MISTAKE-555 attribution correction was read before making literature claims.

Concept board: labelled legal language; used-arc union; ordered strong
components; deletion cuts; named versus fresh re-marking. The anchor is exact
reconstruction, the niche is the sharp five-vertex unused-arc obstruction, and
the wildcard is whether the labelled language is a sufficient continuation
observer. P3 and P4 answer that last question without resolving unnamed
successor-class aggregation.

Let T be a finite tournament on a specified set V, and write

- `L(T)` for the set of vertex permutations `(v1,...,vn)` with every
  `vi -> v(i+1)` an arc of T;
- `U(T)` for the spanning digraph whose arcs occur consecutively in at least
  one member of `L(T)`;
- `R(D)` for directed reachability in D, including length-zero paths.

Singletons are allowed; empty orders can be handled separately and contain no
orientation information. A nontrivial strong tournament has at least three
vertices. Strong components of a tournament have a unique linear order
`C1 -> ... -> Cr`; every cross arc goes forward in that order. A *separating*
vertex means its deletion destroys strong connectivity, not weak connectivity.

Connection contract: `T -> L(T) -> U(T)` retains common vertex labels and directed
consecutive incidence. Grouping paths into their union drops the individual
words syntactically, but under the promise `U=U(T)` P3 reconstructs T, and then
`L(T)` exactly. No extra orientation sidecar is required. At the class level,
allowing one common relabelling gives the corresponding isomorphism result;
independently relabelling each observation is a different quotient.

## 2. The classical separator input

**CITED, not proved newly here.** Arthur H. Busch, *Arc-traceable Tournaments*,
PhD thesis, University of Colorado at Denver, 2005, **Theorem 4.12**, printed
pp.45–47 (PDF page indices 52–54), states the following necessary condition.
For a strong tournament T and an arc `x->y` on no Hamilton path, there is
`z` different from x and y such that T−z has at least four ordered strong
components, x is in the first and y in the last; its second component dominates
z, and z dominates its penultimate component. We read the statement and proof,
and inspected the rendered statement page. [Primary thesis PDF](https://digital.auraria.edu/downloads/ma7yx-pcr25/AA00005831_00001.pdf).

This is **not** an iff characterization of nontraceability: the source explicitly
provides a traceable example with the same necessary structure. The reconstruction
uses only the existence of z with x and y in different ordered components.
The full extra structure is retained in our exact finite tests.

The thesis's Proposition 4.1 (p.39) says every arc incident with a nonseparating
vertex of a strong tournament lies on a Hamilton path. Theorem 4.7 (pp.42–43)
says two internally vertex-disjoint paths from y to x suffice for `x->y` to lie
on a Hamilton path. The proof of 4.12 uses 4.7 and Menger's theorem. These
statements were also read, rather than inferred from the paper's abstract.
They are due to Busch/Jacobson/Reid, not the unrelated Thomassen arc-cycle
results: the later paper is A. H. Busch, M. S. Jacobson and K. B. Reid,
[*On arc-traceable tournaments*, J. Graph Theory 53 (2006), 157–166](https://doi.org/10.1002/jgt.20173).
The later local-tournament paper restates this cut theorem as Theorem 1.10,
but our numbered primary proof dependency is the thesis itself.

We also use the classical Hamilton-cycle fact for strong tournaments
(Camion 1959). For completeness, an elementary proof is given in §7, along
with the Hamilton-path insertion fact, so no uninspected statement about
arbitrary digraphs is imported.

## 3. P1: outer strong components survive in the union

**P1.** The ordered strong components of U(T) are exactly those of T. For each
component C, `U(T)[C] = U(T[C])`. Between consecutive components every cross arc
belongs to U(T); between nonconsecutive components none does.

*Proof.* A Hamilton path cannot leave an earlier strong component and later
return, since all cross arcs go forward. It traverses all components once, in
their fixed order. Conversely, concatenating an arbitrary Hamilton path of each
component gives a Hamilton path of T. This proves the restriction equality.
Each nontrivial strong component has a Hamilton cycle; every edge of this cycle
survives in a Hamilton path by breaking a different cycle edge. Thus U(T)[C]
is strong. Singletons cause no difficulty. A strong component has a Hamilton
path beginning, or ending, at any prescribed vertex, obtained by breaking its
Hamilton cycle at that vertex. Hence every selected cross arc between adjacent
components can be realized as the joining arc of a global path. Nonconsecutive
components cannot be consecutive along a Hamilton path because an intermediate
nonempty component must be traversed between them. ∎

## 4. P2: deleting a vertex preserves reachability inside a strong source

**P2 (deletion lemma).** If T is strong and z is a vertex, then

`R(U(T) - z) = R(T - z)`.

The proof was developed jointly with the root integration lane. The strong-source
hypothesis is essential: for the transitive tournament `0->1->2`, with `0->2`,
the union is just the path; deleting 1 removes its only 0-to-2 route in U.

*Proof.* The inclusion from left to right follows from `U(T) subset T`.
Let `C1 -> ... -> Cr` be the strong components of T−z.

**Case r>=2.** Strong connectivity of T forces `z->a` for some `a in C1`
and `b->z` for some `b in Cr`; otherwise C1 could not be reached from z,
or Cr could not reach z. For each Ci choose a Hamilton cycle when it has
more than one vertex. Concatenations of component Hamilton paths can be
extended to T by either:

1. fixing the C1 path to start at a and prepending z; or
2. fixing the Cr path to end at b and appending z.

Using the first construction, any chosen cycle edge in any `Ci, i>=2`, can
be retained: break a different edge of its Hamilton cycle, while the required
start in C1 stays fixed. Using the second construction retains any chosen
cycle edge in C1, while the required end in Cr stays fixed. Thus every Ci is
strongly connected inside U(T)−z. One prepended full path also supplies an
arc of U(T)−z between each pair of consecutive components. Its condensation
therefore has exactly the same ordered reachability as T−z. No backward
cross-component arc can appear because U(T) is a subgraph of T.

**Case r=1.** When T is a nontrivial strong tournament this case has
`|T−z|>=3`; a two-vertex tournament is not strong. Choose a Hamilton cycle
C of T−z. Fix any edge e of C. For another edge `f=(a,b)` of C, break f to
obtain the Hamilton path starting at b and ending at a. This path retains e.
It extends to a Hamilton path of T by prepending z if `z->b`, or appending z
if `a->z`. If both options fail for every f different from e, then every such
f satisfies `b->z` and `z->a`. The other `|C|-1` edges contain two consecutive
edges f,g because `|C|>=3`. Their shared vertex would both dominate and be
dominated by z, impossible. Hence every e is globally traceable and C lies
in U(T)−z, proving strong connectivity. The one-vertex strong T is trivial
under the empty-deletion convention. ∎

## 5. P3: an explicit inverse of the used-arc map

**P3.** For every finite labelled tournament T, the following deterministic
procedure reconstructs T from `U=U(T)` alone.

1. Compute the ordered strong components of U.
2. Retain all observed arcs. Between different components, orient every
   missing pair from the earlier to the later component.
3. For a missing pair `{x,y}` in the same nontrivial component C, inspect
   `U[C]−z` for vertices `z in C \ {x,y}`. Choose any z for which exactly one
   of `x reaches y` and `y reaches x` holds. Orient the pair in that direction.

Such z always exists, and every qualifying z gives the same answer.

*Proof.* P1 proves step 2 and identifies `U[C]` with `U(T[C])`. In step 3,
let the actual missing arc be `x->y`. The cited separator theorem supplies
z for which x and y belong to the first and last components of T[C]−z.
P2 transfers this asymmetric reachability to U[C]−z, proving existence.
Conversely, any asymmetric x-to-y reachability in U[C]−z is the same in
T[C]−z by P2. The actual orientation cannot be `y->x`, since that edge would
itself give the forbidden reverse reachability. This proves soundness and
consistency, including simultaneous recovery of many missing arcs. ∎

The certificate for a hidden internal orientation is the vertex z, a directed
path in U[C]−z from x to y, and its ordered strong-component cut excluding a
return path. The certificate for a hidden outer orientation is the global SCC
order. This is not ordinary transitive closure: U is itself strong whenever T
is strong, so ordinary reachability alone would make both directions visible.

A conservative implementation bound is `O(n^4)` time and `O(n^3)` space from
an n-by-n supplied U matrix: compute reachability after each vertex deletion
inside each strong component, then inspect those matrices for each pair.
More efficient SCC routines can reduce the bound, but are unnecessary here.
From an explicit list of H Hamilton paths, constructing U costs `O(nH)` first.
Neither that list nor a truth table on all n! permutations is polynomial-sized
in n in general. This is a promise inversion result for a complete used-arc
union, not a recognition algorithm for arbitrary digraphs advertised without
a realizability check.

## 6. P4: continuation consequence and sharp hostile controls

**P4.** If `L(T)=L(T')` on the same vertex set, then `T=T'`. The same holds
if only `U(T)=U(T')` is assumed. Relabelling commutes with reconstruction,
so `U(T) isomorphic U(T')` iff `T isomorphic T'`.

This follows by applying P3. In particular, the complete named legal
re-marking guard vector in the strong-atom packet C5 is injective at **all
orders**, not just the tested orders. Once names and failure states are retained,
no two different labelled current tournaments can have all the same one-step
guards. This does not assert that histograms of successor classes determine a
source, or decide the fresh fixed-port or free re-marking quotients of HYP-9416.
Nor does it prove metagraph connectivity, an unrestricted A000568 recurrence,
or a path-count classification.

The first strong hostile in lexicographic labelled enumeration occurs at n=5,
mask 10, where bit k corresponds to the k-th lexicographic pair `i<j` and bit 1
means `i->j`. Its adjacency rows are `(20,1,2,7,14)` and the missing arc is
`4->1`. The exact rows, all paths, and cut certificate in `normal.json` are the
machine authority for this numerical presentation. No strong tournament of
order at most four has an unused arc, both by the classical four-components
lower bound and by exhaustive checking. The five-vertex example is strong,
so U is strong too: a plain reachability matrix cannot recover the hidden
orientation. Some legal paths place 1 before 4, so an unused arc does not
imply a fixed relative order in every Hamilton path.

Additional controls:

- The transitive 3-tournament has one Hamilton path and a missing long arc;
  P3 recovers it from the ordered singleton components.
- That same directed path, considered as an incomplete digraph, has the same
  Hamilton language as the transitive tournament. Completeness is essential.
- The path `(0,1,2)` belongs to both the transitive tournament and the directed
  triangle. One selected legal path is insufficient.
- Full path words are recoverable from U under its promise, although U does
  not explicitly list them. A count or an independently relabelled histogram
  is not the same input.

## 7. Elementary existence facts used in P1–P2

Every tournament has a Hamilton path: inductively insert a new vertex v into
an existing path. Prepend if v beats its first vertex, append if its last beats
v, and otherwise an adjacent transition `a->v->b` supplies an internal slot.

A strong tournament on at least three vertices has a Hamilton cycle (the
classical Camion fact). It contains some directed cycle. Choose one of maximum
length. An outside vertex with both an in-neighbour and an out-neighbour on
the cycle can be inserted between a consecutive in/out transition, contradicting
maximality. Thus every outside vertex either dominates the whole cycle (class A)
or is dominated by it (class B). If any outside vertices remain, strong
connectivity forces both classes nonempty and a path from B to A. Its first
transition from B to A is an arc `b->a`; replacing any cycle edge `u->v`
by `u->b->a->v` yields a longer cycle. Contradiction. This proves the fact,
including the arbitrary-start/end path consequence used above.

## 8. Reproduction and scope

Run from the repository root:

```
python3 05-knowledge/results/growing_sieve_path_20261026/paths/check_paths.py
python3 -O 05-knowledge/results/growing_sieve_path_20261026/paths/check_paths.py
```

The producer fills Hamilton-path constraint subcubes in the space of **all
labelled tournaments** for each `1<=n<=6`. Independently, a vertex-by-vertex
DFS enumerates every path of every tournament and compares the complete
languages and unions. The checker verifies P1, P2 in each strong component,
P3 on every input, and all four clauses of Busch's separator for every missing
internal arc. No isomorphism, symmetry, strongness, or sampling filter removes
a tournament from the reconstruction census. Explicit `require` calls survive
Python optimization. Exact gate counts and witnesses are in the two matching
JSON receipts; source versions and primary reading receipts are pinned in
`SOURCES.json`.

The proof, not the finite census, supplies all-order reconstruction. The
classical cut theorem remains a necessary condition rather than an iff test
on an arbitrary guessed orientation. The outstanding continuation targets
include the fresh unnamed class observers and efficient production/recognition
of complete observer data; none is silently discharged by the promised inverse.
