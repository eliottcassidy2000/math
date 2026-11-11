# Independent audit of the defect-three escape

**PASS.** The literal constructions and the all-height two-triangle inverse
in [DEFECT_THREE_ESCAPE.md](../reset/DEFECT_THREE_ESCAPE.md) have the stated
scope. In particular, the exclusion of 5993 is not a finite word search.
The finite return at 5995 preserves the same physical marked root and its
original in-neighbour. These results do not give closure of an infinite
interface class or prove Chan--Pak.

## Independent all-height inverse

The audit reconstructs all 64 six-vertex core signatures from actual
tournaments, using nine fresh subset-DP counts per signature. It does not
import the producer's compiler or binary inverse. For each cut it recovers

    F(A,B) = h + 9(g(A)+g(B)-2) + 144AB
               + iA + oB + s(A mod 2) + e(B mod 2),

and checks the coefficients against the receipt. In every row, the linear
and parity coefficients are nonnegative, and the cross-parity coefficient
vanishes. Completeness of this functional form is the inherited marked
core/flank formula; the nine probes independently validate its parameters,
not a replacement proof that finitely many samples determine an arbitrary
function.

For A,B positive and N in {5993,5995}, nonnegativity gives
`AB <= floor((N-h)/144) <= 41`. Direct enumeration of this entire finite
positive-product region therefore covers all positive-coordinate sizes.
On an axis with positive linear coefficient c, the variable is bounded by
`floor((N-h)/c)`. Every zero-slope axis has h=9 and zero parity coefficient,
so its values are exactly `9*g(u)`. Neither target is divisible by 9.
This excludes the whole remaining infinite axis, including every
trailing-one representative, without imposing a word-length cutoff.

This independent rectangle-and-axis calculation reproduces no normalized
5993 solutions and precisely 18 normalized 5995 solutions: the nine listed
cuts, each with (A,B)=(4,9) and (9,4). Arbitrary global source/sink padding
preserves H. If an internal singleton gap separates the two triangles,
the proved first-extension content is 3; such gaps cannot produce either
target. The combination covers every stated two-triangle deletion chart,
including arbitrarily large ones. It excludes no other core types.

## Actual sources, retained names, and return

Fresh counting of the displayed twelve-vertex B gives H=1997,
F_end(0)=2000, and exactly one in-neighbour, vertex 11. Independent
marked-module counting gives the complete coefficient list
`(1997,1134,0,...)` for inflation of vertex 1. The first three inflated
graphs independently verify H=1134*2^m-271, strong connectivity, and
unchanged defect 3. The general indegree-one proof uses the two leading
O letters: after deleting the marked first base vertex, a path starting
at the apex must enter the next base vertex, and then has exactly H(C)
ways through the fixed core. The asserted nonempty final I flank gives
the needed strong-connectivity boundary.

Both retained maps were checked as literal equality of induced adjacency
matrices, rather than equality of counts or isomorphism signatures:

* The earlier thirteen-vertex 1997 source retains ten vertices in B,
  with exact trace `1997,241,619,1997`. Its old marked root becomes B's
  vertex 0, and its old in-neighbour becomes B's vertex 11.
* Starting at B, the two actual births give the displayed fourteen-vertex
  5997 graph. Replaying all fifteen residual arc reversals reproduces
  every stated full-graph count. Reversing only the apex incidence at
  vertex 8 then gives the displayed strong 5995 graph, with every one of
  the fourteen names retained. The peak count is 12215.
* The old root's resulting profile is exactly
  `(H,S,q,F,b)=(5995,339,1897,19759,13764)`. A full endpoint-profile scan
  in both orientation gauges confirms that this particular endpoint has
  no defect-1 or defect-2 root.
* The nine-vertex retained map into the ready 5995 target gives the exact
  trace `5995,107,253,545,1085,2165,5995`. Old vertex 0 maps to the target's
  ready root 13, and old vertex 11 maps to its sole in-neighbour 0. Fresh
  counting gives F_end(13)=5996. Two further actual births give 17987.

The two deletion-and-insertion returns are same-count endpoint returns;
their intermediate graphs have the explicitly lower counts. The initial
5997-to-5995 arc repair likewise has its stated nonneutral intermediate
prices. No step has been silently promoted to H-neutrality.

## Classification dependencies and limits

The inherited all-order ready inverse, branching low-defect classification,
and complete ready-target census are dependencies; this independent audit
does not rerun those larger enumerations. Their use in the note is properly
scoped. The weak-source extension follows from the SCC product: for
0<b<H the root lies in the first SCC, and the suffix count divides both b
and H. At these targets and b<=6 that suffix count is 1. For the exclusion
of b=1,2 specifically, oddness alone forces the suffix count to be 1.
The first SCC cannot be a singleton with positive defect, so the strong
classification applies after discarding transitive padding. This validates
the inference that B attains the minimum positive defect 3 at H=1997,
assuming the cited all-order inverse. Defect zero is a separate boundary.

The singleton-inflation and partial-transport failures retain their stated
restricted universes. Their complete coefficient formulas rule out those
uniform lifts; they do not rule out a different core, different retained
subset, several modules, or a new source-return operation. In particular,
the finite 5995-to-ready return is a positive construction, not a proof of
an infinite repeatable return.

## Reproduction

Run `python -B .../berggren/audit_defect_escape.py` and repeat with
`python -O -B`. The matching `defect_escape_audit.json` records all 64
independently recovered signatures, both complete normalized inverse
answers, the full retained-source traces, and SHA-256 pins of the reviewed
producer files. The 576 signature probes and the literal construction
checks pass in both modes; invariant checks do not depend on assertions.
