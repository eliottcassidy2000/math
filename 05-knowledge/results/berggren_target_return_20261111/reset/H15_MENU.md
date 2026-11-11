# An actual +2 escape through the S4 core

**PROVED / exact finite-core certificates:** the complete all-height
S4→C3 and C3→S4 deletion-chart menus exclude 5993; the S4 menu has six
normalized charts. **FINITE-EXACT:** a named 1997→5993→ready5993 path
retains the same physical root and its sole incoming neighbor. **PROVED:**
one actual infinite family of these retained-subgraph replacements has
explicit positive target drift. **OPEN:** a uniform target-preserving
repeatable return and Chan–Pak.

This continues [DEFECT_THREE_ESCAPE](DEFECT_THREE_ESCAPE.md). It changes
the residual core instead of only transporting singleton gaps between
the two retained triangles. The source is the SAME 12-vertex defect-three
source B used there for the actual +4 branch to 5995.

## Complete next-core inverse

Write S4=(2,4,9,3), with H(S4)=5, and C3=(2,4,1), with H(C3)=3.
For the ordinal joins

    S4→C3=(114,116,121,115,32,64,16),
    C3→S4=(122,124,121,16,32,72,24),

all 128 apex cuts of each core were inverted at target 5993, with
arbitrary transitive left and right flanks. Both answer sets are empty.
These are precisely the H=15 residual intermediate types in the gap
transport, but this inverse also permits every cut and arbitrarily
large flanks. Its scope is these fixed adjacent cores with external
flanks, not arbitrary seven-vertex cores or arbitrary internal gaps.

The inherited exact binary inverse supplies the all-height certificate.
For its positive count coefficients (h,d,L,i,o,s,e), the mixed term
L·A·B bounds the smaller address by floor(sqrt((5993−h)/L)). The other
address is then solved exactly by the binary recurrence. Axes are
handled separately by the inherited trailing-one normalization.
No finite flank-length cutoff is used. All 272 coefficient rows and
inverse answers are recorded in `H15_inverse.json`.

The 16 cuts of S4 alone have exactly six normalized solutions:

| Apex cut | Left address | Right address | Order |
|---:|---:|---:|---:|
| 1 | 3 | 57 | 13 |
| 7 | 57 | 3 | 13 |
| 4 | 0 | 693 | 15 |
| 13 | 693 | 0 | 15 |
| 8 | 0 | 460 | 14 |
| 14 | 460 | 0 | 14 |

All six normalized graphs are strong, verified by independent subset DP.
This counts chart presentations. Source/sink padding may yield larger
weak presentations; it is not counted as another normalized solution.

## A retained-root +2 construction

The first S4 solution is

    T=(4094,4092,4040,8144,8164,8140,
       8064,3840,3584,7168,6144,4096,391).

It has H=5993, root 0, sole incoming neighbor 12, and defect 5. To pass
directly from B (H=1997, defect 3) to T, retain B vertices

    (0,1,2,4,5,6,7,8,10,11)

and map them respectively to T vertices

    (0,1,2,5,3,6,9,10,11,12).

Delete old vertices 3 and 9, then insert target vertices 4,7,8. All
retained arcs agree. The exact count trace is

    1997 → 409 → 1599 → 3067 → 5993 = 3·1997+2.

Ten of twelve source vertices survive. In particular the same physical
root old0 remains root0, and its same sole incoming neighbor old11
becomes target12. This is a concrete defect-3→defect-5 transition.
The previously recorded +4 construction starts from this same B and
returns to a defect-1 representative of 5995. Thus both finite digits
now have actual constructions from one named source.

## Returning 5993 to a two-ready state

The inherited all-order ready inverse supplies exactly four strong
two-ready rooted classes at 5993, of orders 13,12,11,11. Fix the source
T above, its root0, and its sole incoming neighbor12. Exhaustively
testing injective arc-preserving maps on all source subsets that retain
those two labels gives these exact retention optima:

| Target order | Maximum retained source vertices | Source losses |
|---:|---:|---:|
| 13 | 10 | 3 |
| 12 | 9 | 4 |
| 11 | 8 | 5 |
| 11 | 8 | 5 |

These are optima for this fixed source and these four fixed rooted
targets, allowing arbitrary deletions and target births. They are not
claimed for padded weak targets, other S4 source charts, or unrestricted
arc edits. Every failed larger-retention level is recorded. The search
does not enumerate only a preferred embedding or preferred deletions.

For the order-13 target, retain T vertices

    (0,2,3,4,5,6,7,8,11,12)

and map to target labels

    (12,0,2,3,4,7,8,9,10,11).

Insert target vertices 1,5,6. The exact same-count return is

    5993 → 395 → 1027 → 2483 → 5993.

The final root12 is the same physical root throughout, and its sole
incoming neighbor11 is the same physical incoming neighbor throughout.
The target rows and all four optimal maps are in `H15_ready.json`.
Its standard two-birth exit gives 17983=3·5993+4. This gives a finite
+2,+4 segment through the b=3 escape; closure at arbitrary future
targets is not inferred.

## The precise uniformity failure

Inflate B vertex1 and T vertex1 to TT_m, m≥1. The retained map above
extends by mapping the entire transitive module in order. Because all
external incidences are uniform, this is an actual retained-subgraph
replacement for every m. Both deletion charts have first two letters O,
so the indegree-one identity from DEFECT_THREE_ESCAPE gives defects
H(C3)=3 and H(S4)=5 for every m; the last right letter I makes them
strong. Complete module-run coefficients give

    H(B_m)=1134·2^m−271,
    H(T_m)=3684·2^m−1375.

Consequently

    H(T_m)−(3H(B_m)+2)=282(2^m−2).

The exact target law holds only at m=1. Its first failed implication is
that a correct finite named core replacement would have the required
dyadic slope: 3684 differs from the required 3402.

This failure was tested against every vertex placeholder in all six
normalized S4 witnesses: 13+13+15+15+14+14=84 complete coefficient
tables. None has second-run coefficient 3402 with every higher-run
coefficient zero. Hence none can be the target of the proposed
single-TT-module lift of B_m. This is a complete finite template audit,
not an obstruction to multivariate modules, changing templates with m,
or a different source family. In the earlier parameter k=m+1 it
excludes these single-module targets for 1701·2^k−811.

## Reproduction

From the worktree root, run:

    python -B 05-knowledge/results/berggren_target_return_20261111/reset/H15_inverse.py
    python -B 05-knowledge/results/berggren_target_return_20261111/reset/H15_ready.py

Both retain all checks under `-O -B`. The first uses the proved inherited
binary inverse plus independent literal positive counts, complete module
coefficients, named maps, and finite independent controls m=1,...,4.
The second uses the stored complete four-class ready inventory and
exhaustive fixed-root injective adjacency matching; it verifies every
actual deletion/birth trace and the next standard exit by subset DP.
