# The ready replacement crosses a restricted ancestry boundary

**PROVED:** the degree profiles and restricted-birth boundary below, for
every nonempty tournament module M. **VERIFIED:** exact adjacency controls
on six modules, including cyclic and weak modules. Chan–Pak remains OPEN.

## Inherited ancestry and the operation being restricted

The general unique-last-vertex, chronological decoding, and automorphism
mechanism is inherited from
[adaptive ANCESTRY, A1–A3](../../adaptive_metagraph_reformulation_20261103/ANCESTRY.md).
That larger operation adds z beaten only by an **arbitrary** old vertex p
of a strong parent. Nothing below is a new general ancestry theorem.

Here a **strict birth** additionally requires that p have indegree one
before the birth. L and K in
[BINARY_LANGUAGE](../returns/BINARY_LANGUAGE.md) are its two possible
markings; T=KL is two such births. A tournament is called primitive here
only when it has no inverse for this restricted operation. This does not
mean irreducible for arbitrary adaptive births, homogeneous substitution,
or another graph operation.

A strict birth always leaves both p and z at indegree one. If the parent
is strong, every other vertex receives an incoming arc from z and has
indegree at least two. Thus every positive-length strict history from a
strong source ends with exactly two indegree-one vertices p→z.

Conversely, suppose strong U has order greater than three and exactly two
indegree-one vertices p→z. Then U−z is strong. Every old vertex reaches p
along a simple path avoiding z, because entering z already requires p.
The vertex p beats every remaining vertex except its sole predecessor a.
Choose v outside {p,z,a}; such v exists. A simple v-to-p path avoids z and
ends a→p, so p→v followed by its prefix reaches a without z. Hence every
remaining vertex reaches p and is reachable from p. Deleting z is a legal
strict inverse.

Together with the inherited fact that a strong tournament of order
more than three has at most two indegree-one vertices, this gives a scoped
classification: among strong tournaments admitting an indegree-one mark,
the strict primitive cases are C3 and those with exactly one indegree-one
vertex. Strong tournaments with none cannot be marked in this language.
The usual adaptive inverse may still exist in a strict primitive case.

## Exact profiles of the source, raw excursion, and returned source

Use the source and return from [MODULE_RETURN](MODULE_RETURN.md). Write
m=|M|≥1 and let d_M(v) be a module vertex's internal indegree. The source
S(M) is obtained from

    P=(78,12,24,16,35,79,30)

by substituting M at old vertex1 and adjoining z beaten only by old4.
Its complete indegrees are:

| Vertex or module | Indegree in S(M) |
|---|---:|
| old0 | 3 |
| v in M | 5+d_M(v) |
| old2 | m+4 |
| old3 | m+5 |
| old4 | 3 |
| old5 | 2 |
| old6 | 3 |
| z | 1 |

Thus z is the unique indegree-one vertex, and its sole predecessor old4
has indegree3. The high-defect two-arc partner is also strictly primitive:
its sole indegree-one vertex is the same z, whose predecessor old0 has
indegree2.

Two ordinary strict births x then y give the raw successor R(M). The first
leaves z,x at indegree one; the second leaves exactly x,y at indegree one,
with x→y. Old z now has indegree two. This is inside the strict language,
and its newest y has the canonical two-root ancestry signature.

For the returned Q(M), use the funnel names: triangle0,1,2; module M;
zero-tail vertices t1,t2,t3; final tail vertex w; funnel root x; final root y.
Its full degree table is:

| Vertex or module | Indegree in Q(M) |
|---|---:|
| triangle0 | 3 |
| triangle1, triangle2 | 2 each |
| v in M | 4+d_M(v) |
| t1,t2,t3 | m+5,m+6,m+7 |
| w | m+7 |
| x | m+3 |
| y | 1 |

Therefore y is the unique indegree-one vertex and its same physical sole
predecessor x now has indegree m+3. The full-target-preserving replacement
changes the ancestry profile from

    (number of indegree-one vertices, predecessor indegree)
               (2,1)  to  (1,m+3).

It preserves every arc incident to y, yet changes the predecessor's other
incoming arcs. The distinction is structural, not just a scalar defect
change. Both S(M) and Q(M) lie outside the positive-length strict-birth
image, for all internal orientations of M.

They are nevertheless ordinary adaptive births from strong parents:
S(M)=B_old4(P[M]) and Q(M)=B_x(Q(M)−y). The input port degrees are3 and
m+3 respectively. Thus they are **not** irreducible under the unrestricted
adaptive ancestry of the inherited theorem. The replacement changes which
restricted ancestry types are available; it does not leave all adaptive
birth constructions.

## Why this boundary is necessary at the ready targets

Credit the incoming
[two-root low-readiness theorem, R3](../../mutation_collatz_chanpak_20261112/returns/PROOF.md).
For a strong tournament with two indegree-one vertices, the positive
ready root can have defect1 only at H=1+2^j (j≥2), and defect2 only at H=9;
the other root has defect zero. This is an inherited classification,
not a new enumeration here.

For k≥1 the first source count 45·2^k+7 is not1+2^j, because its difference
from1 is divisible by3. The returned count 135·2^k+23 is not1+2^j either:
for k=1 its difference from1 is292=4·73; for k≥2 that difference is2 modulo4
and larger than2. Both counts exceed9.

Consequently **every** strong representative at either count carrying an
indegree-one root of defect1 or2 has a unique indegree-one vertex. It is
strict primitive. This includes the first targets97 and293, and holds at
all transitive heights. A positive-length L/K program cannot end at such
a ready source, even if its initial strong source is chosen differently.
This excludes a restricted operation language, not the existence of the
ready sources: the physical replacement constructs them explicitly.

The uniform module theorem still has defect H(M), so the ready1/2
specialization above must not be extended to arbitrary M merely by its
size. The structural degree tables, however, hold for every nonempty M.

## Checks

```
python -B 05-knowledge/results/cycle_opening_six_return_20261112/compiler/ancestry_boundary.py
```

The receipt checks both complete degree tables, the raw two-root profile,
the high-side predecessor degree, and strongness of all outputs and both
ordinary adaptive parents on TT1, TT2, C3, C4, the H9 five-core, and the
weak module C3 followed by a singleton. These are cheap orientation
checks; the already audited full-H identities are not recomputed here.

## Exact inventory refinement of the restricted grammar

Let u_n count ordinary strong tournaments of order n with exactly one
indegree-one vertex, and let t_n count those with exactly two. The unique
vertex in a u_n object is a distinguished root fixed by every automorphism.
In a t_n object, the two possible ports are individually distinguished
by their edge p→z, so they form two different root orbits. The inherited
unique adaptive inverse prevents any double counting of the next births.
Consequently

    t_4=1,
    t_n=u_(n−1)+2t_(n−1),       n≥5.

The first equation uses the one root orbit of C3. Repeated strict peeling
ends at C3 or a unique-indegree-one core. This is a complete restricted
ancestry grammar, with its as-yet-unclassified primitive inventory retained.
It does not classify arbitrary adaptive irreducibles or all strong
 tournaments with no indegree-one vertex.

Write U(x)=sum_(n≥4) u_n x^n, and let S(x) count all ordinary strong
 tournaments admitting an indegree-one vertex, including C3. Then

    S(x) = (1−x)/(1−2x) · (x^3+U(x)).

Indeed each primitive core has one ordinary zero-birth representative and
2^(c−1) ordinary descendants at c≥1 strict births. For a unique-root core,
all its automorphisms already fix that root and extend over every history;
its descendants retain Aut(core). For C3, the first selected root kills
the cyclic automorphism group, so all positive-depth descendants are rigid.
These statements also follow by restricting the inherited ancestry orbit
formula; they are not new general orbit theorems.

The exact ordinary bank through order8 gives:

| n | u_n: one indegree-one vertex | t_n: two indegree-one vertices |
|---|---:|---:|
| 4 | 0 | 1 |
| 5 | 3 | 2 |
| 6 | 19 | 7 |
| 7 | 170 | 33 |
| 8 | 2209 | 236 |

Thus the recurrence also gives t_9=2209+2·236=2681, without enumerating all
order9 tournaments. This is a proved structural slice, not a new full
A000568 coefficient. The receipt checks all279 two-root strong objects
in the complete7412-object ordinary bank through order8: every claimed
strict inverse is strong, and every displayed recurrence holds.
