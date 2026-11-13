# A finite inverse selector, and a proved sparsity boundary

**PROVED scoped / FINITE-EXACT.** These are actual reusable paired
tournaments. This is a complete selector for the specified finite language,
not a selector for arbitrary tournament Hamilton counts.

## The smallest seed and two reusable moves

At root3 the two strong tournaments

    (2,4,9,3), (2,8,3,5)

have states (5,2,3,6) and (5,1,1,5). They are two different markings of
the same unrooted four-vertex type. They have common (h,E)=(5,9), lie on
F+s=h+q, and differ by −(0,1,2,1). They therefore seed the universal
return of [PG3](PORT_GRAMMAR.md).

In addition to its four-vertex context U, the five-vertex context

    A=(4,5,24,19,3), entry4, exit3

has matrix

    [ 5 −3 −6 15 ]
    [ 1  0 −2  5 ]
    [ 1  0 −2  6 ]
    [ 5 −3 −6 16 ].

It preserves the repair direction d=(0,1,2,1), maps every state into
F+s=h+q, and acts on targets by W below. The four-vertex U reverses d.
Both contexts and the seed are strong; gluing strong tournaments at a
common root remains strong. Thus every word in U and W gives an actual
strong pair, with repair magnitude one and an explicitly named next port.

    U: (h,E) ↦ (h+4E, 2h+7E),     det U=−1;
    W: (h,E) ↦ (2h+9E, 3h+16E),  det W=5.

The first U ray is (h,E)=(5,9),(41,73),(333,593),(2705,4817),... .
W sends (5,9) to (91,159). Order matters: UW gives (727,1295), whereas
WU gives (739,1291). Both use the same number of added vertices.

The original eight-vertex return in `port_compiler.json.gz` factors
literally as U after W. Its outer U occupies vertices(4,5,6,3), in U's
label order, of the eight-vertex context. Deleting the three non-entry
outer vertices leaves precisely the five-vertex W context above.
This checks the geometry, not just a matrix factorization. Its target
matrix is UW=[[14,73],[25,130]], explaining determinant−5 and the factor5
in the next target. U² is the other unimodular matrix in the head below.

For primitive (h,E), U preserves gcd1. W has output gcd either1 or5,
and the latter occurs exactly when h+2E=0 mod5. This follows from det5
and the kernel of W modulo5. Arithmetic primitivity is not unit readiness.

## Complete return-matrix head through eight vertices

Restrict to port contexts that send d to a nonzero multiple of d and
preserve the plane F+s=h+q. From all ordinary context types through8,
all entry vertices and all exit vertices, there are nine distinct full
four-coordinate matrices. Their target actions are:

| Smallest context order | Target matrix | Determinant | Repair multiplier |
|---:|---|---:|---:|
|4|[[1,4],[2,7]]|−1|−1|
|5|[[2,9],[3,16]]|5|1|
|7|[[9,32],[16,57]]|1|1|
|7|[[10,73],[17,122]]|−21|1|
|8|[[14,73],[25,130]]|−5|−1|
|8|[[6,75],[10,141]]|96|2|
|8|[[16,123],[27,218]]|167|5|
|8|[[20,71],[35,124]]|−5|−1|
|8|[[25,244],[46,419]]|−749|−13|

This nine-matrix head is a newly specified universe. It is not the older
“nine-family” obstruction mentioned in the research history. All matrices
have positive target entries and preserve E>=h>0. Their full matrices,
actual contexts and port labels are in `renewal_catalog.json.gz`.

The quotient (h,s,q,F)→(h,E) is sufficient for composition INSIDE this
language: on the preserved plane its kernel is precisely the repair
direction, and every rule preserves that direction nontrivially. It
forgets s, q, defect and unit readiness. Those predicates cannot be read
from this quotient and are not used as entry conditions here.

## Exact decreasing membership

For any matrix N=[[a,b],[c,d]] in the table, with determinant D≠0, a
candidate target (H,E) has only one numerical predecessor under that rule:

    ((dH−bE)/D, (−cH+aE)/D).

Reject nonintegral predecessors, nonpositive predecessors, or those not
strictly decreasing H. Recursively accept exactly when a branch reaches
(5,9). Every accepted word reconstructs an actual pair by PG1. Conversely
the last letter of every construction word supplies one of these branches.
Thus this is an iff membership test for the language, with termination.
It is not an assertion that every integer pair satisfying a congruence is
realizable outside the language.

Starting at the H5 seed, forward closure through H<=100,000,000 contains
exactly1,346 target pairs, with1,346 different H values. The independent
decreasing selector accepts every one and agrees on1,000 additional
deterministic pseudorandom controls. The longest word in this bounded
language has eight letters, counting U² as a separate available letter.
The attained global target53,084,231 is absent from this entire language.
The exclusion is complete for these words at every order, since H grows
strictly and can never return from above the cutoff.

## Why these nine rules cannot finish Chan–Pak

For the nine rules the H growth factors on E>=h are at least

    L=(5,11,41,83,87,81,139,91,269).

For s=3/4, use k=(3,6,16,27,28,27,40,29,66). Each k_i^4<=L_i^3, so

    sum L_i^(−3/4) <= sum 1/k_i
      = 3602617/4823280 < 3/4.

For a word w with H(w)<=X, its growth product L(w)<=X/5. Hence
1 <= (X/5)^(3/4)L(w)^(−3/4). Summing over all finite words gives

    #{words with H<=X} < 4 (X/5)^(3/4).

Distinct attained counts are no more numerous than words. Thus this
entire recursively entered nine-rule language has density zero among the
positive integers. Any finite initial seed set gives the same exponent,
with a changed constant. This is an elementary weighted-word-count bound,
analogous to a Kraft sum in prefix coding; it is not a universality claim
about cellular automata or Turing machines.

The bound excludes only the stated language. Arbitrary larger contexts,
other repair slopes, or extra operations inserted between returns are
not covered. Neutral padding AFTER a completed word changes no H; padding
used as a new endpoint preparation before another rule must be analyzed
separately. This is the specific missing resource: more branching relative
to target growth, or a genuinely unbounded selector, not more iterations
of the same nine matrices.

There is also a direct arithmetic precedent for this viewpoint:
[Lagarias's overview, section2](https://arxiv.org/html/2111.02635v1) frames
backward Collatz as closure under guarded affine maps and discusses
Klarner-type integer sets generated by such maps. Our two-register grammar
has a literal tournament constructor and an independently proved growth
bound; the analogy does not transfer a convergence or density conclusion
from a different affine system.
