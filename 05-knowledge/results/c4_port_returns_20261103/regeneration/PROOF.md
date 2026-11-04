# Two missing digits cannot regenerate inside the whole attachment chart

**Status: PROVED from FINITE-EXACT exhaustive classifications and the CITED
Moon--Busch strong-order minimum.** This excludes specified all-clock
signatures, not Hamiltonian-path counts in arbitrary tournaments. In
particular, it does not claim that 91, 97, 217 or 223 are forbidden counts.
The global `{7,21}` spectrum conjecture and reusable two-digit grammar remain
OPEN.

## Inheritance and the changed target

The nearest positive mechanism is the literal common family in
[`TWO_DIGIT_ENTRY.md`](../../multiport_defect_grammar_20261103/transport/TWO_DIGIT_ENTRY.md):
both digits start from `H0(t)=-5+7*2^t`, with actual outputs

```
H2(t)=-13+21*2^t,       H4(t)=-11+21*2^t.
```

The C2 macro has exact auxiliary defect transport but resets the new
auxiliary defects to zero; the inherited all-marking ready obstruction
prevents assuming the digit prepares its own next use. See
[`construction/RESULT.md`](../../multiport_defect_grammar_20261103/construction/RESULT.md).
The corrected near miss is the older unit-response dial: its scalar unit
difference does not preserve the required joint cover. See
[`recovery/NOTE.md`](../../zero_ready_component_returns_20261103/recovery/NOTE.md).
The less-used sidecar here is the whole-graph strong-order minimum: it turns
a small-clock specialization into a global bound on the suffix size.

Our four live objects are the literal marked source, the seven-coordinate
path-cover signature, arbitrary suffix isomorphism classes, and the
small-clock strong-order bound. The map from a literal source to a suffix
isomorphism class preserves the target signature existence question but
forgets physical names, birth order, and native guards. That loss is harmless
for an exclusion; a positive class would still require a source-preserving
construction.

## Exact convention and strongness

Let B be a nonempty finite tournament and W a subset of its vertices. In
`M(TT_t;B,W)`, t>=1, the prefix dominates B, the new root r dominates the
prefix and B\W, and W dominates r. Let H(empty)=1. With ordered slots
allowing empty paths, define

```
h = H(B),
p = sum_S H(B[S]) H(B[V\S]),
f = sum_S sum_(u in W) E_u(B[S]) H(B[V\S]),
j = sum_S sum_(u,v in W) E_u(B[S]) S_v(B[V\S]).
```

Here E/S count paths ending/starting at the named vertex and are zero if
that vertex is absent; these forced endpoint slots are therefore nonempty.
The inherited attachment formula is

```
H(M(TT_t;B,W)) = c + 2^t f,       c=h-j.
```

The full attachment is strong iff W meets the terminal SCC of B. Indeed,
every vertex of B reaches a marked vertex in that SCC, then r, then the
prefix, which reaches all of B. Conversely, if the terminal SCC is unmarked,
it cannot leave B or reach r. The computation tests this exact condition;
it does not require B itself to be strong.

## Finite classification

`catalog.cpp` recursively extends every representative of order n-1 by
every one-vertex cut. Every order-n tournament has such a deletion, so this
construction is complete. Isomorphism keys first sort by outdegree and then
minimize the adjacency bit word over every permutation within equal-degree
blocks. All isomorphisms preserve those blocks. Equal keys are equivalent
to isomorphism, so the quotient is exact. No conjectured metagraph recurrence
or assumed list of cores enters the construction.

The observed representative counts at orders 1 through 8 are

```
1, 1, 2, 4, 12, 56, 456, 6880.
```

For each representative, all nonempty W are tested, retaining precisely
those meeting the terminal SCC. Endpoint subset DP gives the exact path
counts; disjoint-subset convolution gives f and j. There are no floating
point operations. The supplementary `small_bank.cpp` independently builds
the order-6 bank from all 32768 labelled tournaments before taking all
56*64 births to cover order7.

The complete queries are:

| Suffix orders | f | c | Hits | Receipt |
|---|---:|---|---:|---|
|1--7|63|-35, -29|0|`catalog7_f63.json`|
|1--8|126|-35, -29|0|`catalog8_f126.json`|
|1--7|63|-7|0|`catalog7_f63_c-7.json`|

Each catalog exports all representatives, rows, H, and ordered two-slot
count p. The first two queries concern missing words24 and44. The third is
the exact suffix target needed by the smallest surviving H=5 prefix escape;
its relation to arbitrary fixed prefixes is proved in
[`prefix/PROOF.md`](../prefix/PROOF.md).

## Why the finite no-hits exclude every suffix order

**CITED:** Arthur H. Busch, *A Note on the Number of Hamiltonian Paths in
Strong Tournaments*, EJC13 (2006), N3,
[Theorem1 and Corollary1, pp3--4](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v13i1n3/pdf/).
The exact minimum for strong order N is

```
m(3k)=3*5^(k-1), m(3k+1)=5^k, m(3k+2)=9*5^(k-1), k>=1.
```

In particular m(10)=125, m(11)=225, and the minimum is increasing. Repository
route: **THM-4744**,
[`interval-chains-attain-busch-and-are-the-extremal-strong-tournaments-through-order-twelve`](../../../../01-canon/theorems/THM-4744-interval-chains-attain-busch-and-are-the-extremal-strong-tournaments-through-order-twelve.md),
with the prior primary-source audit in
[`core_reduction/THEORY.md`](../../session133_fractran_spectrum_refinement_20261101/core_reduction/THEORY.md).

For f63,c=-35 or-29, specialize the valid attachment to t=1. Its count is
91 or97, its order is |B|+2, and it is strong. Thus its order is at most9,
so |B|<=7. The complete order7 query has no solution. This proves:

**No finite marked suffix (B,W) in this attachment chart realizes either
missing all-clock signature `(f,c)=(63,-35)` or `(63,-29)`.**

Consuming one clock means an output prefix `TT_(t-1)`, with t>=2. Coefficient
comparison now requires f126 and unchanged c. At output clock1, H is217 or
223. By m(11)=225, the full order is at most10 and |B|<=8. The complete
order8 query is empty. Thus **consuming exactly one clock also cannot realize
either missing all-clock return in this chart**.

For f63,c=-7, the same specialization gives H119<125 and |B|<=7. Its empty
query is likewise an all-order signature exclusion.

These conclusions remain valid if a template is intended only for large
clock values: its fixed (B,W) still defines a strong attachment at clock1.
They do not exclude clock-dependent suffixes whose signatures change, larger
clock consumption, nonhomogeneous prefix interactions, or a different
attachment chart.

## Literal-source controls and retained information

The actual five-vertex suffixes are

```
B2=[30,28,8,16,4], W2={0,1,2};
B4=[30,28,24,16,0], W4={0,2,4}.
```

`local_small.json` tests every single new-vertex cut and every final W.
It finds nine labelled22 controls from B2 and three42 controls from B4,
but no24 or44. `two_small.json` tests exactly two successive births and
every final W: eighteen22 controls from B2, no24, and no42/44 from B4 in
that exact two-birth menu. All old suffix names and arcs are retained in
these local menus. Changing W changes old arcs incident to the outer root;
a numeric hit alone does not certify a native schedule. In particular,
the zero two-birth42 count does not invalidate the one-birth42 controls.

An all-clock preparation preserving both H and the old root's endpoint
two-cover response preserves h,f,j separately. Since B2 has h3 and B4 has
h1, inserting a strong C4 as an ordinal component cannot provide that
particular neutral preparation: its Hamilton count contributes a factor5.
The whole-chart exclusions above are stronger and allow arbitrary suffix
SCCs and arbitrary port changes.

For the prefix sidecar, the catalog gives:

| Strong H | Number of classes through7 | p values | Odd parts of p |
|---:|---:|---|---|
|5|1|32|1|
|29|6|320, 580|5, 145|
|35|2|676|169|

These are finite exact outputs here; the SCC and order-bound argument
extending their use to arbitrary prefixes belongs to `prefix/PROOF.md`.

## Reproduction and next structural target

From the repository root:

```
python 05-knowledge/results/c4_port_returns_20261103/regeneration/search_small.py
python 05-knowledge/results/c4_port_returns_20261103/regeneration/catalog.py --order 7 --slope 63
python 05-knowledge/results/c4_port_returns_20261103/regeneration/catalog.py --order 8 --slope 126
python 05-knowledge/results/c4_port_returns_20261103/regeneration/catalog.py --order 7 --slope 63 --intercept -7
python 05-knowledge/results/c4_port_returns_20261103/regeneration/check.py
python -O 05-knowledge/results/c4_port_returns_20261103/regeneration/check.py
```

The inherited strong-order bound makes the current exclusions exhaustive;
enlarging the finite C4 bank at the same clock cannot repair them. A useful
next construction must change the coefficient balance, make prefix/suffix
interaction nonhomogeneous, or provide an explicit new chart with preserved
named ports. The positive two-digit entry and zero-digit C2 amplifier remain
available; no acceptance or physical identity is inferred from equal H.
