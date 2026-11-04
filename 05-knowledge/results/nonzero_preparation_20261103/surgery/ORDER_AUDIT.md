# Independent audit of the retained-order C2 obstruction

**PASS — full analytic and source audit, complete normal/optimized replays,
and an independent ordinary-core and order11 reconstruction.** This audits
[ORDER_OBSTRUCTION.md](ORDER_OBSTRUCTION.md). Its finite classifications and
the cited strong-order minimum imply the stated all-order obstruction; they
do not classify arbitrary tournament operations or the global H spectrum.

## 1. Complete universes, with all scalar cases retained

I recovered [READY §§1–3](../../two_entry_funnels_20261103/READY.md) and
[TWO_ENTRY §§1–3](../../two_entry_funnels_20261103/TWO_ENTRY.md), including their
actual SCC inverse, the unique-start entry-set law, and the exceptional
indegree-two language. Converse carries these end-gap classifications to
the start-gap auxiliaries used here while preserving H and strongness.

For scalar gap2, the universe consists of every indegree-one newborn on a
strong two-start rooted parent, plus five exceptional rooted classes at each
relevant order. The proof and source retain both parts. They do not assume
the newer joint-preparation theorem to discard the exceptional classes.
For gap1 the universe is exactly an indegree-one newborn on every strong
one-start parent. This is why the order11 calculation has987 classes.

In a two-start parent, deletion of the root gives a strong first core and
only singleton later components. Entry is either a single two-start vertex,
or two unique-start vertices in a unique-cycle core. In both cases the
entry set is proper; an empty tail is legal, and a nonempty tail is legal
exactly when its last bit is1. The source enumerates every one of these
choices. At target order n, the core order m and tail length obey
`length=n-m-2`. This exact size identity fixes the finite universe; no
arbitrary depth or port cutoff is present.

The generating series independently give

```
U=x^3/(1-3x+x^2),
P=x^3(1-x^2)/(1-3x+x^2),
V=x^4(1-x)(1-x^2)/(1-3x+x^2)^2.
```

I independently expanded them before reading the producer receipt. At
orders14,15,16 the degree-one ready counts are44,394;127,171;361,595.
The five exceptional cases give44,399;127,176;361,600 scalar rooted classes.
These are class counts, not labelled graph counts or distinct H values.

## 2. Canonicalization and independent ordinary-core coverage

The one-start unranking follows a unique rooted Hamilton path. Appending
the forced tail and prepending the new root makes the labels consecutive
on the unique directed Hamilton cycle, including the empty-tail boundary.
An ordinary isomorphism must carry that unique directed cycle to itself,
so taking the minimum over cyclic rotations is the complete ordinary
quotient. No reflection or converse quotient is silently imposed.

Every unique-cycle class has a unique-start vertex, so forgetting the root
of the complete one-start inventory loses no ordinary core. At orders at
least4 these cores are rigid; every unordered pair from the proved entry
set is therefore needed. At order3 the three pairs form exactly one orbit,
so replacing the prototype's three rotations by one representative is the
correct correction from258,902 to258,900 core-entry rows.

For two-start rooted objects, the two actual starting Hamilton paths give
the complete rooted canonicalization: any rooted isomorphism carries that
two-element path set to itself. A tournament has no nontrivial involution,
so its rooted action on those two paths is trivial. The source transports
both actual words through every funnel and records the minimum adjacency
code. The words remain directed by construction: new root, chosen core
path, forced tail. This is an actual graph/path quotient, not a quotient
by a seven-coordinate response or by H.

I separately generated all ordinary unique-cycle cores through order14
using the classical maximum-outdegree/unique-newest insertion grammar,
recovered from the cited-input section of the inherited
[operation-closed synthesis](../../operation_closed_synthesis_20261029/THEORY.md#ii-a-classical-grammar-with-a-new-finite-observer-experiment).
This route starts at the cyclic triple and inserts a newborn immediately
after its allowed parent on the cycle. It uses direct bit rotations,
rejects duplicate canonical children, and compares the entire resulting
sets with the producer's different root-deletion unranking. Every order
agrees. At order14 there are17,711 ordinary cores and39,603 marked entry-pair
classes. Complete pair-set digests at every order are recorded in
[ORDER_AUDIT.json](ORDER_AUDIT.json).

## 3. Exact evaluation and the independent order11 fibre

I read every line of both Python producers and the C++ counter, as well as
the frozen ready-transfer dependency. The C++ input is each literal strong
core with its actual one/two-entry set. It validates tournament adjacency,
strongness, entry types and actual starting budget2. Its forward subset
endpoint recurrence and the recurrence on the converse compute precisely
`(h,p2,e,F,G,b,C)`. The graph is the source of these numbers; the generator
never supplies guessed responses.

The largest core has14 vertices. Each one-path table entry is at most14!;
a spanning ordered two-path cover is a permutation with a cut, so its total
is at most15!. Entry-restricted responses are subsets of these covers.
Thus every nonnegative product and partial sum fits uint64, before the
post-computation checks. The C++ always-active rejection controls cover
invalid adjacency and the wrong actual starting budget.

The all-tail formula is the proved FUNNELS F3–F4 expression already checked
in the first surgery audit. Every legal word of the exact required length
is evaluated. The Python producer checks34 complete core-response rows by
a separate arbitrary-integer implementation, then constructs each attaining
minimum graph and directly counts its H and gap. The five exceptional values
are evaluated by their inherited all-size formulas. No selected sample
replaces a class in the minimum computation.

For the especially important order11 absence I used a separate full path:
from the independently generated ordinary order10 cores, mark every actual
unique-start vertex, producing exactly987 distinct rooted objects. A fresh
backward endpoint DP checks all unique-start sets and computes every newborn
value `h+F-1`. The entire histogram matches the producer, not merely its
missing1429 fibre:

```
SHA256 0b1ab5c2a805f91008fab5787a25d9106e7e61eadff05e73a64b5b2596f21537
minimum423, maximum9637, no H1429.
```

This matters because the rank has values on both sides of1429; a lower-bound
argument alone would not prove the empty fibre. The already audited actual
order12 H1429 target supplies the positive retained-order boundary.

## 4. Accepted minima and the all-order consequence

The complete scalar minima are4841,9545,18905 at orders14,15,16. Each is
attained by the recorded C3/two-entry/tail construction, and its auxiliary
has the retained root as sole exit. The exceptional scalar minima are
5119,10239,20479, so they cannot lower these three minima. In particular no
strong scalar gap2-ready target of these orders has H4289.

The cited primary Busch theorem was already directly recovered in this
session's induced-extension proof. The relevant exact strong-order minima
are `m(3q)=3*5^(q-1)`, `m(3q+1)=5^q`, and
`m(3q+2)=9*5^(q-1)`. They increase strictly, with `m(17)=5625>4289`.
The complete three-rank calculation therefore pays every possible order
at least14. This proves an all-order statement about the declared property,
not an extrapolation from three sampled ranks.

In particular the sharp global minimum over **all** strong scalar gap2-ready
tournaments of order at least14 is4841, attained at order14. The other two
explicit ranks have larger minima, and every order at least17 has
H>=5625>4841 by the cited bound. This corollary requires no further search
or monotonicity assertion about the ready-class minima themselves.

For any eleven-name strong H1429 source, an H-neutral strong preparation
for the +2 C2 digit cannot have order11, by the complete gap1 fibre. It must
retain at least12 names. C2 adds two permanent names, giving an H4289 child
with at least14 names. Any next H-neutral strong preparation retaining all
those names is prohibited from having scalar gap2 by the theorem just
proved. Therefore the digit sequence(+2,+4) is impossible in this grammar,
even with arbitrary old-arc changes, arbitrary choices of auxiliary/root,
and arbitrary changes in the root-response defect.

The argument does not require exact F-tripling and so is stronger than a
failure of the joint guard. It still concerns H-neutral preparations followed
by the specified two-permanent-birth C2 operation. The actual weak prefix
branch in CHILDREN is a sharp category boundary. Direct non-C2 operations,
nonneutral preparations, and deletion of permanent old names remain outside
the exclusion. A global Hamilton-path spectrum classification is not implied.

## 5. Reproduction and exact pins

Run from the session worktree:

```
python3 -B 04-computation/nonzero_preparation_20261103/surgery_orders.py --output /tmp/nonzero_orders_audit.json
python3 -B -O 04-computation/nonzero_preparation_20261103/surgery_orders.py --output /tmp/nonzero_orders_audit_optimized.json
python3 -B 04-computation/nonzero_preparation_20261103/surgery_order11.py --output /tmp/nonzero_order11_audit.json
python3 -B -O 04-computation/nonzero_preparation_20261103/surgery_order11.py --output /tmp/nonzero_order11_audit_optimized.json
python3 -B 04-computation/nonzero_preparation_20261103/surgery_orders_audit.py > /tmp/nonzero_orders_inventory_audit.json
```

The full producer modes reproduce the frozen receipts and stdout byte for
byte:5,564,616 rank gates and4,580 order11 gates. The rank replay processes
all258,900 input rows; its full input digest and full core-response digest
also match. The independent75,804-check inventory/histogram reconstruction
reproduces ORDER_AUDIT.json. Its source is4,095 bytes, SHA256
`3cb28ae09ddcbdb365afd74bcbc4c0f76e6eed9584d0795304f52aea0e5e3da7`;
its2,476-byte receipt has SHA256
`b2af222da835c60c8b6c477e8f28a5dc40b6f13556f85117eebf6d997cb35988`.

|Producer artifact|Bytes|SHA256|
|---|---:|---|
|surgery_orders.py|10212|831f652ed3370ba0ce599592793be1a4395b8c63bc5958644067fb6cdaca040f|
|surgery_orders.cpp|2163|f9bc4f7d47c85236669c73ae6edda3252bdd941ce731a4c51c1c8082364568e0|
|surgery_order11.py|2764|f8564a63192e417d23deea1005180e65fdadb41447ea893c2ec62fce571f397e|
|orders.json|7470|3be3db0232f596fd768616d7d2928cb9a21ee1580cd76f648a5122083702a17e|
|orders.out|225|fbb1c1d299b51bc7fb691256928706d9f19c8c2a98349caab267bf552cf1294b|
|order11.json|2045|e337679fbac124ea947ff67b58acd8c7e8ca29943c198615357425801a578dbb|
|order11.out|152|90f79a38c64567810ec324dc6d8cc6ea0d94d9cec94796eb3f809a288a34142e|
|ORDER_OBSTRUCTION.md, audited prepromotion version|12491|d712f56044ccc42e395b556327e17214b07c7f76a54a2acae11efd284f038a3c|

All ORDER_MANIFEST pins match. No correction remains. Status promotion and
integration belong to the parent; producer mathematical files were preserved.
