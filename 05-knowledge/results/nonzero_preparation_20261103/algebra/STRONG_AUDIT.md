# Independent audit: the all-strong stabilized native compiler

**VERIFIED / ANALYTIC AND EXACT-SOURCE AUDIT PASS.** This audits the complete
argument in [STRONG_COMPILER.md](STRONG_COMPILER.md), including its sufficient-order
strong-target criterion. The parent producer owns status promotion. No claim of
an unrestricted Hamilton spectrum or an unbounded two-branch target family is
made here.

The standalone verifier was constructed for this audit rather than imported
from the primary proof's inherited compiler. Its bit-matrix native interpreter
is checked against a second adjacency-set interpreter. The parent has separately
read and replayed this verifier. The all-order conclusion follows from the proof
below; the finite bank supplies controls rather than an extrapolation.

## 1. Contract and inherited scope

The input and output are strong tournaments T,U on the same n>=3 **physical
names**. An arbitrary number of old marks may be attached. A native move protects
the consecutive arcs of a current directed Hamilton word and reverses every
other arc. Births have explicitly supplied cuts; only temporary names may be
deleted. No order-preserving or fixed native-colour restriction is imposed.

The conclusion is actual adjacency transport through strong graphs, with at
most three simultaneous temporary names. It does not assert intermediate H or
F invariance. In particular it is distinct from fixed-order black-core
connectivity, and it cannot create a prescribed scalar target that does not
already exist.

The closest inherited objects were recovered directly:

* [METAGRAPH §5](../../adaptive_metagraph_reformulation_20261103/METAGRAPH.md#5-one-temporary-vertex-compiles-any-old-arc-reversal)
  proves the generic stabilized one-pair compiler without a universal strongness
  promise.
* [NATIVE_SIX §4](../../multiport_grammar_20261103/component/NATIVE_SIX.md#4-the-inherited-legal-compiler-instantiated-on-every-discrepancy)
  supplies the pair-block word and a family-specific strong instantiation.
* The endpoint path fact cited in
  [construction/RESULT §4](../../multiport_defect_grammar_20261103/construction/RESULT.md#4-why-the-natural-auxiliary-choices-do-not-rearm-the-digit)
  applies in converse to obtain a Hamilton path starting at any prescribed
  vertex of a strong tournament. The displayed anchor exists.

A targeted current-tree search for the universal all-strong/scaffold assertion
recovered those generic and family-specific routes, not an earlier proof of
this precise statement. This is bounded inheritance recovery, **not** an
external priority claim. The C3 native-complement hostile remains applicable:
a legal complement of a strong graph need not itself be strong.

## 2. The guarded middle and its reverse

Write P=(f,...,u,v,...,l), with u not first and v not last. Let z initially beat
only u. The word Pz is directed because l!=u. After its native complement,
the word Q obtained by replacing the edge u→v by u→z→v is directed.
All other successive old edges of Q were protected by Pz.

The symmetric difference of the two protected-edge sets on the old names is
exactly {u,v}. Thus the second complement changes only that old pair. On helper
incidences the final outgoing cut is exactly {v,l}; the three names u,v,l are
distinct. The source of each incidence reversal was checked independently,
not inferred from an observer equality.

The first native state contains the Hamilton cycle Pz: its path edges remain
protected, and z→f is obtained by reversing f→z. The latter was unprotected
because f differs from both u and l. This is the complete middle-state
strongness certificate. The cut {u} at one endpoint and {v,l} at the other are
both nonempty and proper in the applications, so strong old endpoints give
strong helper extensions.

Traversing Q and then Pz backwards is valid because native complementation
protects its own guard and is an involution on that guarded state. The reverse
birth cut is {v,l}. It has the same middle Hamilton cycle. The verifier checks
both actual command directions with a separately implemented set interpreter.
No old-pair reversal is silently accepted as a primitive native operation.

## 3. Scaffold, editing, and cleanup

Choose w in the strong old graph. The first birth a loses only to w and beats
all other old names. The second birth b loses to every old name and beats a.
Each extension has a nonempty incoming and outgoing cut against a strong base,
so both are strong.

For an old Hamilton path R starting at w, P=(b,a,R) is legal in the desired
scaffold graph. The pair (a,w) is strictly interior. The reverse guarded gadget
therefore changes w→a to a→w while every intermediate is strong. The resulting
relation is exactly

    a → old core → b → a.

It makes the entire graph strong for **every nonempty old core**, including a
weak one. This explicitly pays the loss of strongness that arbitrary individual
old-pair reversals could cause inside the core.

The direct pair-block construction is also valid for a weak core. Starting
with a Hamilton path on the other old names, insert the pair before the first
vertex beaten by either endpoint, putting an endpoint that beats this vertex
last in the block. All preceding vertices beat both endpoints. If there is no
such vertex, append the pair. One of its internal orientations is therefore
legal, and the selected pair is its only possible disagreement with the current
core. Extending this word by a first and b last makes the chosen pair strictly
interior. The forward or inverse guarded gadget edits exactly that pair and
leaves all scaffold incidences unchanged.

Once the induced old relation is U, use an old Hamilton path beginning at w
and the forward gadget at (a,w) to restore w→a. Then delete b, then a. Every
remaining graph is strong. Both superficially tempting shortcuts fail:

* deleting b before restoring w→a leaves a universal source a;
* after restoration, deleting a before b leaves a universal sink b.

These are literal failed controls in the receipt. They explain why both the
opening/closing gadget and the deletion order are essential to this proof.

## 4. Budgets, names, and the strong-target criterion

There are d old discrepancies and exactly two scaffold-changing gadgets.
Each gadget uses two native complements and one helper birth/deletion. Adding
the two scaffold names gives 2d+4 complements and d+4 births and deletions.
Only a,b and the current gadget helper coexist, so the peak is n+3. A helper
slot can be reused after deletion; every original name keeps its index and
survives every command. The proof gives an explicit finite construction, not
a shortest program. In the identity case the empty program is also available.

The sufficient-order criterion in primary §6 is correct. If a strong ready
witness has m>=n vertices, first add m−n permanent vertices with proper
nonempty cuts, for example sole exit at the retained root. The enlarged source
stays strong. Relabel the target positions onto all these physical names,
keeping the marked root, and apply the compiler. Conversely, any such final
strong ready output must itself have at least n names and be a witness at the
required H. No transitive-prefix padding is used to manufacture a strong target
of adequate order. This is a criterion for reachability **given a witness**,
not a proof that the needed witnesses exist throughout a ternary orbit.

Every guard and every adjacency matrix is available, so all requested
observers can be recounted. The example deliberately has changing H and F.
The coloured endpoint grammar is mentioned only for operations within its
separately proved scope; this theorem does not claim that a bounded observer
quotient determines arbitrary compiled old-arc edits.

## 5. Exact finite universe and independent paths

The verifier enumerates all labelled tournaments on 3,4,5 vertices, with no
isomorphism quotient or discarded weak classes. Their strong counts are
2,24,544. It checks:

| Universe | Exact size |
| --- | ---: |
| Every old pair inside a scaffold over every labelled 3–5 vertex core | 10,648 |
| Setup and cleanup at every root of every strong labelled 3–5 vertex core | 2,822 |
| Every ordered strong source/target pair and every root at orders 3 and 4 | 2,316 |

The first size is 8·3+64·6+1024·10. The second is
2·3+24·4+544·5. The third is 2²·3+24²·4. These formulas independently pay the
coverage and the meaning of the labelled universe.

Across these controls and one displayed trace there are 27,860 gadgets,
including 13,770 inverse gadgets. Each checks its legal full guards, exact
single-pair restriction, reverse helper cut, Hamilton middle cycle, all endpoint
and deletion strongness predicates, and agreement of the bit-matrix and
adjacency-set native implementations. The test gates use explicit exceptions,
so optimized Python does not erase them.

The displayed source (8,1,3,6) and target (2,4,9,3), with retained root 0,
have d=4. The receipt records every graph and guard, together with independently
recounted endpoint-path observers at every stage. Both endpoints have H=5,
but their F values are 7 and 9; intermediate H values also vary. This is a
positive adjacency test and a hostile to count-invariance extrapolation.
Both wrong cleanup orders are retained as actual weak graphs.

## 6. Frozen reproducibility and verdict

Source: [strong_compiler.py](../../../../04-computation/nonzero_preparation_20261103/strong_compiler.py),
11,052 bytes, SHA256
`71c64d0caf21f22d454a7589dc2e05e87d66ab6a52378e13700e16ed0c362326`.

Receipt: [strong_compiler.json](../../../../04-computation/nonzero_preparation_20261103/strong_compiler.json),
10,168 bytes, SHA256
`b914205ecba7e9d805ece970b1bda3dc2ee5a67a059054f2add1d0ccf7803fff`.

Both normal and optimized runs reproduce the frozen receipt byte for byte:
**424,665 exact gates**. From the repository root:

```bash
python3 04-computation/nonzero_preparation_20261103/strong_compiler.py --output /tmp/strong_compiler.normal.json
python3 -O 04-computation/nonzero_preparation_20261103/strong_compiler.py --output /tmp/strong_compiler.optimized.json
cmp /tmp/strong_compiler.normal.json /tmp/strong_compiler.optimized.json
cmp /tmp/strong_compiler.normal.json 04-computation/nonzero_preparation_20261103/strong_compiler.json
```

The audited pre-promotion primary was 10,660 bytes, SHA256
`3e49173fba9a8b2c63bff0e8b4fb5b8bc64c69638405ad7433322da3eb9cfa24`.
Its pending-control footer may be replaced by these receipts and its status
promoted by the parent; this pin identifies the complete mathematical text
reviewed, not a requirement that promotion retain an obsolete footer.

**PASS:** analytic all-order native legality and strongness, physical-name
retention, budgets, sufficient-order target criterion, source coverage and
normal/optimized replay. No mathematical repair was required. Minimal helper
count, restricted native-colour languages, fixed-order connectivity, global
Hamilton spectrum coverage and unbounded two-branch target selection are not
consequences of this theorem.
