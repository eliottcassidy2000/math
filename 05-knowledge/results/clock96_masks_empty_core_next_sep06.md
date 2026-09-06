# Joint mod-eight shadows exclude both clock-96 signatures

**PROVED ANALYTICALLY + FINITE-EXACT + INDEPENDENTLY AUDITED.** The
independently addressed padded masks have union at most93 labels for
`g=(1,1,4,4,6,12)` and at most95 for `g=(1,3,4,4,6,12)`. Neither can
cover `Z/96Z`. These are upper bounds, not claimed exact union maxima.

Root's concurrent clock90 construction is independently checked here as
an actual full tail cover above a body-safe phase. Its full row is safe
at another time. The boundary is sharp for uniform body-phase lifting,
not an LRC(14) counterexample. No theorem ID or Git mutation is made.

## 1. Inheritance and the labelled model

The closest proved mechanism is the balanced CRT/Hunter bound in
[the hereditary gcd sieve](lrc14_recursive_gcd_empty_core_next_sep06.md).
Its complete compiler retains exactly these two signatures at the
largest seven-body clock96; their original-clock tree bounds are97 and98.
The new step retains a **common excluded residue shadow** before
minimizing the pairwise overlaps.

Relevant inherited boundaries are
[THM-3387, exact labelled sheet-cover atlas and endpoint repair](../../01-canon/theorems/THM-3387-exact-cyclic-sheet-cover-atlas-and-q2-gcd-graph.md),
[THM-3414, fixed-zero six-owner classification](../../01-canon/theorems/THM-3414-fixed-zero-six-owner-base-classification.md),
and the [clock-nine body-safe cover hostile](gcd_nine_audit_empty_core_next_sep06.md).
THM-3414 permits some covers at96, which is divisible by16 and24; it
neither imposes these signatures nor transfers its fixed centre to a
mobile phase. The corrected near miss is treating independent pair
minima as simultaneously attainable. The least-used sidecar is their
shared high-class residue support. The live board is: labelled masks,
effective orders, pair credit, quotient shadows, conditional trees,
and physical phase addresses. No external priority is asserted.

For `q|c`, the allowed padded mask is

```text
M(q;a,b)={j mod c: j mod q in {a,a+b,...,a+(k(q)-1)b} mod q},
k(q)=ceil(q/7),                    gcd(b,q)=1.             (1)
```

It has `(c/q)k(q)` labels. Actual strict-danger masks are contained in
such masks, but the starts and units in (1) may be chosen independently.
All intersections below keep the same original labels. No affine
normalization, phase transport, or divisor absorption is used: absorption
at another clock repackages the body and can change the safe phase.

## 2. One common complement forces two new edges

Name the masks `X,Y,A,B,D,E`:

| Mask | X | Y, first | Y, second | A | B | D | E |
|---|---:|---:|---:|---:|---:|---:|---:|
| Order | 96 | 96 | 32 | 24 | 24 | 16 | 8 |
| Block length | 14 | 14 | 5 | 4 | 4 | 3 | 2 |
| Capacity | 14 | 14 | 15 | 16 | 16 | 18 | 24 |

Capacity totals are102 and103. For every choice of starts and units:

- `|X intersect E|>=2`: the fourteen-term unit block meets every class
  modulo8, while E consists of two whole classes. The first signature
  also has `|Y intersect E|>=2`.
- `|X intersect D|>=1`: modulo16, their shadows occupy fourteen and
  three distinct classes.
- The second signature has `|Y intersect A|,|Y intersect B|>=1`:
  their shadows modulo8 have five and four distinct classes. Each
  compatible residue pair modulo32 and24 has one preimage modulo96.
- `|E intersect D|` is a multiple of6, and each of
  `|E intersect A|,|E intersect B|` is a multiple of4. This follows
  directly from their periods16 and24. A nonempty intersection pays
  at least the respective quantum.

Suppose now that E is disjoint from D,A,B. Their modulo8 shadows all
lie in the **same six-class complement** of E. The D shadow has three
classes and each A/B shadow has four, all distinct because the steps
are units. Pigeonhole gives a common class for D with each A/B.
Each compatible order16/order24 residue pair has
`96/lcm(16,24)=2` labelled preimages. Therefore simultaneously

```text
E avoids D,A,B  ==>  |D intersect A|>=2,
                     |D intersect B|>=2.                 (2)
```

The unconditional pair minima for these two edges are zero. Keeping
one shared forbidden shadow is exactly what strengthens them.

## 3. Explicit tree certificates

For a tree T on the six masks, Hunter's inequality is

```text
|union_i M_i| <= sum_i |M_i| - sum_(ij in T)|M_i intersect M_j|. (3)
```

At a label contained in `m>=1` masks, the induced tree has at most
`m-1` edges; unoccupied labels contribute zero. Summation proves (3),
including arbitrary higher intersections.
Zero-credit edges can connect components.

For the first signature:

| Condition | Five tree edges | Guaranteed credit | Union bound |
|---|---|---:|---:|
| E meets D | `ED,EX,EY,EA,EB` | `6+2+2+0+0=10` | 92 |
| E meets A | `EA,EX,EY,XD,AB` | `4+2+2+1+0=9` | 93 |
| E meets B | `EB,EX,EY,XD,AB` | `4+2+2+1+0=9` | 93 |
| E avoids D,A,B | `EX,EY,XD,DA,DB` | `2+2+1+2+2=9` | 93 |

For the second signature:

| Condition | Five tree edges | Guaranteed credit | Union bound |
|---|---|---:|---:|
| E meets D | `ED,EX,YA,YB,EY` | `6+2+1+1+0=10` | 93 |
| E meets A | `EA,EX,XD,YA,YB` | `4+2+1+1+1=9` | 94 |
| E meets B | `EB,EX,XD,YA,YB` | `4+2+1+1+1=9` | 94 |
| E avoids D,A,B | `EX,XD,DA,DB,YA` | `2+1+2+2+1=8` | 95 |

Every listed graph is a tree, and at least one condition holds. Thus
both unions are strictly smaller than96, proving cover nonexistence
for the padded relaxation and consequently for actual masks.

Removing the two clock96 words from the complete inherited sieve
lowers its largest possible seven-body gcd to90. The
[joint-shadow master](lrc14_joint_shadow_empty_core_next_sep06.md)
also integrates the independent clock32 exclusion and the refined
profile compiler. Lower-runner LRC and primitivity remain prerequisites
for the resulting restrictions on arbitrary full thirteen-speed rows.

## 4. Actual clock90 boundary, credited to root

At c90, let `g=(2,2,2,3,3,6)`. The following unit-one blocks partition
all labelled sheets:

| g | q | Block modulo q |
|---:|---:|---|
| 2 | 45 | `8,...,14` |
| 2 | 45 | `23,...,29` |
| 2 | 45 | `38,...,44` |
| 3 | 30 | `3,...,7` |
| 3 | 30 | `18,...,22` |
| 6 | 15 | `0,1,2` |

Modulo15, the last mask covers classes0--2; the two order30 masks
partition all lifts of classes3--7; the three order45 masks partition
all lifts of classes8--14. Total capacity is exactly90.

Root supplied the actual realization

```text
C={1,...,7},   body=90C,
T=(542,55082,25292,211773,30513,51126),
y=126/1009,    x_j=(y+j)/90.
```

Each `gcd(90,T_i)` is the displayed g_i and each effective unit is one.
The strict physical masks

```text
{j: ||T_i(126+1009j)/(90*1009)||<1/14}
```

are exactly the six blocks above in the same order. Also
`||v y||>1/14` for every `v in C`, so the body is safe at the covered
phase. The full primitive, distinct thirteen-speed row `90C union T`
is nevertheless safe at `x=1/1260`.

This is a sharp boundary for demanding a surviving lift above **every**
body-safe phase. It is not an unsafe full row or a sharp counterexample
gcd for LRC(14). Changing only the last block to `{1,2,3}` destroys this
particular cover while keeping all orders and capacities, a control for
the importance of shared phase addresses.

## 5. Verification and scope

The source-to-target map projects the six actual labelled masks to
`Z/8Z`, keeps E's occupied classes, and lifts conditional intersections
back by CRT. It preserves the full-cover predicate and restores the
common-complement information lost by individual pair minima. General
body ratios and physical phase addresses remain absent; the c90
control supplies them only for that named row. The decisive operation
is the alternative: a coarse mask either pays a large intersection,
or its avoidance forces new edges in one small complement.

The standalone [verifier](../../04-computation/clock96_masks_empty_core_next_sep06.py)
uses standard-library integers and imports no repository producer.
Its complete local universe contains all distinct affine masks at
orders8,16,24,32,96, with counts `16,64,96,256,1536`. It checks all150,016
pair cases used above and all9,216 relevant disjoint triples `(E,D,A)`.
All trees are checked as graphs, and all63 nonempty membership words
check the labelwise inequality (3). This audits the written certificate;
it is not a six-mask census or its completeness argument.

The c90 control is constructed independently both as affine blocks and
as literal strict rational danger inequalities. Full labelled equality,
body safety, full-row safety elsewhere, and the shifted-block hostile
are checked. No floating-point signs or optimization-dependent assertions
are used.

```text
python3 -B 04-computation/clock96_masks_empty_core_next_sep06.py
python3 -B -O 04-computation/clock96_masks_empty_core_next_sep06.py
```

Both modes pass **191,773 explicit gates** with byte-identical
[output](clock96_masks_empty_core_next_sep06.out).

```text
mask atlas bfe4220ccb9f7fd43f1ecd74a054252a402220c20712eebaab77be42578ab602
semantic   73edd07034335984d6bf0a1c1a9c78eb51a0d6c17b4602383c65269291440a31
source     b928341ccf3e2e5c636ea83513a498330b4440d8a90e828deef690d1cd548a40
output     345f8d93f7147b054911111784b8dd167ec4adaebcc02e253d2b5065d862652c
```

The last two hashes use raw LF bytes. **Independent analytic audits:
PASS.** Root and certificate_audit separately reconstructed all baseline
credits, the6/4 intersection quanta, the conditional2/2 CRT credit and
the spanning trees. Root also read the full source and checked the actual
c90 row. The sibling's [master acceptance](clock32_audit_empty_core_next_sep06.md)
checks the refined signature universe. No fixed-zero theorem or generic
phase assumption substitutes for the certificate above.
