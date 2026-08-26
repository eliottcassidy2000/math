---
id: THM-4207
title: "Two-newcomer sharp depth transition, base-surplus composition, and variable-pool chart number"
status: >
  PROVED ANALYTIC LEMMAS + VERIFIED-EXACT + INDEPENDENTLY AUDITED FIXED-POOL
  RESULT. For the displayed thirty-label pool and newcomers (50,51), the
  joint repair decks satisfy tau(E_8)=8 but tau(E_9)>9. Consequently every
  one of the C(30,9)=14,307,150 eleven-label bodies B union {50,51} is
  1/14-safe in the Haar-mass sense mu(G)>=4/63; THM-4150 turns each into a
  genuine LRC(14) family after common scaling, doubling, and adjoining any
  distinct positive odd-tail pair. This is a sharp depth-eight to depth-nine
  transition for this fixed repair-deck certificate, not an arbitrary
  two-newcomer theorem or full LRC(14). A general t-newcomer Bonferroni core
  and a general chart-number lemma are proved analytically.
source: codex-lrc14-two-newcomer-20260826
depends_on:
  - THM-4150-safe-set-haar-measure-universal-odd-tail-lrc14-transfer
  - THM-4191-complete-full-pool-newcomer-haar-transfer
related:
  - THM-4174-six-deletion-completion-of-divisor-complete-newcomer-haar-transfer
  - THM-4188-all-newcomer-zero-original-anchor-hierarchy-and-resonance-filtration
  - THM-4203-fixed-pool-seventeen-body-depth-eight-haar-completion
  - THM-4214-two-newcomer-pascal-complete-eleven-body-haar-charts
script: 04-computation/lrc14_two_newcomer_depth_eight_obstruction_thm4207.cpp
output: 05-knowledge/results/lrc14_two_newcomer_depth_eight_obstruction_thm4207.out
independent_audit_script: 04-computation/lrc14_two_newcomer_depth_eight_obstruction_independent_audit_thm4207.cpp
independent_audit_output: 05-knowledge/results/lrc14_two_newcomer_depth_eight_obstruction_independent_audit_thm4207.out
depth_nine_script: 04-computation/lrc14_two_newcomer_depth_nine_closure_thm4207.cpp
depth_nine_output: 05-knowledge/results/lrc14_two_newcomer_depth_nine_closure_thm4207.out
depth_nine_independent_audit_script: 04-computation/lrc14_two_newcomer_depth_nine_closure_independent_audit_thm4207.cpp
depth_nine_independent_audit_output: 05-knowledge/results/lrc14_two_newcomer_depth_nine_closure_independent_audit_thm4207.out
composition_reproducer_script: 04-computation/lrc14_two_newcomer_base_surplus_and_variable_pool_chart_thm4207.py
composition_reproducer_output: 05-knowledge/results/lrc14_two_newcomer_base_surplus_and_variable_pool_chart_thm4207.out
script_sha256: 325642dfd01d3fe12ce79e8f6e5d03128321af22808009d32fd117515a77dc61
output_sha256: 11891ac9166c269563995aa27d2799005221d5201a5f35714790cdf02da3ba16
independent_audit_script_sha256: 1a9f59950b3c73ad9519d1ff9ed59df294d40879dfcb67e108d73e8ea467b965
independent_audit_output_sha256: d8961fbd46a897a3ec5136938b4a3097a1d258cde25df407f516befe15c74fd9
depth_nine_script_sha256: c436c704395997b4e0a274c6a67e2263d76bbec51fba35855fd4c1af95ea2a69
depth_nine_output_sha256: 248bfc0ab786a8e90a4083e0685c75eca5133a911492da1d26c3639cf5812a69
depth_nine_independent_audit_script_sha256: 981e378756dcdd563cb15015ea6e59e2e725cfec838c9680a7bb76114135522e
depth_nine_independent_audit_output_sha256: 0e2a7e349358f14a136eae452d28eb485631d8fd233bc4bd5e344d8cb3980cef
composition_reproducer_script_sha256: 6c51acd7259cfd6629755c6450ce72b2f214fd22570af4d4c7cbcd10eddd3cf6
composition_reproducer_output_sha256: ec2ca5376e177c11152e2ed87b043e6fecc9f280721c492268f67fb748409241
hash_basis: raw LF bytes
primary_audit: >
  PASS. Literal joint-wall integration gives the complete depth-six through
  depth-eight staircase and exact depth-eight transversal number. A second
  exact pass enumerates all fourteen million nine-bodies; 436 cover E8, and
  complement-local zeta transforms produce 35 explicit lawful E9 separators
  which resolve all exceptions. GCC O0/O2 outputs match, and Clang UBSan
  reports no diagnostics.
independent_audit: >
  ACCEPT. The depth-eight audit uses a separate fixed-cell / safe-comb
  interval-overlay implementation. The standalone depth-nine audit instead
  uses an endpoint-event overlay and combinadic reverse-incidence masses; it
  rebuilds the full E8 and E9 decks, again finds exactly 436 E8 covers, and
  verifies the frozen 35-edge E9 certificate. GCC and Clang outputs match,
  and Clang UBSan reports no diagnostics.
---

# THM-4207 -- two-newcomer sharp depth transition, base-surplus composition, and variable-pool chart number

**PROVED ANALYTIC LEMMAS + VERIFIED-EXACT + INDEPENDENTLY AUDITED FIXED-POOL
RESULT; LRC(14) REMAINS OPEN.**

## 1. Statement and inheritance pass

Retain the fixed pool

```text
P={8,10,15,16,20,30,40,42,60,63,80,84,85,88,95,
   120,126,132,143,145,168,170,176,190,193,240,252,
   264,286,290}.                                        (1)
```

For a finite positive label set `S`, put

```text
G_S={y in R/Z:min_(s in S)||sy||>=1/14},
alpha=4/63.                                              (2)
```

For distinct newcomers `q,r notin P`, define the joint deletion deck

```text
E_d(q,r)={R in binom(P,d):
          mu(G_((P\R) union {q,r}))>=alpha}.             (3)
```

> **Fixed-pool sharp-transition theorem.** At `(q,r)=(50,51)`,
>
> ```text
> tau(E_8(50,51))=8,
> tau(E_9(50,51))>9.                                    (4)
> ```
>
> Consequently, for every `B in binom(P,9)`,
>
> ```text
> mu(G_(B union {50,51}))>=4/63.                        (5)
> ```
>
> Hence, for every positive integer `c` and every two distinct positive odd
> integers `a,b`, some `x in R/Z` satisfies
>
> ```text
> min_(v in 2c(B union {50,51}) union {a,b})||vx||
>   >=1/14.                                             (5a)
> ```

Thus `(5)` is a finite exact family of

```text
binom(30,9)=14,307,150                                  (6)
```

safe eleven-label body cores. Haar invariance gives the same mass bound for
`c(B union {50,51})`, and THM-4150 proves `(5a)`. Thus these are genuine
LRC(14) families with thirteen nonzero relative speeds. Equation `(4)` is
sharp in deletion depth for this global repair-deck certificate: depth eight
has a nine-body transversal, whereas depth nine has none. It does **not** say
that an E8-transversal body is unsafe, and it does not provide arbitrary
entry into the fixed-pair family.

The closest proved mechanism is THM-4191's complete full-pool one-newcomer
repair transfer. The canonical hostile is now the comparable pair `(50,51)`,
where ordinary one-newcomer deck intersection fails already at depth four and
the joint depth-eight deck has an eight-cover. The corrected near miss is
again exact: certificate failure is not body danger. The least-used sidecar is
the base mass retained before adjoining newcomers; Section 6 turns it into a
general Bonferroni composition coordinate.

The live concept board was

```text
joint deletion deck | transversals | base surplus | marginal deck overlap
replacement-pool chart | complement-local zeta oracle. (7)
```

The depth-eight hostile and the base-surplus deficit forced the pull to the
next deletion boundary. There, the exceptional E8 transversals form a tiny
structured residual class that the complement-local E9 oracle closes exactly.

## 2. Why `tau(E_9)>9` proves every body

If `B in binom(P,9)` is not a transversal of `E_9(50,51)`, choose a disjoint
`R in E_9(50,51)`. Then

```text
B union {50,51} subset (P\R) union {50,51},
G_((P\R) union {50,51}) subset G_(B union {50,51}).     (8)
```

Safe-set monotonicity and `(3)` prove `(5)`. Conversely, `tau(E_9)>9` means
that every nine-set is not a transversal. A hypothetical smaller transversal
could be extended inside `P` to a nine-transversal, so exhausting nine-sets
also excludes all smaller transversals.

The exact connection contract is

```text
source:       lawful joint repair R subset P
target:       nine-body B subset P with newcomers fixed
map:          disjoint incidence R intersect B=empty
preserved:    both newcomers, labelled repair, exact Haar threshold
destroyed:    safe phase and wall address after integration
sidecar:      full failure-mask atom mass
decisive test: tau(E_9(50,51))>9.                       (9)
```

## 3. Exact depth-eight obstruction

The primary literal joint-wall implementation and the independent
fixed-cell/comb implementation agree on the full staircase

| depth `d` | `|E_d(50,51)|` | `tau(E_d(50,51))` |
|---:|---:|---:|
| `0,1,2,3,4,5` | `0` | not used |
| `6` | `39` | `2` |
| `7` | `10,114` | `5` |
| `8` | `311,544` | `8` |

There are zero threshold equalities in every row. At depth eight, an exact
minimum transversal is

```text
C={16,88,95,143,168,193,240,290}.                     (10)
```

It meets all `311,544` edges. For the lower bound, the primary audit exhausts

```text
binom(30,7)=2,035,800                                  (11)
```

seven-sets and finds no cover. Its closest seven-set

```text
{95,143,168,193,240,252,290}                           (12)
```

misses the lawful edge

```text
{8,42,88,132,145,170,176,264},
63N-4D=1,334,721,427,452>0.                            (13)
```

The independent implementation uses a different edge ordering and performs
`42,832,658` incidence checks, versus `42,630,041` in the primary path. It
again proves no seven-cover and verifies `(10)`. Hence the equality in `(4)`
is exact, not a greedy upper bound.

## 4. Exact depth-nine closure

The primary closure enumerates every member of `binom(P,9)` against the full
depth-eight deck. Exactly

```text
436 of 14,307,150 nine-sets are E8 transversals.        (14)
```

For each exceptional target `B`, the 21-dimensional complement `P\B` is
compressed. The exact failure-atom masses are zeta-transformed on that
complement, so every nine-subset `R subset P\B` receives its exact value

```text
Delta(R)=63D*mu(G_((P\R) union {50,51}))-4D.           (15)
```

Here `D=91,205,797,082,400`. Greedy set-cover compression of the lawful local
repairs produces 35 explicit E9 separator edges. Their complete labelled list,
positive deltas, and resolution counts are frozen in the primary output. The
smallest certificate margin is still strictly positive:

```text
R_min={15,30,60,63,85,126,170,176,193},
Delta(R_min)=1,448,234,172.                             (16)
```

The certificate separates all 436 exceptional bodies.

It remains to close the other `14,306,714` bodies without silently replacing
E8 by E9. If `B` is not an E8 transversal, take `R in E8` disjoint from `B`.
Since

```text
|P\(B union R)|=30-9-8=13,                             (17)
```

choose `x` there. Deleting one more label enlarges the safe set, so

```text
R union {x} in E9,
(R union {x}) intersect B=empty.                       (18)
```

The 35-edge certificate handles the remaining 436. Thus every nine-set misses
an E9 edge, proving `tau(E9)>9` and `(5)`.

The standalone independent depth-nine path instead starts all 32 comb states
at the cyclic danger interval, sweeps endpoint events, and accumulates each
failure atom into fixed-cardinality supersets through combinadic arrays. It
rebuilds both complete decks and additionally finds

```text
|E9(50,51)|=3,159,764,  threshold equalities=0.        (18a)
```

It recursively enumerates the same `14,307,150` bodies in a different order,
finds the same 436 E8 covers, recomputes every one of the 35 E9 masses twice
(direct atom scan and ranked array), and checks every separation. Its ordered
incidence total is `39,296,089,046`, versus `1,469,767,402` in the primary
hash-ordered scan. The full E9 census is supplementary; the proof needs only
the 35 frozen lawful edges plus deletion monotonicity.

## 5. Certificate failure is not body failure

Extend the depth-eight transversal `(10)` to

```text
B_0={8,16,88,95,143,168,193,240,290}.                 (19)
```

This nine-body meets every E8 edge, so depth eight cannot certify it. Direct
joint-wall integration nevertheless gives

```text
mu(G_(B_0 union {50,51}))
  =5,181,339,311,468 / 30,401,932,360,800,
63N-4D=204,816,647,179,284>0.                          (20)
```

The independent path uses the three-times larger ambient denominator and
gets the same normalized mass and three-times the displayed delta. Equation
`(20)` is the canonical hostile to the implication “deck certificate fails,
therefore the body is unsafe.”

## 6. General `t`-newcomer base-surplus composition

Let `t>=1`, let `U` be any measurable subset of a probability space, and let
`A_1,...,A_t` be measurable. Applying the union bound to the complements
inside `U` gives

```text
mu(U intersect A_1 intersect ... intersect A_t)
 >= sum_i mu(U intersect A_i) - (t-1)mu(U).             (21)
```

Indeed,

```text
U\(intersection_i A_i)=union_i(U\A_i),                (22)
```

and the union bound proves `(21)`. For a base repair `R`, take

```text
U_R=G_(P\R),  A_i=G_(q_i),
delta_R=mu(U_R)-alpha,
sigma_i=mu(U_R intersect G_(q_i))-alpha.               (23)
```

Then

```text
mu(U_R intersect intersection_i G_(q_i))-alpha
 >= sum_i sigma_i-(t-1)delta_R.                        (24)
```

Thus the exact sufficient composition condition is

```text
sum_i sigma_i >= (t-1)delta_R.                         (25)
```

The information-theoretically sharp lower bound from only
`u=mu(U)` and the marginals `a_i=mu(U intersect A_i)` is

```text
max(0, sum_i a_i-(t-1)u).                              (26)
```

On a nonatomic space this is attainable: when the total complement mass is
at most `u`, make the complements disjoint; when it is at least `u`, arrange
complements of the prescribed masses to cover `U`. This sharpness is only for
arbitrary measurable sets with the listed marginal data. It is **not** a claim
that runner combs realize every extremizer.

For the first depth-four row in Section 7, the exact base coordinate is

```text
mu(U_R)=4,128,114,491 / 47,256,889,680,                (27)
```

and the pair lower bound from `(21)` is

```text
11,216,461 / 222,071,850,
lower bound-alpha=-320,371 / 24,674,650<0.             (28)
```

So ordinary one-newcomer membership has discarded precisely the base-surplus
coordinate needed by `(25)`.

## 7. Exact depth-four marginal-intersection hostile

Write `E_d(q)` for the corresponding one-newcomer deck. The complete
`E_4(50)` consists of the following four repairs. Every one also belongs to
`E_4(51)`, but none belongs to `E_4(50,51)`:

| `R` | `mu(G_((P\R) union {50}))-alpha` | `mu(G_((P\R) union {51}))-alpha` | joint margin |
|---|---:|---:|---:|
| `{88,95,176,193}` | `69892919/236284448400` | `1096763/103633530` | `-382639/37011975` |
| `{88,145,176,193}` | `543443/4073869800` | `3377207/328172845` | `-35566763/3281728450` |
| `{88,145,193,290}` | `751147/2715913200` | `158243/13693680` | `-81007249/8147739600` |
| `{145,168,193,290}` | `7237/35735700` | `1154897/116396280` | `-43760897/4073869800` |

For example, on the first row the three masses themselves are

```text
m_50    =15,072,080,119 / 236,284,448,400,
m_51    =23,030,009 / 310,900,590,
m_50,51 =5,901,983 / 111,035,925.                      (29)
```

Therefore

```text
E_4(50) intersect E_4(51) not_subset E_4(50,51),
E_4(50,51)=empty.                                      (30)
```

This is a minimal-depth hostile for the naive operation “intersect the two
lawful marginal decks.” It does not contradict the base-surplus lemma: the
first row's exact lower bound is negative in `(28)`.

## 8. One-label replacement pools and the chart number

Fix `p in P` and `q notin P`, and form the one-label replacement pool

```text
P_p^q=(P\{p}) union {q}.                               (31)
```

THM-4191 plus safe-set monotonicity proves that every ten-subset of `P_p^q`
is safe. If `T` does not contain `q`, apply THM-4191 to `T union {q}` and
delete `q`. If `T=C union {q}` with `|C|=9`, extend `C` to a ten-subset
`B subset P`, apply THM-4191 to `B union {q}`, and delete `B\C`. Thus

```text
mu(G_T)>=alpha for every T in binom(P_p^q,10).          (32)
```

This is the full, and only, safety transfer asserted for the replacement
pool.

The underlying chart count has a general form. Let `|P|=n`, let `k<n`, and
define the omitted-label chart

```text
X_p=binom(P\{p},k).                                    (33)
```

For `O subset P`, the charts `{X_p:p in O}` cover `binom(P,k)` if and only if
no k-set contains `O`, equivalently

```text
|O|>=k+1.                                              (34)
```

If `|O|<=k`, extend `O` to a k-set; that set lies in none of the selected
charts. If `|O|>=k+1`, every k-set omits some member of `O`. Hence the exact
minimum chart number is `k+1`. More precisely, `d<=k` selected charts leave

```text
binom(n-d,k-d)                                         (35)
```

uncovered k-sets. At `n=30,k=9`, the chart number is exactly `10`; nine
charts leave one target and ten leave none.

The scope firewall is essential. The chart map says that a nine-set `C`
omitting `p` is represented inside `P_p^q` as the ten-set `C union {q}`. It
preserves membership and the individual safe conclusion `(32)`, but destroys
the fixed-pool repair deck and the resonance data on which THM-4191 depends.
THM-4191 cannot be reapplied to `P_p^q` to adjoin a second newcomer `r`.
Therefore `(34)` is a combinatorial coverage/atlas fact, **not** a
two-newcomer safety transfer. The pair safety in `(5)` comes only from the
separate exact E9 certificate in Sections 2--4.

## 9. Exact audit and replay

The primary depth-eight program constructs one joint wall lattice, labels
each open cell by its failed fixed-pool mask, and applies submask zeta sums.
The independent program constructs the fixed pool's cells first, intersects
explicit safe-comb interval lists for 50 and 51, and overlays that intersection
on the fixed cells. Their geometries, traversal orders, and depth-eight
no-seven-cover audits are distinct.

The primary depth-nine path generates its 35 separators with complement-local
zeta transforms. The standalone independent depth-nine path hardcodes only
that labelled 35-edge certificate; its endpoint-event geometry, combinadic
reverse-incidence masses, recursive body generator, and full E8/E9 censuses
import neither maintained base implementation. Thus it does not trust the
primary mass oracle, its zeta transform, or its nine-body enumeration order.

Replay from the repository root:

```bash
mkdir -p /tmp/thm4207-replay

g++ -std=c++20 -O2 -Wall -Wextra -Werror -pedantic \
  04-computation/lrc14_two_newcomer_depth_eight_obstruction_thm4207.cpp \
  -o /tmp/thm4207-replay/d8-primary
/tmp/thm4207-replay/d8-primary

g++ -std=c++20 -O2 -Wall -Wextra -Werror -pedantic \
  04-computation/lrc14_two_newcomer_depth_eight_obstruction_independent_audit_thm4207.cpp \
  -o /tmp/thm4207-replay/d8-independent
/tmp/thm4207-replay/d8-independent

g++ -std=c++20 -O2 -Wall -Wextra -Werror -pedantic \
  04-computation/lrc14_two_newcomer_depth_nine_closure_thm4207.cpp \
  -o /tmp/thm4207-replay/d9-primary
/tmp/thm4207-replay/d9-primary

g++ -std=c++20 -O2 -Wall -Wextra -Werror -pedantic \
  04-computation/lrc14_two_newcomer_depth_nine_closure_independent_audit_thm4207.cpp \
  -o /tmp/thm4207-replay/d9-independent
/tmp/thm4207-replay/d9-independent

python3 \
  04-computation/lrc14_two_newcomer_base_surplus_and_variable_pool_chart_thm4207.py
```

Compare each stdout stream with its matching frozen output. Both C++ base
implementations also byte-match under GCC `-O0` and `-O2`. All four C++
programs pass Clang with

```text
-O1 -g -Wall -Wextra -Werror -pedantic
-fsanitize=undefined -fno-sanitize-recover=all
```

and emit the same semantic outputs with no sanitizer diagnostics. All hashes
in the frontmatter are SHA-256 of raw LF bytes.

## 10. Strict scope

This theorem proves exactly:

1. the analytic measurable-set inequalities `(21)--(26)`;
2. the exact chart lemma `(33)--(35)` and the one-replacement ten-skeleton
   consequence `(32)` of THM-4191;
3. the finite exact `(50,51)` joint-deck transition `(4)` and the
   `14,307,150` safe bodies `(5)`;
4. THM-4150's genuine fixed-pair LRC(14) family `(5a)` for every common
   positive scale and every distinct positive odd-tail pair.

It does **not** prove arbitrary pair safety, a variable-pool version of
THM-4191, safety after a second use of a replacement chart, arbitrary entry
into this family, safety for an arbitrary thirteen-speed set, or full
LRC(14).
