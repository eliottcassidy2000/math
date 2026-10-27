# Independent audit: Session75 dyadic pairing and native debt

**Verdict: PASS for producer revision 2.** The elementary proofs are sound in
their stated scopes, and independent exact reconstruction confirms the native
examples and all reported debt statistics. One source-record serialization bug
was found, reported, corrected by the producer, and independently regression
checked. The native fixed-low-bit target remains **OPEN**.

Audit target: `/tmp/session75-dyadic`, frozen revision 2. The self-checklist
named `AUDIT.md` in that package was not treated as an independent audit.
This audit was authored independently by `/root/audit_dyadic_pairing_debt`.
No producer code was imported or executed for the primary audit. No shared
repository files were edited, and no publication or canon promotion is implied.

## 1. Binding and independence

Repository source pin:
`96e35642110104cc547516a63f37eebcf930f707`.
Every declared source was retrieved with targeted `git show COMMIT:PATH` and
compared byte-for-byte to the producer extraction; SHA-256 values also agree.
In addition, the frontier was independently bound to its original path
`04-computation/experiments/first_arrival_20261018/V_closed_Y1000.txt`, not only
the later paid-root copy. Both are identical.

`SOURCE_BINDINGS.json` records these bindings. The audited producer manifest
has SHA-256
`ad8d60b20884230f7bd2180e22115b16f970d4cd0aa284a1149d45cd51f321e8`.
Its complete declared file hashes were checked. A local frozen copy is in
`producer_revision2/`; reproducibility does not depend on its former `/tmp`
location or access to the repository.

The independent verifier implements the source charge with multiplication in
`Z[T]/(T^2+T+1)`, literally summing the entire B-chain. This differs from the
producer's charge lookup tables. Three separate real rows provide another
arithmetic check. Native edges are reconstructed from exact integer divisions,
including every inverse braid and the corrected affine clock identity.
Every source includes its complete B-chain and every exterior native edge in
`reconstructed/independent_sources_N*.json`.

The frontier check establishes 627 distinct states, forward closure, and strict
rank decrease away from fixed state 1. It does not assume universal entry:
every tested finite source is traced until its actual first entry; a cycle or
unresolved trace raises an exception rather than dropping the source.

## 2. Proof audit

### Perfect matchings and blocks — PASS

For labelled even `n>=2`, a fixed unordered pair occurs in a uniform perfect
matching with probability `1/(n-1)`. Therefore

`E M = C+(P-C)/(n-1) = [P+(n-2)C]/(n-1)`.

This retains repeated values and zero-valued labels. For uniform blocks of
size `b|n`, replacing that probability by `(b-1)/(n-1)` gives
`E E_b=[(b-1)P+(n-b)C]/(n-1)`. The cases `b=1,n` are correct. A shared uniform
permutation supplies actual compatible nested levels. Independent level
partitions do not supply compatibility.

The expected energy path is monotone between `P` and `C`, but convexity of the
positive part gives `E D >= D(E energies)`, not the reverse. The note retains
this direction and does not promote a mean ledger to an actual tree.

### Four-child insertion and sharp excess bounds — PASS

For any nonnegative endpoint and middle energies,

`(P-M)+ + (M-C)+ = (P-C)+ + dist(M,[min(P,C),max(P,C)])`.

Consequently inserting a level never reduces debt. The fixed-endpoint
stability estimate in the note follows from the total-variation identity,
with the stated constant 1 on interior perturbations.

The three four-child matching energies obey

`P/2 <= M_i <= 2C`, and `M_1+M_2+M_3=P+2C`.

The lower bound applies to the two pair sums; the upper bound applies
separately within the pairs. Feasibility already implies `P<=4C`.
The deterministic and expected excess bounds are correct:

`min excess <= min(P/2, max(C-P/2,0))`,

`mean excess <= (P+C)/3` for `P<=C`,
`<=2(2C-P)/3` for `C<=P<=2C`, and `=0` for `P>=2C`.

For the expected bound, the feasible middle-energy triple is contained in the
simplex with vertices given by permutations of `(2C,P/2,P/2)`. A convex sum of
interval distances is bounded by its values on those vertices. For the
deterministic bound, unless an energy is already in the endpoint interval,
the mean identity forces a below-interval energy; its lower bound gives the
claimed estimate. In the last regime all energies lie in the interval.

The sharp family `x1=x2=(p+a)/4`, `x3=x4=(p-a)/4` realizes the vertex triple
and attains both estimates, including the zero and boundary cases. Every
admissible real endpoint pair `C>=P/4` is realizable by norms of `p,a`.

The exact Eisenstein witness `(2-omega,2-omega,omega,omega)` has
`P=C=16`, middle energies `32,8,8`, and debts `16,8,8`. Thus averaging the
linear energies erases debt only from the averaged path: no actual pairing
has zero debt. Its eight labelled unit-source realization has original source
norm 8, not child-bin energy 16. The note correctly distinguishes them.

### Global selected-tree theorem — PASS, changed-tree scope

At a bottom-up stage with `m` bins and energy `E`, the minimizing perfect
matching produces energy `E' <= [P+(m-2)E]/(m-1)`. The debt edge in the
root-to-leaf orientation is `(E'-E)+`, hence

`(E'-E)+ <= (P-E)+/(m-1) <= P/m`,

using the Hilbert-space inequality `E>=P/m`. Summing over
`m=n,n/2,...,2` yields `(1-1/n)P`. The induction direction and geometric sum
are correct. Actual unions of predecessor blocks certify compatibility and
retain each labelled source once. The `n=1` case has no edges. Equal charges
give equality and establish sharpness.

Terminal coordinates must remain orthogonal. Pairing a vector supported in
terminal A with one in terminal B creates a block with two terminal
coordinates; it does not introduce their cross term. The independent tests
include this hostile case: charges `+1` and `-1` in different terminals have
`P=2`, not `P=0`. Alternatively, separately selected terminal trees can be
padded to a common power of two and their whole-level debts summed from above.
The integer `n` in the displayed bound is the power-of-two leaf count,
**including padded empty bins**, not an unpadded source count.

The theorem changes the grouping according to actual charges. Native source
weights, clocks, ports, identities and terminal metadata can all be retained,
but the tree is not asserted to be the native low-bit tree or a native dynamic
reordering. Its right side `P` is terminal direct-sum energy, not original
source norm. It therefore does not establish `D_native<=S_N/2` or an upper
bound on terminal energy from `S_N`.

The alternating-unit fixed-low-bit example is also correct: for `n=2^h`,
`E1=0`, `E_(2^k)=n^2/2^k` for `k>=1`, and
`D=n^2/2-n` while `S=n`. Its alternative opposite-sign pairing has zero debt.
This is a generic labelled-source obstruction, not a native `A_N` example.

## 3. Native source audit

The exact finite universe is all positive odd `s<4N` with `s!=5 mod8`.
Zero-charge labels are retained. All 768, 1536, and 6144 records at
`N=512,1024,4096` respectively agree field-by-field and byte-for-byte with the
producer's canonical JSON: source label, address, original charge, signed
charge, first-entry terminal, exterior port, complete clock, and native depth.

For every declared cutoff, the independent audit computes both raw bin
energy differences and a separate direct enumeration of every same-terminal
labelled pair. Its owner is exactly
`q=2^v2(m_s-m_t)`, where `m_s=(s-1)/2`, equivalently
`2^v2((s-t)/2)`. Each pair is charged exactly once; no cross-terminal pair is
included. Parent-bin polarization provides a third equality check.

Every residue bin exists conceptually. Occupied zero-sum bins and genuinely
empty bins are distinguished in the receipt. Empty bins are arithmetic zeros,
not dropped source labels. At partial cutoffs, `S_N` retains all sources,
including those not yet arrived, and is distinct from `S_arrived`.

### Recovered 433 versus 79 — PASS

At N512 complete arrival (maximum native depth 42):

| Statistic | Exact value |
|---|---:|
| Original labels | 768 |
| Full original source norm | 1244 |
| Terminal direct-sum energy | 1313 |
| Positive labelled-pair contributions | 433 |
| Positive parent-bin interactions | 293 |
| Positive part after each whole-level sum | 79 |

Raw energies are
`1313,1274,1268,1261,1255,1242,1252,1247,1244,1244,1244`.
Their differences are `39,6,7,6,13,-10,5,3,0,0`; taking positive parts only
after these whole-level sums gives 79. All three positive quantities are
different, with the correct ordering. Numerical recovery does not identify
the historical computation cited by the user; that provenance remains unknown.

### Three-source native witnesses — PASS

For `(N,q,u,r)=(512,1,435,0)`, the four modulus-4 children contain exactly
`[1833], [435], [], [1375]` and have charges
`omega,omega,0,-1-omega`. The full ledger contains no additional contributor
to this terminal. The sources have `(port,clock,depth)` respectively
`(4643,10,5)`, `(0,0,0)`, and `(4643,8,4)`.
The endpoints are `P=C=3`, middles `5,2,2`, and debts `2,1,1`.

For `(N,q,u,r)=(4096,4,443,0)`, the four modulus-16 children contain exactly
`[12609], [13993], [9457], []` and have charges
`omega,-1-omega,-1-omega,0`. All three exterior ports are 9457; clocks/depths
are `8/2`, `24/12`, `6/1`. The endpoints are again 3, middles `2,2,5`.
The witness is present at depth 12 and persists at later cutoffs because
there are no additional source labels in this parent at complete arrival.

The native middle always uses `(0,2)|(1,3)`, not `(0,1)|(2,3)`. Three nonzero
children are minimal: with at most two, any middle is `P` or `C` and hence
creates no excess. Neither local witness implies a global lower bound, since
other parents can cancel its level contribution.

At N512 and whole-level `q=16`, the independently recovered endpoints are
`P=1255,C=1252` and same-pattern global middles `1251,1242,1266`, as claimed.
This rules out zero excess for that three-member correlated ensemble, not for
all independent parentwise choices. The note's stated qualification is needed
and present.

### Marked observer pair — PASS after serialization repair

At N512:

| source | A=v | terminal | port | clock | depth | m |
|---:|---|---:|---:|---:|---:|---:|
|799|`-2-omega`|799|0|0|0|399|
|1065|`omega`|799|1065|2|1|532|

The exact native edge is `1065 -> 799`, clock 2. Their address difference is
133, so their pair belongs to `q=1`. Norm polarization is 0, but marked rows
`c=2a-b` are `-3,-1`, giving cross term `2*(-3)*(-1)=6`.

**Defect found in revision 1:** producer variable `r`, initially source 799's
record, was reused for a parent residue. Its final receipt serialized
`records=[0, source1065]`. This did not invalidate the prior calculation or
the complete source ledger, but did corrupt the promised marked-pair sidecar.
The producer renamed the records and residue variable, added explicit final
JSON-round-trip regression checks, and regenerated both receipts and manifest.
Revision 2 now has both full authentic records. The independent verifier
compares the entire serialized observer object, so this exact regression
would fail rather than silently pass. The producer's superseded manifest
preserves the correction lineage.

## 4. Finite coverage and results

Independent algebra tests cover:

- All `9^4=6561` ordered Eisenstein four-tuples from `{-1,0,1}^2`.
- 49 rational sharp-family controls, including endpoint and zero cases.
- All perfect matchings of explicit four-coordinate, two-terminal samples at
  arities 2,4,6,8,10 (1,3,15,105,945 matchings).
- All `8!=40320` permutations for a compatible hierarchy: mean energies
  `0,32/7,48/7,8`, debt of mean zero, mean actual debt `8/7`.
- All 256 eight-source sign words for the selected-tree construction;
  direct-sum terminal, empty, padding, singleton, and coherent sharpness cases.
- Literal fixed-low-bit bins through 512 leaves, and all 1331 nonnegative
  integer endpoint/middle triples in `[0,10]^3` for the insertion formula.

The main native comparisons are:

| N | depth cutoff | arrived | S_N | terminal energy | whole-level D | parent-positive | pair-positive |
|---:|---:|---:|---:|---:|---:|---:|---:|
|512|42|768|1244|1313|79|293|433|
|1024|43|1536|2485|2803|335|1619|4288|
|4096|1|1093|9949|1850|66|436|493|
|4096|3|2115|9949|3692|144|2185|4593|
|4096|10|4361|9949|7051|86|6978|30964|
|4096|67|6144|9949|9629|460|11419|86664|

Additionally, the audit reconstructs energies for **every integer native depth
cutoff** from zero through complete arrival: 43+44+68 = 155 cutoffs. Maximum
debts are 90 at N512 depth15; 400 at N1024 depth10; and 494 at N4096 depth36.
This strengthens the finite audit coverage only. It is not an all-scale proof,
and no untested clock/depth rectangular cutoff family is claimed.

## 5. Reproduction and decision

From this audit directory, with only Python >=3.9 standard library and bash:

```bash
bash run_audit.sh
```

Equivalent direct checks:

```bash
python3 independent_verify.py --producer producer_revision2 > replay.json
python3 -O independent_verify.py --producer producer_revision2 > replay.optimized.json
cmp replay.json replay.optimized.json
cmp replay.json receipt.compared.json
```

To regenerate the detailed literal B-chains and full first-entry routes:

```bash
python3 independent_verify.py --producer producer_revision2 --write-sources reconstructed > replay.json
```

All mathematical checks use explicit exceptions, not Python `assert`.
Normal and optimized runs pass and have identical output. The independent
comparison receipt SHA-256 is
`5a7177ac3f205a7a58682d4fcb6f867e32f74bee09a23b54f6d19203f749fdbe`.

**Accept** the scoped algebra, sharp obstruction, changed-tree theorem, native
finite examples, recovered 433/293/79 configuration, and corrected marked-pair
sidecar. **Do not infer** the native fixed-low-bit estimate, universal arrival,
historical computation identity, or a source-norm bound from terminal-energy
control. No unresolved mathematical or reproducibility blocker remains for
this coherent audit milestone.
