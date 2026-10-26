# Independent audit: native alignment rank and named path pairing

2026-10-26. **ACCEPTED, PROVED in the stated scopes below.** Exact controls:
**FINITE-EXACT**, independent verifier, normal and `-O` byte-identical.
This initial audit is frozen for Session70 checkpoint 1. Later graph-prose
and graceful-guard review belongs to a separate supplement; it must not
silently replace this record.

Base input: `fb872272125510c317b8114dd3f697b433c595f9`.
No repository or producer files were edited. Sparse-excluded inherited
files were read from pinned git objects. The coordination processor could
not run because `.machine-id` is absent; no identity was fabricated.

## 1. Audited native report and provenance

Accepted report:
`session70-vitali/REPORT.md`, SHA-256
`da137bfa8ad50378a29f963b87e7b7950b8878c4cfea79590c06783c0d4d4792`.

Accepted checker:
`session70-vitali/check_alignment_moments.py`, SHA-256
`e3836182474bf9a336c148f73ab1fb24640bef79cd31b077545592da3e734d88`.

The producer checker was separately replayed normally and under `-O`;
both match its 24,741-check receipt SHA-256
`ada554eb657e00356f91c676e16995cd141630ebac3c47b9ff4e8ce8b71f6327`.
Our independent checker imports no producer code and gives a second path.

Dependencies, not new results:

- C3/C4 of `05-knowledge/results/strong_atom_continuation_vitali_20261026/README.md`
  give the seed and the `2^k` native paired-rooted fibre at `n=7k+2`.
- THM-840, `01-canon/theorems/THM-840-hamming-five-continuation-congruence-boundary.md`,
  §1a gives the exact general linear-sidecar rank `dim N(ker R)`.
- THM-2237, `01-canon/theorems/THM-2237-truncated-boolean-moment-interval-and-parity-top-atom-majorants.md`,
  §3 already proves abstract proper-marginal parity ambiguity.

The new contribution is actual simultaneous-guard realization inside the
native tournament fibre and its resulting exact rank hierarchy, plus source
rigidity. It is not priority for Boolean inversion or parity laws.

## 2. Native guard and orientation audit

Independently decoding free masks 1056 and 4160 gives exactly

```
X=(2,12,41,17,39,75,31), Y=(2,20,9,51,37,71,31).
```

The seed guard is the old-vertex order `(0,1,2,3,5,6,4)`.
It fails in X and succeeds in Y. The old-to-new isomorphisms are

```
s=(0,1,4,2,3,5,6): X -> Y,
t=(0,2,4,1,5,3,6): tau(X) -> tau(Y).
```

All 7! permutations independently show that each of these maps is unique
and both X and tau(X) are rigid. Their mismatch is
`t^-1 s=(1 3)(4 5)`, which is not an automorphism of X.

For bitword b, build the ordered sum of k seed blocks and append u,v with
`S -> u -> v -> S`. Uniform inter-block arcs make the altered exit of the
local guard harmless. Its last selected-block exit goes to u, and u goes
to v. Thus every subset A defines one actual permutation `P_A`, with

```
legal(P_A,T_b) = product_(i in A) b_i.
```

This is simultaneous one-step legality, not multiplication of sequential
conditional probabilities. Re-marking changes coordinates, so the latter
claim would require separate transport and is not inferred.

The source and target isomorphisms both fix each block's original first
and last vertex. They therefore preserve every inter-block hinge in the
complemented target; both rooted endpoint types really are equal. The outer
C3 makes each source strong. First/last degrees are `7k-5,7k`; their sum is
not `n-1` for integer k>=1, proving the literal presentations black.

For rigidity, u has the unique degree 1, v the unique degree 7k, and all
S-vertices have degrees in `[2,7k-1]`. Fixing u,v leaves the ordered SCCs of
S, each a rigid seed. Thus the source is rigid for every k. Growing alignment
memory persists without source symmetries and is not an equivariant-selector
obstruction in disguise.

## 3. Rank, probabilities, and boundary audit

With rows A and columns B=supp(b), the matrix is `Z[A,B]=1_(A subset B)`.
Numerical subset-mask order makes Z upper triangular with diagonal one.
Its inverse has entries `(-1)^(|A|-|B|)` when B is a subset of A. This is an
integer inverse, so rank `2^k` is valid over **every field**, including
characteristic two.

For `0<=t<=k`, retain all rows of size at most t, including the empty row.
They are independent, so `rank R_t=d_t=sum_(j=0)^t binom(k,j)`.
The full guard map N is invertible. Hence

```
minimum repairing linear-sidecar rank = dim N(ker R_t) = 2^k-d_t.
```

Appending the omitted rows attains the bound. The lower bound is exactly
THM-840's inherited mechanism applied to this native map.

Required scope distinctions, all respected by the accepted report:

- A single deterministic presentation b needs only k singleton guard bits;
  it does not need `2^k` deterministic bits.
- The exponential dimension is for linear representations recovering every
  guard, or arbitrary signed/distribution weights. It is not a bound on an
  arbitrary nonlinear, unbounded-precision encoding.
- For real probability laws, normalization removes one affine dimension.
  Since R_t includes mass, the missing affine dimension remains `2^k-d_t`.
  Strictly positive laws allow sufficiently small kernel perturbations.
- The normalized even/odd laws are real/rational probability laws, not
  characteristic-two probability objects. Every proper marginal agrees;
  the full native guard gap is exactly `2^(1-k)` and tends to zero.
- For fixed t, k=t+1 gives a native obstruction to all t-local statistics;
  it does not imply blindness of arbitrary global named permutations.
- k=0 is a formal one-state constant-map boundary, not a pair of parity
  laws or a strong C3 construction at order two.

## 4. Stable graph checks included in this initial verifier

These independently verified mathematical statements are stable. The final
graph report is evolving and is not frozen as a prose artifact by this
initial audit; the supplement will record its exact reviewed hash.

### 4.1 Named endpoint matching

For r nontrivial selected paths with 2r distinct named exposed endpoints,
mandatory retained edges, and arbitrary fresh internally-disjoint connectors,
the endpoint perfect matching is a sufficient closed signature for forest
legality and final spanning-path acceptance. Joining endpoints from distinct
pairs merges them and pairs their unused mates. Joining partners creates a
cycle. Distinct pairings P,Q are separated by a pair in Q but not P: insert
a length-two fresh connector, producing a cycle in Q and a merge in P.
Connect P's remaining r-2 components in a chain. The final positive object
is a spanning path, while Q retains its isolated cycle. Therefore the exact
behavioral state count is `(2r-1)!!` and the fixed-length bit lower bound is
its base-two logarithm rounded up.

This is unrestricted named gluing, **not** a claim about realization of
every pairing or separator inside a fixed square-sum graph Q_n. It also
does not reconstruct hidden path interiors or weighted counts. r=1 is the
single-pair case. Singleton paths require endpoint-occurrence coincidence
data, and are excluded from the distinct-endpoint theorem. Hidden closed
cycles must either be rejected by the forest invariant or separately marked;
an empty boundary is not by itself a certificate of an empty graph. If a
cycle is allowed only terminally, closure is accepted only when it spans
the required whole object. These boundary conventions cannot be omitted.

### 4.2 Petersen ports and the K3,3 mark refinement

Contracting the five Petersen spokes gives K5. Retain the outer/inner
identity of each incident edge. A quotient edge subset lifts to a spanning
2-factor iff each contracted vertex has equal outer and inner selected
degrees in `{1,2}`. The spoke is selected exactly for common degree one.
This follows directly from `d_O+x=d_I+x=2`, `x in {0,1}`; the lift is unique.

The bare K5 has 12 Hamilton cycles, none with balanced outer/inner ports
(an odd alternating cycle is impossible). This alone only excludes
one-run-per-cluster lifts. The full local rule allows degree-four quotient
vertices, representing two runs in a contracted cluster. It yields exactly
the six native Petersen 2-factors, all two disjoint 5-cycles. Global
connectedness remains a separate necessary test.

Delete outer vertex 0 and suppress the actual paths
`2-1-6`, `3-4-9`, `7-5-8`. The quotient is K3,3 with parts
`{2,8,9}` and `{3,6,7}`. A quotient Hamilton cycle lifts to a spanning cycle
of Petersen-minus-0 iff it uses all three marked suppressed edges. Exactly
two of K3,3's six Hamilton cycles do; their expansions are exactly the two
Hamilton cycles of Petersen-minus-0. The remaining four skip at least one
original subdivision vertex. Thus bare quotient Hamiltonicity is not a
spanning lift certificate, while the marked-edge carrier repairs the loss.

These are graph statements. No tournament orientation or square-sum
arithmetic compatibility is manufactured from their common vocabulary.

## 5. Reproduction and evidence

```
python3 -B audit_bridge_bounds.py > audit.normal.json
python3 -O -B audit_bridge_bounds.py > audit.optimized.json
cmp audit.normal.json audit.optimized.json
```

**PASS: 78,653 explicit checks in each mode; byte-identical receipts.**
Checks use exceptions, not removable assertions. Universe:

- all seed permutations, independently decoded seeds;
- all `2^k` native states and all `2^k` subset guards for k=1,...,6,
  totaling 5,460 actual guard tests;
- the exact integral inverse for each matrix and every truncation rank
  modulo 2, 3, and 5; every proper normalized marginal;
- all endpoint matchings through r=5 (1,3,15,105,945 states), and all
  11,136 ordered distinct-pair final contexts through r=4, including every
  positive intermediate forest prefix;
- all 1,024 K5 edge subsets, all native Petersen 2-factors, all K3,3 cycles,
  and all spanning Petersen-minus-0 cycles;
- an intentionally false bare-K5 lifting claim, rejected in both modes.

No claim is made here to audit the square-sum or graceful-tree external
literature, their full inherited computation corpus, or all-order tournament
guard reconstruction outside the explicit Boolean subfamily.

The independent verifier SHA-256 is
`b33318c27b063344a283721d3773d162603f68d3a255e24f0374f55f817d9946`.
Both independent receipts have SHA-256
`256d1988e29f0091c5c55f868acb65c6a8cb8294c8384aa674e28e50c35ec8ea`.
`SOURCES-INITIAL.json` pins inherited bytes and reviewed producer hashes;
`SHA256SUMS-INITIAL` freezes this initial artifact set.
