# Independent audit of transposition rigidity and all cumulative signed-cycle gaps

**Verdict: PASS, including the necessary complete finite base.** The
[incoming proof candidate](overnight_hexagon_sep05_d7_d8_gap.md), inspected
at worktree commit `545a20c2b5eb`, has a correct analytic Hamilton-layer
argument for every order at least nine, a correct fixed-layer transfer, and
a correct cumulative consequence for every `D>=7`, `n>=D+1`. Its equality
classifications are complete in the stated ranges. No mathematical repair
was needed.

This audit independently replays **all** `2^21` root-gauged K8 switching
classes using packed cycle-sign bitsets and Gray-code updates, with no Walsh
transform or repository mathematics import. It recovers the full required
base, including its zero and equality sets. The `2^28` K9 computation was
not rerun and is not needed by the accepted proof. This report is an audit
artifact, not an automatic promotion of a reserved theorem or navigation
surface.

## 1. Exact accepted scope and dependencies

Let `c_k(H)` count negative **unoriented simple** `k`-cycles of a signing of
`K_n`, considered modulo cut switching. Put

```text
e_(n,k)=(n-2)!/(n-k)!,
S_D(H)=sum_(k=3)^(D+1)c_k(H),
A_(n,D)=sum_(k=3)^(D+1)e_(n,k).
```

Write `B` for the balanced switching class and `A` for the antibalanced
class. The accepted conclusions are:

| Scope | Minimum / equality |
|---|---|
| Hamilton layer `n>=9` | Minimum outside its zero classes is `(n-2)!`, by analytic proof |
| Hamilton layer `n=8` | Minimum outside `{B,A}` is `720`, by the complete finite base |
| Every fixed `k>=8`, all `n>=k` | Minimum outside `Z_k` is `e_(n,k)` |
| Every `D>=7`, all `n>=D+1` | Minimum of `S_D` outside `B` is `A_(n,D)` |

Here `Z_k={B}` for odd `k` and `Z_k={B,A}` for even `k`; these are exactly
the zero sets in the accepted fixed-layer range. Odd-layer equality consists
exactly of single-negative-edge switching classes. Even-layer equality
consists exactly of these classes and their global negatives. For `n>=8`
these are respectively `binom(n,2)` and `2 binom(n,2)` distinct labelled
classes. Cumulative equality consists exactly of the `binom(n,2)`
single-edge classes, one relabelling orbit.

The all-order Hamilton proof is elementary and has **no finite-base
dependency**. The fixed-layer result additionally needs the K8 Hamilton
base. The cumulative result additionally inherits the all-order `D=6`
inequality and its single-edge equality classification from
[THM-4416 / even-graph-cumulative-d5-d6-spectral-gap](../../01-canon/theorems/THM-4416-even-graph-cumulative-d5-d6-spectral-gap.md).
The deletion rigidity is also re-proved below; its prior source is
[THM-4083 / even-graph-cumulative-d3-d4-spectral-gap](../../01-canon/theorems/THM-4083-even-graph-cumulative-d3-d4-spectral-gap.md).
Both source theorems were read with their stated ranges and exceptions.

Only the spectral interpretation needs
[THM-4078 / even-graph-triangle-quotient-spectrum-and-boolean-noncommutation](../../01-canon/theorems/THM-4078-even-graph-triangle-quotient-spectrum-and-boolean-noncommutation.md).
For its multiplicity-weighted cycle operator, the unnormalized Laplacian
gap is `2A_(n,D)`, with labelled multiplicity `binom(n,2)` and quotient
multiplicity one. This is not a Booleanized adjacency or absolute mixing
gap statement.

The inheritance hostile is `A`, which has zero negative even cycles but
nonzero odd layers. MISTAKE-496 requires excluding `B` explicitly from every
nontrivial minimum; MISTAKE-495 requires retaining the canonical all-cycle
operator and its multiplicities. Both corrections were read and are
respected. The short Hamilton hostile is a negative C5 plus a positive apex
on K6, whose Hamilton count is `20<24`. No extension of the accepted
single-layer theorem to every shorter length is inferred.

## 2. Independent analytic Hamilton audit

For a transposition `tau=(u v)`, the edgewise xor `H+tau H` is exactly the
negative two-star graph `K_(2,r)`, where

```text
r=#{x outside {u,v}: H(ux)!=H(vx)},
s=n-2-r.
```

The edge `uv` is unchanged, as are all edges outside the two stars. A cycle
is negative in the xor precisely when its original and transposed parity
bits disagree. Consequently

```text
c_n(H+tau H)<=2c_n(H).                                (1)
```

This is a Hamming triangle inequality on the **complete labelled
Hamilton-cycle set**. It does not assume that `H` is fixed by the
transposition. Switching invariance permits choosing a root gauge first;
the same argument then applies to every transposition in that gauge.

### Two independent counts of the difference

Set `m=n-2`, `q=rs`. Delete `u,v` from a Hamilton cycle. On a fixed cycle
of the other `m` vertices, let `X` count gaps crossing between its R and S
parts. Inserting both vertices into one gap gives `2X` negative extensions.
Putting them in distinct gaps gives `2X(m-X)`: one chosen gap must cross
and one must not, and the assignments of `u,v` are distinct. Thus the total
negative extensions are `2X(m+1-X)`.

For a uniform labelled Hamilton cycle on these `m` vertices,

```text
E X=2q/(m-1),
E X^2=4q(q-1)/[(m-1)(m-2)].                           (2)
```

One direct derivation of the second identity counts ordered pairs of
crossing edges. There are `q(m-2)` adjacent ordered pairs, each present with
probability `2/[(m-1)(m-2)]`, and `q(q-m+1)` disjoint ordered pairs, each
present with probability `4/[(m-1)(m-2)]`. Adding `E X` gives (2).
Multiplying the extension mean by `(m-1)!/2` yields

```text
T_n(r)=2q[(n-2)(n-3)-2q](n-5)!.                       (3)
```

I separately checked the producer's positional derivation. Adjacent and
distance-two positions of `u,v` each contribute `2rs(n-4)!`; the remaining
positions contribute
`2rs[(r-1)(r-2)+(s-1)(s-2)](n-5)!`. Shared neighbours cancel at distance
two. These three contributions simplify to (3), with the reversal factor
included exactly once.

### Strict threshold and its precise boundary

For `2<=r<=n-4`, one has
`2(n-4)<=q<=(n-2)^2/4`. The polynomial
`F(q)=q[(n-2)(n-3)-2q]` is concave. Subtracting
`(n-2)(n-3)(n-4)` at the two real endpoints gives

```text
(n-4)(n^2-13n+38),
(n-2)(n-4)(n^2-12n+28)/8.
```

At `n=9+t`, the quadratic factors are respectively `t^2+5t+2` and
`t^2+6t+1`, strictly positive for every `t>=0`. Hence
`T_n(r)>2(n-2)!` throughout the interior range. The use of the real upper
endpoint is valid: it enlarges the possible q interval, and concavity still
places its minimum at an endpoint.

Combining with (1), if `c_n(H)<=(n-2)!`, every row-disagreement count is in
`{0,1,n-3,n-2}`. At `n=8,r=3`, however, `T_8(3)=1296<1440`; the analytic
threshold genuinely fails there. This is a method boundary, not a
counterexample to the K8 minimum.

### Exhaustive classification forced by the threshold

Root-gauge all edges at vertex zero positive and let `G` be the negative
graph on the other `m=n-1` vertices. Comparing each row with the root gives
degrees in `{0,1,m-2,m-1}`. Partition the vertices into low and high degree
sets of sizes `a,b`. Their number of crossing edges lies between
`b(a-1)` and `a`. If both sizes are at least two, this implies
`(a-1)(b-1)<=1`, hence `a=b=2`, impossible for `m>=8`.

If all vertices are low, `G` is a matching; two distinct matching edges
produce an excluded disagreement count two, leaving at most one edge.
Complementation handles the all-high case. If precisely one vertex is
high, it is adjacent to every other vertex or misses exactly one. In the
second case the missed vertex is isolated, and comparing the high vertex
with one neighbour gives `m-3`, also excluded. Thus only a full star
survives. Complementation handles exactly one low vertex.

The six resulting graph families are exactly empty, complete, single edge,
complement of a single edge, full star, and complement of a full star.
The last four represent all globally signed single-edge switching classes,
including the edges incident to the gauge root. The complete root-gauged
graph represents `A`, not an additional class.

For even `n`, `B,A` have weight zero and both signed single-edge families
have weight `e=(n-2)!`. For odd `n`, their weights are respectively

```text
B:0, single edge:e, A:(n-1)e/2,
globally negative single edge:(n-3)e/2.
```

The latter two exceed `e` for `n>=9`. This finishes the all-order Hamilton
minimum and every equality/zero case without a finite search.

## 3. Fixed-layer transfer audit

Fix `k>=3`, and suppose a full base on `K_(n0)`, with
`n0>=max(k,4)`, proves the claimed `Z_k` exclusions, lower bound, and
equality set. For the induction step `n>n0`, call a deletion exceptional
when it lies in `Z_k`.

For even `k`, a balanced and an antibalanced deletion cannot coexist:
their overlap is a `K_(n-2)` with a triangle, whose sign would have to be
both positive and negative. This compares invariant triangle signs, not
incompatible switching gauges. Two balanced deletions force every negative
triangle to contain their two deleted vertices. Four-set parity forces all
such triangles to have the same sign, and nonbalancedness makes it negative;
triangle signs therefore identify a single-edge switching class. Two
antibalanced deletions give its global negative by applying the argument to
`-H`.

Outside these equality classes, at least `n-1` deletions are nonexceptional.
Each `k`-cycle survives `n-k` deletions, so

```text
c_k(H) >= (n-1)e_(n-1,k)/(n-k)
        = ((n-1)/(n-2))e_(n,k) > e_(n,k).
```

Division is legal because `n>n0>=k`. This proves the transfer, the exact
zero sets, and strictness away from the equality classes. The base equality
classification remains necessary at `n=n0`; it is not replaced by the
induction. Applying this to each analytic Hamilton base `k>=9`, plus the
independently replayed K8 base below, proves every fixed layer `k>=8`.

## 4. Cumulative consequence and antibalanced exception

For `H` outside `{B,A}`, add the inherited THM-4416 `D=6` bound to the
accepted fixed-layer bounds for `8<=k<=D+1`. This gives
`S_D(H)>=A_(n,D)` for every `D>=7`, `n>=D+1`. Equality forces equality in
the D6 summand and hence a single-edge class; each such class attains every
summand.

For `A`, pair each odd length `k-1` with the next even length `k`.
Its odd cycle count is `n!/[2(k-1)(n-k+1)!]`, and its even count is zero.
The pair pays `e_(n,k-1)+e_(n,k)` exactly when

```text
n(n-1)-2(k-1)(n-k+2)>=0.
```

At `n=k+t` this is `(k-1)(k-4)+t+t^2`, nonnegative for every even
`k>=4`, `t>=0`. The first pair `k=4` is strictly positive since `n>=8`.
If the last length is odd, it separately exceeds its edge value, with
ratio `n(n-1)/[2(D+1)]>1`. Thus `A` is strictly above the cumulative
minimum throughout the claimed range.

The candidate's separate expanded S7 excess was also checked by exact
polynomial arithmetic. At `n=8+t`, its coefficients in increasing powers
are

```text
(1652,244667/105,2723/2,1328/3,189/2,73/5,3/2,1/14).
```

All are positive. This is a valid alternate D7 check; the odd/even pairing
already proves the needed antibalanced statement for every D simultaneously.

## 5. Independently reproduced finite evidence

The required K8 base is now independent of the incoming evidence ledger.
The new C++ companion enumerates all 2,520 unoriented Hamilton cycles,
stores for each of the 21 nongauge edges a bitset of its incident cycles,
then visits every root-gauged signing in Gray-code order. A signing step
flips one edge, so xor with its incidence bitset updates **every negative
cycle indicator exactly**. Popcounts give the Hamilton weight. No transform,
relabelling filter, frustration filter, or connectivity filter is used.

The edge order is colexicographic, different from the candidate's order.
The output spectrum is independently remapped to the candidate's canonical
order only for the final digest comparison. Minimum, zero classes, and exact
equality masks are checked independently of that digest.

```text
complete switching classes: 2,097,152
Hamilton cycles: 2,520
zero classes: exactly B,A
nonzero minimum: 720
labelled equality classes: 56, exactly both signed single-edge families
complete spectrum FNV64: e2dfba14125e7983
spectrum storage: 4,194,304 bytes
cycle incidence storage: 6,720 bytes
explicit gates: 2,099,699
```

The new Python companion separately checks:

- Every two-star parameter `0<=r<=n-2` by literal cycle parity at `n=6,...,10`.
- Both strict-threshold polynomial identities and the full S7 excess exactly.
- First and second crossing-gap moments on 81 fixed-cycle color universes.
- The complete disagreement classification on all 1,024 graphs of order five
  and 32,768 graphs of order six.
- All 1,024 root-gauged K6 switching classes and all 15,360 vertex
  transpositions, including the typed zero-deletion conclusions.
- The K6 negative-C5-apex hostile and the K8 strict-threshold hostile.

It passes 32,182 explicit gates. These finite checks corroborate the
all-order analytic derivations; only the full K8 Hamilton enumeration is a
new finite premise needed by the theorem.

## 6. Evidence boundaries and filing notes

The incoming K9 Walsh and path-attachment source/output files were read.
Their displayed full per-character agreement is consistent, including the
optional reference-spectrum check and the common keyed digest. The final
K9 run was **not independently replayed here**, so this audit does not claim
a second execution of those `2^28` characters. It is redundant for the
accepted mathematical consequence.

At the inspected commit `545a20c2b5eb`, the incoming candidate still had
provisional status and a reproduction footer saying hashes and outputs
were to be frozen. Its owner subsequently resolved this filing issue in
incoming commit `4d1ad2a39` and promoted
[THM-4427 / all-cumulative-signed-cycle-gaps-by-transposition-rigidity](../../01-canon/theorems/THM-4427-all-cumulative-signed-cycle-gaps-by-transposition-rigidity.md).
The promoted theorem was read and has the same audited mathematical scope.
This report supplies additional distinct frozen evidence for its K8 premise;
it neither creates nor claims ownership of that theorem namespace.

No Hamilton equality statement below order eight, complete cycle-profile
cone, near-minimizer classification, Boolean spectrum, tournament
inequality, or LRC(14) result follows from this audit.

The live concept board was: switching gauge; typed zero deletions;
transposition difference; Hamming weight; Hamilton insertion; cumulative
parity pairing. The completed connection contract is:

| Field | Content |
|---|---|
| Source | Signed switching class with its labelled cycle-parity vector |
| Target | Explicit two-star count and rigid equality classes |
| Map | Difference with a transposed signing; then root-gauge classification |
| Preserved predicate | Cycle xor and the Hamming triangle inequality |
| Lost information | Which original negative cycles cancel in the difference |
| Sidecar | Both original weights, disagreement counts, and parity of the layer |
| Decisive hostile | K8 two-star threshold; K6 negative C5 with positive apex |

The META-PATTERNS used are the quotient-loss discipline and testing the first
failed order before generalization. No new pattern is promoted from this
audit alone.

## Reproduction and frozen companions

From the repository root, for example:

```text
python 04-computation/overnight2_20260906_signed_cycle_audit.py
python -O 04-computation/overnight2_20260906_signed_cycle_audit.py
g++ -O2 -std=c++17 -mpopcnt 04-computation/overnight2_20260906_signed_cycle_audit_k8.cpp -o C:/w/overnight2_20260906_signed_cycle_audit_k8.exe
C:/w/overnight2_20260906_signed_cycle_audit_k8.exe
```

Python normal and optimized streams byte-match. The C++ source was also
built with `-DNDEBUG`; its stream is identical. Every gate is explicit and
remains active in those optimized modes. Frozen text files use LF line
endings; native Windows C++ stdout may require CRLF-to-LF normalization when
comparing file hashes.

Companions:
[Python controls](../../04-computation/overnight2_20260906_signed_cycle_audit.py),
[Python output](overnight2_20260906_signed_cycle_audit.out),
[K8 bitset verifier](../../04-computation/overnight2_20260906_signed_cycle_audit_k8.cpp),
[K8 output](overnight2_20260906_signed_cycle_audit_k8.out).

```text
SHA-256, LF bytes:
Python source 4d39ab3ae0cabc726b24f2ddb35f3483654df29b3b91edc251ee066465eaa39d
Python output d784dd19b7c869e84814d8b9bcd5afc3473a9d7da47034bcbb2579e7114a77b1
C++ source aedf2a78bb7045c1a663206ecc1c9279cd14d148bc4ef882775f0c409ae5558f
C++ output 3d674902b9e8a4e690501762bfd2b33064c3b04e12d25d8e3a287cb4eba509c1
```
