---
id: THM-4083
title: "Even-graph cumulative D3/D4 spectral gap"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY VERIFIED-EXACT for D=3 and D=4.
  THM-4416 separately closes D=5,6; the cumulative statement for D>=7 remains OPEN.
source: codex-frontier-synthesis-creative-20260825e / cumulative signed-dual niche
audit: >
  PASS. The proof-local path directly gauges and exhausts all 33,864 labelled
  switching classes through n=7, performs 38,606,835 direct cycle-parity
  gates, audits the balanced-deletion trichotomy and four-set rigidity, and
  certifies n=8 by the exact deletion bridge. The independent path generates
  cycles from ordered tuples and applies exact integer Walsh transforms to
  2,131,016 characters through n=8, using 1,511 cycle vectors and 66,813,528
  butterfly gates. Both paths freeze the n=4,...,8 D3/D4 minima and labelled
  equality multiplicities. Normal and optimized outputs byte-match; both
  scripts have zero assert nodes and zero floating literals.
depends_on:
  - THM-4078-even-graph-triangle-quotient-spectrum-and-boolean-noncommutation
related:
  - THM-4073-even-graph-diameter-layer-exact-cycle-distance
  - THM-4069-even-graph-basis-dependence-and-canonical-cycle-envelope
script: 04-computation/even_graph_cumulative_d3_d4_gap_thm4083.py
output: 05-knowledge/results/even_graph_cumulative_d3_d4_gap_thm4083.out
script_sha256: f43a38eb1c4f549e6884ceb159e01fef14c94b97253c278920e20972a0c22508
output_sha256: d8daaf5bd804dd7d8ee3b9a27e84e17d5c4dfddeffebbacb5660500ff07c672d
independent_audit_script: 04-computation/even_graph_cumulative_d3_d4_gap_thm4083_independent_audit.py
independent_audit_output: 05-knowledge/results/even_graph_cumulative_d3_d4_gap_thm4083_independent_audit.out
independent_audit_script_sha256: 0cec330209dc2c1e77d31a96c1e9faa7ac985d8e3ec0540c5df1ea33b2dc3372
independent_audit_output_sha256: 6b99f3e7e1b66bd088f6a36ca49bbb2640a98aa5601229f4442d1aee154f1773
hash_basis: raw LF bytes
---

# THM-4083 -- even-graph cumulative D3/D4 spectral gap

**PROVED + VERIFIED-EXACT + INDEPENDENTLY VERIFIED-EXACT.** The signed
dual from THM-4078 reduces a cumulative diameter-layer gap to a negative-cycle
extremum. For `D=3` and `D=4`, balanced vertex deletion gives enough rigidity
to solve that extremum for every admissible `n`. The analogous statement for
`D=5,6` is now closed separately by
[THM-4416](THM-4416-even-graph-cumulative-d5-d6-spectral-gap.md).
The remaining `D>=7` boundary is explicit in Section 8.

## 1. Typed signed-dual reduction

Let `C_(n,k)` be the set of unoriented simple `k`-cycles of `K_n`, and put

```text
N_(n,k)=|C_(n,k)|=n!/(2k(n-k)!),
c_k^-(H)=#{C in C_(n,k): |E(C) intersect H| is odd}.          (1)
```

Here `H` is a signed complete graph, represented by its negative-edge set,
and `[H]` is its labelled switching class. For `2<=D<=n-1`, define

```text
S_D(H)=sum_(k=3)^(D+1)c_k^-(H),
Q_(n,D)=sum_(k=3)^(D+1)N_(n,k),
M_(<=D+1)=sum_(k=3)^(D+1)M_k.                                (2)
```

On the full labelled Eulerian Cayley lift, the character belonging to `[H]`
has cycle-layer eigenvalue

```text
lambda_k(H)=sum_(C in C_(n,k))(-1)^|E(C) intersect H|
           =N_(n,k)-2c_k^-(H).                               (3)
```

The same value labels the corresponding orbit-sum eigenvector on the
unlabelled quotient of THM-4078. Therefore

```text
lambda_(<=D+1)(H)=Q_(n,D)-2S_D(H),
Laplacian eigenvalue=Q_(n,D)-lambda_(<=D+1)(H)=2S_D(H).       (4)
```

The balanced switching class is the trivial character and has `S_D=0`.
If one fixed edge is negative, exactly the cycles through that edge are
negative, so

```text
e_(n,k)=(n-2)!/(n-k)!,
A_(n,D)=sum_(k=3)^(D+1)e_(n,k).                              (5)
```

## 2. Exact D3 and D4 extrema

For every `n>=4`,

```text
boxed: min_([H] nontrivial) S_3(H)
      =(n-2)+(n-2)(n-3)=(n-2)^2.                            (6)
```

For `n>=5`, equality in `(6)` consists exactly of the single-negative-edge
switching classes. At `n=4`, the antibalanced class also ties.

For every `n>=5`,

```text
boxed: min_([H] nontrivial) S_4(H)
      =(n-2)^2+(n-2)(n-3)(n-4)
      =(n-2)(n^2-6n+10).                                    (7)
```

Equality in `(7)` consists exactly of the single-negative-edge switching
classes.

## 3. Balanced-deletion trichotomy

Let `T(H)` be the set of negative triangles. A complete signing is balanced
exactly when `T(H)` is empty, and triangle signs determine its switching
class: two signings with the same triangle signs have balanced edgewise
ratio. Every four-set contains an even number of members of `T(H)`, because
the product of its four triangle signs contains each edge sign twice.

Call `v` **balanced-deleting** when `H-v` is balanced, and write

```text
B(H)={v: H-v is balanced},             b(H)=|B(H)|.          (8)
```

The following classification holds for every nontrivial signed `K_n`,
`n>=4`.

### Case I: at least two balanced deletions

Suppose `u,v in B(H)`. Every negative triangle must contain both `u` and
`v`, so it has the form `uvx`. Four-set parity on `{u,v,x,y}` gives

```text
1_(uvx in T(H))=1_(uvy in T(H)).                             (9)
```

Nontriviality makes the common value one. Thus

```text
T(H)={uvx: x notin {u,v}},                                  (10)
```

which is the triangle set of a single negative edge `uv`. Hence `[H]` is
that single-edge switching class, and `(5)` gives all of its cycle counts.

### Case II: exactly one balanced deletion

Suppose `B(H)={v}`. All negative triangles contain `v`. On
`W=V(K_n)\{v}`, define a graph `R` by

```text
xy in E(R) iff vxy in T(H).                                  (11)
```

Four-set parity on `{v,x,y,z}`, together with the positivity of `xyz`, says
that every triangle of `R` has even edge parity. Fixing `r in W` and putting
`p_x=1_(rx in E(R))` gives

```text
1_(xy in E(R))=p_x+p_y mod 2.                                (12)
```

Thus `R=K_(X,Y)` is a cut. Both parts have size at least two: if, say,
`X={x}`, then every negative triangle contains both `v` and `x`, making `x`
a second balanced deletion.

Switching gives a representative whose only negative edges are `vx` for
`x in X`. A negative `k`-cycle must contain `v` and have one neighbour of
`v` in each part. If `r=|X|`, `s=|Y|`, then `r+s=n-1` and

```text
c_k^-(H)=rs (n-3)!/(n-k)!,
S_D(H)=rs/(n-2) A_(n,D).                                    (13)
```

Since

```text
rs-(n-2)=rs-(r+s-1)=(r-1)(s-1)>0,                           (14)
```

Case II is strictly above the single-edge value for every admissible `D`.

### Case III: no balanced deletion

If `b(H)=0`, every `H-v` is nontrivial. This is exactly the class to which
the deletion averages in Sections 5 and 6 apply.

## 4. Four-set rigidity for D3 equality

On a fixed four-set, the number of its negative triangles is `0`, `2`, or
`4`. A direct local switching check gives respectively `0`, `2`, or `0`
negative Hamiltonian four-cycles. Consequently

```text
c_4^-(H)=2 #{U in C(V,4): |T(H) intersect C(U,3)|=2}.         (15)
```

In particular, if `c_4^-(H)=0`, all four triangle signs on every four-set
are equal. Any two triangles sharing an edge lie in one four-set, and the
graph of triples joined when they share two vertices is connected. Hence all
triangle signs are equal. A nontrivial such signing is antibalanced, so

```text
c_3^-(H)=C(n,3).                                             (16)
```

For `n>=5`,

```text
C(n,3)-(n-2)^2=(n-2)(n-3)(n-4)/6>0.                         (17)
```

This closes the equality branch that the deletion average alone does not
separate.

## 5. Proof of the D3 formula

Write `F(H)=c_3^-(H)+c_4^-(H)`. At `n=4`, four-set parity gives the three
possibilities

```text
(c_3^-,c_4^-)=(0,0), (2,2), (4,0).                           (18)
```

The middle vector is the single-edge class and occurs for each of the six
labelled edges; the last is the antibalanced class. Thus the nontrivial
minimum is four, with the claimed tie.

Proceed by induction for `n>=5`. Cases I and II are settled by `(5)` and
`(14)`. In Case III, every deletion is nontrivial, so the induction
hypothesis and the fact that a `k`-cycle survives exactly `n-k` deletions give

```text
sum_v F(H-v)
 =(n-3)c_3^-(H)+(n-4)c_4^-(H)
 =(n-4)F(H)+c_3^-(H)
 >=n(n-3)^2.                                                 (19)
```

If `c_3^-(H)>=(n-2)^2`, then `(6)` follows. Otherwise,

```text
(n-4)F(H)>n(n-3)^2-(n-2)^2
 =(n-4)(n-2)^2+(n-3)(n-4),                                  (20)
```

so the inequality is strict. If equality in `(6)` were to occur in the
first branch of Case III, it would force
`c_3^-(H)=(n-2)^2` and `c_4^-(H)=0`; equations `(16)`--`(17)` rule this out.
Thus Case III is also strict at equality. This proves `(6)` and its complete
equality classification.

## 6. Proof of the D4 formula

Put

```text
G(H)=c_3^-(H)+c_4^-(H)+c_5^-(H),
b_n=(n-2)(n^2-6n+10).                                       (21)
```

The proof-local exact census supplies the bases

```text
n                   5    6    7
min G(H)           15   40   85
labelled minimizers 10   15   21,                            (22)
```

and checks that every minimizer is a single-edge class. These are exhaustive
root-gauge universes of sizes `64`, `1,024`, and `32,768`; the independent
audit also exhausts `n=8`.

For completeness, `n=8` follows analytically from `(22)`. Cases I and II are
already settled. In Case III, if `G(H)<=b_8=156`, deletion gives

```text
680=8b_7
 <=sum_v G(H-v)
 =4G(H)+c_3^-(H)-c_5^-(H)
 <=4G(H)+C(8,3)
 <=680.                                                       (23)
```

If `G(H)<156`, the last bound is at most `676`, already impossible. At
equality, every inequality in `(23)` is an equality, forcing
`c_3^-(H)=56` and `c_5^-(H)=0`. All triangles are then negative, so `H` is
antibalanced and every odd cycle is negative; in particular
`c_5^-(H)=N_(8,5)=672`, a contradiction. Thus the `n=8` minimum is `156`,
attained only by the `28` single-edge classes.

Now let `n>=9` and suppose Case III holds. Induction gives

```text
sum_v G(H-v)
 =(n-3)c_3^-(H)+(n-4)c_4^-(H)+(n-5)c_5^-(H)
 =(n-4)G(H)+c_3^-(H)-c_5^-(H)
 >=n b_(n-1).                                                (24)
```

If `G(H)<=b_n`, then the left side of the last inequality is at most
`(n-4)b_n+C(n,3)`. But exact algebra gives

```text
n b_(n-1)-[(n-4)b_n+C(n,3)]
 =5(n-8)(n-4)(n-3)/6>0,                                     (25)
```

a contradiction. Cases I and II finish the induction, and their strictness
properties give the equality classification in `(7)`.

## 7. Spectral gaps and typed multiplicities

Equations `(4)`, `(6)`, and `(7)` give the unnormalized Laplacian gaps

```text
gap_(D=3)=2(n-2)^2,                         n>=4,
gap_(D=4)=2(n-2)(n^2-6n+10),                n>=5.             (26)
```

The exact degrees are

```text
Q_(n,3)=C(n,3)+3C(n,4),
Q_(n,4)=C(n,3)+3C(n,4)+12C(n,5).                             (27)
```

For the normalized nonlazy operator `P_D=M_(<=D+1)/Q_(n,D)` and its lazy
version `P_D^lazy=(I+P_D)/2`, the algebraic gaps are therefore

```text
gamma_(n,3)=2(n-2)^2/Q_(n,3),
gamma_(n,3)^lazy=(n-2)^2/Q_(n,3),

gamma_(n,4)=2(n-2)(n^2-6n+10)/Q_(n,4),
gamma_(n,4)^lazy=(n-2)(n^2-6n+10)/Q_(n,4).                  (28)
```

Multiplicity must be typed by the ambient space.

- In the **full labelled Fourier lift**, the D3 gap has multiplicity `7` at
  `n=4` and `C(n,2)` for `n>=5`. The D4 gap has multiplicity `C(n,2)`.
- In the **unlabelled switching quotient** of THM-4078, all single-edge
  classes form one relabeling orbit, so their orbit sum contributes
  multiplicity one. At `(n,D)=(4,3)`, the antibalanced orbit is distinct, so
  the quotient gap has multiplicity two. In every other case proved here it
  has multiplicity one.

Thus `C(n,2)` is a labelled-character count, not a quotient multiplicity.

## 8. Exact audits, OPEN boundary, and loss ledger

Replay the new proof-local and independent paths from the repository root:

```bash
python3 -B 04-computation/even_graph_cumulative_d3_d4_gap_thm4083.py
python3 -B -O 04-computation/even_graph_cumulative_d3_d4_gap_thm4083.py
python3 -B 04-computation/even_graph_cumulative_d3_d4_gap_thm4083_independent_audit.py
python3 -B -O 04-computation/even_graph_cumulative_d3_d4_gap_thm4083_independent_audit.py
```

The primary audit uses direct cycle parity, not Fourier transformation. It
checks every nontrivial root-gauged signing through `n=7`, the complete
balanced-deletion classification `(9)`--`(14)` for every cycle length, the
local identity `(15)`, all equality masks, and the `n=8` bridge `(23)`. The
independent audit does not inspect balanced deletions or invoke induction; it
generates cycles by a different ordered-tuple route and recomputes every D3
and D4 prefix on all `2^21` characters at `n=8` using integer Walsh
transforms. THM-4078's older finite-exact Walsh sidecar supplies a third
consistent census, but is not the sole base certificate here.

For `D>=7` and `n>=D+1`, the statement

```text
min_([H] nontrivial) S_D(H)=A_(n,D)                           (29)
```

remains **OPEN**. Cases I and II reduce its value assertion to the weak
inequality

```text
S_D(H)>=A_(n,D) for every H with b(H)=0.                      (30)
```

The predicted equality classification--only the single-edge classes--is the
strict version of `(30)`. Thus a `b(H)=0` signing below `A_(n,D)` would refute
the value, while one exactly at `A_(n,D)` would preserve the value but refute
its expected rigidity.

The exact deletion identity in that unresolved class is

```text
sum_v S_D(H-v)=sum_(k=3)^(D+1)(n-k)c_k^-(H).                 (31)
```

It underweights the longest layers and gives Hamiltonian cycles coefficient
zero. A termwise repair is false: an antibalanced signing has no negative
even cycles. THM-4416 repairs D5/D6 with joint local-profile inequalities and
exact subset multiplicities; larger D still needs such a sidecar. The
finite-exact box through `n=8` does not prove `(29)` for arbitrary `D`.

Finally, the source object is a signed switching character; the map

```text
(c_3^-,...,c_(D+1)^-) -> S_D -> 2S_D                         (32)
```

preserves the cumulative weighted eigenvalue but destroys the separate
length counts and, away from the equality classification, the signing. The
full negative-cycle vector and `B(H)` are the required sidecars. At
`(n,D)=(4,3)`, even the minimizing scalar does not identify the single-edge
orbit because the antibalanced orbit ties.

MISTAKE-495 is respected: `(3)`--`(4)` use the canonical all-simple-cycle
layers. They make no claim that an arbitrary fundamental-cycle generator set
or its edge quotient is basis independent.
