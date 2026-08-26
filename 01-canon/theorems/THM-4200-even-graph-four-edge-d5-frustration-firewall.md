---
id: THM-4200
title: "Even-graph four-edge D5 frustration firewall"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY VERIFIED-EXACT. For every n>=6,
  every signed K_n represented by four negative edges has cumulative odd
  cycle count S_5 at least the single-edge value A_n. Equality occurs only
  for a K_(1,4) negative star at n=6, and switching at its center leaves one
  negative edge. Consequently every switching class that is neither balanced
  nor single-edge and lies at or below A_n has frustration index at least five
  and no balanced vertex deletion.
  The complete D=5 single-edge conjecture remains OPEN.
source: open-frontiers-incoming-20260826b / recovered four-edge support lane
depends_on:
  - THM-4084-even-graph-matching-character-profile-and-d5-firewall
related:
  - THM-4083-even-graph-cumulative-d3-d4-spectral-gap
  - THM-4078-even-graph-triangle-quotient-spectrum-and-boolean-noncommutation
script: 04-computation/even_graph_four_edge_d5_frustration_firewall_thm4200.py
output: 05-knowledge/results/even_graph_four_edge_d5_frustration_firewall_thm4200.out
independent_audit_script: 04-computation/even_graph_four_edge_d5_frustration_firewall_thm4200_independent_audit.py
independent_audit_output: 05-knowledge/results/even_graph_four_edge_d5_frustration_firewall_thm4200_independent_audit.out
script_sha256: cca54f8c845191171d4d5fa689bea38e04cd018138730729e5326af00e16c817
output_sha256: 73e8f29c4f3bac713206f752f47b9c0802ead4780bbb86a0ca5732a04e4b293a
independent_audit_script_sha256: e24f2fc5630454dc69da183dd0c772ad11b142143d60b2f5dc10ee02362fe183
independent_audit_output_sha256: d0fc9c09fff5706302e707c2436291072ec06dc45f26650ef7ac72fbdd0e9e85
hash_basis: raw LF bytes
---

# THM-4200 -- the four-edge D5 frustration firewall

**PROVED + VERIFIED-EXACT + INDEPENDENTLY VERIFIED-EXACT.** THM-4084
excluded every class that is neither balanced nor single-edge and has
frustration index at most three
from the still-open cumulative `D=5` single-edge conjecture. The next support
layer has only eleven graph types. Exact cycle-containment inversion closes
all eleven at once and raises that firewall from four to five.

The equality boundary is load-bearing. A four-edge star does tie at the
smallest order, but it is another representative of the single-edge class,
not a new minimizing switching class.

## 1. Statement and inherited notation

For a signed complete graph represented by its negative-edge set `H subset
E(K_n)`, let

```text
C_(n,k)  = set of unoriented simple k-cycles of K_n,
c_k^-(H) = #{C in C_(n,k): |E(C) intersect H| is odd},
S_5(H)   = sum_(k=3)^6 c_k^-(H).                         (1)
```

Cycle parity is unchanged by switching because every cycle crosses a cut an
even number of times. Thus every quantity in `(1)` depends only on the
labelled switching class `[H]`. Retain THM-4084's single-edge value

```text
A_n=(n-2)(n^3-11n^2+41n-50).                           (2)
```

> **Four-edge theorem.** For every `n>=6` and every four-edge negative set
> `H subset E(K_n)`,
>
> ```text
> S_5(H)>=A_n.                                         (3)
> ```
>
> Equality holds if and only if `n=6` and the four edges form `K_(1,4)`.
> In that equality case `[H]` is a single-edge switching class.

Consequently, if `[H]` is neither balanced nor a single-edge class and

```text
S_5(H)<=A_n,                                           (4)
```

then

```text
frustration([H])>=5,             b(H)=0.               (5)
```

Here `b(H)` is THM-4083's number of balanced vertex deletions. The second
condition in `(5)` was already proved in THM-4083/4084; the new content is
the strict frustration-four exclusion.

## 2. Exact cycle-containment inversion

For a nonempty edge set `F`, write

```text
N_k(F)=#{C in C_(n,k):F subset E(C)}.                  (6)
```

If `q=|E(C) intersect H|`, the elementary parity identity

```text
1_(q odd)=sum_(j=1)^q (-2)^(j-1) binom(q,j)
```

gives, after summing over cycles and exchanging the finite sums,

```text
c_k^-(H)
 =sum_(emptyset != F subset H)(-2)^(|F|-1) N_k(F).     (7)
```

The containment term has a complete graph-theoretic form. Let `j=|F|`.

- If some vertex of `F` has degree at least three, then `N_k(F)=0`.
- If `F` has a cycle component, then `N_k(F)=1` precisely when `F` itself is
  one `k`-cycle, and it is zero otherwise. A cycle component already uses both
  incident places at each of its vertices, so no larger simple cycle can
  contain another component or leave that cycle.
- Otherwise `F` is a linear forest. If it has `c` path components and
  `v=j+c` vertices, then

  ```text
  N_k(F)=2^(c-1)(k-j-1)! binom(n-v,k-v).                (8)
  ```

For `(8)`, choose the `k-v` extra vertices, contract each prescribed path to
one block, cyclically order the resulting `k-j` objects, and orient the `c`
path blocks. Cyclic order modulo reversal supplies `(k-j-1)!/2`, while the
orientations supply `2^c`. The degenerate contracted cases with one or two
objects give the same formula by direct endpoint closure.

Equations `(7)--(8)` are an exact all-layer profile for every four-edge
support. No asymptotic estimate or finite-order extrapolation is used in the
primary proof.

## 3. The eleven supports and their strict gap polynomials

Ignoring isolated vertices, a connected four-edge graph is one of the three
four-edge trees `K_(1,4)`, the fork, and `P_5`, or one of the two unicyclic
graphs `C_4` and the paw. For disconnected graphs, the component edge
partitions `3+1`, `2+2`, `2+1+1`, and `1+1+1+1` give the other six types

```text
K_(1,3)+K_2, K_3+K_2, P_4+K_2,
2P_3, P_3+2K_2, 4K_2.                                  (9)
```

This is a complete eleven-type classification. Substitute each type into
`(7)--(8)` and subtract `(2)`. In the table, `n=n_0+t` and `t>=0`.

| negative support | `n_0` | `S_5(H)-A_n` |
|:---|---:|:---|
|`K_(1,4)`|6|`3t^4+21t^3+51t^2+48t`|
|paw|6|`3t^4+23t^3+61t^2+70t+32`|
|fork|6|`3t^4+25t^3+63t^2+76t+40`|
|`K_(1,3)+K_2`|6|`3t^4+27t^3+57t^2+66t+36`|
|`K_3+K_2`|6|`3t^4+27t^3+57t^2+66t+40`|
|`C_4`|6|`3t^4+25t^3+71t^2+92t+48`|
|`P_5`|6|`3t^4+27t^3+65t^2+74t+36`|
|`P_4+K_2`|6|`3t^4+29t^3+63t^2+56t+36`|
|`2P_3`|6|`3t^4+29t^3+59t^2+72t+32`|
|`P_3+2K_2`|7|`3t^4+43t^3+172t^2+257t+177`|
|`4K_2`|8|`3t^4+57t^3+333t^2+732t+596`|

Every coefficient is nonnegative. Every row has positive constant term
except `K_(1,4)`, and that row is strict for `t>0`. Hence `(3)` holds and its
only possible equality is the displayed star at `n=6`.

At `n=6`, the center cut of this `K_(1,4)` has five edges. Symmetric
difference with the four negative star edges leaves exactly the edge from the
center to the sixth vertex. Thus switching at the center turns the equality
row into a one-edge representative. Its frustration index is one, proving the
claimed equality classification at the level of switching classes.

## 4. Frustration-four exclusion

Suppose a switching class has frustration index four. Choose a minimum
representative `H` with exactly four negative edges. It is one of the eleven
types in Section 3, so `(3)` applies. Equality would force the `n=6` star,
whose class has frustration one. Therefore every frustration-four class
satisfies

```text
S_5(H)>A_n.                                             (10)
```

THM-4084 already treats the balanced, single-edge, frustration-two, and
frustration-three strata. Combining it with `(10)` proves the frustration
bound in `(5)`. THM-4083's balanced-deletion trichotomy proves the simultaneous
condition `b(H)=0` for any additional equality or counterexample class.

In THM-4078's canonical all-simple-cycle operator, the corresponding
Laplacian eigenvalue is `2S_5(H)`. Thus no frustration-four orbit can lower
the proposed D5 spectral gap `2A_n` or enlarge its equality space.

## 5. Exact and independent audits

The primary symbolic path enumerates every one of the

```text
binom(28,4)=20,475
```

four-edge subsets of `K_8`, recovers exactly the eleven component/degree
signatures in `(9)`, applies `(7)--(8)` in exact SymPy arithmetic, and checks
the nonnegative shifted polynomials in Section 3. As literal controls, it
independently materializes every simple cycle of lengths three through six at
the first three admissible orders for each type and matches every symbolic
`c_k^-` value.

The independent path does not use SymPy, parity inclusion--exclusion, or the
linear-forest formula. It hard-codes one representative of each of the eleven
types, enumerates simple cycles literally, and evaluates five consecutive gap
values. The gap has degree at most four because every negatively counted cycle
contains at least one fixed support edge and `k<=6`. Five values therefore
determine it. Forward differences give its exact Newton expansion

```text
sum_(j=0)^4 a_j binom(t,j),                             (11)
```

with every `a_j>=0`; only the star has `a_0=0`, and its other coefficients
are positive. A sixth literal value in every row is a hostile interpolation
control. The independent route also re-enumerates the complete `K_8` type
inventory and checks the center switch directly.

Normal and optimized runs of both paths byte-match their frozen outputs.

## 6. Scope, connection contract, and replay

```text
source:       a minimum negative-edge support in one switching class,
target:       its length-3/4/5/6 cycle character and cumulative D5 value,
map:          support-subset containment -> parity Mobius sum (7),
preserved:    every c_k^-, S_5, switching class, and spectral character,
lost by S_5:  support type, separate cycle lengths, and frustration index,
sidecar:      one minimum support graph plus the full (c_3^-,...,c_6^-) row,
hostile:      K_(1,4) at n=6 ties but switches to one negative edge,
decisive test: all eleven support types with exact nonnegative gap polynomials.
```

This theorem treats only the cumulative `D=5` all-simple-cycle envelope. By
MISTAKE-495 and THM-4069, it does not apply to a basis-dependent fundamental-
cycle quotient. It excludes frustration index four; it does not control the
frustration-at-least-five, no-balanced-deletion core. Hence the all-order D5
single-edge conjecture and its spectral consequence remain **OPEN**.

Primary replay:

```bash
python3 -B 04-computation/even_graph_four_edge_d5_frustration_firewall_thm4200.py
python3 -O -B 04-computation/even_graph_four_edge_d5_frustration_firewall_thm4200.py
```

Independent replay:

```bash
python3 -B 04-computation/even_graph_four_edge_d5_frustration_firewall_thm4200_independent_audit.py
python3 -O -B 04-computation/even_graph_four_edge_d5_frustration_firewall_thm4200_independent_audit.py
```

Both streams match the outputs and SHA-256 values in the front matter.
**QED.**
