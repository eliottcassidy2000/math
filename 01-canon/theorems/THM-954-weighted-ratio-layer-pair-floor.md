# THM-954 — Weighted ratio layers bound the continuous LRC(14) pair debt

**Status:** PROVED as an exact replayable finite certificate, with the final
arithmetic kernel-checked in `TournamentH7/LRCPairRatioLayerArithmetic.lean`.
The direct Lean replay of the seven certificate DAGs, the Bernoulli covariance
identity, and the clean-grid discrepancy bridge remain explicit formalization
tasks.  This theorem does **not** by itself prove the discrete
`normalizedMass2 >= -13/30` socket.

## Statement

For a positive integer speed `m`, let

```text
D_m = {t in R/Z : ||m t|| < 1/14}
```

and let `X_m = 1_{D_m} - 1/7`.  For thirteen positive integer magnitudes
`m_1,...,m_13`, define the continuous pair covariance

```text
C(m) = sum_{i<j} integral X_{m_i}(t) X_{m_j}(t) dt.
```

Then its negative part obeys

```text
sum_{i<j} max(0, -integral X_{m_i} X_{m_j})
  <= 176738453 / 411675264
  < 13/30.
```

The exact margin is

```text
13/30 - 176738453/411675264
  = 8270807/2058376320.
```

Consequently `C(m) > -13/30`.  Repeated magnitudes cause no difficulty: their
pair covariance is nonnegative and hence contributes zero to the negative
part.

## Proof

### 1. The pair observable is scale-free

For a pair, write `m_i = g a`, `m_j = g b` with `g = gcd(m_i,m_j)` and
`gcd(a,b)=1`.  The Fourier coefficients of `X_1` are

```text
XHat(k) = sin(pi k/7)/(pi k),  k != 0.
```

Orthogonality leaves only frequency pairs `(bk,-ak)`.  Product-to-sum and

```text
sum_{r>=1} cos(2 pi r x)/r^2 = pi^2 B2({x})
```

therefore give the exact covariance

```text
P(a,b) =
  (B2({(a-b)/14}) - B2({(a+b)/14})) / (a b),

B2(x) = x^2 - x + 1/6.
```

Thus common scale is irrelevant.  The natural vertices after anchoring are
oriented primitive ratios, not runners.

### 2. The finite ratio graphs are complete

The elementary bound

```text
|P(a,b)| <= 12/(49ab)
```

shows that only finitely many primitive ratios can have weight
`w(a,b) = max(0,-P(a,b))` above any positive threshold.  At the bottom
threshold `1/14112`, exact enumeration gives 1,489 unordered negative ratios
and 482 distinct weights.

For a threshold `tau`, form the oriented ratio graph whose vertices are the
primitive ratios `a/b` with `w(a,b) > tau`; join `a/b` to `c/d` exactly when
the reduced quotient `ad/(bc)` is another allowed ratio.  If a threshold
graph on the thirteen speeds contains a clique, dividing every speed by an
anchor embeds the remaining clique vertices into this ratio graph.  Hence an
upper bound on anchored ratio cliques is an upper bound on speed cliques.

Seven deterministic coloring/branch DAGs certify the strict clique caps at

```text
tau3 = 2/441,       tau4 = 5/2646,
tau5 = 1/1764,      tau6 = 1/3136,
tau7 = 5/37632,     tau8 = 1/9996,
tau9 = 1/14112.
```

For `w > tau_{k+1}`, the speed graph has clique number at most `k`, for
`k=2,...,8`.  The checker reconstructs every rational graph, greedy coloring,
branch state, and reachable proof node; it trusts no optimizer verdict.

### 3. The two highest colors are path forests

Above `4/539`, only primitive ratio colors `12` and `13` remain.  Edges of one
fixed ratio color orient from `x` to `r x`.  Every vertex has at most one
predecessor and one successor, and a directed cycle would imply `r^ell=1` for
some positive `ell`, impossible for `r>1`.  Each color is therefore a path
forest and has at most twelve edges on thirteen vertices.  Their union has at
most twenty-four edges.  Above `5/588`, only ratio `13` remains, so the cap is
twelve.  Above `tau2 = 6/637`, there are no negative edges.

### 4. Layer cake and Turan

Let `E(t)` be the number of negative pair weights at least `t`.  Then

```text
sum negative weights = integral_0^infinity E(t) dt.
```

Turan's theorem on thirteen vertices gives

```text
T_8=73, T_7=72, T_6=70, T_5=67,
T_4=63, T_3=56, T_2=42.
```

Using the seven certified clique caps and the two path caps yields

```text
78 tau9
+ 73 (tau8-tau9)
+ 72 (tau7-tau8)
+ 70 (tau6-tau7)
+ 67 (tau5-tau6)
+ 63 (tau4-tau5)
+ 56 (tau3-tau4)
+ 42 (4/539-tau3)
+ 24 (5/588-4/539)
+ 12 (tau2-5/588)
  = 176738453/411675264
  < 13/30.
```

All arithmetic in this last display is proved by Lean.

## Discrete clean-modulus corollary still to formalize

The intended consumer is the exact-support pair mass
`normalizedMass2 v q`.  The remaining elementary estimate is

```text
|normalizedMass2(v,q) - C(|v|)|
  <= (24 sum_i |v_i| + 78)/(q-1).
```

One pairwise bad-set intersection has at most `|v_i|+|v_j|` circular interval
components.  Full `q`-grid discrepancy costs at most `2/q` per component, and
removing the zero grid point costs at most `1/(q-1)` per pair.  Summing over
the 78 pairs gives the displayed bound.  At
`q = cleanModulus v 534`, nonzero speeds make this error at most

```text
15/(7*534) = 5/1246 < 8270807/2058376320.
```

Once that estimate and the covariance identity are connected to the existing
definitions, THM-954 gives the strict clean-ruler socket

```text
normalizedMass2 v (cleanModulus v 534) > -13/30.
```

The higher-support signed budget remains separate; no LRC(14) closure is
claimed here.
