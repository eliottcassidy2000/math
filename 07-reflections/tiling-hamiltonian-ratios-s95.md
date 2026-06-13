# Tiling-Hamiltonian Ratios S95

**Session:** codex-2026-05-29-S95
**Script:** `04-computation/tiling_hamiltonian_ratios_s95.py`
**Stored output:** `05-knowledge/results/tiling_hamiltonian_ratios_s95.out`

## Vocabulary

Fix the explorer's base Hamiltonian path

```text
P0 = 0 -> 1 -> ... -> n-1.
```

The tournament-tiling-explorer is the **Szele slice**: the set of tournaments
that contain this one ordered path. It has

```text
m = C(n-1,2)
2^m tilings.
```

For an isomorphism class `C`, define:

- **path-fiber** `F(C)`: number of explorer tilings landing in `C`.
- **tiling ratio** `rho(C) = F(C)/H(C)`.
- **automorphic tax** `tau(C)=|Aut(C)|`: the symmetry factor compressing
  Hamiltonian paths into explorer tilings.
- **rigid class**: `|Aut(C)|=1`; then every Hamiltonian path contributes one
  explorer tiling.
- **taxed class**: `|Aut(C)|>1`; then Hamiltonian paths collapse in orbits.

## Main Theorem

For every tournament isomorphism class `C`,

```text
F(C) = H(C) / |Aut(C)|.
```

Equivalently,

```text
rho(C) = F(C)/H(C) = 1/|Aut(C)|.
```

### Proof

Count pairs `(L,P)` where `L` is a labeled copy of a representative of `C` and
`P` is a Hamiltonian path in `L`.

There are `n!/|Aut(C)|` labeled copies and `H(C)` Hamiltonian paths in each, so
the total number of pairs is

```text
(n!/|Aut(C)|) H(C).
```

The symmetric group is transitive on the `n!` possible ordered vertex paths.
Therefore each ordered path appears in exactly the same number of pairs. The
explorer fixes one path `P0`, so the number of labeled copies containing `P0` is

```text
((n!/|Aut(C)|) H(C)) / n! = H(C)/|Aut(C)|.
```

But labeled copies containing `P0` are exactly the explorer tilings in class `C`.

## Consequences

1. **Total slice identity**

   ```text
   Σ_C H(C)/|Aut(C)| = 2^C(n-1,2).
   ```

   This is the exact meaning of "Hamiltonian paths are tilings" in the explorer:
   the fixed-path slice partitions into path-fibers.

2. **Ratio spectrum**

   The only possible tiling ratios are reciprocals of automorphism group orders.
   Exact data through `n=7`:

   ```text
   n=3: {1, 1/3}
   n=4: {1, 1/3}
   n=5: {1, 1/3, 1/5}
   n=6: {1, 1/3, 1/5, 1/9}
   n=7: {1, 1/3, 1/5, 1/7, 1/9, 1/21}
   ```

3. **Complement-merged explorer nodes**

   If the UI merges complement/transpose pairs, the merged path-fiber is

   ```text
   F_merged(C) = H(C)/|Aut(C)|        if C is self-complementary,
               = 2H(C)/|Aut(C)|      otherwise.
   ```

   Complement does not change `H` or `|Aut|`; it only doubles the fiber when the
   complement class is distinct.

4. **Generic collapse disappears**

   For rigid classes, `F=H`. The mismatch between tiling count and Hamiltonian
   path count is entirely a symmetry phenomenon.

5. **Automorphic defect**

   Define

   ```text
   D_n = Σ_C (H(C)-F(C)) = Σ_C H(C)(1 - 1/|Aut(C)|).
   ```

   This measures how much Hamiltonian-path mass is lost to symmetry when passing
   into the fixed Szele slice.

## Exact Sequences

For `n=3..7`, exact enumeration gives:

```text
classes:
  2, 4, 12, 56, 456

ΣF = total explorer tilings = 2^C(n-1,2):
  2, 8, 64, 1024, 32768

ΣH over isomorphism classes:
  4, 12, 92, 1224, 35180

ΣF·H = path-fiber H moment = Σ H(C)^2/|Aut(C)|:
  4, 32, 632, 29696, 3251200

ΣF² = ordered same-class tiling pairs = Σ H(C)^2/|Aut(C)|²:
  2, 28, 540, 27800, 3186296

ΣC(F,2) = unordered same-class tiling pairs:
  0, 10, 238, 13388, 1576764

automorphic defect Σ(H-F):
  2, 4, 28, 200, 2412
```

These pair sequences look like natural OEIS candidates. They are not ordinary
tournament-class counts; they are second moments of the explorer quotient.

## Cube-Edge Laws

Let `W_ij` be the number of directed single-tile flips from class `i` to class
`j` in the fixed-path tiling cube. Then:

```text
Σ_j W_ij = m F_i
W_ij = W_ji
```

The first identity is just "each tiling has `m` tile flips." The second holds
because every cube edge can be traversed in both directions.

Therefore the quotient flip chain is reversible with stationary distribution

```text
π_i ∝ F_i = H_i/|Aut_i|.
```

This is the exact Markov-chain meaning of the explorer: the random tile-flip
walk samples isomorphism classes by Hamiltonian-path mass, taxed by symmetry.

Exact edge sequences for `n=3..7`:

```text
silent cube edges:
  0, 5, 20, 264, 5504

cross-class cube edges:
  1, 7, 172, 4856, 240256

edge tiling energy Σ_edges F_i F_j:
  1, 27, 12296, 3832396, 2482098524

edge H energy Σ_edges H_i H_j:
  3, 71, 18612, 4393168, 2585540208
```

The silent edge sequence is the "neutrality" sequence of the Szele slice. The
energy sequences measure how much the cube quotient preferentially connects
large path-fibers.

## Structural Pattern

The explorer is not merely a visualization of tournaments. It is a quotient of
a Hamiltonian-path incidence design:

```text
labeled tournament copies  x  Hamiltonian paths
             |
             v
fixed-path Szele slice / isomorphism
```

Every visual feature of the explorer must respect this incidence design:

- node area should scale with `H/|Aut|`, not just `H`;
- complement merging either preserves or doubles this mass;
- single-tile flip weights have row sum `mH/|Aut|`;
- stationary mass of the quotient walk is exactly `H/|Aut|`;
- all deviations from `F=H` are automorphism shadows.

## New Hypotheses

### H1. Pair-Moment Asymptotics

The same-class pair moment

```text
K_n = Σ_C F(C)^2
```

should be asymptotic to the second moment of `H` over random labeled tournaments,
because almost all classes become rigid. More precisely, the automorphic tax
should contribute a lower-order correction controlled by the same odd-cycle
cycle-index terms that appear in A000568.

### H2. Edge-Energy Concentration

The ratio

```text
E_F / E_H = (Σ_edges F_i F_j)/(Σ_edges H_i H_j)
```

should tend to `1` as `n` grows, because almost all adjacent classes are rigid.
Observed values:

```text
n=3: 0.3333
n=4: 0.3803
n=5: 0.6606
n=6: 0.8723
n=7: 0.9600
```

This is a quantitative form of "tilings and H become the same coordinate."

### H3. Automorphic Denominator Spectrum

The observed tiling-ratio denominators are exactly observed tournament
automorphism orders in the fixed-path slice. The first appearances

```text
3 at n=3, 5 at n=5, 9 at n=6, 7 and 21 at n=7
```

should be explainable by the minimum vertex support of cyclic automorphism
actions. This gives a new way to read automorphism-group taxonomy from the
explorer without inspecting permutations directly.

## Working Slogan

Hamiltonian paths are the mass; automorphisms are the tax; tilings are the
taxed mass visible in one fixed Szele slice.
