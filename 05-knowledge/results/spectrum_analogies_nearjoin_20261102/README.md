# Spectrum analogies: the near-join sweep, count spectra, and strong realisers

klein-2026-11-02. Computations behind the reflection [`partners-windows-and-bands-analogies-for-the-completeness-of-the-hamiltonian-path-spectrum-klein-20261102.md`](../../../07-reflections/partners-windows-and-bands-analogies-for-the-completeness-of-the-hamiltonian-path-spectrum-klein-20261102.md). **Status: PROVED elementary (the near-join formula `H(T′) = h(T) + 2^k B(T, t₀)` and `B(T′, z) = H(T′)`, verified exhaustively on small orders) + FINITE-EXACT (count spectra through order 9; near-join coverage with bases of order `≤ 8`) + HYP-9453 (open: the inductive sweep, the density of `B`-sets, the two-strong-realisers law). Nothing about the Chan–Pak conjecture is closed.**

## 1. The near-join family ([`nearjoin.py`](scripts/nearjoin.py), [`nearjoin8.py`](scripts/nearjoin8.py))

`T′ = T → M_k → z` with the single arc `z → t₀` reversed (`M_k` a transitive chain of `k ≥ 1` vertices). Then `H(T′) = h(T) + 2^k B(T, t₀)` with `B(T, t₀) = Σ_{S ⊆ V(T) − t₀} h(T[S]) · h_{t₀}(T[S^c ∪ t₀])` (checked on 170 cases of orders `≤ 5`, `k ≤ 2`, zero mismatches), `B(T, t₀) ≥ h(T)` with equality at a sink (the join `T → E_{k+2}`), and `B(T′, z) = H(T′)` (1,260 cases, orders `≤ 6`, `k ≤ 3`); `B(T′, t₀)` is not a function of `(h, B, k)`.

Coverage of odd numbers by `{h + 2^k B : k ≥ 1}` over all `(T, t₀)`:

| base orders | distinct `(h, B)` pairs | odd numbers missed below the first large gap | first large gap |
|---|---|---|---|
| `≤ 7` | 2,378 | `1, 7, 11, 13, 21, 47, 67` | 1,439 |
| `≤ 8` | 43,015 | `1, 7, 11, 13, 21, 47, 67` | 22,229 |

With bases of order `≤ 8` every odd number in `(X/2, X]` is covered for `X ≤ 20,000` and 96.8 % for `X = 50,000` ([`nearjoin8.out`](data/nearjoin8.out); the pairs are in [`nearjoin_pairs_le8.json`](data/nearjoin_pairs_le8.json)). The `B`-sets per `h` (`h ≤ 31`, bases `≤ 7`) have 14–42 members and always contain `h`.

## 2. Count spectra ([`spectra.py`](scripts/spectra.py), [output](data/spectra.out))

Directed Hamiltonian-cycle counts: no value below the maximum is missing at orders 3–6 and 8; order 7 misses `16, 18–23` of 24, order 9 misses a band from 215 of 249. Paths from a fixed start vertex (orders `≤ 8`): first gap 14 at order 6, then only top bands. 3-cycle counts: every value from 0 to the maximum at every order `≤ 9`. Bipartite tournaments on `≤ 7` vertices: `H ∈ {0, …, 8, 12}`, first missing 9. **None of these spectra forbids 4 with 11 or 13.** Strong realisers through order 9: the strong singletons are exactly `3, 5, 11, 13` and twenty-five band values above 2,881; every attained odd `h ∈ [15, 2881]` has at least two strong realisers; 47 has two, 67 six.

## 3. Reproduction

```bash
cd 05-knowledge/results/spectrum_analogies_nearjoin_20261102
cc -O2 -o scripts/hampath scripts/hampath.c
python3 scripts/spectra.py        # about 2 min (gentourng through order 9)
python3 scripts/nearjoin.py       # about 2 min (bases <= 7)
python3 scripts/nearjoin8.py      # about 25 min (bases <= 8)
shasum -a 256 -c SHA256SUMS
```

Receipts for `spectra.py` and `nearjoin.py` are byte-identical under `python3` and `python3 -O`.
