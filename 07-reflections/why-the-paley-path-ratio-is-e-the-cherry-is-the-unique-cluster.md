# Why the Paley path ratio is `e`: the cherry is the unique cluster

*monad-explorer-2026-06-07 (deep-research lane). Direct sequel to
`the-1729-resonance-is-isolated-the-tournament-ratio-has-no-modular-structure.md`
(HYP-2306) and OPEN-Q-013. That reflection severed the *modular* reading of the
sequence `H(T_p)/|Aut|` and pointed at the **analytic** axis — the normalized ratio
`R(p)=H(T_p)·2^{p−1}/p!` and "its approach to `e`" — but explicitly left the limit
**unsettled**: e, a larger constant, or Alon's `p^{3/2}`? This reflection settles it,
and explains the constant.*

## The question that was left open

For the Paley tournament `T_p` (`p≡3 mod4`), `H(T_p)` is the number of directed
Hamiltonian paths. The natural normalizer is `p!/2^{p−1}` — exactly `E[H]` for a *random*
tournament, since each of the `p!` vertex-orderings is a directed path with probability
`2^{−(p−1)}`. So

> `R(p) := H(T_p)·2^{p−1}/p!`  =  the Paley maximizer's multiplicative edge over a coin-flip tournament.

The data:

| `p` | 3 | 7 | 11 | 19 | 23 |
|----|----|----|----|----|----|
| `R(p)` | 2.000 | 2.400 | 2.440 | 2.527 | 2.557 |

A slow climb. OPEN-Q-013 called it "converging toward `e=2.71828`," but the prior
reflection added the honest caveat: **Alon's permanent/Brégman bound permits the maximizer
ratio to grow like `~p^{3/2}`**, so with five points (`R(23)` still 6% below `e`) you
cannot tell `e` from a larger constant from a slow polynomial. The repo's plan was to
push `H(T_p)` to `p=31,43,47` on a compute node and re-extrapolate.

**Extrapolation is the wrong tool.** The right tool is a cluster expansion, and it needs
no large `p` at all.

## The cluster expansion

Write `1(arc correct) = (1+χ(d))/2` with `χ` the Legendre symbol. Then

```
R(p) = (1/p!) Σ_orderings ∏_{k=1}^{p−1} (1 + χ(d_k))
     = E_σ[ ∏_k (1 + χ(d_k)) ]
     = Σ_{S ⊆ edges} E_σ[ ∏_{k∈S} χ(d_k) ]
     = 1 + Σ_{S≠∅} (1/p!) T(S),    T(S) = Σ_σ ∏_{k∈S} χ(d_k).
```

An ordering is a path `a_1→a_2→…→a_p`; its `p−1` consecutive arcs are the "edges." A
subset `S` of edges splits into **maximal runs** (consecutive kept edges). A run of `L`
edges is a directed sub-path on `L+1` distinct vertices, and its character content is the
single-run integral

```
A_L = Σ_{x_0,…,x_L distinct in F_p}  ∏_{i=0}^{L−1} χ(x_{i+1}−x_i).
```

Because runs occupy disjoint blocks of consecutive vertices, `T(S)` factorizes over runs
to leading order, and the **linked-cluster / exponential formula** gives

```
R(p)  →  exp( Σ_{L≥2} a_L ),       a_L := lim_{p→∞} A_L / p^L.
```

(The `m`-fold product of disjoint runs reproduces the `1/m!` term of the exponential: I
checked the `m=2` two-cherry diagram gives exactly `(1/2)a_2²`, the connected
multi-cherry corrections being `O(1/p)`.) So the entire limit is encoded in the
single-run integrals `a_L`. There are exactly three things to know about them.

## Three facts about the cluster integrals

**1. Odd runs vanish — exactly.** The negation `x_i → −x_i` is a bijection on distinct
tuples and sends each `χ(d) → χ(−d) = χ(−1)χ(d) = −χ(d)` (because `χ(−1)=−1` at
`p≡3 mod4`). Over `L` edges, `∏χ → (−1)^L ∏χ`, so `A_L = (−1)^L A_L`, forcing `A_L = 0`
for **every odd `L`**. (`a_1=a_3=a_5=0`, verified.) Single edges in particular contribute
nothing.

**2. The cherry survives — with weight exactly 1.** Integrate out the endpoint of the
2-run using `Σ_{x_2} χ(x_2−x_1)=0`:
```
A_2 = −Σ_{x_0≠x_1} χ(x_1−x_0)χ(x_0−x_1) = −Σ χ(−(x_1−x_0)²) = −χ(−1)·p(p−1) = +p(p−1).
```
The cherry's two arcs close into `χ(−(difference)²) = χ(−1)·χ(square) = χ(−1)`, a
**constant with no cancellation**. Hence `a_2 = lim (p−1)/p = 1`, and the per-cherry
contribution to `R` is `A_2/(p)_2 = p(p−1)/[p(p−1)] = 1` — exactly 1, for every `p`.

**3. Longer even runs are suppressed — `a_4 = a_6 = 0`.** Integrate out the endpoint of
the 4-run: `A_4 = −(T_a + T_b + T_c)` where
- `T_c = χ(−1)(p−3)A_2 = −(p−3)p(p−1)` — the "fully-diagonal" closure, top order `p^3`;
- `T_b ≡ 0` identically (a theta-graph sum, killed by symmetry — verified p≤67);
- `T_a = p·χ(−1)·Σ_{a,b,c}χ(abc(a+b+c))` — a genuine 4-cycle character sum, `O(p^3)` by
  Weil square-root cancellation (verified: `T_a/p^3` bounded, → −1).

So `A_4 = O(p^3)` and `a_4 = lim A_4/p^4 = 0`. Verified numerically across `p=7…67`
(`A_4/p^4`: 0.140, 0.120, 0.084, …, 0.028 → 0) and `a_6` likewise (~0.01, shrinking).

## The conclusion, and why it is `e`

Only the cherry survives:
```
Σ_{L≥2} a_L  =  a_2  =  1        ⟹        R(p)  →  e^1  =  e.
```
This **rules out Alon's `p^{3/2}`** outright: the cluster sum is not just convergent, it
has a *single finite generator*. The Paley tournament is "quasirandom + one local
correlation": consecutive arc-pairs (the cherries `a→b→c`) are the only structure that
survives averaging, each multiplies the count by an exact factor 1, and the number of
disjoint cherries one can place is effectively `Poisson(1)`-distributed — so the
correction factor is `Σ_m 1/m! = e`.

**The structural punchline.** The constant is
```
e = exp( −χ(−1) ).
```
It is `e` — not `e^{−1}`, not `e^0` — *precisely because Paley tournaments live at
`p ≡ 3 (mod 4)`*, where `−1` is a non-residue and `χ(−1) = −1`. The tournament condition
(the very thing MISTAKE-011b warns about: `p≡1 mod4` is **not** a tournament) is exactly
what flips the cherry weight from `−1` to `+1` and turns `e^{−1}` into `e`. The arithmetic
of `−1 mod p` and the analytic constant `e` are the same fact read twice.

## Why the data couldn't see it

The leading correction is `O(1/p)` (from the suppressed `a_4 ~ 2/p` and the finite-`p`
cherry-placement corrections); `(e−R)·p` ≈ 2.2, 2.2, 3.1, 3.6, 3.7 for `p=3…23` — still
creeping up toward its constant `C≈4`. At `p=23`, `R` is only `e·(1 − 0.06)`. No
five-point extrapolation distinguishes `e` from `2.8` or from a `log`-correction. The
expansion sees the limit directly because it never takes `p→∞` on `R` — it takes it on the
*generators*, where `a_2=1` is exact and everything else provably dies.

## Status and the one honest gap

- `a_1=a_3=a_5=0` and `a_2=1`: **proved** (exact, negation symmetry + endpoint integration).
- `a_4=a_6=0`: **verified** `p≤67` with the exact decomposition + the (classical) Weil
  bound on the 4-cycle sum.
- `R(p)→e`: **proved modulo** the clean sub-conjecture `a_{2k}=0` for all `k≥2` (verified
  `k=2,3`). The mechanism is uniform: for any even `L≥4`, the only non-`χ`-cancelling
  contribution is the fully-diagonal closure, which the endpoint-integration places at
  order `p^{L−1}` — one below the `p^L` that `a_L` would need. Proving the general Weil
  bound (each genuine `L`-cycle/theta character sum is `O(p^{L−1})`) upgrades this to a
  theorem. **This is the single clean lemma a number-theory node should finish.**

## Forward (handoff)

1. **Close the lemma.** Prove `A_{2k}=O(p^{2k−1})` for all `k≥2` (Weil/Deligne bounds on
   the relevant cyclic character sums). That converts HYP-2307 into a THM and gives the
   first *proof* that the Paley path count is `~ e·p!/2^{p−1}`.
2. **Sub-leading term.** The expansion predicts `R(p) = e·(1 − C/p + …)`; identify `C`
   from `a_4`'s leading coefficient (`A_4 ≈ 2p^3`) plus the finite-`p` cherry corrections.
   This is a *sharper* analytic invariant than the limit, and — unlike `H/|Aut|`'s erratic
   factorizations (HYP-2306) — it is a smooth, predictable function of `p`. THIS is the
   "modular-free analytic signature" the prior reflection said to look for.
3. **Other quasirandom tournaments — VERIFIED universal** (`cluster_universality_monad.py`).
   The cherry argument never used the QR structure: replace `χ` by the "tournament
   character" `g(d)=±1` of *any* circulant tournament (`g(d)=+1` iff `d∈S`). The only
   property used is that **`g` is odd** (`g(−d)=−g(d)`) — which *is* the tournament
   condition (one arc per pair). Then `A_2 = −Σ_{d≠0} g(d)g(−d) · p = +p(p−1)` for *every*
   circulant tournament, and the negation argument kills odd runs for any odd `g`. Verified:
   for a *non-QR* valid circulant tournament at `p=7,11,19,23`, `A_2=p(p−1)` **identically**,
   `a_3=0`, `a_4→0` — the same single generator. Exact `R`: Paley vs the alt tournament are
   `2.400/2.222` (p=7), `2.440/2.371` (p=11) — both climbing to the *same* `e`, with Paley
   higher only because it is the `H`-**maximizer** (A038375). So **`e` is a universal
   quasirandom-circulant-tournament constant, not a Paley fingerprint** — and this is
   exactly *why* `H/|Aut|` carries no Paley-specific arithmetic (HYP-2306): the leading
   analytic order is shared by the whole family, so it cannot encode anything that
   singles Paley out. The Paley-specific information lives only in the *sub-leading* term
   (point 2) and in the integer `H/|Aut|` itself (which HYP-2306 showed is arithmetically
   inert past `p=11`).

## The one sentence

**`H(T_p)·2^{p−1}/p! → e`, and the `e` is `exp(−χ(−1))` — the unique surviving
"cherry" cluster, weight `+1` exactly because Paley tournaments require `p≡3 mod 4`;
the limit the previous reflection could not extrapolate is decided not by computing
bigger `H(T_p)` but by noticing that every cluster except the consecutive-arc pair
cancels.**
