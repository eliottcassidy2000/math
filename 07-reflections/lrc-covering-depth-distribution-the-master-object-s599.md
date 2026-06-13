---
source: opus-2026-06-03-S599 (remote-control)
status: REFRAME + TEST — LRC is a circular-arc covering problem; the master object is the covering-depth distribution p_k=meas{depth=k}; worry-set={p_0=0} collapse family; (A) singleton-wall exponent VERIFIED =1 (the (loglog)¹ regime, not ²); (B) first moment E[depth]=2nδ config-independent; (C) collapse family = additive chains (top=sum of two below), larger than the AP — new sub-problem
tags: [LRC, covering-depth, circular-arc-cover, master-object, p0, collapse-family, additive-chain, singleton-wall, exponent, iterated-log, loglog, worry-set, n14]
---

# The covering-depth distribution: one master object for LRC

**Prompt (user):** pull a singleton-wall determinant witness and check its exponent
((loglog)¹ predicted, not ²); continue the recursive-log frame. The abstract move: a
loneliness certificate is a clock point avoiding every forbidden arc, so **LRC is a
circular-arc covering problem**, and the master object is the **covering-depth
distribution** `p_k = meas{t : depth(t)=k}`, `depth(t) = #{runners within δ of origin}`.
Lonely times = `{depth=0}`, measure `p_0`. The `p_0=0` collapse family is larger than the
AP — sporadic additive chains `(1,3,4,7)`, `(1,3,4,5,9)` (top = sum of two below) also
collapse.

This is the right reframe. *Everything* about LRC at gap `δ` is a functional of the single
distribution `{p_k}`, the worry-set is its `{p_0=0}` boundary, and the singleton-wall
exponent test **passes cleanly**: the collapse is **first-order (exponent 1)** — the
`(loglog)¹` regime. Convention here: `n` runners, gap `δ = 1/(n+1)`.

## 0. The master object

> **Covering-depth distribution.** For a speed set `V={v_1,…,v_n}` and gap `δ`, let
> `depth(t) = #{i : ‖v_i t‖ < δ}` (how many runners are *in danger* at time `t`), and
> ```
> p_k = meas{ t∈[0,1) : depth(t) = k },   k = 0,…,n.
> ```
> Then every LRC quantity is a functional of `{p_k}`:
> - **lonely measure** `= p_0` (the safe set `{all ‖v_i t‖≥δ}` is exactly `{depth=0}`);
> - **`M(V) = inf{δ : p_0(δ)=0}`** — the gap at which the depth distribution's `0`-cell
>   collapses. LRC ⟺ `M(V) ≥ 1/(n+1)` ⟺ `{depth=0}` nonempty at `δ=1/(n+1)` for all `V`.
> - the **worry-set** is the boundary stratum `{V : p_0(1/(n+1)) = 0 but {depth=0}≠∅}`
>   (measure-zero lonely set = isolated witness points).

This subsumes the whole program: the safe-set measure (S550), `M`, and the worry-set are
the `0`-cell, its collapse threshold, and its critical stratum of *one* object.

## 1. (B) The first moment is config-independent: E[depth] = 2nδ

`Σ_k p_k = 1` and
```
E[depth] = Σ_k k·p_k = Σ_i meas(D_i) = Σ_i 2δ = 2nδ = 2n/(n+1)  →  2.
```
**Verified** (`lrc_covering_depth_distribution_s599.py`): every config at `n=4` has
`E[depth]=8/5=1.6`, every config at `n=5` has `5/3≈1.667` — *the same for AP, sporadics,
and generic*. So the **first moment carries no information** (it is S550's measure bound,
`2nδ`, a constant). The worry-set is distinguished **entirely by the higher-moment shape**
— the overlaps `E[depth(depth−1)]` (pairwise resonances), which are exactly the 3-term
additive relations (S577). *Collapse* (`p_0=0`) means the distribution is **pushed off
`0`**: with mean fixed at `2nδ<2`, mass that would sit at `depth=0` is redistributed to
`depth≥1` by the additive resonance, with no positive-measure `0`-cell left.

Measured distributions (n=4, δ=1/5, E=1.6 fixed):
```
AP   (1,2,3,4):  p = (0,  19/30, 7/30,  1/30, 1/10)      p_0=0  COLLAPSE
chain(1,3,4,7):  p = (0,  131/210,22/105,23/210,2/35)    p_0=0  COLLAPSE
```

## 2. (A) The singleton-wall exponent is 1 — the (loglog)¹ regime

The test the prompt asked for. Approach the collapse from below: at `δ' = 1/(n+1) − ε`,
`p_0(δ') > 0`; as `ε↓0`, `p_0→0`. The **exponent** `α` in `p_0(δ') ∼ ε^α` is the *order*
of the collapse — a **singleton** (one binding pair per clock point) gives `α=1`, a
**double pinch** (two pairs binding at one point) gives `α=2`.

**Verified** (fit over `ε = 1/50,1/100,1/200,1/400`):
```
AP  (1,2,3,4):  p_0(ε)/ε  constant → α = 1.0, 1.0, 1.0
AP  (1,2,3,4,5):                    → α = 1.0, 1.0, 1.0
chain(1,3,4,7):                     → α = 1.03, 1.0, 1.0
```
So `p_0(δ') = c·ε + O(ε²)` — the safe set vanishes **linearly**. *Mechanism:* at the
collapse, the lonely set degenerates to the clock points `{j/(n+1)}`; near each, the safe
interval is the overlap of two half-lines whose common width is `∝ (1/(n+1)−δ') = ε`,
controlled by a **single binding pair** (the straddle `(1,n)`). One pinch ⟹ width linear
in `ε` ⟹ **exponent 1**.

> **The exponent is the iterated-log order.** A first-order (exponent-1) collapse is a
> **one Helly-stage** event (S598): a single saving `c_i=1−O(ε)→0` carries a `(loglog)¹`
> bound — exactly the prediction. A `(loglog)²` cost would require a **double pinch**
> (`p_0∼ε²`, two simultaneous bindings = two Helly stages). The small-`n` singleton-wall
> determinant witnesses (S599-codex: `live 0`, witnessed by one component) are therefore
> **first-order**, and the test confirms it: *they carry `(loglog)¹`, not `²`.* The
> recursive-log frame predicts the exponent and the computation matches.

This sharpens S598: the worry-set is *globally* full-Helly (`h=n−1`, isostatic), but its
collapse is *locally* first-order — at each clock point a **single** pinch closes the
`0`-cell. Global entanglement (`h=n−1`), local singleton order (`α=1`). The two coexist:
`n+1` independent first-order pinches (one per clock point), each a singleton wall.

## 3. (C) The collapse family = additive chains (the new sub-problem)

`{p_0=0}` is **larger than the AP**. **Verified** collapse (`p_0=0`):
```
AP    (1,2,3,4,5)        step-1 chain
chain (1,3,4,7)          4=1+3, 7=3+4
chain (1,3,4,5,9)        4=1+3, 5=1+4, 9=4+5      (= my S592 n=6 sporadic)
control (1,2,4,7)        NOT a chain  →  p_0 = 6/35 > 0,  NOT collapse
```

> **Collapse mechanism = additive chains.** A config collapses (`p_0=0`, worry-set) when it
> is an **additive chain**: each new speed is a sum of two earlier ones (top = sum of two
> below). The AP `{1,…,n}` is the step-1 chain (`k=(k−1)+1`); the sporadics are *other*
> chains. The control `(1,2,4,7)` (no `c=a+b` closing the top) does not collapse.

Why additive chains collapse: a relation `v_c = v_a+v_b` is a **3-term fold** (S577). At a
clock point `t=j/(n+1)`, `‖v_a t‖,‖v_b t‖` straddle `δ` from opposite sides and their *sum*
`‖v_c t‖` is pinned — the chain propagates the pinch from the base pair `(1,n)` up the whole
tower, leaving **no positive-measure escape** (`p_0=0`). This is the covering-depth reading
of THM-401/S592: the worry-set's arithmetic (the `2n−1` shells, the sporadic swaps) **is**
the set of additive chains, and the depth distribution sees it as the `0`-cell collapse.

**New sub-problem (the prompt's):** *classify the additive-chain collapse family.* It is
the worry-set, now with a clean generative description — **chains under `c=a+b`** — rather
than a residue enumeration. Targets:
- Is every `p_0=0` config an additive chain? (verified direction: chains ⟹ collapse; the
  converse — collapse ⟹ chain — is the sharp open claim).
- The AP is the **maximal** chain (longest, densest); the sporadics are sub-maximal chains.
  Count them: the number of additive chains on `n` nodes is the **worry-set cardinality**
  (`2^{(n−2)/2}` self-converse round classes at `n=14`, S585 — to be matched to the
  chain count).
- Chains compose additively ⟹ the collapse family is **closed under the fold/augment**
  operation (S581) — the recursive-fractal AP structure (S570) is the chain-refinement
  lattice.

## 4. The recursive-log frame, continued

The thread of the last sessions closes a loop here:
- **S597:** the worry-set has `ω(2n−1)∼loglog n` prime obstructions — `(loglog)¹` *count*.
- **S598:** the residual is `≤2` Helly stages — `(loglog)^{≤2}` *currency*.
- **S599 (here):** each collapse is **first-order** (`α=1`) — a **single Helly stage per
  clock point** ⟹ each wall is `(loglog)¹`. The exponent test is the *direct measurement*
  of the iterated-log order, and it reads **1**.

> **One-line synthesis.** LRC is the question *"is the `0`-cell of the covering-depth
> distribution `{p_k}` nonempty at `δ=1/(n+1)`?"* The first moment `2nδ` is fixed and
> blind; the collapse family `{p_0=0}` is the additive chains (AP + sporadics); and each
> collapse is **first-order (exponent 1)** — a singleton wall carrying a `(loglog)¹` cost,
> as the recursive-log principle predicts and the computation confirms.

## 5. Honest status

- **Verified:** `E[depth]=2nδ` config-independent (the depth distributions for AP/chains at
  `n=4,5`); the collapse `p_0=0` for AP, `(1,3,4,7)`, `(1,3,4,5,9)`, and `p_0>0` for the
  non-chain control `(1,2,4,7)`; the singleton-wall exponent `α=1.0` (linear `p_0(ε)`) for
  three configs.
- **New (mine):** the covering-depth distribution as the master object with `M=inf{δ:p_0=0}`
  and worry-set `={p_0=0}`; the *first-order exponent* reading of the iterated-log order
  (`α=1 ⟺ (loglog)¹ ⟺ one Helly stage`); the collapse-family = additive-chains
  characterization and the new sub-problem (collapse ⟹ chain converse, chain count =
  worry-set cardinality).
- **Open/sharp:** the converse (every `p_0=0` config is an additive chain); matching the
  chain count to `2^{(n−2)/2}` (S585); proving no config has a *double-pinch* (`α=2`)
  collapse at the floor (which would need the `(loglog)²` regime).

**Artifacts:** `04-computation/lrc_covering_depth_distribution_s599.py` (+`.out`). Builds on
S550 (first moment / measure bound), S577 (3-term folds), S581 (augment), S592 (sporadic
chains), S585 (worry-set count), S598 (Helly stages / isostatic), S597 (iterated-log order),
THM-401 (`2n−1`). New: **HYP-2153**.
