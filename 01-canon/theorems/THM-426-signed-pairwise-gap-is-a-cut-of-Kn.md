# THM-426 — The signed pairwise LRC gap depends on the sign gauge only through a CUT of K_{n−1}

**Status:** PROVED (exact, elementary). Exhaustively consistent n=3..6 (`…census_monad_s2.py`,
`…floor_conjecture_monad_s2.py`); bounded follow-up S711 extends the floor search to
`n=6, B<=10` and `n=7, B<=8`.
**Source:** monad-explorer-2026-06-06-S2, monad-explorer-2026-06-12-S711.
**Convention** (repo canon, as THM-420/THM-369): `n` runners total, stationary observer speed `0`,
movers `v_1<⋯<v_{n−1}` distinct positive integers, gap threshold `1/n`,
`‖x‖` = distance to nearest integer. Observer loneliness `M_obs(S)=max_t min_i ‖v_i t‖`;
LRC(n) ⟺ `M_obs ≥ 1/n`.

## The two gauges and their content

The **sign gauge** maps `v_i ↦ ε_i v_i`, `ε ∈ {±1}^{n−1}`.

> **Prop 0 (observer-blindness — exact restatement of T1/HYP-2286).** `M_obs` is sign-invariant:
> `‖(ε_i v_i) t‖ = ‖v_i t‖` for every `t`, since `‖−x‖=‖x‖`.
> **Exhaustive certificate (NEW, exact arithmetic):** over ALL `2^{n−1}` sign patterns and ALL
> gcd-1 speed sets at n=3..6 (B≤9/9/8/7), `M_obs` is constant — **0 violations**
> (`signed_lrc_pairwise_gap_census_monad_s2.out`). So the sign gauge carries **no** observer content.

The gauge IS visible on the **pairwise** structure (s674: "observer-blind, pair-visible"). Define the
**signed pairwise gap** over the `C(n−1,2)` moving pairs:
```
   G_pair(ε,S) = max_t  min_{i<j} ‖(ε_i v_i − ε_j v_j) t‖,
   Gstar(S)    = max_{ε} G_pair(ε,S).
```

## The cut theorem

> **THM-426.** `G_pair(ε,S)` depends on `ε` only through the bipartition
> `(A,B)`, `A={i:ε_i=+1}`, `B={i:ε_i=−1}`. Concretely, the unordered relative-speed magnitude of
> pair `{i,j}` is
> ```
>     |v_i − v_j|   if i,j on the SAME side (within A or within B),
>     |v_i + v_j|   if i,j on OPPOSITE sides (across the cut).
> ```
> Hence the reachable relative-speed multisets are exactly the **cuts of the complete graph
> K_{n−1}**, and
> ```
>     Gstar(S) = max over bipartitions (A,B) of  max_t min_{pairs} ‖w_{ij} t‖,
>       w_{ij} = v_i−v_j (same side),  v_i+v_j (across).
> ```
> Modulo the global flip `A↔B` (which negates every `w_{ij}`, leaving `‖·‖` fixed) there are
> `2^{n−2}` distinct cuts.

**Proof.** `ε_i v_i − ε_j v_j = ε_i(v_i − (ε_j/ε_i)v_j)`. If `ε_i=ε_j` this is `±(v_i−v_j)`; if
`ε_i≠ε_j` it is `±(v_i+v_j)`. Since `‖x‖=‖−x‖`, only same-side-vs-across matters, which is precisely
the cut `(A,B)`. Global flip swaps `A,B` and negates all differences/sums. ∎

So "enumerate sign-reversal patterns" = "enumerate cuts of `K_{n−1}`": same-side pairs keep their
**difference**, cross-cut pairs become a **sum**. The all-`+` (empty cut) pattern is the
**all-differences** set; a maximal cut converts the most pairs to sums.

## Corollary: the n-grid cut witness

> **Corollary.** Fix a cut `ε` and let `W_ε={ε_i v_i−ε_j v_j}` be its relative-speed multiset.
> Then the time `t=1/n` witnesses `G_pair(ε,S) >= 1/n` iff **no** element of `W_ε` is divisible by
> `n`. Equivalently:
> - same-side pairs must have DISTINCT residues mod `n`;
> - across-cut pairs must NOT be additive inverses mod `n`.
>
> Hence any cut with no `n`-multiple relative speed is an immediate certificate for the pairwise
> floor `Gstar(S) >= 1/n`; if every cut has some `n`-multiple, the floor cannot be witnessed on the
> `n`-grid.

**Proof.** At `t=1/n`, each pair contributes `‖w/n‖`. This is `0` exactly when `n|w`, and otherwise
it is at least `1/n`. Apply THM-426's cut dictionary:
same-side pairs contribute differences, across-cut pairs contribute sums. ∎

## Empirical consequences (n=3..6, exhaustive, exact)

- **The all-differences pattern is usually NOT optimal.** A nontrivial cut strictly raises `G_pair`
  in 51/79 (n=4), 52/69 (n=5), 15/21 (n=6) of gcd-1 sets — converting some small differences (e.g.
  the difference `1` of consecutive speeds) into large sums lifts mutual loneliness.
  (`…census_monad_s2.out`.)
- **Why "maximize the smallest relative speed" is the WRONG heuristic (C1, REFUTED).** The cut
  maximizing `min_{i<j}|w_{ij}|` need NOT maximize `G_pair`: at `V=(1,2,5,6)` the `G`-optimal cut has
  a relative speed `1` yet beats a cut with min-rel-speed `3` (`3/13>2/9`). A single small clock is
  not fatal — the maximin is global. (40/79, 30/69 failures; `…optimal_pattern_law_monad_s2.out`.)

## Tie to the shell-partner / synchronization program (THM-425, THM-420)

A **shell-partner** pair (`v_i+v_j ≡ 0 mod q`) is realized as a *sum* exactly when `{i,j}` is placed
**across** the cut. By THM-425 (synchronization) the across-cut sum-clock then has
`‖(v_i+v_j)·k/q‖ = 0` on the entire `q`-grid `{k/q}` — a synchronized zero pair-clock. For the
canonical shell modulus `q = 2n−1`, putting a shell-partner across the cut **blanks the
(2n−1)-grid** for that pair, forcing the pairwise witness time off that grid. So the cut picture
unifies the sign gauge (S674/HYP-2286), the binding-pair = synchronized-shell-partner statement
(THM-425/S1), and the shell-partner witness lemma (THM-420/S700): **a sign pattern is a cut, and a
shell-partner is the pair the cut sends across.**

## The pairwise floor is STRICTLY below the observer floor (HYP-2293 REFUTED)

One might hope the sign gauge always lifts the pairwise gap to the observer floor:
`Gstar(S) ≥ 1/n` (**HYP-2293**). It does NOT.

> **Counterexample (n=6).** `V=(2,3,4,6,8)` (gcd 1) has `Gstar = 3/19 = 0.1579 < 1/6 = 0.1667`.
> Verified two ways: exact maximin over all 16 cuts, and an independent float grid `N=5·10⁶`
> (`signed_lrc_counterexample_verify_monad_s2.out`). The best cut is `A={2,3,4}, B={6,8}`, giving
> relative-speed multiset `{1,1,2,2,8,9,10,10,11,12}`.
>
> **Bounded follow-up (S711).** The `n=6` bounded minimum drops further to `2/13` over `B<=10`,
> attained by `V=(1,4,8,9,10)`, `(2,3,7,9,10)`, `(2,4,5,9,10)`, and the first `n=7` floor failure
> appears already at `V=(1,2,4,6,7,8)` with `Gstar = 4/29 < 1/7` over `B<=8`
> (`signed_lrc_floor_residue_obstruction_s711.out`).

So the signed pairwise gap is **not** governed by `n` alone. The obstruction is **not** imprimitivity
(all 16 cuts here yield `gcd=1` relative-speed sets): it is an **unavoidable cluster of small relative
speeds**. The triple `{2,3,4}` forces differences `{1,1,2}`, and no single cut can send enough of them
across to remove the cluster — moving a vertex across the cut trades its small differences for sums but
creates new small differences elsewhere. The small clocks `‖t‖,‖2t‖` then cap the maximin at `3/19`.

At the `1/n` witness level, the new reduction shows a second layer: all current bounded failures lie in
the **no-n-cut** class, where every cut has some same-side equal residue or across-cut opposite residue,
so the `n`-grid cannot witness the floor at all. But this residue quotient is not the whole story:
already at `n=4,5` there are many no-n-cut sets with `Gstar > 1/n`, so the quotient preserves the
divisibility obstruction and loses the off-grid metric slack.

The true pairwise floor therefore lies **strictly between** the naive pair-count floor
`1/(C(n−1,2)+1)` (here `1/11`) and the observer floor `1/n` (here `1/6`): `1/11 < 3/19 < 1/6`.
Determining `inf_S Gstar(S)` as a function of `n` is the open question (T764). This sharply
distinguishes the signed *pairwise* LRC from the observer LRC: the sign gauge is a cut, and a cut is
not strong enough to lift every set to `1/n`.
