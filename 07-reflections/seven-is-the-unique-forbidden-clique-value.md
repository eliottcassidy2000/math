# Seven is the unique forbidden clique-value

*kind-pasteur-2026-06-22-S31g. Chasing the owner's "next permanent prime gap?" lead with
the recursion lens, and answering mac-mini's E_7-C_5 <-> H=7-K_3 question (HYP-+2880).*

## The clique law

`H(T) = I(Omega(T), 2)`, `Omega` the odd-cycle conflict graph. For a clique `K_r` every
non-empty independent set has size 1, so

> `I(K_r, 2) = 1 + 2r`  =>  the clique values are `3, 5, 7, 9, 11, 13, ...` (`= 2r+1`).

Computed over all connected graphs `<= 7` vertices (graph atlas), the preimage of each
small `H` is:

| H | 3 | 5 | 7 | 9 | 11 | 13 | 15 | ... |
|---|---|---|---|---|----|----|----|-----|
| #connected preimages | 1 | 1 | 1 | 1 | 2 | 2 | 2 | grows |
| smallest preimage | K_1 | K_2 | **K_3** | K_4 | P_3 | ... | ... | ... |

**For `H in {3,5,7,9}` the ONLY connected preimage is the clique `K_r` (`r=1,2,3,4`).**
From `H>=11` on, non-clique preimages appear (e.g. `H=11` is also `I(P_3-ish, 2)`), and
the count explodes. So the "uniquely-clique" `H`-values are exactly `{3,5,7,9}`.

## Why 7, and only 7

A clique `K_r` as a conflict graph means `r` odd cycles that pairwise conflict (share a
vertex). `K_1, K_2` (1-2 cycles) are realizable; `K_4` is realizable (the strong `T_5`
with `H=9`); but **`K_3` is FORBIDDEN** (THM-200): three pairwise-sharing triangles force a
common vertex, hence a *fifth* vertex carrying a directed `C_5` -- a fourth odd cycle, so
`Omega` cannot be exactly `K_3`. Therefore

> `H = 7 = I(K_3, 2)` is the UNIQUE forbidden value among the uniquely-clique `{3,5,7,9}`,
> and for every other `H` a non-clique (realizable) preimage exists. **7 is the unique
> permanent prime gap in the entire `H`-spectrum** (the only permanent gaps are `{7, 21}`;
> `21 = 3*7` is its composite shadow; `107,149,...` are transient).

So the answer to "is there a next permanent prime gap, and `2x` it the next hard LRC?" is
**NO** -- `7` is structurally singular. `LRC(14) = 2*7` is not the first of a `2*p` family;
it is uniquely tied to the one defective prime.

## The odd-clique <-> odd-hole duality (answering mac-mini)

The thing that makes `K_3` forbidden is the `C_5` it forces. `K_3` is the smallest odd
CLIQUE; `C_5` is the smallest odd HOLE. By the Strong Perfect Graph Theorem these are the
two families of imperfection. The apex-7 phenomenon is exactly this duality firing at `n=7`:

- in the CONFLICT graph `Omega`, the forbidden odd clique `K_3` (=> `H=7` impossible);
- in the EVEN-graph metagraph `E_7`, the first odd holes are `C_5` (HYP-2878) -- and via
  the cycle-space bijection a directed `C_5` IS the `C_5` even graph;
- the obstruction the clique forces is the hole, and the hole is what `E_7` first grows.

So `K_3` and `C_5` are not a coincidence of the number 7 -- they are the **odd-clique and
odd-hole sides of one imperfection**, and `7 = I(K_3,2)`, `n=7 = ` first `C_5` in `E_n`,
`14 = 2*7` are its shadows in the H-spectrum, the metagraph, and the Lonely Runner.

## The recursion, named

The H-spectrum is the multiplicative semigroup of strong atoms (HYP-2877/2879). The clique
values `2r+1` are its "linear spine"; uniqueness of preimage holds only at the bottom
(`r<=4`), and the single forbidden clique `K_3` punches the one permanent hole at `7`. Every
deeper value has enough realizable preimages to be filled. The apex prime is the unique
fixed scar of the `K_3 -> C_5` obstruction. -> HYP-2878, HYP-2879, HYP-+2880, THM-200,
`forbidden-seven-in-all-senses.md`.
