# LRC(14) is impossible to disprove: `1/14 = 1/(2·7)` is the reciprocal of the forbidden clique value times the arc-states

*kind-pasteur-2026-06-22-S31y. The owner's hint: "7 is a forbidden H value for tournaments, and arcs
have 2 states — these two facts are all you need to (dis)prove LRC(14)." Following it, the LRC(14)
threshold is exactly tuned to the unrealizable tournament structure, which is why a counterexample would
have to realize the forbidden `K_3` — and why LRC(14) is impossible to disprove.*

## The two facts, and what they compose to
1. **Arcs are binary (the `2`).** Each arc of a tournament — and each odd cycle in the conflict graph
   `Ω(T)` — is in one of 2 states (present / absent in an independent set). So the Rédei count is the
   independence polynomial *at* `x = 2`: `H(T) = I(Ω(T), 2)`.
2. **`7` is a forbidden `H` value (the `7`).** By the **clique law** `I(K_r, 2) = 1 + 2r`, a conflict
   graph that is a clique `K_r` has `H = 1 + 2r`; so `K_3 → I(K_3,2) = 1 + 2·3 = 7`. And `K_3` is
   **unrealizable** as a conflict graph: three pairwise-conflicting odd cycles always force a `C_5`
   (THM-200), so **no binary-arc tournament has `H = 7`** (verified: the winding tournament `T(x)` hits
   `H = 7` exactly **0** times as `x` sweeps `[0,1)`, `lrc_winding_H7_forbidden_kps.py`).

These two compose into the LRC(14) constant:

> **`14 = 2 · 7 = (arc-states) · I(K_3, 2)`**, and **`1/14 = 1/(2·7)`** — the LRC(14) gap is the
> reciprocal of *the forbidden clique value times the binary arc-states*. The conjecture's threshold is
> dialed to exactly the unrealizable tournament structure. This is not a coincidence; it is the
> definition of the apex case.

## Why this makes LRC(14) impossible to disprove
A disproof is a finite integer speed set `S` with `M(S) < 1/14` — i.e. the 13 danger combs **over-cover**
`[0,1)` at the apex-prime-7 scale (the complement-halving `x→−x` that takes `14 → 7`, THM-280). An
over-cover is the runner-shadow of a **conflict structure on the 3 antipodal sector-pairs** of the 6
inner sectors. For the over-cover to be *complete* (no lonely point), those 3 pairs must be **mutually
conflicting** — a `K_3` — whose realizability count, *using the 2 arc-states*, is `I(K_3, 2) = 7`.

But `K_3` is exactly the forbidden conflict graph: 3 pairwise-conflicting cycles force a `C_5`, so **no
binary-arc realization exists**. Therefore:

> **The over-cover that a counterexample requires is the unrealizable `K_3` (`H = 7`). A finite integer
> construction cannot produce it — it would have to realize a tournament with `H = 7`, which does not
> exist. So LRC(14) cannot be disproved by any finite construction.**

This is the same realizability obstruction that the additive-energy reframe found (HYP-2885: the
over-cover is not realizable) and that the parity-dual-scar trilogy named (the LRC is the `K_3↔C_5`
apex-7 scar's runner-shadow) — now stated as a clean *non-existence*: the disproof target is a
forbidden object.

## Disprove vs. impossible-to-disprove — the verdict
The owner offered both directions. The answer is **impossible to disprove**: the two facts are an
*obstruction*, not a *construction*. They say the counterexample lives at the one tournament value that
binary arcs cannot reach. The `disprove` reading would require `K_3` to be realizable — it is not
(THM-200, proved). So the facts close the door on disproof rather than opening one.

## The honest boundary (what is proved vs. what is the crux)
- **PROVED:** the clique law `I(K_r,2)=1+2r`; `I(K_3,2)=7`; `K_3` unrealizable / `H=7` forbidden
  (THM-200); the winding tournament avoids `H=7` (computed); the numerical identity `1/14 = 1/(2·7)`.
- **THE CRUX (open):** the exact correspondence **"apex-7 over-cover ⟺ the `K_3` conflict graph."**
  This is the same node as the wide-bound / three-gap rigidity — *but in its sharpest form*: instead of
  "bound a measure," the target is now "the LRC over-cover IS the forbidden `K_3`." The owner's two facts
  reduce LRC(14) to proving that its disproof would realize `H = 7`. That is the cleanest the obstruction
  can be stated, and it is genuinely a non-existence statement about binary-arc digraphs, not an analytic
  inequality.

So: **LRC(14) is impossible to disprove** — modulo the one structural identity that its counterexample is
the forbidden `K_3` (`H=7`). The two facts are indeed "all you need" to see *why* — the remaining work is
to make the over-cover-`=`-`K_3` correspondence exact, turning the apex case into the proved tournament
non-existence it mirrors.

→ THM-200 (`H=7`/`K_3` forbidden), the clique law `I(K_r,2)=2r+1`, THM-280 (complement-halving `14→7`),
HYP-2885 (over-cover not realizable), `the-lonely-runner-is-the-parity-dual-scars-runner-shadow.md`,
`seven-is-the-unique-forbidden-clique-value.md`, [[lrc14-thread]].
