> **REFUTED / HISTORICAL — read MISTAKE-231 first.** There is no single
> repo-wide entropy invariant here. The displayed fiber size is only a
> pointwise Hartley ambiguity after fixing a finite universe and observable;
> the CF digit mean and moment depth are different quantities. Already at
> `n=5`, the regular tournament has score-fiber size `1`, not the claimed
> maximum, while `(1,2,2,2,3)` has a three-class non-cospectral fiber. Retain
> the small exact tables and elementary score-distribution calculation as
> separate facts. This reflection proves no LRC statement.

# Arithmetic entropy is a repo-wide invariant: the rigid extremum is the zero-entropy point

*boxeph-2026-07-21-S218. Owner: extend the arithmetic-entropy idea (S217) and apply it to as many repo
pieces as possible. Builds on S217 (class number = arithmetic entropy) and its correction (MISTAKE-225 — the
LRC *gate* is not a binary form; the general principle stands), S213/S214 (tournament chirality / rank-11
vertex), S206 (Fibonacci foil), S211 (detection depth), kps (reconstruction wall at 7). Verified in
`04-computation/arithmetic_entropy_across_the_repo_boxeph_S218.py`.*

## One definition

> **`H_arith(X ∣ L) = log₂ |{X' : L(X') = L(X)}|`** — the number of bits of a *global* object `X` that are
> **hidden** from a chosen *local* invariant `L(X)`. **Zero** entropy = local determines global (**rigid**);
> **positive** entropy = a genuinely hidden global object.

S217 found one instance (class number vs Legendre). The point of this session: **this is a repo-wide
invariant**, and its two extremes organize every recent thread — the *rigid* extremum (the AP / transitive
tournament / Heegner `h=1`) is the **zero-entropy** point, and the *difficulty* of each problem is exactly
its **positive-entropy hidden object**.

## Four instances, verified

**(1) Binary forms ∣ genus — refining S217.** The hidden entropy is not the whole class number but its
**deep (within-genus) part**, because the **genus is locally detectable** (by congruences — the genus
characters). Verified `h(D) = (#genera) × (forms/genus)`:

| disc | `h` | genera (local) | deep (hidden) | `H_genus` + `H_deep` |
|---|---|---|---|---|
| `−3,−7,−11` | 1 | 1 | 1 | `0 + 0` — Heegner, rigid |
| `−15,−20,−24,−84` | 2–4 | =`h` | 1 | `H + 0` — **pure genus** (congruence-visible, no hidden bits) |
| `−23` | 3 | 1 | 3 | `0 + 1.58` — **pure deep** (invisible to any congruence) |
| `−47` | 5 | 1 | 5 | `0 + 2.32` — pure deep |

So the *truly* hidden arithmetic entropy is `H_deep = log₂(h / #genera)` — the odd class group, the part no
local (mod-anything) test can see (the Hilbert class field). `−23,−47` are pure hidden information; `−15` has
none (all its `h=2` is genus, i.e. congruence-detectable).

**(2) Tournaments ∣ score sequence.** `H = log₂(#iso classes sharing a score sequence)`. Verified: the
**transitive tournament** (score `0,1,…,n−1`) is the **unique** realization (Landau) — `H=0`, *rigid*, the
AP/nullcone/rank-11 vertex (S214). The largest fibers sit at the **near-regular** scores (`n=5`: score
`(1,2,2,2,3)` has 3 classes, `H=1.58`) — the hidden reconstruction entropy, and exactly where kps's
"reconstruction wall at 7" (cospectral, local-invariants-fail) lives.

**(3) Reals ∣ continued-fraction prefix.** The CF is the optimal code for a real; the partial-quotient sizes
are the information. Verified: **golden** `[0;1,1,…]` has geometric mean `1 ≪` Khinchin `2.685` — the
*worst-approximable*, minimal CF-entropy, the **LRC foil** (S206); the extremal **`t*=14/183=[0;13,14]`** has
geometric mean `13.5 ≫` Khinchin — *well-approximable*, maximal CF-entropy, the **extremal**.

**(4) Nullcone ∣ moment depth.** The **detection depth** is the moment depth at which local data pins the
global object — a *certificate entropy*. Verified-by-structure: **LRC(14)** has a bounded danger alphabet
`X∈{0,…,13}`, so its Bonferroni/moment certificate **terminates at finite depth `~5`** (klein-S389/THM-671)
— *finite* arithmetic entropy; **GMC/NC2** has *unbounded* radial degree, so **no finite depth** determines
the nullcone (Watson/Laplace determinacy, S211) — *infinite* arithmetic entropy.

## The dual entropies (a structural bonus)

Two entropies act on a tournament in *opposite* directions (verified): the **score-distribution** entropy
(how spread the scores are) is **maximal** for the transitive (`0..n−1`, `log₂ n`) and zero for the regular
(delta); the **reconstruction** entropy (`H_arith`) is **zero** for the transitive (unique) and **maximal**
for the regular (huge fiber). So the AP/transitive is simultaneously **maximum order (spread) and zero
hidden information** — *that is what rigidity is*: an object so ordered its local invariant pins it exactly.
The regular/Paley pole is the reverse: minimal spread, maximal hidden info.

## The unification

```
   H_arith(X | L) = log2 |{X' : L(X')=L(X)}|   (hidden global bits beyond a local invariant)

   thread          local invariant L      RIGID (H=0)              HIDDEN entropy (H>0)
   ─────────       ─────────────────      ───────────              ────────────────────
   binary forms    genus / Legendre       Heegner -3,-7,-11        deep class group (-23,-47)
   tournaments     score sequence         transitive (=AP vertex)  near-regular cospectral fiber
   reals (LRC)     CF prefix              golden (foil, min)       t* tail (extremal, max)
   nullcone        moments to depth d     — (depth where H=0)      LRC finite / GMC infinite
```

**Every rigid extremum in the repo is a zero-arithmetic-entropy point** — the AP, the transitive tournament,
the Heegner `h=1` form, the reify-ladder nullcone vertex (S214) — an object whose *local* invariant already
*determines* it. And **the difficulty of every thread is its positive-entropy hidden object** — the deep
class group, the cospectral fiber, the CF tail, the deep moment. Rigidity is *why* the extremal is unique;
hidden entropy is *where the proof still has to go*. (This is the honest, gate-independent core surviving
S217's correction: it is an organizing information-theoretic lens across the repo's objects, not a claim
about codex's polygonal LRC gate.)

## Scope

The genus/class-group split, Landau's score-sequence rigidity, the Khinchin/Lévy CF facts, and the detection
depths are classical/verified. The contribution is the **unification** — one local-to-global information
deficit `H_arith` realized across binary forms (genus-refined S217), tournament reconstruction, Diophantine
approximation, and moment nullcones — and the observation that **the repo's rigid extrema are exactly its
zero-entropy points, its difficulties exactly its hidden-entropy objects.** An organizing invariant and a map
of where the hidden information lives, not a new proof step.

Links: HYP-8875, HYP-8870, THM-587, THM-1750,
[[class-number-is-arithmetic-entropy-hidden-binary-forms-and-why-7-is-rigid-boxeph-S217]],
[[the-rank11-ap-core-is-the-achiral-vertex-where-the-rank-or-euler-frontier-meets-boxeph-S214]].
