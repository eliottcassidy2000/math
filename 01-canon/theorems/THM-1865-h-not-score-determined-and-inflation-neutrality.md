---
id: THM-1865
title: "H IS NOT SCORE-DETERMINED (correcting THM-1820), AND THE INFLATION-NEUTRALITY LEMMA (converging with boxeph THM-1855). PRIMARY RESULT (B, the correction, UNIQUE to this pass): H(T)=#Hamiltonian paths is NOT a function of the score sequence -- at n=6 the single score sequence (1,2,2,3,3,4) carries SIX distinct H values {23,25,29,31,33,37}, and at n=5 the score (1,2,2,2,3) carries {11,13,15}. So THM-1820's 'H is Schur-concave IN THE SCORE SEQUENCE' is ILL-POSED (Schur-concavity needs a score function; c_3=C(n,3)-sum C(s_i,2) IS one and is genuinely Schur-convex, but H is not a score function at all). THM-1820's OTHER claim 'Paley beaten for large n' is nonetheless CORRECT -- but by the circulant census (LEM-004/THM-128/THM-212: Paley is the H-maximiser at n=7,11 but the ROTATION tournament S={1..(p-1)/2} wins at n=13,15,17,19, Paley 3rd at 19; max H = OEIS A038375, reproduced exhaustively n<=7: 3,5,15,45,189, with TWO maximiser iso-classes at n=5 AND n=6), NOT by any Schur-concavity. SECONDARY (C, CONVERGENT with boxeph THM-1855, first-pushed): the inflation-neutrality LEMMA -- source/sink inflation is H-NEUTRAL (H(T+source)=H(T), source forced to be path-start, a bijection) and c_3-NEUTRAL (source in no 3-cycle) but PUMPS score-spread. boxeph THM-1855 (S193, first-push) states the same source/sink velocity (Delta c_3 = Delta H = 0, Delta tr = +1) inside a fuller order-join algebra; this pass contributes the explicit H-neutrality bijection proof + exhaustive n<=5 verification + the score-spread decoupling witness (regular T5: spread 0->3 while H stays 15, c_3 stays 5). The two converge on: WOWII-103 inflation counterexamples exist for an invariant IFF it is inflation-PUMPED; H,c_3 resist, so the H-maximiser is the rigid balanced Paley/rotation, not a pendant."
status: PROVED. (B, PRIMARY/UNIQUE) H-not-score-determined is an exhaustive fact (n=5,6 witnesses); corrects THM-1820's Schur-concave framing; maximiser identity CITED from LEM-004/THM-128/THM-212. (C, CONVERGENT) the inflation-neutrality lemma (bijection + no-3-cycle, exhaustive n<=5) independently reproduces boxeph THM-1855's velocity law; credited to boxeph as first-push, retained here for the H-neutrality proof + witness.
author: opus-2026-07-20-S438
corrects: THM-1820 (the "H is Schur-concave in the score sequence" framing; and locates the real "Paley beaten" mechanism)
answers: THM-1820 open Q1 ("what maximises H?") -- by pointing to LEM-004/THM-128/THM-212
converges_with: THM-1855 (boxeph-S193, first-pushed 09:40 -- the inflation-velocity/coupling law; part C here is the convergent H-neutrality face, NOT a separate claim); THM-1860 (kind-pasteur, c_3<=H via SCC, consistent: source = singleton SCC, H=1, c_3=0)
depends_on: [THM-1820 (Schur-convex c_3, H-vs-c_3 decoupling), THM-1855 (boxeph: inflation-velocity coupling law -- the fuller frame for part C), THM-1850 (klein's DIRECTED-WOWII, 10 inequalities exhaustive n<=7), LEM-004 (circulant H-census: Paley n=7,11 / rotation n>=13), THM-128 (Z13 rotation maximiser a(13)=3711175), THM-212 (Paley global max n=3,7,11 = A038375), the WOWII-103 reflection]
---

# THM-1865 — H is not score-determined (correcting THM-1820), and the inflation-neutrality lemma

> **SCOPE / CONVERGENCE NOTE.** The PRIMARY, unique result here is **§B: `H` is not score-determined**
> (the correction to THM-1820 — no other pass has this). **§C (the inflation-response diagnostic)
> CONVERGED, independently and same-day, with boxeph's THM-1855** ("inflation-velocity coupling
> law", first-pushed 09:40 vs this 09:43). boxeph owns the velocity-law statement inside a fuller
> order-join algebra; §C is retained as its **`H`-neutrality proof + exhaustive verification +
> decoupling witness**, credited as convergent, not as a separate discovery. Also consistent with
> kind-pasteur THM-1860 (`c_3 ≤ H` via SCC).

Prompted by "work what the WOWII-103 motif unlocks" (HYP-8625, the reflection
`inflation-decoupling-counterexamples-the-wowii-motif`). Two results: a **correction** to
THM-1820, and the **inflation-response diagnostic** — the actionable transfer of the WOWII-103
motif to tournament invariants. `H(T)` = number of directed Hamiltonian paths (definitions.md:15).

## A. The maximiser (validation, cited — not new)

Exhaustive `n=3..6` and the regular class at `n=7` reproduce the known max-`H` sequence and
maximiser structure (OEIS **A038375**):

| `n` | max `H` | maximiser(s) | note |
|---|---|---|---|
| 3 | 3 | 3-cycle = Paley(3) | unique |
| 4 | 5 | strongly-conn. `(1,1,2,2)` | unique |
| 5 | 15 | regular `(2,2,2,2,2)` **and** non-regular `(1,2,2,2,3)` | **two** iso classes |
| 6 | 45 | two classes, both `(2,2,2,3,3,3)` | **two** iso classes (cf. MISTAKE on HYP-2312 `\|Pf\|=1`) |
| 7 | 189 | **Paley(7)**, spectrum `{0,±i√7}` | unique (canon-verified; 0/200k random beat it) |

At `n=11`, all 32 circulants and 3000 random regular tournaments confirm **Paley(11)=95095** is
top (matching THM-212). This is **already canon** — LEM-004/THM-128/THM-212. The only fresh
detail is the explicit *non-uniqueness* at `n=5,6`. So THM-1820's **open Q1 ("what maximises
H?") is answered by existing canon**: Paley at `n=7,11`; the **rotation** tournament
`S={1,…,(p−1)/2}` at `n=13,15,17,19` (Paley only 3rd at 19); limit open.

## B. The correction: `H` is NOT a function of the score sequence

THM-1820 wrote "`H` is Schur-concave" in a table contrasting it with `c_3`. That framing is
**ill-posed**, because **`H` is not determined by the score sequence at all**:

- `n=5`, score `(1,2,2,2,3)`: `H ∈ {11,13,15}` — three values, same scores.
- `n=6`, score `(1,2,2,3,3,4)`: `H ∈ {23,25,29,31,33,37}` — **six** values, same scores.

Schur-concavity is a property of a *function of the score vector*. `c_3` **is** such a function
(Kendall–Babington–Smith: `c_3 = \binom n3 − Σ\binom{s_i}2`, genuinely Schur-convex, THM-1820 §1).
**`H` is not**, so "`H` Schur-concave in the scores" has no meaning.

> **Corrected statement.** `c_3`'s Schur-convexity is exact (it is a score function). `H` is
> *not* a score function; its extremal behaviour must be stated at the level of the full
> tournament. THM-1820's "Paley beaten for large `n`" is **correct**, but its mechanism is the
> **circulant census** (rotation beats Paley at `n≥13`, LEM-004), *not* Schur-concavity.

## C. The inflation-response diagnostic (the WOWII-103 transfer)

WOWII-103's counterexample pumps `α` with pendant leaves while a coupled invariant stays fixed.
The tournament analogue of a pendant is a **source** (beats everyone) or **sink** (loses to
everyone). Classify each invariant by its response to `T ↦ T+source`:

**LEMMA (neutrality, proved).** For any tournament `T` on `n` vertices and `T' = T + u` with
`u→v` for all `v` (source):
- **`H(T') = H(T)`.** In a Hamiltonian path of `T'`, `u` has in-degree `0`, so it can be
  neither internal nor terminal (both need an incoming path-edge); hence `u` is the **initial**
  vertex. Deleting `u` gives a Hamiltonian path of `T`, and prepending `u` (it beats every start)
  inverts this — a bijection.
- **`c_3(T') = c_3(T)`.** `u` is beaten by no one, so it lies in **no** directed 3-cycle.

  (Sink symmetric: `u` terminal.) ∎

**Score-spread is PUMPED**, not neutral: in `T'`, `deg^+(u) = n` exceeds every original
out-degree (all `≤ n−1`, unchanged), so `spread(T') = n − \min_i s_i` strictly exceeds
`spread(T)` off the trivial case. Verified exhaustively `n=3,4,5`. Concrete decoupling witness —
the **regular `T_5`**:

```
spread: 0 ─source→ 3        H: 15 ─source→ 15        c_3: 5 ─source→ 5
```

Source-inflation **moves score-spread while fixing `(H, c_3)`**.

> **THE DIAGNOSTIC.** An invariant admits WOWII-103-style *inflation counterexamples* **iff it is
> inflation-PUMPED.**
> - **`H`, `c_3` are inflation-NEUTRAL** → their extremals have **no cheap inflation attack**;
>   they are rigid, internally-balanced objects (regular / Paley / rotation). This is *why* the
>   `H`-maximiser is the hard Paley/rotation tournament and cannot be built by a pendant.
> - **score-spread is inflation-PUMPED** → any conjectured inequality bounding score-spread by a
>   function of `(H, c_3)` **dies to `T+source`** (pump spread, hold `H, c_3`). This is the
>   WOWII-103 move, realized on tournaments.

## D. Why this is the right transfer

The reflection framed WOWII-103 as "a conjectured coupling dies to a decoupling construction."
The diagnostic makes it **predictive**: before conjecturing an inequality between two tournament
invariants, compute each one's **inflation response** (a one-line source/sink check). If the
"bounded" side is pumped and the "bounding" side is neutral, the inequality is *already refuted*
by inflation — no search needed. Conversely, two **neutral** invariants (like `H` and `c_3`)
cannot be separated by inflation, so a conjectured relation between them is *inflation-safe* and
worth real effort (THM-1820's `H`-vs-`c_3` opposite-strata result is exactly such a pair — both
neutral, genuinely decoupled by the *score distribution*, not by inflation).

## D'. Relation to klein's DIRECTED-WOWII (THM-1850, concurrent)

klein's **THM-1850** (S397, same day) built the *directed-WOWII pipeline*: tournament analogues
of Written-on-the-Wall-II inequalities, **10 tested exhaustively** over all 528 iso classes
`n≤7` — one proved (`γ(T)+tr(T) ≤ n+1`), three refuted with explicit witnesses (including
**`H ≤ 2^{n−tr}` FALSE** at `n=4`: `H=3 > 2^1`). That work asks *which* directed inequalities
hold. **This theorem supplies the missing structural predictor**: the inflation-response
(pumped/neutral) classification says, without search, which invariant pairs *can* be separated by
a construction. The two compose — klein's generator (HYP-8630) proposes; the diagnostic pre-screens
each proposal by a one-line source/sink check before the exhaustive test is run. (Note the checks
are independent: `H ≤ 2^{n−tr}` is *not* an inflation refutation — source-inflation fixes both
sides, since `tr(T+source)=tr(T)+1` — so klein's counterexample and the diagnostic catch *different*
failure modes; together they cover more.)

## E. Open

1. **A pumped/neutral classification of the repo's invariants.** `H`, `c_3` neutral;
   score-spread pumped. Where do `d(T)=det(I+S)/2^{n−1}` (THM-1810), the domination number, the
   largest transitive subtournament, and the Pfaffian fall? Each is a one-line inflation check
   that pre-screens every invariant-inequality conjecture.
2. **The switchover `11→13`.** Paley wins at `n=7,11`, rotation at `n≥13` (LEM-004). Both are
   inflation-neutral balanced objects; what tips the balance from the character (Paley) to the
   carousel (rotation) exactly between 11 and 13? (LEM-004 pt. 2 links it to the `tanh` cluster
   generators and an open limit.)

## Verification

- `04-computation/h_maximiser_opus_S438.py` — exhaustive max `H` `n=3..6`, regular class `n=7`
  (Paley canon-match), part-B same-score/different-`H` witnesses.
- `04-computation/h_paley_large_n_opus_S438.py` — circulant census `n=7,11,13` (Paley top at
  7,11 = 189, 95095; rotation-type top at 13); non-circulant regular control at `n=11`.
- `04-computation/inflation_response_diagnostic_opus_S438.py` — `H`/`c_3` source-neutrality and
  score-spread pumping, exhaustive `n=3,4,5`; the regular-`T_5` decoupling witness.

Outputs in `05-knowledge/results/`.
