# The AP-base case of the loose-branch rigidity: an elementary proof (ladder + CRT free-rider)

*mac-mini-2026-07-05-S59. Owner: work the EXACT mathematics of the LRC(14) crux before formalizing.
Building directly on klein-S140 (HYP-4151): the loose branch is the 12-runner AP-rigidity, the gap
`(1/13, 2/25)` is a Farey window, `r=1` is proved, and the open half is the "equal-spacing rigidity uniform
in `r`." This note proves that open rigidity for the **AP-base** sub-case — every 12-family that contains an
LRC-tight 11-subfamily — with a clean, elementary argument. It is the near-AP regime, exactly where the
danger of a gap-violator is greatest.*

## The statement

> **Theorem (AP-base rigidity).** Let `V = c·{1,…,11} ∪ {X}` with `c ≥ 1`, `X` a positive integer not in
> `c·{1,…,11}`, and `gcd(V) = 1` (equivalently `gcd(c, X) = 1`). Then
> **`M(V) ∈ {1/13} ∪ [2/25, ∞)`** — `V` is never in the gap `(1/13, 2/25)`.
> Moreover `M(V) = 1/13` only for `V = {1,…,12}` (`c=1, X=12`).

Verified exactly over all primitive `c·{1,…,11}∪{X}` with `c ≤ 12` (`lrc14_apbase_rigidity_macmini_S59.py`):
zero in the gap. This handles klein's rigidity whenever *some* 11-subfamily of the 12-set is LRC(12)-tight
(a dilated AP) — the case that is closest to the extremal AP.

## The proof

The 11-runner AP `c·{1,…,11}` is LRC(12)-tight: `M(c·{1,…,11}) = 1/12`, attained at every
`t_j = j/(12c)` with `gcd(j,12)=1`. *(At `t_j`, runner `ck` sits at `kj/12`; `gcd(j,12)=1` makes
`{kj mod 12 : k=1..11}` a permutation of `{1,…,11}`, so no runner hits `0` and the minimum distance is
exactly `1/12`.)* There are such witnesses for `j ≡ 1,5,7,11 (mod 12)`.

Add the twelfth runner `X`. At a witness `t_j`,
`min_{v∈V} ‖v t_j‖ = min( 1/12,\ ‖X·j/(12c)‖ )`.
`X` is **safe** at `t_j` (contributes `≥ 1/12`) iff `‖X j/(12c)‖ ≥ 1/12`, i.e. iff
`X j \bmod 12c \notin (-c, c)` (the width-`2c` "danger band" around `0`). If `X` is safe at *even one*
witness, then `M(V) ≥ 1/12 > 2/25`, done. So assume `X` is unsafe at every witness and derive the ladder.

**Case `c = 1`.** Witnesses are `t_j = j/12`, `j ∈ {1,5,7,11}`. `X` unsafe means `Xj \bmod 12 \in \{0\}`
(band `(-1,1)`), i.e. `12 \mid Xj`; taking `j=1`, this forces `12 \mid X`. So the only unsafe case is
`X = 12k`, giving the **ladder** `V = \{1,…,11,12k\}`, whose value is `M = k/(12k+1)` (klein-S140): `1/13`
at `k=1` (the AP), and `k/(12k+1) ≥ 2/25 ⇔ k ≥ 2`. Every other `X` is safe ⇒ `M ≥ 1/12`. Either way
`M(V) ∈ {1/13} ∪ [2/25, ∞)`.

**Case `c ≥ 2` (the CRT free-rider).** Consider the `c` witnesses with `j ≡ 1 (mod 12)`:
`j = 1, 13, 25, …, 1+12(c-1)`. Their residues are, with `e := X \bmod 12c`,
```
   X·(1+12t) ≡ e + 12·(tX \bmod c)   (mod 12c),   t = 0,…,c-1.
```
Since `gcd(X,c)=1` (primitivity), `t ↦ tX \bmod c` is a bijection of `{0,…,c-1}`, so this set equals
`\{ e + 12s \bmod 12c : s = 0,…,c-1 \}` — an arithmetic progression of step `12`, spanning `[e, e+12(c-1)]`.
These `c` values are spaced `12` apart; the danger band `(-c, c) \bmod 12c` is a single arc of width `2c`, so
it can contain **at most `⌊c/6⌋ + 1`** of them. For every `c ≥ 2`, `⌊c/6⌋ + 1 < c` — so **at least one of the
`c` witnesses is safe.** (Primitivity is what forbids the degenerate `e = 0`, i.e. `12c \mid X`, which would
put a runner at `0` and is excluded by `gcd(c,X)=1`.) Hence `M(V) ≥ 1/12 > 2/25`. ∎

So for `c ≥ 2` there is **no ladder at all**: `M(V) = 1/12` for *every* admissible `X` (confirmed —
`2·{1,…,11}∪{X}` gives exactly `1/12` for all odd `X`). The dilation `c ≥ 2` "spreads the killer's phase"
across the AP's many witnesses, and primitivity guarantees a hole it cannot plug — the same CRT free-rider
mechanism as the covering-min work (mac-mini S46/S47), now doing load-bearing work in the *gap* proof.

## Where this sits, and what remains

- **What it closes.** klein's open rigidity for every 12-set with an LRC-tight (dilated-AP) 11-subfamily —
  the near-extremal regime. It also makes explicit the handler the fleet has been routing tight-subfamily
  cases to ("the S47 CRT case"): here it is, as a proof.
- **The clean reduction.** A gap-violator `V` (if any) would, on removing *any* runner, give an 11-subfamily
  `B'` with `M(B') ≥ 1/12`. If `M(B') = 1/12` then `B'` is a dilated 11-AP (LRC(12) rigidity) and the theorem
  above rules `V` out. Therefore **a gap-violator must have *every* 11-subfamily strictly loose**
  (`M(B') > 1/12`). That is the entire remaining open case.
- **The recursion it suggests.** "Every 11-subfamily strictly loose" invites the same argument one level
  down: an 11-runner family is either a dilated AP (`M=1/12`) or loose (`M ≥ 2/23`, the 11-runner gap,
  klein-S126). The general loose branch is thus a **descent through the runner count**, with the AP-base
  case (this note) as the rung and the "all-subfamilies-loose" case as the step. Whether the descent
  terminates cleanly — or needs a genuinely global (non-inductive) input at the "all loose" case — is the
  sharp open question. (Note the base rung `r=1` is klein's LRC(13)-sandwich; this note is the `c`-dilation
  companion at the family level.)

## Links

- Builds on: klein-S140 HYP-4151 (loose branch = AP-rigidity, Farey window, r=1), mac-mini S46/S47 (the CRT
  free-rider / offset-forcer for the covering-min — the same mechanism), the ladder `k/(12k+1)` (klein-S140,
  mac-mini-S38 Ostrowski), klein-S126 (11-runner gap `2/23`, the one-lower analogue for the recursion).
- Open: the "every 11-subfamily strictly loose" case of the 12-runner gap-emptiness (the loose-branch
  residual after this note). HYP-4152. Script: `lrc14_apbase_rigidity_macmini_S59.py`.
