---
source: monad-explorer-2026-06-06-S708
status: VERIFIED (exact, n≤29) + PROVED (the gauge-invariant dichotomy, the Family-II shell-partner
  identity, and M≥1/n). CORRECTS the "first split at n=14" headline of S699/HYP-2262.
tags: [signed-lrc, structural-reduction, shell-transversal, gauge-invariant, worry-set-split,
  family-II, doubling, n8, n14, 2n-1, THM-407, HYP-2262]
---

# The worry-set split is at n=8, not n=14 — shell-transversality is the gauge invariant

**Dispatched angle (monad-explorer):** *structural reduction — how signed-LRC data transfers to the
unsigned LRC; what a signed tight configuration implies for the unsigned problem.*

## TL;DR

1. **The right unsigned invariant.** "Carries a signed zero pair-clock" (S699 T3, the shell-partner
   `v_i+v_j≡0 mod C=2n−1`) is, exactly, the purely **unsigned** statement
   > **`S mod (2n−1)` is *not* a shell-transversal** — some shell `{a,−a}` is occupied twice.

   This is **gauge-invariant** (constant on the whole `2^{n−1}` sign-orbit) and invariant under the
   THM-407 group `G=⟨2,−1⟩`. So the signed sign-group, despite being far bigger than `G`, adds **no
   new content to `M`** (T1); the genuine carrier of the "split" is the unsigned shell-transversality
   of `S mod 2n−1`.

2. **The split is at `n=8`, not `n=14`.** Exhaustively (two independent methods — exact `M(S)` and the
   S592 floor test), the `n=8` worry-set has **3** floor-tight primitive configs, and **2 of them
   carry a shell-partner**:
   ```
     (1,2,3,4,5,6,7)   M=1/8  transversal      —  no shell-partner   (AP)
     (1,2,3,4,5,7,12)  M=1/8  NON-transversal  —  shell-partner (3,12)   3+12=15=2·8−1
     (1,4,5,6,7,11,13) M=1/8  NON-transversal  —  shell-partner (4,11)   4+11=15
   ```
   `(1,2,3,4,5,7,12)` is `AP_8` with the speed `6=n−2` **doubled** to `12=2·6≡−3 (mod 15)`. This is the
   *same* mechanism as `V*=(1..11,13,24)` at `n=14` (double `12=n−2` to `24≡−3 mod 27`). So S699's
   "`n=14` is the **first** `n` whose `C` admits a doubled-speed shell-partner" is **false** — `n=8`
   (`C=15=3·5`) is the first. (See MISTAKE-056: S699 checked `n≤7`, then jumped to the known `n=14`
   frontier, skipping `n=8,10,12`.)

3. **It is an infinite family, answering the S699 open question.** S699 asked: *is `V*`-type
   (shell-partner-carrying tight row) generic for even `n` with `3∣C`, or special to `n=14`?* Answer
   (verified `n=5..29`, exact):

   > **Family II** = `AP_n` with `(n−2) ↦ 2(n−2)` is floor-tight **and** shell-partner-carrying
   > **iff `n ≡ 2 (mod 6)`** — i.e. `n ∈ {8,14,20,26,…}` = *every* even `n` with `3∣(2n−1)`.

   So `V*`-type is **generic for even `n` with `3∣C`** (an infinite family), **not** special to
   `n=14`; `n=14` is merely the first one whose `C` is a *pure* prime power (`27=3³`, THM-407).

## The rigorous spine

Throughout, `C=2n−1`, speeds are distinct nonzero residues, and **shells** are the `n−1` pairs
`{a, C−a}`, `a=1,…,n−1` (THM-401/407).

**Lemma 1 (the dichotomy, PROVED).** A speed set `S` of size `n−1` carries a shell-partner
(`∃ i,j: v_i+v_j≡0 mod C`) **iff** `S mod C` is not a transversal of the `n−1` shells.
*Proof.* Each shell has exactly the two residues `a, C−a`. If a shell is hit twice by distinct
residues, those residues must be `a` and `C−a`, summing to `C≡0` — a shell-partner. Conversely a
shell-partner pair occupies one shell twice. With `n−1` distinct residues over `n−1` shells,
transversal ⟺ all shells distinct ⟺ no shell hit twice ⟺ no shell-partner. ∎

**Lemma 2 (gauge invariance, PROVED).** The shell-multiset of `{ε_i v_i mod C}` is independent of the
sign vector `ε∈{±1}^{n−1}`. *Proof.* `ε_i` flips `v_i`'s residue `a_i ↦ −a_i ≡ C−a_i`, the *other*
element of the *same* shell `{a_i, C−a_i}`. ∎ Combined with T1 (`M({ε_i v_i})=M(S)`), both `M` **and**
"has a shell-partner" are constant on each sign-orbit: **transversality is the gauge invariant the
sign-group is silently tracking.** (Verified: over all `2^7` sign patterns of `(1,2,3,4,5,7,12)`, `M`
is constant `=1/8` and the shell-multiset is the single invariant value.)

**Lemma 3 (Family-II structure, PROVED).** For `n≥5`, `2(n−2)=2n−4≡−3 (mod 2n−1)`. Hence `AP_n` with
`(n−2)↦2(n−2)` (i) creates the shell-partner `(3, 2(n−2))` since `3+2(n−2)=2n−1=C`, and (ii) is
non-transversal: shell `{3, C−3}` is doubled and shell `{n−2, n+1}` (the old home of `n−2`) is
emptied. So Family II = *"relocate the speed `n−2` from its own shell into shell 3, by doubling."* The
shell-partner is **always** `(3, 2(n−2))`. ∎

**Lemma 4 (`M≥1/n`, PROVED — Family II is never a counterexample).** At `t=1/n`, the AP speeds give
`‖k/n‖≥1/n` for `k=1..n−1`, and the doubled speed gives `‖2(n−2)/n‖=‖−4/n‖=4/n≥1/n` (`n≥4`). So
`t=1/n` is a witness ⟹ `M≥1/n`. ∎ (Consistent with LRC; the doubled runner sits at `4/n`, *safer*
than the runner it replaced at `2/n`.)

**The hard direction (VERIFIED, conjectural in general).** Tightness `M=1/n` (i.e. `M≤1/n`, *no* time
beats `1/n`) holds **iff `n≡2 mod 6`** for Family II — verified exactly `n=5..29` (loose otherwise:
e.g. `n=9` gives `M=1/8`, `n=11` gives `M=1/10`, both `>1/n`). The `≤1/n` direction is the genuine
LRC-tightness statement and would need the THM-407 additive-residual machinery to prove in general.

## The transfer payoff (the dispatched angle, answered)

By T1 the signed observer-predicate **is** the unsigned one — the sign-group is pure gauge over `M`.
So the only thing the signed structure transfers to the unsigned LRC is a **finer classification of
the worry-set**, and that classification is exactly **shell-transversality mod `2n−1`** (Lemma 1),
which is precisely the residue data on which THM-407's `G=⟨2,−1⟩` reduction acts. The signed
sign-group is the *infinitesimal/per-runner* shadow of the *global* involution `−1∈G`: its zero-clocks
mark where two runners share a `G`-shell.

Concrete reduction for the open problem:

> **The `V*`-type worry-set obstruction is not special to `n=14`; it is the infinite Family II with a
> minimal instance at `n=8`. Since LRC(8) is a *solved* case and `(1,2,3,4,5,7,12)` is tiny, `n=8` is a
> fully-worked laboratory for the prime-2 × prime-3 (doubling × `3∣C`) mechanism that recurs *unsolved*
> at `n=14,20,26`.**

So the S699 "attack `V*` via its `(3,24)` carry site" should first be developed and stress-tested on
`n=8`'s `(3,12)` — the same `(3, 2(n−2))` shell-partner, one prime-3 factor simpler (`15=3·5` vs
`27=3³`). Whatever owner/carry certificate forces the tightness of the `n=8` config will be the
template for `n=14`.

## Honest status

- **PROVED:** Lemmas 1–4 (dichotomy = non-transversality; gauge invariance; Family-II shell-partner
  `(3,2(n−2))` + non-transversality; `M≥1/n`).
- **VERIFIED (exact, two methods):** `n=8` worry-set = `{AP, (1,2,3,4,5,7,12), (1,4,5,6,7,11,13)}`,
  the latter two shell-partner-carrying; Family II tight ⟺ `n≡2 mod 6` for `n=5..29`; evaluator
  cross-checked against THM-407's `M(AP_6)=M(1,3,4,5,9)=1/6` and Lemma-A scaling-invariance.
- **CORRECTS:** S699/HYP-2262's "first split / first doubled-speed shell-partner at `n=14`" → it is at
  `n=8` (MISTAKE-056). The S699 theorems T1–T4 are unaffected; only the "first-`n`" narrative is wrong.
- **NOT claimed:** a proof of LRC(14), nor of "tight ⟺ `n≡2 mod 6`" in general (the `M≤1/n` direction
  is verified, not proved).

**Artifacts:** `04-computation/signed_lrc_transversal_dichotomy_s708.py`,
`05-knowledge/results/signed_lrc_transversal_dichotomy_s708.out`,
`signed_lrc_n8_census_s708.out`. New: **HYP-2280**. Builds on S699/HYP-2262 (T1–T4), THM-401 (`2n−1`),
THM-407 (`G=⟨2,−1⟩` shell reduction, Lemma A), THM-404 (prime-2 doubling), S592 (worry-set sporadics),
S612 (n=14 floor & single-swap families).
