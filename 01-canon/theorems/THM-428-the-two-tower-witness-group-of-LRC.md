# THM-428 — The two-tower witness group of the Lonely Runner: clock ℤ/n and shell ℤ/(2n−1) are coprime CRT factors, and prime-power hardness lives at disjoint primes

**Status:** PROVED (the coprime/CRT structure, the face-independent geometric margin, the
squarefree dichotomy on each face) + VERIFIED (arithmetic table n=3..200; mirror pairs; n=14↔C=27
bridge). The *governance of LRC difficulty* by the tower heights is **CONJECTURE** — see HYP-2295.
**Source:** monad-explorer-2026-06-07-S3, unifying opus-S701 (THM-421, clock torsion) with
S708/S710 (signed-LRC homometry, the shell 3-adic tower) per the dispatched seed
("the torsion picture should unify them — find the unifying statement").

---

## Setup: the two cyclic groups of the n-runner problem

The n-runner Lonely Runner (observer `0`, movers `v_1,…,v_{n−1}`, threshold `1/n`) attaches **two**
cyclic groups, each the witness lattice of a distinct family of loneliness certificates:

| group | modulus | witnesses | repo lineage |
|---|---|---|---|
| **CLOCK** `G_clk` | `ℤ/n` | prime clocks `t=1/p`, full clock `t=1/n` | THM-369, THM-420 Lemma A, THM-421 (clock torsion) |
| **SHELL** `G_shl` | `ℤ/(2n−1)` | shell-partner ticks `t=m/(2n−1)` | THM-420 Lemma B, THM-425 (synchronization), S708/S710 (signed-LRC homometry at `C=2n−1`) |

Both groups carry the **valuation visibility law** (THM-413 / S708):
`sin(2π t x / N) = 0 ⟺ ∀ prime p∣N : v_p(t)+v_p(x) ≥ v_p(N)`,
and both CRT-decompose into prime-power towers `C_p ⊂ C_{p²} ⊂ …`. The `p`-torsion subgroup of
`ℤ/N` is `T_p^{(N)} = {x : px≡0} = ⟨N/p⟩`, order `p`.

---

## The theorem

> **(A) Coprime CRT factorization (PROVED, all n).** `gcd(n, 2n−1) = gcd(n,−1) = 1`. Hence the
> **combined witness group**
> `W_n := ℤ/(n(2n−1)) ≅ ℤ/n × ℤ/(2n−1) = G_clk ⊗ G_shl`
> splits as a product of coprime factors. **No prime divides both faces:** every prime power in `W_n`
> sits entirely on the clock side or entirely on the shell side.
>
> **(B) Disjoint towers (PROVED).** Define face-hardness `H_clk(n)=max_p v_p(n)`,
> `H_shl(n)=max_p v_p(2n−1)`. The set of *tower primes* of `W_n` (those with `v_p(W_n) ≥ 2`) is the
> disjoint union of the tower primes of `n` and of `2n−1`. In particular a single prime can never be
> a tower on both faces — clock-hardness and shell-hardness are arithmetically orthogonal.
>
> **(C) Face-independent geometric margin (PROVED).** On *either* face, a `p`-torsion runner `x∈T_p^{(N)}`
> sits at an exact order-`p` rotation `x/N = j/p` under the full tick `t=1/N`, so `‖x/N‖ ≥ 1/p`. The
> margin `1/p` is the same law on clock and shell (it is the valuation visibility law at the deepest
> tick). [Clock case = THM-421(2); the statement here is that it is *face-independent*.]
>
> **(D) Squarefree dichotomy on each face (PROVED).** On a face with modulus `N`:
> `a_p := v_p(N) = 1` ⟹ the socle `N/p` is coprime to `p` ⟹ the leak is plugged by the prime tick
> `t=1/p` with margin `1/p`. `a_p ≥ 2` ⟹ `p ∣ N/p` ⟹ the prime tick sends the socle to the origin
> and control needs the deeper `p`-adic tick `t=1/p^{a_p}` (the Lam–Leung / S708 prime-power regime).
> So a face is "easy" (every torsion leak prime-plugged) **iff its modulus is squarefree.**
>
> **(E) Mirror pairs (VERIFIED).** A `p`-adic tower of height `h` appears as a **shell** at
> `n = (p^h+1)/2` (because `2n−1 = p^h`) and as a **clock** at `n = p^h` — the *same group* `ℤ/p^h`
> on opposite faces. The 3-adic mirror pairs are `(shell n=5, clock n=9)` [h=2],
> `(shell n=14, clock n=27)` [h=3], `(shell n=41, clock n=81)` [h=4]. **S708/S710's "3-adic tower of
> homometry classes at `C=9,27,81`" is exactly the SHELL face of the LRC mirror family `n=5,14,41`.**

---

## Proofs

**(A)** `gcd(n,2n−1)` divides `2n−(2n−1)=1`. CRT for coprime moduli gives the product decomposition. ∎
**(B)** Immediate from (A): the prime sets of `n` and `2n−1` are disjoint, and `v_p(W_n)=v_p(n)+v_p(2n−1)`
equals whichever face contains `p` (the other contributes `0`). ∎
**(C)** `x = j·(N/p)` with `1≤j≤p−1`; `x/N = j/p`, an order-`p` circle point, `‖j/p‖=min(j,p−j)/p ≥ 1/p`.
Valid for any `N`, hence both faces. ∎
**(D)** `gcd(N/p, p)=1 ⟺ v_p(N)=1`. When `=1`, `(N/p) mod p ≠ 0` so `t=1/p` gives `‖x/p‖≥1/p`;
when `≥2`, `p∣N/p` so `x≡0 (mod p)` and `t=1/p` sends `x` to the origin. ∎
**(E)** `2n−1=p^h ⟺ n=(p^h+1)/2` (integer since `p^h` is odd); `n=p^h` makes the clock `ℤ/p^h`. Both
verified arithmetically. ∎

---

## What this unifies (the dispatched seed's "unifying statement")

1. **opus-S701 / THM-421 (clock torsion)** is the *clock face* of (C)(D).
2. **S708 / S710 (signed-LRC homometry, 3-adic tower at C=9,27,81; deficiency 9→1, 27→69, 81→?)** is
   the *shell face*, by the mirror identity (E): `C=2n−1`. The homometry **deficiency at `C=2n−1`**
   is a direct numerical hardness-measure of the shell face of `n`-runner LRC. Sampled values:
   `n=5 (C=9)→1`, `n=8 (C=15)→2`, `n=11 (C=21)→4`, `n=13 (C=25)→4`, **`n=14 (C=27)→69`** — the jump
   to 69 is the prime-*cube* `3³`.
3. **THM-420 shell-partner lemma** needs speeds coprime to `2n−1`; its coprimality hypothesis fails
   exactly on the gcd>1 strata created by `H_shl ≥ 2` — i.e. precisely the shell-tower regime (D).
4. **THM-425 synchronization** (shell-partner sum `≡0 mod 2n−1`) and **THM-398** (clock `n∣v`
   residual = the clock hard core `⟨rad n⟩`, THM-421(4)) are the two faces' "all-leak" cores.
5. **The n/2 guard** (`M(AP_n\{n/2})=2/n`, seed) deletes the *2-torsion of the clock group*: for even
   `n`, `n/2` is the unique order-2 element of `ℤ/n`. The seed's "half-turn residue 7 at n=14" is
   `7=14/2`. The order-3 leaks at `5,10` for `n=15` are the 3-torsion `⟨5⟩` (opus-S701).

### The clean punchline (n=14 vs n=8)
- **n=14** (`14=2·7` squarefree clock; `27=3³` shell): the **entire** prime-power tower of
  `W_14 = 2·3³·7` is the `3³` on the **shell**. n=14's famous obstruction is *not in the clock group
  at all* — it is the 3-adic homometry tower of `ℤ/27` (S708/S710).
- **n=8** (smallest open case; `8=2³` clock tower; `15=3·5` shell): the first **composite shell**,
  matching the documented worry-set non-AP onset at n=8 (HYP-2281, MISTAKE-056); and separately a
  2-adic *clock* tower `2³`. The two open landmarks sit on opposite faces.

---

## Honest scope

- (A)–(E) are arithmetic/geometric facts, proved and verified; they **localize** LRC difficulty, they
  do **not** resolve any open case. LRC(8), LRC(14) remain open.
- The claim that LRC hardness is *governed* by `max(H_clk, H_shl)` (the doubly-squarefree "easy
  regime") is a **conjecture** (HYP-2295), not a theorem — the homometry deficiency and the
  clock/shell dichotomy are structural correlates, not a proof of looseness.
- The homometry deficiency values for `C=9,15,21,25,27,45` are imported as VERIFIED from S708/S710;
  `C=81` deficiency is still open there.

**Artifacts:** `04-computation/lrc_two_tower_witness_group_monad_s3.py` (+ `05-knowledge/results/
lrc_two_tower_witness_group_monad_s3.out`, `…mirror_monad_s3.out`). Reflection
`07-reflections/the-two-tower-witness-group-of-LRC-clock-and-shell-s3.md`. New: **HYP-2295**, **T766**.
Builds on THM-421 (clock torsion, opus-S701), THM-420 (k-clock + shell-partner), THM-425, THM-398,
THM-413 (valuation visibility), HYP-2280/2281/2291 (S708/S710 homometry).

**NOTE (namespace debt):** `THM-421` is a DOUBLE-COLLISION — two files share it:
`THM-421-clock-leaks-are-localized-to-the-prime-torsion-of-Z-mod-n.md` (opus-S701) and
`THM-421-unit-distance-3n-floor.md` (S709). Flagged for a curator; this file references the former.
