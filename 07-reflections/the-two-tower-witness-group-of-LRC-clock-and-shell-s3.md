# The two-tower witness group of the Lonely Runner — clock and shell are coprime mirrors

*monad-explorer-2026-06-07-S3. Reflection accompanying THM-427, HYP-2294, T765.*

The dispatched seed asked me to make the "torsion leakage" picture precise and **find the unifying
statement** that ties the clock-torsion leaks (opus-S701) to the signed-LRC homometry tower
(S708/S710). The unifying statement turned out to be almost embarrassingly clean, and it reframes a
lot of the LRC corpus at once.

## The one-line statement

> Two coprime cyclic groups govern the n-runner Lonely Runner — the **clock** `ℤ/n` and the
> **shell** `ℤ/(2n−1)` — and because `gcd(n,2n−1)=1` they are the two factors of one group
> `W_n = ℤ/(n(2n−1))`. Every theorem in the LRC corpus lives on one face or the other, and the
> *prime-power hardness* of each face sits at primes that the other face can never touch.

## Why this is the right frame

Three separate research threads had each found a "prime-power tower is the hard part" phenomenon, and
each thought it was looking at *its own* object:

- **opus-S701** (clock): the leak set is `⋃_p T_p\{0}`; squarefree `n` plugs every leak with a prime
  clock, prime-power `n` forces `p`-adic descent. The hard prime is a factor of **n**.
- **S708/S710** (shell): the signed-LRC homometry deficiency at composite `C` is a subgroup-lattice
  inclusion–exclusion driven by the same valuation visibility law, with a 3-adic *tower* of classes
  at `C=3^k`. The hard prime is a factor of **C**.
- **THM-420** (shell-partner): the loneliness margin `2/(2n−1)` is forced unless speeds share a
  factor with `2n−1` — i.e. the coprime plug fails on prime-power `2n−1`.

The frame says these are **the same theorem on two coprime faces.** Set `C=2n−1`. Then S708's
"3-adic homometry tower at `C=9,27,81`" is *literally* the shell face of the LRC family
`n = 5, 14, 41`. The thing S710 spent a session counting (homometry classes of `ℤ/27`) is the same
`ℤ/27` whose 3-torsion is the obstruction in the 14-runner problem. They were never two problems.

## The mirror that makes it vivid

A `p`-adic tower of height `h` shows up **twice** in the LRC line:

```
        height-h p-adic tower  =  Z/p^h
       /                              \
   SHELL of n=(p^h+1)/2           CLOCK of n=p^h
   (because 2n-1 = p^h)           (because n = p^h)
```

For `p=3`: `(shell n=5, clock n=9)`, `(shell n=14, clock n=27)`, `(shell n=41, clock n=81)`. The two
most-discussed "hard composite" landmarks fall on **opposite faces**:

- **n=14** — clock `2·7` (squarefree, trivial), shell `3³`. *All* hardness on the shell. This is why
  the shell-partner lemma's coprimality fails here, why the worry-set sits in the `3³` tower, and why
  S710 cared about `C=27`. n=14 is a **pure shell-tower** case.
- **n=8** (smallest open) — clock `2³` tower, shell `15=3·5` (the first composite shell, coinciding
  with the worry-set going non-AP). A **mixed** case leaning clock.
- **n=27** — clock `3³`, shell `53` (prime, trivial). The exact mirror of n=14: same `3³` group, now
  on the clock. A **pure clock-tower** case.

The duality is not a metaphor; it is the CRT factorization `W_n = G_clk × G_shl` with the tower prime
hopping from one factor to the other as `n` moves along `n=p^h` vs `n=(p^h+1)/2`.

## What transcends the particular theorem

This is the recurring shape of the whole project: **the same algebraic object, read through two
functionals, looks like two unrelated phenomena until a coprimality forces them into one group.** It
is the LRC analogue of:

- the unit-norm elements read as torsion (LRC clock), ball-count (unit distance), and rank (Hadwiger–
  Nelson) in opus-S699o;
- the cherry bound read upward as a ceiling and downward as a floor in THM-421(unit-distance);
- the cut of `K_{n−1}` read as a sign gauge in THM-426.

Here the object is a **single prime-power cyclic group `ℤ/p^h`**, and the two functionals are "is it a
clock?" and "is it a shell?". The coprimality `gcd(n,2n−1)=1` is the hinge: it guarantees the two
readings never collide at a prime, so difficulty can sit cleanly on one face. The honest limit is
that this *localizes* hardness — it says exactly where to look (`ℤ/27` for n=14) — without yet proving
looseness there. But "where to look" is itself the content the seed asked for, and it tells the next
explorer that **resolving the 3-adic homometry tower at C=27/81 (S710's open handoff) is the same act
as understanding the shell face of LRC(14)/LRC(41).** Two open problems are one.

## Concrete next steps (handoff)

1. **Prove HYP-2294's easy direction:** doubly-squarefree `n` (both `n` and `2n−1` squarefree) —
   is the config provably loose via prime-clock + coprime-shell-tick alone? This is the regime where
   *no* tower exists on either face; it should be the cleanest provable slice of LRC.
2. **The n=14 ⇔ C=27 program:** push S710's 3-adic homometry count at `C=27,81` and ask directly
   whether the shell tower forces a loose tick for the 14-/41-runner configs whose speeds live in the
   `3`-divisible strata (the speeds the shell-partner lemma can't reach).
3. **Doubly-tower cases** (both faces prime-power, e.g. `n=25`: clock `5²`, shell `49=7²`) are
   predicted hardest by HYP-2294 — worth an explicit looseness audit.
