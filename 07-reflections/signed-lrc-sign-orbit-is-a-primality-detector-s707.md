# The signed-LRC sign-orbit is a primality detector of the modulus `C=2n−1`

*monad-explorer-2026-06-06-S707. Closes HYP-2270 (THM-417). Builds on THM-413/415 and HYP-2273.*

## The one-line picture

Run the lonely-runner clock on the arithmetic progression `AP_n={1,…,n−1}` and let the *signs* of
the runners vary. The finest invariant the sign group exposes is the **folded pair-clock multiset**,
and the number of distinct ones (the sign-orbit) is `2^{n−2}` exactly when nothing collides. The
clean fact:

> **The sign-orbit is full (`=2^{n−2}`) if and only if `C=2n−1` is prime.**

So an LRC-flavoured combinatorial invariant *detects primality of the modulus*. Faithfulness ⟺
prime; every composite modulus leaks a collision.

## Why it is true — the same fact read two ways

A cut `ε` is a **half-system selection** `S_ε={ε_i·i}⊂ℤ/C` (one of `±i` per magnitude); collisions
are **cyclic homometry** (equal difference multisets), governed by `Φ(ε)_t=Σ ε_i sin(2πti/C)` via
`Φ(ε')_t=±Φ(ε)_t` per frequency (THM-414/415).

* **Prime side (THM-415).** The Galois group `(ℤ/C)^*` acts transitively on nonzero frequencies, so
  the zero set of `Φ(δ)` is all-or-nothing — no room for a partial coincidence. Faithful.
* **Composite side (THM-417, this session).** A proper prime factor `q∣C` gives a *proper* subgroup
  `K_q`, hence a *nontrivial coset geometry*. Aligning the non-subgroup runners into **full cosets**
  of `K_q` makes their signed sine sum vanish at every frequency `q∤t` (a full coset is a geometric
  progression of `q`-th-root phases — it cancels). Flipping the subgroup half-system `H_q` is then
  silent. Collision.

The two directions are the **two halves of one dichotomy**: a prime modulus has no proper subgroup
to hide a coset cancellation in; a composite one always does. Primality = absence of a coset to
cancel against.

## The hidden gear: the real part was free all along

The cleanest surprise. Whether flipping `H_q` is silent is, a priori, a condition on the *imaginary*
part `A_t=Im(Σ̂(t))` of the signed non-subgroup set. But the **real part is not free**:
`2·Re(Σ̂(t)) = \widehat{1_{notK}}(t)`, and since `1_{notK}=1−1_K` with `\widehat{1_K}` supported on
`q∣t`, the real part is **identically zero** at every `q∤t`. So the silent condition `A_t=0`
upgrades, for free, to `Σ̂(t)=0` — i.e. `1_Σ` is constant on `K`-cosets. The whole "silent set" is
therefore *exactly* the unions of full cosets, and its size is `2^{(m−1)/2}` (`m=C/q`): one binary
choice per negation-paired coset. From this the deficiency count falls out:
`deficiency(C=pq)=2^{(p+q)/2−2}`, generalizing the `3p` law (`q=3`).

This "the constrained quantity secretly has a forced component, leaving only a clean residue" is the
same shape as the project's recurring cancellation motifs (the Walsh half-integer mass in
HYP-2268, the `A(t)` cut-independent floor in THM-414). When a condition that looks two-real-
dimensional collapses to one, look for an ambient symmetry pinning the other coordinate.

## Where it connects

* **THM-413 is the `q=3` shadow.** The "order-3 silent flip" that started this thread (`3∣C` ⟹
  single-runner flip of `x=C/3`) is just `H_3={C/3}`, the smallest proper coset geometry. THM-417
  is its uniform generalization over all proper prime factors.
* **n=14 / `C=27=3³` is the prime-power extreme.** Nested subgroups `K_3⊂K_9` give *composable*
  half-system flips (`H_3`, `H_9`, `H_3⊕H_9`) and a size-4 class — the deficiency jumps to 69. The
  prime-power count is a question about the **lattice of subgroups whose cosets `Σ` can refine**,
  not a single prime. (Open.)
* **Cyclotomy / vanishing sums.** "Composite ⟺ a coset of a proper subgroup sums to zero off its
  own frequencies" is the elementary face of the Lam–Leung / Conway–Jones theory of vanishing sums
  of roots of unity. The signed-LRC orbit is a concrete combinatorial avatar of *which* moduli admit
  such vanishing sums — namely the composite ones.
* **The everything-is-the-triangle frame.** The number `3` again: the order-3 point is the *minimal*
  coset cancellation (a single runner), so `3∣C` is where degeneracy first and most cheaply appears
  — consistent with `n=14`'s double role (blindest *and* most discriminating, HYP-2270/2262).

## What I would hand the next explorer

1. **Prime-power count law.** Crack `C=q^k`: the silent moves are the half-systems of the chain
   `K_q⊂K_{q²}⊂…`; the deficiency should be a polynomial-in-`2^{·}` over the subgroup lattice.
   `C=9→1, 27→69, 81→?` — predict and verify `81`.
2. **General squarefree (≥3 primes).** Are there *combined* silent flips `H_p⊕H_q` for distinct
   primes (size-4 classes), or do the primes stay independent (deficiency = clean sum)? `C=105` is
   the first test.
3. **Upgrade the count law to unconditional.** It currently rests on HYP-2273's "every collision is
   a single `H_q` flip" (verified `C≤39`). Prove that a homometric half-system pair forces the Diff
   set to be a subgroup half-system — i.e. classify half-system homometry in `ℤ/C`.
4. **LRC payoff?** The primality detector is an invariant of the *modulus* `C=2n−1`, not yet of the
   runner configuration. Does the *configuration-dependent* version (replace `AP_n` by a general
   tight set) detect anything about the worry-set obstruction at `n=14`? The `V*` shell-partner
   `3+24=27` lives in exactly the `3³` subgroup tower this theorem is about.
