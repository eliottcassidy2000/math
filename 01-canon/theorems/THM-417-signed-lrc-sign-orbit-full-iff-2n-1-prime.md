# THM-417 — Signed-LRC: the AP_n sign-orbit is full (`=2^{n−2}`) **if and only if** `C=2n−1` is prime

**Status:** PROVED (both directions). The forward direction (`prime ⟹ full`) is THM-415. This
theorem supplies the **converse** (`composite ⟹ a collision exists`) by an explicit construction,
**closing HYP-2270** entirely. Verified constructively for all 63 composite `C∈[5,219]`
(`signed_lrc_converse_proof_s707d.py`), and against exhaustive brute-force orbit counts for all
`C≤39` (`signed_lrc_collision_converse_s707.py`, matching THM-413/THM-415/HYP-2273 data).
**Source:** monad-explorer-2026-06-06-S707. Builds on THM-413 (order-3 silent flip — the special
case `q=3`), THM-415 (prime ⟹ full, the homometry/Galois argument), THM-401/403 (shells, modulus
`C=2n−1`), HYP-2262 (signed-LRC theory: sign = cut, folded clock-multiset), HYP-2273 (homometry
reframe; the `H_q` half-system flips).

---

## Statement

Setup (THM-415, HYP-2262). Runners `V={1,…,n−1}`, modulus `C=2n−1`. A **cut** `ε∈{±1}^{n−1}`
(up to the global swap `ε↦−ε`, so `2^{n−2}` of them) sends runner `i` to the point
`u_i = ε_i·i ∈ ℤ/C`. Since `{0}∪{±1,…,±(n−1)} = ℤ/C`, the point set `S_ε={ε_i i}` is a
**half-system selection** (one of `{i,C−i}` per magnitude). The **folded clock-multiset** is the
multiset of circular distances `ρ(u_i−u_j)`; two cuts **collide** when these multisets agree. The
**sign-orbit** is the number of distinct folded clock-multisets among the `2^{n−2}` cuts.

> **THEOREM.** The `AP_n` sign-orbit `= 2^{n−2}` (no nontrivial collision) **⟺ `C=2n−1` is prime**.

Equivalently: a nontrivial collision exists **iff** `C` is composite.

---

## Decomposition recalled (THM-414/THM-415)

With `ζ=e^{2πi/C}` and `f̂_ε(t)=Σ_i ζ^{t ε_i i} = A(t)+i·Φ(ε)_t`, where `A(t)=Σ_i cos(2πti/C)` is
cut-independent and `Φ(ε)_t=Σ_{i=1}^{n−1} ε_i sin(2πti/C)` is the signed sine sum, one has
`|f̂_ε(t)|²=A(t)²+Φ(ε)_t²`. Hence the

> **Collision criterion.** `ε,ε'` collide ⟺ `Φ(ε)_t² = Φ(ε')_t²` for every `t=1,…,(C−1)/2`,
> i.e. `Φ(ε')_t = ±Φ(ε)_t` with an **independent sign per frequency `t`**.

---

## Forward direction (prime ⟹ full) — THM-415

Proved in THM-415 via the Galois action on the zero set `Z={t:Φ(δ)_t=0}` (`δ=ε−ε'`): for prime
`C`, `Z` is `(ℤ/C)^*`-stable, forcing `Z=∅` or `Z=ℤ/C`, both of which collapse a putative
collision to the trivial swap. Not reproved here.

---

## Converse (composite ⟹ collision exists) — the construction

Let `C` be composite and fix **any** prime `q ∣ C` with `q<C`; put `m=C/q>1`. Both `q,m` are odd
(`C` is odd). Let
```
   K  = ⟨m⟩ = { jm mod C : j=0,…,q−1 }        (the unique order-q subgroup of ℤ/C),
   H_q = { m, 2m, …, ((q−1)/2)·m }            (its half-system; |H_q|=(q−1)/2 ≥ 1).
```
Two elements lie in the same `K`-coset iff they are congruent mod `m`; cosets ↔ residues `r∈ℤ/m`,
and `r` ↔ `m−r` are negation-paired. As `m` is odd, the only self-paired coset is `K` itself
(`r=0`).

**The base cut.** Define `ε` by aligning every non-`K` runner into a *canonical* coset:
```
   ε_i = +1  if (i mod m) ∈ {0,1,…,(m−1)/2},      ε_i = −1  otherwise.            (def)
```
Set `ε' := ε` with the signs on `H_q` negated (`ε'_x=−ε_x` for `x∈H_q`, else `ε'_x=ε_x`).

**Claim 1 — the signed non-`K` points form full cosets.** For a runner `i` with `i≡r (mod m)`,
`ε_i i ≡ r (mod m)` if `r≤(m−1)/2` and `ε_i i ≡ −i ≡ m−r (mod m)` if `r>(m−1)/2`. So for each
canonical residue `r∈{1,…,(m−1)/2}`, the runners in the two paired cosets `{r,m−r}` contribute,
after signing, exactly the `q` points of the full coset `r+K = {r, r+m,…, r+(q−1)m}` (one point per
magnitude of the half-system, all routed into `r+K`). Hence
```
   { ε_i·i : i∉K } = ⊍_{r=1}^{(m−1)/2} (r+K)      (each a complete coset of K).
```

**Claim 2 — `A_t:=Σ_{i∉K} ε_i sin(2πti/C)=0` whenever `q∤t`.** By Claim 1,
`A_t = Σ_{r=1}^{(m−1)/2} Σ_{x∈ r+K} sin(2πtx/C)`. For a fixed coset,
```
   Σ_{x∈ r+K} sin(2πtx/C) = Im[ ζ^{tr} Σ_{j=0}^{q−1} ζ^{t j m} ]
                          = Im[ ζ^{tr} · Σ_{j=0}^{q−1} e^{2πi t j/q} ] .
```
The inner geometric sum is `q` if `q∣t` and `0` if `q∤t`. Thus `A_t=0` for every `q∤t`. ∎(Claim 2)

**Claim 3 — `B_t:=Σ_{x∈H_q} ε_x sin(2πtx/C)=0` whenever `q∣t`.** Each `x=jm∈H_q` gives
`sin(2πt·jm/C)=sin(2πtj/q)=0` when `q∣t`. ∎(Claim 3)

**Conclusion.** Split `Φ(ε)_t = A_t + B_t` and `Φ(ε')_t = A_t − B_t` (negating `H_q` flips the sign
of the `B`-part only). Then for every `t`:
* if `q∣t`: `B_t=0` (Claim 3), so `Φ(ε')_t = A_t = Φ(ε)_t`;
* if `q∤t`: `A_t=0` (Claim 2), so `Φ(ε)_t = B_t`, `Φ(ε')_t = −B_t`, and `Φ(ε')_t² = Φ(ε)_t²`.

By the collision criterion, `ε` and `ε'` collide. They are distinct modulo the global swap:
`ε'=ε` is impossible (`H_q≠∅`), and `ε'=−ε` would require `H_q` to be **all** `n−1` runners, i.e.
`(q−1)/2=(C−1)/2`, i.e. `q=C`, contradicting `q<C`. Hence the sign-orbit is `<2^{n−2}`. ∎

**Per-frequency reflection (observed, S707c).** For the constructed pair the per-frequency sign
`s_t = Φ(ε')_t/Φ(ε)_t` equals exactly `+1` if `q∣t` and `−1` if `q∤t` — i.e. flipping the
half-system `H_q` *realizes the subgroup reflection* `σ_q` on the frequency line. This is the
runner-space mechanism behind HYP-2273's "frequency-localized" half-system flips.

---

## The special cases recover prior theorems

* `q=3`: `H_3={C/3}` is a single runner (the order-3 torsion point), and the move is the
  single silent flip of `x=C/3`. This is exactly **THM-413** (`3∣C ⟹` degenerate via the order-3
  point), now seen as the `q=3` instance of one uniform construction.
* The construction uses **any** proper prime factor, so it certifies a collision for every
  composite `C` (squarefree or prime-power) — e.g. `q=5` for `C=25`, `q=7` for `C=49`,
  `q=11` for `C=121`.

---

## Corollary (complete classification, closes HYP-2270)

> For the arithmetic progression `AP_n={1,…,n−1}`:
> **sign-orbit `= 2^{n−2}` ⟺ `2n−1` prime; sign-orbit `< 2^{n−2}` ⟺ `2n−1` composite.**

The signed-LRC pair-clock invariant is therefore a *primality detector* of the modulus `C=2n−1`:
its faithfulness is equivalent to `C` being prime.

---

## The exact silent set — full-coset characterization (PROVED)

The converse construction is one point of a structure that is exactly describable. For a fixed
prime `q∣C` (`m=C/q`, `K=K_q` the order-`q` subgroup), let `Σ_ε = {ε_i·i : i∉K}` be the signed
non-subgroup point set of a cut `ε`, and `Σ̂_ε(t)=Σ_{x∈Σ_ε} ζ^{tx}`.

> **LEMMA 1 (the real part is forced).** For every half-system selection `Σ` of the non-`K` runners
> and every `t` with `q∤t` (`t≢0`), `Re(Σ̂(t)) = 0`.

*Proof.* `Σ ⊍ (−Σ)` is the full symmetric set `notK = ℤ/C∖K∖{0}` of non-subgroup elements, so
`2·Re(Σ̂(t)) = \widehat{1_{notK}}(t)`. As `1_{notK}=1_{ℤ/C}−1_K` and `\widehat{1_K}(t)=q·[q∣t]`
(geometric sum over `K=⟨m⟩`), `\widehat{1_{notK}}(t) = C·[t≡0] − q·[q∣t] = 0` for `t≢0, q∤t`. ∎

Since the silent condition is exactly `A_t = Im(Σ̂(t)) = 0` for `q∤t`, Lemma 1 upgrades it to
`Σ̂(t)=0` for all `q∤t`, i.e. `1_{Σ}` is Fourier-supported on `q∣t`, i.e. **`1_Σ` is constant on
`K`-cosets**. Hence:

> **CHARACTERIZATION (PROVED).** Flipping `H_q` is silent for `ε` **⟺** `Σ_ε` is a **union of full
> `K_q`-cosets**. (Verified independently: `#{silent} = #{full-coset}` exactly for all composite
> `C≤39`, `signed_lrc_count_law_s707e.py`; Lemma 1's `Re≡0` confirmed to `10⁻¹³` for all `C≤77`,
> `signed_lrc_re_zero_s707f.py`.)

> **LEMMA 2 (count).** The number of half-system selections `Σ` that are unions of full `K_q`-cosets
> is `2^{(m−1)/2}` — one independent binary choice (`r+K` vs `(m−r)+K`) for each of the `(m−1)/2`
> negation-paired nonzero cosets. (Matches every enumerated value.)

## Count law for the deficiency

Combining: the cuts for which `flip H_q` is silent number `2^{(m−1)/2}` (choice of `Σ`) `× 2^{(q−1)/2}`
(free signs on `H_q`); modulo the pairing `ε↔ε⊕H_q` and the global swap, the number of `H_q`-flip
collision **pairs** is `2^{(m−1)/2 + (q−1)/2 − 2}`, `m=C/q`. Summing over the prime factors (valid
when every collision class has size 2 and comes from a single subgroup — the squarefree regime,
HYP-2273 verified `C≤39`):

> **COUNT LAW (squarefree `C=pq`, PROVED modulo "every collision is a single `H_q` flip").**
> deficiency `= 2^{(p−1)/2+(q−1)/2−2} + 2^{(q−1)/2+(p−1)/2−2} = 2^{(p+q)/2 − 2}.`

Verified: `15→4, 21→8, 35→16, 33→32, 39→64`, and predicts `55→64, 65→128, 51→256, 57→512, …`
(`signed_lrc_re_zero_s707f.py`). This *generalizes* the established `C=3p` law `2^{(p−1)/2}`
(HYP-2273; set `q=3`).

**Prime powers need the full subgroup lattice.** `C=27=3³` has deficiency `69` (a size-4 class plus
size-2 classes) because the nested subgroups `K_3⊂K_9` produce *combined* half-system flips
(`H_3`, `H_9`, `H_3⊕H_9`); the single-prime sum gives only `8`. The uniform count over all
factorizations remains **open** (it is a question about the lattice of `Σ̂`-vanishing sets, i.e.
which subgroups' cosets `Σ` simultaneously refines) — see HYP-2270/HYP-2273.

---

## Artifacts

`04-computation/signed_lrc_collision_converse_s707.py` (flip-sets = subgroup half-systems),
`signed_lrc_explicit_pairs_s707c.py` (per-frequency reflection `s_t=+1⟺q∣t`),
`signed_lrc_converse_proof_s707d.py` (constructive verification, all composite `C≤219`),
and outputs in `05-knowledge/results/`. Reflection:
`07-reflections/signed-lrc-sign-orbit-is-a-primality-detector-s707.md`.
