# The converse — digit 0 implies a near-join count: the dichotomy, the composite right source, and what the target may choose

Owner: *"Try to prove the converse: digit 0 implies a near-join count"*; *"let the target choose the core or admit a
composite right source or something even more creative to prove chan-pak."*

THM-5263 (I3): the near-join `T → M_k → z` (`k ≥ 1`, `z → t_0 ∈ T`, everything else `→ z`) has `H = h(T) + 2^k B_{t_0}(T)`
and `z` is a digit-0 vertex. So every near-join value is the count of a tournament with a digit-0 vertex (strong when
`t_0` is in the initial component of `T`). The converse: *if a strong `S` has a digit-0 vertex `z`, is `H(S)` a near-join
value?* The near-join family of HYP-9453 has `k ≥ 1` throughout (spectrum_analogies_nearjoin_20261102 §1); with `k = 0`
allowed the converse would be the tautology (I3).

## 1. The dichotomy (PROVED)

**Theorem 6.** Let `S` be strong with digit-0 vertex `z` and port `t_0`; let `C_1 → … → C_m` (`m ≥ 2`) be the strong
components of `S − z`, `D = C_m`, `A = C_1 ∪ … ∪ C_{m−1}` (nonempty, `t_0 ∈ C_1 ⊆ A` by Corollary 4 of OUTDEG.md). Then

    H(S) = h(A)·h(D) + B_{t_0}(A)·P_D(2),      P_D(2) = Σ_{W ⊆ V(D)} h(D[W]) · h(D[V(D)∖W]).

- **(N) `D` is a single vertex** (the base `S − z` has a sink). Let `k ≥ 1` be the number of trailing single-vertex
  components and `T = C_1 ∪ … ∪ C_{m−k}` (if every component is a single vertex, `T = C_1 = {t_0}` and `k = m − 1`). Then
  `S = T → M_k → z` with `z → t_0`: **literally a near-join**, `H(S) = h(T) + 2^k B_{t_0}(T)`.
- **(C) `D` is strong of order ≥ 3** (the base has no sink): `H(S) = h(A)h(D) + B_{t_0}(A)P_D(2)` — the dead-end append
  over a **composite right source** `D`.

*Proof.* `H(S) = h(S − z) + B_{t_0}(S − z)` (THM-5263 (I2), Theorem 3). `S − z = A → D` is an ordinal sum with `t_0 ∈ A`;
THM-5260 (B) gives `h(A → D) = h(A)h(D)` and `B_{t_0}(A → D) = B_{t_0}(A)·P_D(2)`. In (N), `D = M_1` has `P = 2`; the
trailing singletons form `M_k`, `h(T → M_k) = h(T)` and `B_{t_0}(T → M_k) = 2^k B_{t_0}(T)` (THM-5260 (B), transitive
suffix). ∎

Replay ([receipts/converse_descent.out](receipts/converse_descent.out) (3)): over all 341 strong digit-0 pairs `(S, z)` of
orders 4–8 the identity has 0 failures; 279 pairs are in case (N), 62 in case (C) (P-3 confirmed). By Theorem 5 each
fertile strong tournament of order 4–8 carries exactly one pair.

So **the converse is PROVED in case (N)** — the near-join sits inside `S` itself, with `k ≥ 1`. The content is case (C).

## 2. Case (C): finite-exact, and the exact equivalence

[receipts/converse_census.out](receipts/converse_census.out) (A): the 62 case-(C) pairs carry 41 distinct tuples
`(H, order, h(A), B_{t_0}(A), h(D), P_D(2))`, with counts from 15 to 477, and **every one of them is a near-join value**
`h(T) + 2^k B_t(T)`, `k ≥ 1`, from a base `T` of order ≤ 7 (least base order 3 for five of the 62, 4 for three, 5 for 17,
6 for 35, 7 for two; none needs order 8). P-6: first clause confirmed, second clause (that one would need a base of
order 8) refuted.

Beyond order 8. By Theorem 3, `H(S)` for strong digit-0 `S` is a dead-end-append count; by (I3) every near-join value is
one. So the converse at all orders says exactly

    {counts of strong tournaments with a digit-0 vertex} = {near-join values}.

HYP-9453's census (bases ≤ 8) makes every odd number ≤ 22,227 other than 1, 7, 11, 13, 21, 47, 67 a near-join value;
THM-5263 (B) shows the seven have no strong digit-0 realiser (7 and 21 are no counts at all, 1 has no strong realiser);
and `H` is always odd (Rédei). Hence **for every strong tournament `S`, of any order, with `H(S) ≤ 22,227`: digit 0 ⟹
near-join count.** What stays open is the structural statement for counts above 22,227 in case (C): is
`h(A)h(D) + B_{t_0}(A)P_D(2)` always `h(T) + 2^k B_t(T)` with `k ≥ 1`? The census bases `T` for the 41 tuples bear no
visible relation to `(A, D)` (least base orders 3–7 for tournaments `S` of order 8), and no rewriting of a strong right
source as a transitive suffix was found. OPEN; recorded under HYP-9453.

## 3. The composite right source is a linear reader

`P_D(2) = B_{t_0}({t_0} → D)` is the number of ordered pairs of paths, either possibly empty, that cover `D`. So the right
source reads the rooted core `(h, B) = (h(A), B_{t_0}(A))` through the linear form `v = h(D)·h + P_D(2)·B`; the near-join is
the transitive reader `(1, 2^k)`. Strong readers of orders 3–5: `(3, 12)`, `(5, 32)`, `(9, 84)` (two non-isomorphic
order-5 sources), `(11, 100)`, `(13, 104)`, `(15, 108)`, `(15, 120)`; the sources of orders 6–7 occurring in the 62 pairs add
`(15, 216)`, `(17, 224)`, `(19, 256)`, `(23, 272)`, `(23, 296)`, `(25, 304)`, `(27, 288)`, `(29, 320)`, `(31, 328)`,
`(33, 336)`, `(33, 360)`, `(37, 352)`, `(37, 376)`, `(41, 416)`, `(43, 424)`, `(45, 432)`.

`P_D(2)/h(D)` is a power of two for only 3 of the 396 strong `D` of orders 3–7 — `(3, 12)`, `(13, 104)`, `(15, 120)` — and
there the composite value is `h(D)·(h + 2^j B)`, a count times a near-join value (P-2 confirmed). All the others are new
linear readers; those with `gcd(h(D), P_D(2)) = 1` — `(5, 32)`, `(11, 100)`, `(17, 224)`, `(19, 256)`, `(23, 272)`,
`(23, 296)`, `(25, 304)`, `(29, 320)`, `(31, 328)`, `(37, 352)`, `(37, 376)`, `(41, 416)`, `(43, 424)` — run through every
residue class modulo `P_D(2)` as `h` varies. That is the divisor-breaking the previous session's READERS.md could not get
from homogeneous readers, at the price the owner names: an inhomogeneous vertex (the digit vertex `z`).

The smallest core, `A = {t_0}` (`h = B = 1`), gives `H(t_0 → D → z → t_0) = h(D) + P_D(2)` — a strong `D` with a two-arc
handle: 15, 37, 93, 111, 117, 123, 135 at orders 5–7 and 231, 241, 275, 295, 315, 319, 329, 349, 359, 369, 389, 393, 413,
457, 467, 477 at order 8.

## 4. The target chooses the core ([receipts/converse_census.out](receipts/converse_census.out) (B))

Cores `A` of order ≤ 6 (every root `t`), readers `M_1, …, M_8` and the strong tournaments of order ≤ 5. The odd `v ≤ 125`
with no representation `v = h(A)h(D) + B_t(A)P_D(2)` are exactly `1, 7, 11, 13, 21, 47, 67` (P-10 (a) confirmed, the
dead-end-append table of THM-5263 again). **No odd `v ≤ 125` needs a composite right source** (P-10 (b) refuted):
composite representations exist only at `15, 27, 37, 45, 51, 69, 75, 81, 87, 93, 99, 111, 117, 123` — all `≡ 3 (mod 6)`
except `37 = 5 + 32` — and always beside near-join ones. Below 125 the target never has to choose a composite source, and
the seven exceptions are untouched by composite readers, as THM-5263 (B)'s argument predicts: a prime composite value has a
strong realiser with a digit-0 vertex, and 11, 13, 47, 67 have none.

## 5. Two more things the family can do (PROVED)

**(E) Closure under left multiplication by counts.** For any tournament `R` and near-join `T′ = T → M_k → z`:
`R → T′` is the dead-end append of `z` at `t_0` to `R → T → M_k`, i.e. the near-join with base `R → T`, and THM-5260 (B)
(`B_{t_0}(R → T) = h(R)B_{t_0}(T)` for `t_0` in the second summand) gives `H(R → T′) = h(R)·H(T′)`. So
`{counts} · {near-join values} ⊆ {near-join values}` (P-11: 400 cases, 0 failures,
[receipts/closure_check.out](receipts/closure_check.out)). With THM-4094's reduction of Chan–Pak to primes and the `7p`
lane, the near-join conjecture (every odd number except the seven is a near-join value) reduces — modulo Chan–Pak for the
cofactor, so not as a proof strategy — to the primes `≥ 23` other than 47, 67 and the lanes `7p, 11p, 13p, 21p, 47p, 67p`.

**(S) The sink trick.** For any count `v = H(T)`, the padded sink `s` of `T → s` has `B_s = v` (MISTAKE-776), so the
near-join `(T → s) → M_k → z` with `z → s` has count `(2^k + 1)·v`: `3v, 5v, 9v, 17v, …` are near-join values for every
count `v` (the port `s` is not a dead end, `e_s(T → s) = v`; the near-join formula does not need one). This is the `3q`
lane of THM-5256 once more and adds nothing below 22,227, but it is the cheapest proof that the near-join family is closed
under the Fermat multipliers `2^k + 1`.

## 6. Chan–Pak

Nothing is closed. The converse and THM-5263 identify the strong digit-0 counts with the near-join values below 22,228
and make the two families coincide at every order iff case (C) always has a near-join rewriting. Chan–Pak for the
near-join family is HYP-9453 (open: the density of `B`-sets). The composite readers enlarge the set of linear forms but,
below 125, add no value; the closure (E) and the sink trick (S) are the only new exact tools, and both act by
multiplication, which Chan–Pak already has from ordinal sums.
