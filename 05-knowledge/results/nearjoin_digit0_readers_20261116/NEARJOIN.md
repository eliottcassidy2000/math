# Why the counts with no strong digit-0 vertex are the near-join exceptions

Session klein-2026-11-16-nearjoin. Scripts `scripts/nearjoin_check.py`, `scripts/digit0_structure.py`; receipts
`nearjoin_check.out`, `digit0_structure.out`; the strong census `strongB_3..9.out` of the previous packet. The claim being
chased is the corrected one (MISTAKE-776): 11, 13, 47, 67 are the counts below 125 none of whose **strong** realisers has a
vertex `u` with `B_u = H` (a padded sink always has `B = H`), and these are HYP-9453's near-join exceptions other than
the holes and 1.

## 1. Three identities (PROVED; `nearjoin_check.out`)

For a tournament `S`, a vertex `u`, `B_u = Σ_{A ⊆ V∖u} h(S[A]) s_u(S[V∖A])` (THM-4865) and the digit `σ_u = B_u − H(S)`:

- **(I1) Backward junctions.** `B_u` counts the ordered two-path covers `(P_1, P_2)` of `S` with `P_2` starting at `u`
  and `P_1` possibly empty; cutting a Hamiltonian path before `u` shows `H` counts those with `P_1` empty or
  `last(P_1) → u`. Hence `σ_u` = the number of covers with `P_1` nonempty and `u → last(P_1)` ("backward junctions").
  Digit 0 at `u` means no such cover exists. (3,610 rooted tournaments of orders 3–7, no failure.)
- **(I2) Sole exits.** If `u = z` has out-degree 1 with out-neighbour `t_0`, then `σ_z = e_{t_0}(S − z)`, the number
  of Hamiltonian paths of `S − z` ending at the port: `z` is a digit-0 vertex iff no Hamiltonian path of `S − z` ends
  at `t_0`. (Paths from `z` go `z → t_0 → …`; the one-exit law THM-4862 (a) gives `H(S) = h′ − e′ + B′` with
  `h′ = H(S − z)`, `e′ = e_{t_0}(S − z)`, `B′ = B_{t_0}(S − z)`, and `B_z = h′ + B′`; 380 cases, no failure.)
- **(I3) The near-join.** `T′ = T → M_k → z` with `z → t_0` (spectrum_analogies_nearjoin_20261102) is the rooted
  append of `z` at `t_0` to the base `T → M_k`, where no Hamiltonian path ends at `t_0`; so `H(T′) = h + 2^k B_{t_0}(T)`
  (the suffix law of THM-5260 (B): `B_{t_0}(T → M_k) = 2^k B_{t_0}(T)`) and **`z` is a digit-0 vertex of `T′`**
  (1,260 cases `(T, t_0, k)`, orders ≤ 6, `k ≤ 3`, no failure — the packet's own check, re-verified). When `T` is
  strong, `T′` is strong; in general `z`'s strong component carries the full count exactly when the part of `T` not
  reachable from `t_0` is transitive (all 69 failures among the 1,260 have a non-strong base).

## 2. Digit-0 vertices of strong tournaments are sole exits (FINITE-EXACT, orders 3–8)

Over all 50,782 rooted strong tournaments of orders 3–8, every one of the 344 digit-0 vertices has out-degree exactly
1 (`digit0_structure.out`); the first-order consequence of (I1) — no Hamiltonian path of `S − x` starts at `u` for
any out-neighbour `x` — holds at 2,618 further vertices that are not digit-0, so digit 0 is strictly stronger. With
(I2) this gives the working characterisation: **a strong tournament has a digit-0 vertex iff it is the rooted append
of a vertex `z` at a port `t_0` of a base `B′ = S − z` with `e_{t_0}(B′) = 0`**, and then `H(S) = h(B′) + B_{t_0}(B′)`.
(The out-degree-1 statement is proved here only by exhaustion; (I1) gives the necessary conditions
`N^+(u) ∩ C_m = ∅` for the terminal component `C_m` of `S − u`, and the same for every `S − P_2`.)

## 3. The theorem (THM-5263)

- **(A) Every prime near-join value has a strong realiser with a digit-0 vertex** (its own `z`), and more generally
  every count of the "dead-end append" family `{h(B′) + B_{t_0}(B′) : e_{t_0}(B′) = 0}` with the append strong.
  The near-join (`k ≥ 1`) is the sub-family with base `T → M_k`.
- **(B) 11, 13, 47, 67 are never near-join values, from any base of any order**, and never dead-end appends with a
  strong result: every strong realiser of these primes has order ≤ 9 (strong floor `f(10) = 125`) and the exhaustive
  census finds no digit-0 vertex in any of them. With the holes 7, 21 and the trivial 1 this makes HYP-9453's seven
  exceptions **permanent**: the near-join family misses exactly `{1, 7, 11, 13, 21, 47, 67}` below 125 whatever the
  base order, and the hypothesis is precisely the statement that it misses nothing else.
- **(C) The dead-end append family already from bases of order ≤ 7 realises every other odd number below 125**
  (`digit0_structure.out`: missed exactly `1, 7, 11, 13, 21, 47, 67`), while the near-join proper (`k ≥ 1`) from the
  same bases misses 21 numbers; so below 125 the counts with a strong digit-0 vertex are exactly the odd numbers other
  than the seven, and the two lists coincide for a reason, not by accident: **a count has a strong digit-0 vertex iff
  it is a dead-end append, and the near-join is the dead-end append whose base ends in a transitive chain.**

## 4. Why the converse is not "the near-join of `S − u`"

Of the 55 strong `(S, u)` with a digit-0 vertex at orders 4–7, only 23 have `H(S) = h(S − u) + 2^k B_t(S − u)` for
some port `t` and `k ≥ 1` (pre-registered P-3 confirmed): the base of the near-join is not `S − u` in general. What
holds exactly is (I2): `H(S) = h(S − u) + B_{t_0}(S − u)` with `t_0` the sole exit — the `k = 0` member of the family,
which the near-join census never enumerated. Reading HYP-9453 and HYP-9454 together: the near-join needs the base to
end in a transitive chain so that `e_{t_0} = 0` for free; the dead-end append asks only for `e_{t_0}(B′) = 0`.

## 5. The field with one element

The repository's `F_1`-tournament is the linear order (CONCEPT-MAP row 29): one Hamiltonian path, and its sink is the
universal digit-0 vertex (`B_s = H` for every sink, which is MISTAKE-776). A strong digit-0 vertex is an *internal*
sink: a sole exit whose port is a dead end for the rest. The `F_1`-point itself (one vertex, `h = 1`, `B = 1`) has the
near-join ray `1 + 2^k = 3, 5, 9, 17, …` — THM-5126's line `L_1` is the orbit of the field with one element under the
near-join — and the two holes are the first points of the next two lines `L_3`, `L_5` (COLLATZ_ROWS). In this reading
the seven exceptions are: the `F_1`-point (1), the two holes, and the four primes whose strong realisers contain no
internal sink.
