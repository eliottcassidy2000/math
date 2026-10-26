# Pépin's test is the divisibility of a Collatz landing point: the Fermat line of F_{2^n} lands on L_{2^n}, F_n is prime iff F_n | L_{2^n}, and the landing points obey L' = L² + (L−1)²

**kind-pasteur-2026-10-26-S5.** Owner dispatch: creative arguments that `{3,5,17,257,65537}` is comprehensive;
`{1,2,4,8,16}` and `{0,1,2,3,4}` as groups of five against the repository's other fives; Pépin's theorem reworked
for Collatz, geometrically and graph-theoretically; the constructible orders `3,4,5,6,8,10,12,15,16,17,20,24,30,32,34,…`.

**Status.**
- **INHERITED, first-pushed (klein, 2026-10-25):** the packet
  [fermat_five_pepin_collatz_20261025](fermat_five_pepin_collatz_20261025/README.md) and
  [THM-4793](../../01-canon/theorems/THM-4793-the-fermat-line-under-collatz-valuation-two-descent-exact-merge-of-consecutive-lines-and-pepin-as-an-affine-cycle-type.md)
  already contain Pépin in affine-cycle form, the Fermat line, the exact merges, the Sierpiński identity, the Fermat
  Paley graphs, the ladder `{0..4} → {1,2,4,8,16} → F_i` with its tournament floor, and a typed survey of the repository's
  other fives (verdict: coincidences of cardinality). This note does not restate them as new.
- **AUDIT (this session):** THM-4793's four statements are CONFIRMED by independent code (checks A1–A8 of
  [`fermat_landing_points_pepin_kps_20261026.py`](../../04-computation/experiments/fermat_landing_points_pepin_kps_20261026.py));
  an `audit:` line was added to the theorem file.
- **PROVED, new here (Theorem L, with proof):** the landing point `L_m = (3^{2^{m-1}} + 1)/2` of the Fermat line of
  `F_m` satisfies `L_{m+1} = 2L_m^2 - 2L_m + 1 = L_m^2 + (L_m - 1)^2`, and for `n >= 1`, **`F_n` is prime iff `F_n` divides
  `L_{2^n}`**. The packet's dictionary line "endpoint `(3^{2^{n-1}}+1)/2` always, no primality content" is correct at
  equal index and false one tower level up: the endpoint of the line of `F_{2^n}` carries exactly Pépin's test for `F_n`.
  Euler-criterion generalization and its hostiles (Euler pseudoprimes base 3, including `L_4 = 3281` itself) are proved
  and listed.
- **FINITE-EXACT:** 28 checks; `omega(Paley(257)) = 7` by two independent searches; primality of `L_2..L_10`.
- **CITED:** Boklan–Conway, "Expect at most one billionth of a new Fermat prime!", Math. Intelligencer 38 (2016),
  arXiv:1605.01371 (probability of a sixth Fermat prime below `10^{-9}`, heuristic); `F_5..F_32` composite, `F_33` of
  unknown character (status pages, 2025).
- **OPEN:** comprehensiveness of the five. Nothing here proves it; the "tricky" arguments are reformulations with their
  losses recorded, as in the packet's §4.

## 0. Inheritance, board, and what was already ruled

- Closest proved mechanisms: THM-4793 (run law, Fermat line, merges, affine Pépin), THM-4556/4605 (the Mersenne line),
  THM-4748/4749 (`F_i = H(NT_{2^i+2})`), THM-4644/871 (Fermat primes are the all-rigid circulant orders), the Sierpiński
  identity (`gilbreath_fermat_platonic_20260926`, packet §2), and the Collatz-sheet alternation modulo `m | 2^t + 1` of
  [sheets_root27_rank_untouchables_20261025](sheets_root27_rank_untouchables_20261025.md) Theorem A (Fermat moduli are
  the case `t = 2^n`, where negation is `2^n` halvings).
- Canonical hostiles: `F_5` (the line is prime-blind at equal index); the Euler pseudoprimes base 3 below `10^5`
  (`703, 1891, 3281, 8911, 12403, 16531, 44287, 63139, 79003, 97567`): composite numbers that pass the landing-point
  divisibility of Theorem L.2, so the divisibility is a primality test only for Proth-shaped moduli.
- Corrected near miss: the packet's own dictionary (§3) and conclusion that the Collatz dictionary carries no primality
  information; repaired here by changing index, not by changing the dictionary.
- Least-used sidecar: the landing points as a sequence in their own right (base-3 generalized Fermat numbers).
- Board: the base-2 tower `F_n = 2^{2^n}+1` / the base-3 tower `L_m = (3^{2^{m-1}}+1)/2` / Pépin's squaring chain /
  Syracuse valuation-two runs / Proth shape / Euler criterion.

## 1. Theorem L: Pépin as divisibility of a landing point

Setting. `S(x) = (3x+1)/2^{v_2(3x+1)}` on odd `x`; `F_n = 2^{2^n}+1`; by THM-4793 (Fermat line), for `m >= 2` the
Syracuse orbit of `F_m` runs through `3^k 2^{2^m - 2k} + 1` for `k < 2^{m-1}` and lands after `2^{m-1}` steps on the odd
number `L_m := (3^{2^{m-1}} + 1)/2`. Put `L_1 := 2` (the formal value `(3+1)/2`; the actual orbit of `F_1 = 5` lands on 1).

**Theorem L.**
1. (Recursion.) `L_{m+1} = 2L_m^2 - 2L_m + 1 = L_m^2 + (L_m - 1)^2` for all `m >= 1`. So each landing point is a sum of
   two consecutive squares (a centered square number): `5 = 2^2 + 1^2`, `41 = 5^2 + 4^2`, `3281 = 41^2 + 40^2`,
   `21523361 = 3281^2 + 3280^2`. Geometrically `L_{m+1}` is the number of lattice points of the taxicab ball
   `|x| + |y| <= L_m - 1`.
2. (Pépin.) For `n >= 1`: `F_n` is prime iff `F_n | L_{2^n}`, i.e. iff `F_n` divides the Syracuse landing point of the
   Fermat line of `F_{2^n}`. Equivalently: iterate `L -> L^2 + (L-1)^2` from `L = 2` for `2^n - 1` steps modulo `F_n` and
   test for `0`. Instances: `F_1 = 5 | L_2 = 5`; `F_2 = 17 | L_4 = 3281 = 17·193`; `F_3 | L_8`; `F_4 | L_{16}`;
   `F_5 ∤ L_{32}` (residue `5162152`).
3. (Euler criterion in Collatz form.) For every odd `N >= 3`, with `M = (N-1)/2`: the orbit of `4^M + 1` passes through
   `3^j 4^{M-j} + 1` for `j < M` and lands on `(3^M + 1)/2^{v_2(3^M+1)}`, and
   `N` divides that landing point iff `3^{(N-1)/2} ≡ -1 (mod N)`.
   For `N = F_n` this is Pépin and hence a primality test; for general `N` it is Euler's criterion with residue `-1`,
   passed by every prime `N ≡ ±5 (mod 12)` and by the Euler pseudoprimes base 3.
4. (The landing points always pass.) For `m >= 2`, `v_2(L_m - 1) = m` and `3^{(L_m - 1)/2} ≡ -1 (mod L_m)`: every
   landing point is prime or an Euler pseudoprime base 3. `L_2, L_3, L_5, L_6` are prime, `L_7` is a strong probable
   prime, `L_4 = 17·193`, `L_8, L_9, L_{10}` are composite (Miller–Rabin, exact for the composites).

*Proof.* (1) `3^{2^m} = (3^{2^{m-1}})^2 = (2L_m - 1)^2`, so `L_{m+1} = ((2L_m-1)^2 + 1)/2 = 2L_m^2 - 2L_m + 1`.
The lattice count is the standard centered-square identity. (2) `F_{2^n} = 2^{2^{2^n}} + 1 = 4^M + 1` with
`M = 2^{2^n - 1} = (F_n - 1)/2`, so by the Fermat line `L_{2^n} = (3^M + 1)/2`; Pépin's theorem says `F_n` is prime iff
`3^M ≡ -1 (mod F_n)`, i.e. iff `F_n | 3^M + 1`, i.e. iff `F_n | L_{2^n}` because `F_n` is odd. The recursion form is (1).
(3) The staircase `S(3^j 4^{M-j} + 1) = 3^{j+1} 4^{M-j-1} + 1` for `j < M - 1` is the identity `3(3^j 4^{M-j}+1) + 1 =
4(3^{j+1} 4^{M-j-1} + 1)`, and the last step lands on `(3^M+1)/2^{v}`; `N` odd divides it iff `N | 3^M + 1`. (4) By the
lifting-the-exponent lemma `v_2(3^{2^{m-1}} - 1) = m + 1`, so `(L_m - 1)/2 = (3^{2^{m-1}} - 1)/4 = 2^{m-1}·(odd)`, and
`3^{(L_m-1)/2} = (3^{2^{m-1}})^{odd} ≡ (-1)^{odd} = -1 (mod L_m)` since `L_m | 3^{2^{m-1}} + 1`. `QED`

**Reading.** The packet separated a multiplicative tower (`3^{2^k} mod F_n`) from an additive one
(`3^k 2^{2^n-2k} + 1`) and concluded that only the base 3 and the power-of-two certificate are shared. Theorem L says the
two towers are the same object read at different indices: the Pépin iterates `3^{2^k}` are the landing points
`2L_{k+1} - 1`, and Pépin's final comparison with `-1` is the divisibility of a landing point. "Pépin reworked for
Collatz" is therefore exact: **the primality of `F_n` is decided by whether the Syracuse orbit of `F_{2^n}` lands, after
`(F_n - 1)/2` forced steps, on a multiple of `F_n`.** The loss ledger: the test needs the Proth shape `F_n = 1·2^{2^n} + 1`;
for other moduli (3) is only Euler's criterion, and the hostiles of §0 show composite passers.

## 2. The fives, with the floors this note adds

| floor | five | meaning |
|---|---|---|
| index | `0,1,2,3,4` | `n`; also `Z/5` |
| exponent | `1,2,4,8,16` | `2^n`; `<2>` in `(Z/31)^*` (packet); also the first half of `<2>` in `(Z/17)^*`, whose second half `16,15,13,9` is its negatives (`16 = -1`): the same five sit inside the Mersenne prime `2^5 - 1` as a closed group and inside the Fermat prime `2^4 + 1` as a half-clock |
| halving clock | `2,4,8,16,32` | `ord_{F_n}(2) = 2^{n+1}` with `2^{2^n} = -1`, prime or not (the sheet-alternation shift of Theorem A, previous session, is `t = 2^n`) |
| squarings | `0,1,3,7,15` | `2^n - 1` squarings in Pépin's test = recursion steps from `L_1 = 2` to `L_{2^n}`: the Mersenne numbers |
| cyclotomic index | `2,4,8,16,32` | `F_n = Phi_{2^{n+1}}(2)`; among `k <= 64` with `Phi_k(2)` prime, the powers of two are exactly `2,4,8,16,32` (and `64` is not: `F_5`) |
| prime | `3,5,17,257,65537` | product `2^{32}-1`; Pascal rows `1,2,4,8,16` (packet) |
| landing | `L_2..L_6 = 5, 41, 3281, 21523361, 926510094425921` | the base-3 tower; `F_n | L_{2^n}` iff prime |

The owner's two sets are thus the index and exponent floors of one tower whose prime floor is the Fermat five and whose
dual base-3 floor is the landing-point sequence. On the repository's other fives the packet's verdict stands (coincidences
of cardinality); this note checked no further identification and claims none. The Paley floor: `Paley(5) = C_5`,
`Paley(17)` has clique number 3 (packet), and `Paley(257)` has clique number **7** (two independent exact searches), so the
Fermat Paley graphs give the Ramsey bounds `R(3,3) > 5`, `R(4,4) > 17`, `R(8,8) > 257`, the last far below the published
lower bounds for `R(8,8)`. Recorded as data, not as a mechanism.

## 3. Comprehensiveness: the arguments, typed

- **Exact reformulations** (packet §4.2–4.6 plus this note): no `n >= 5` with `F_n | L_{2^n}`; no 33rd Sierpiński row
  type among odd constructible orders; the von Staudt–Clausen denominators of `B_{2^r}` constant for `r >= 16`; the
  affine clock `x -> 3x+1` on `Z/F_n` never a single `(F_n-1)`-cycle for `n >= 5`; no sixth Fermat Paley graph; no
  near-transitive path count `2^k + 1` prime beyond `k = 16`; no prime `p > 65537` with `p - 1` a power of two, i.e. no
  prime beyond the five at which every quadratic non-residue is a primitive root.
- **Heuristic (CITED):** Boklan–Conway refine Hardy–Wright: with the Euler–Lucas sieve `p ≡ 1 (mod 2^{n+2})` on the
  prime factors, the expected number of Fermat primes beyond `F_4` is below `10^{-9}`. Hostile (packet): the same
  heuristic predicts infinitely many Mersenne primes; it proves nothing in either direction.
- **Why no tricky argument closes it:** `F_n = Phi_{2^{n+1}}(2)` is a cyclotomic value, so no algebraic factorization
  exists (Aurifeuille's split `2^{4k+2} + 1 = (2^{2k+1} - 2^{k+1} + 1)(2^{2k+1} + 2^{k+1} + 1)` needs an exponent
  `≡ 2 (mod 4)` and touches only `F_1 = 5 = 1·5`); every prime factor is `≡ 1 (mod 2^{n+2})`, so small-prime congruences
  are silent; and Theorem L.4 shows the dual tower always passes the base-3 test, so no base-3 argument can be turned
  against the base-2 tower. A proof would have to control `3^{2^k} mod F_n` for all `n`, which is Pépin's statement
  itself. Verdict: OPEN, with the exact content being the reformulations above.

## 4. Open and next

- Is the landing-point form of Pépin in the literature? The identity `L_{m+1} = 2L_m^2 - 2L_m + 1` is the standard
  squaring chain written in `(s+1)/2` coordinates; the Collatz reading (landing points of Fermat lines) was not found in
  the packet's survey. No priority claim is made.
- The base-3 tower: which `L_m` are prime? `m = 2,3,5,6` prime, `7` strong probable prime, `4,8,9,10` composite; compare
  OEIS for `(3^{2^k}+1)/2` before claiming anything.
- Theorem A of the previous session gives the alternation set `{m odd : m | 2^t + 1}`; Fermat numbers are its primitive
  members. Whether the alternation shift carries any Pépin content at composite `F_n` is answered negatively here
  (`ord_{F_n}(2) = 2^{n+1}` prime or not).

## 5. Reproduction

```bash
python3 04-computation/experiments/fermat_landing_points_pepin_kps_20261026.py
```

Standard library; about two minutes (the Paley(257) searches dominate); 28 checks; output beside the script.
