# Pépin's test reworked for Collatz, and what a proof that {3, 5, 17, 257, 65537} is complete would need

2026-10-25, klein. **Status.** The completeness of the five known Fermat primes is OPEN (CITED:
`F_n` is composite for `5 <= n <= 32`; every prime factor of `F_n`, `n >= 2`, is `≡ 1 (mod 2^{n+2})`,
Euler–Lucas). Everything proved below is elementary; the "creative arguments" are reformulations and
hostile probes, labelled EXACT, HEURISTIC or OPEN. No proof of completeness is claimed.

## 1. Pépin's test as a permutation statement (EXACT)

Pépin: for `n >= 1`, `F_n = 2^{2^n} + 1` is prime iff `3^{(F_n - 1)/2} ≡ -1 (mod F_n)`. The proof is
that `F_n - 1 = 2^{2^n}` is a power of two, so `3` is a primitive root iff `3^{(F_n-1)/2} ≡ -1`, and
`(3/F_n) = (F_n/3) = (2/3) = -1` by reciprocity since `F_n ≡ 2 (mod 3)`. Equivalently:

**Theorem P (affine form).** For `n >= 1`, `F_n` is prime iff the affine Collatz branch
`x -> 3x + 1` on `Z/F_n` has exactly two cycles: the fixed point `x = -1/2 = 2^{2^n - 1}` and one
cycle of length `2^{2^n}`.
*Proof.* `x -> 3x+1` is conjugate to `x -> 3x` by the translation `x -> x + 1/2` (`2` is invertible
mod `F_n`). Multiplication by `3` fixes `0` and permutes the units; it is a single cycle on the units
iff `3` generates `(Z/F_n)^x`, which for a power-of-two group order is Pépin's condition; if `F_n` is
composite the nonunits add further cycles and the units split. QED.
Data: `n = 1..4` give one unit cycle (orders `4, 16, 256, 65536`); `n = 5` has `ord(3) = 3·2^29·…`
far below `2^32`, many cycles (receipt `fermat-normal.json`).

So Pépin reads: *the "odd branch" of Collatz, reduced modulo a Fermat number, is transitive on
nonzero residues exactly when the number is prime.* The halving branch cannot be reduced modulo an
odd `N` (division by `2^v` depends on the integer, not the residue); that is the exact point where
the Collatz map stops being a permutation and the analogy must change coordinates.

## 2. The Fermat line under Collatz (EXACT; the dual of the Mersenne line)

Write `v(x) = v_2(3x+1)` for odd `x`, `S(x) = (3x+1)/2^{v(x)}` (Syracuse).

**Lemma (constant-valuation-two runs).** For odd `x` and `L >= 0`: the first `L` Syracuse steps of `x`
all have valuation exactly `2` iff `x ≡ 1 (mod 2^{2L+1})`, and then `S^k(x) = 3^k (x-1)/4^k + 1`.
(This is the valuation-two case of THM-4601's residual two-branch count `J = floor((j-1)/2)` for `x = 1 + 2^j u`; the closed form is the only addition.)
*Proof.* `v(x) = 2` iff `3x+1 ≡ 4 (mod 8)` iff `x ≡ 1 (mod 8)`; if `x = 1 + 2^{2L+1} s` then
`S(x) = (3x+1)/4 = 1 + 3·2^{2L-1} s ≡ 1 (mod 2^{2L-1})`; induct. The closed form is the same
induction. QED. (Verified for all odd `x < 20000`, `L <= 19`.)

**Theorem F (Fermat line).** For `n >= 2`, the Syracuse orbit of `F_n = 2^{2^n} + 1` begins with
`2^{n-1} - 1` steps of valuation `2`, then one step of valuation `3`, and arrives at
`(3^{2^{n-1}} + 1)/2`, which is odd:

    F_n -> 3·2^{2^n - 2} + 1 -> 9·2^{2^n - 4} + 1 -> … -> 4·3^{2^{n-1} - 1} + 1 -> (3^{2^{n-1}} + 1)/2.

*Proof.* `F_n ≡ 1 (mod 2^{2^n})` gives `L = 2^{n-1} - 1` by the lemma (`2L + 1 = 2^n - 1 <= 2^n`).
At the end `x = 4·3^{L} + 1`, `3x + 1 = 4(3^{L+1} + 1)` and `3^{L+1} = 3^{2^{n-1}} ≡ 1 (mod 8)`, so the
valuation is exactly `3` and the next value is `(3^{2^{n-1}}+1)/2`, odd because `3^{2^{n-1}} + 1 ≡ 2
(mod 4)`. QED. Receipt: `F_2`: `17 -> 13 -> 5`; `F_3`: `257 -> 193 -> 145 -> 109 -> 41`;
`F_4` and `F_5` as predicted (`fermat-normal.json`).

The survey found no prior analysis of `2^m + 1` or of Fermat numbers as a Collatz family. This is the **dual** of THM-4556's Mersenne line: `2^p - 1 ≡ -1 (mod 2^p)` runs `p - 1` steps of
valuation `1` (growth by `3/2` each, "debt"), while `2^{2^n} + 1 ≡ 1 (mod 2^{2^n})` runs
`2^{n-1} - 1` steps of valuation `2` (decay by `3/4` each, "credit"). Valuations `1` and `2` bracket
the Collatz critical slope `log_2 3 = 1.585` (the barrier of the Antihydra note): the Mersenne line is
the longest pure ascent available at height `2^p`, the Fermat line the longest pure descent available
at height `2^{2^n}` — the two extremal constant-valuation words, selected by the residues `-1` and
`+1` modulo the power of two. The general `2^m + 1` has `floor((m-1)/2)` valuation-two steps
(receipt, `m = 2..40`).

## 3. The dictionary (what transfers, what does not)

| Pépin side (multiplicative, modulo `F_n`) | Collatz side (additive, on integers) |
|---|---|
| `3^{2^k} mod F_n`, `k = 0..2^n - 1`: a tower of squarings | `3^k·2^{2^n - 2k} + 1`, `k = 0..2^{n-1}-1`: a run of "times 3, over 4" |
| endpoint `-1 mod F_n` iff prime | endpoint `(3^{2^{n-1}}+1)/2` always (no primality content) |
| group order `2^{2^n}` is a power of two | residue `1 mod 2^{2^n}` is the run certificate |
| base `3` is a quadratic non-residue (`F_n ≡ 2 mod 3`) | multiplier `3`; `F_n ≡ 2 (mod 3)` means `F_n` is not a Collatz leaf (multiples of 3 have no odd predecessor) |
| transitivity of `x -> 3x+1` on `Z/F_n \ {-1/2}` | the halving branch destroys the permutation structure |

The transfer is exact at the level of the two towers and of the role of `3`. At **equal index** the
line carries no primality information: it is identical for prime and composite `F_n` (both `F_4` and
`F_5` follow Theorem F). But the landing points themselves do: with `L_m = (3^{2^{m-1}}+1)/2` the
landing point of `F_m`, one has `L_{m+1} = L_m^2 + (L_m - 1)^2` (since `2L^2 - 2L + 1 = (t^2+1)/2` for
`L = (t+1)/2`), and `F_n` is prime iff `F_n | L_{2^n}`, which is Pépin's residue `-1` read as a
divisibility of a Collatz landing point one tower level up. This addendum and an independent audit of
Theorems F, M, P are kind-pasteur-2026-10-26-S5's
(`05-knowledge/results/fermat_landing_points_pepin_20261026.md`; THM-4793's `audit` field). The honest
boundary of "Pépin for Collatz" is therefore: prime-blind at equal index, Pépin-complete along the
landing-point tower.

## 4. Creative arguments for completeness, and their hostile probes

1. **Heuristic density (HEURISTIC).** If `F_n` were "random" odd of its size, `P(F_n prime) ≈
   2/ln F_n ≈ 2.9/2^n`; the Euler–Lucas sieve (`p ≡ 1 mod 2^{n+2}`) improves this to about
   `c·n/2^n`? — no: it removes small primes, raising the chance by a factor `O(n)`; either way
   `sum_{n >= 33} O(n)/2^n` is about `10^{-8}`. Expected further Fermat primes: essentially zero.
   Hostile: the same heuristic predicts infinitely many Mersenne primes, and it is a heuristic in
   both cases; it proves nothing.
2. **Pépin as a single-cycle condition (EXACT reformulation).** Completeness says: for every
   `n >= 5` the permutation `x -> 3x+1` of `Z/F_n` has at least three cycles. A proof would need a
   second cycle for every `n`, i.e. an element of `(Z/F_n)^x` outside `<3>` for all `n >= 5` — a
   uniform non-generation statement, which is exactly as hard as the original.
3. **Constructible polygons (EXACT reformulation).** By Gauss–Wantzel, completeness says the odd
   constructible polygon orders are exactly the 31 products of distinct elements of
   `{3,5,17,257,65537}`, i.e. the rows `0..31` of Pascal's triangle modulo 2 read in binary
   (already in `gilbreath_fermat_platonic_20260926.md` §1–2, audited; `fermat-normal.json` re-checks this identity: `prod_{i in S}(2^{2^i}+1) = sum_{T ⊆ S} 2^{sum_T 2^i}`
   is the binary number whose ones sit at the subsets of the bits of `S`, Lucas' theorem). A new Fermat
   prime `F_n` would add a 33rd row-type with a one at bit `2^n`. Geometrically: completeness is
   the claim that the only regular polygons with constructible vertices and odd order are the 31
   binary Sierpiński rows below `2^{32}`.
4. **Paley graphs (EXACT reformulation).** `F_n` prime ⟹ the Paley graph `P(F_n)` is a
   self-complementary strongly regular graph on `F_n` vertices with spectrum `{(F_n-1)/2,
   (-1 ± sqrt F_n)/2}`; `P(5) = C_5` and `P(17)` (clique number 3, the Ramsey `R(4,4)` witness) are
   its two smallest members. Completeness would say the Fermat Paley family has exactly five members.
   No graph-theoretic obstruction is known to a sixth; the spectrum is constructible by Gauss sums
   for every prime `≡ 1 (mod 4)`.
5. **Bernoulli denominators (EXACT reformulation, inherited).** By von Staudt–Clausen,
   `denom(B_{2^r}) = 2·prod_{2^i <= r} F_i` over the Fermat primes, so completeness is the statement
   that these denominators are constant, `2(2^{32} - 1)`, for all `r >= 16`
   (`forest_20260926_bernoulli.md`). A sixth Fermat prime would appear in `denom(B_{2^{32}})`.
6. **Tournaments (EXACT reformulation, inherited).** `F_i = H(NT_{2^i+2})` (THM-4749 via Moon's
   `H(NT_{k+2}) = 2^k + 1`), and Fermat primes are exactly the primes where every circulant
   tournament is rigid (THM-4644). Completeness says the near-transitive path counts `2^k + 1` are
   prime only for `k in {1,2,4,8,16}`, and that `3,5,17,257,65537` are the only "all-rigid" primes.
   Neither reading supplies a proof mechanism; both are faithful restatements.
7. **The "tricky" route and why it fails (hostile).** One might hope to show `3^{2^{2^n - 1}} ≠ -1
   (mod F_n)` for `n >= 5` by a congruence on the tower `3^{2^k}` modulo small primes dividing
   `F_n`: but `F_n` has no small prime factors (Euler–Lucas), and the tower modulo `F_n` has no
   known structure beyond Pépin's own statement. The Collatz dictionary does not help: Theorem F
   is prime-blind.

Conclusion: no argument here, creative or otherwise, proves completeness; the exact content is
Theorem P, Theorem F, the Sierpiński identity, and the Paley reformulation, each a correct
restatement with its loss recorded.

## 2a. The line pairs merge by the `4x+1` rule (EXACT)

Let `m = 2L + 1`, `L >= 2`. After their common run of `L` valuation-two steps, `2^m + 1` and
`2^{m+1} + 1` land on `a = 2·3^L + 1` and `b = 4·3^L + 1`.

**Theorem M.** If `L` is even, the two Syracuse orbits first coincide after exactly `L + 2` odd
steps; if `L ≡ 3 (mod 4)`, after exactly `L + 3`; if `L ≡ 1 (mod 4)`, not within `L + 3` steps.
*Proof.* `a ≡ 3 (mod 4)` gives `S(a) = 3^{L+1} + 2 =: p`. For `L` even, `3^{L+1} ≡ 3 (mod 8)`, so
`3b + 1 = 4(3^{L+1} + 1)` has valuation `4` and `S(b) = (3^{L+1}+1)/4 =: q`, whence `p = 4q + 1`
and `3p + 1 = 4(3q + 1)`: `S(p) = S(q)`, the classical rule that `x` and `4x+1` share a successor.
Before that step all values differ (the run values carry different powers of two, and `a != b`,
`p != q`). For `L` odd, `3^{L+1} ≡ 1 (mod 8)`, `S(b) = (3^{L+1}+1)/2 =: q'` and `p = 2q' + 1`; one
more step gives `S(p) = (3^{L+2}+7)/2` and, when `L ≡ 3 (mod 4)` (so `3^{L+2} ≡ 11 (mod 16)`),
`S(q') = (3^{L+2}+5)/8` with `S(p) = 4 S(q') + 1`, merging at the next step; when `L ≡ 1 (mod 4)`
the valuation of `3q'+1` is at least four and no `4x+1` relation holds at that depth. QED.
Receipt `merge-normal.json` confirms the exact indices for `L = 2..19`.

THM-4556 (ii) is the Mersenne twin: `M_{2k}` and `M_{2k-1}` merge at equal time. Theorem M is the Fermat-side statement with the exact merge depth. It is also an explicit instance of THM-4581's pair-chain coalescence on a named line: the two
landing points are `3^L c + 1` with `c in {2, 4}`, and the merge time is `L + 2` or `L + 3`, not
merely "almost surely finite". The equal total stopping times observed in blocks of `m` (for
example `32` for `m = 15..19` and `56` for `m = 21..26`) are these merges plus further coincidences
of the merged tails.
