# Two sheets modulo odd m, the shape-free rank on the root-27 ray, and the untouchable share: what is exact, what is Goldbach, what is golden

**kind-pasteur-2026-10-25-S4.** Owner dispatch: merge the prior work on modular tables mod odd `n` alternating
between the `3n-1` and `3n+1` graphs; pursue a source-preserving rank for arbitrarily changing tree shapes on the
root-27 orbit, with exact carries and guards and hostile tests; use the Frobenius-defect and clone-merger machinery
where useful; consider `{2,3,5,7}` split as `{3,5}` (15) and `{2,7}` (14); look for or rule out a second odd
untouchable number, or find the exact share of untouchable numbers.

**Status.**
- **PROVED (elementary, complete proofs below):** Theorem A (sheet tables modulo odd `m`: negation conjugacy, the
  halving-shift alternation exactly when `-1` is a power of `2` mod `m`, and its failure mod 7, 15), Theorem B (the
  rank `rho_z(w) = |w| - v_2(Q_z(w))` of a ternary-tree word depends on the shape only through the first disagreement
  with the source's carry word, so no shape change creates a certificate; and the Haar share of halting odd
  sources is `sqrt(5) - 2`), Theorem C (an odd `m >= 9` is touchable whenever `m - 1` is a sum of two distinct primes,
  so a second odd untouchable number is a distinct-Goldbach exception).
- **FINITE-EXACT:** three scripts, 49 checks in all, every count an exact integer:
  [`sheet_tables_mod_odd_kps_20261025.py`](../../04-computation/experiments/sheet_tables_mod_odd_kps_20261025.py) (7),
  [`antihydra_tree_rank_kps_20261025.py`](../../04-computation/experiments/antihydra_tree_rank_kps_20261025.py) (30),
  [`untouchable_census_goldbach_dependence_kps_20261025.py`](../../04-computation/experiments/untouchable_census_goldbach_dependence_kps_20261025.py) (12);
  `.out` files beside them. The untouchable census is complete to `10^6` by a proved case split (no sieve factor),
  reproduces A057978 at `10, ..., 10^6` and the first forty terms of A005114, and agrees with an independent brute
  force to `3000`.
- **CITED:** Erdős 1973 (positive lower density of untouchables), te Riele 1976 (`> 0.0324`), Banks–Luca (`> 1/48`),
  Chen–Zhao 2011 (`>= 0.06`), Pollack–Pomerance 2016 (heuristic density about `0.17`, Trans. AMS Ser. B 3), OEIS
  A057978 (counts to `10^12`, density `0.1711` there), Oliveira e Silva–Herzog–Pardi 2014 (Goldbach verified to
  `4·10^18`; record minimal Goldbach prime `9781`, from the authors' page).
- **OPEN, not touched:** the root-27 inequality (Antihydra nonhalting), Goldbach, the existence of the untouchable
  density. **NOT CLAIMED:** any Collatz statement; any Frobenius-defect or clone-merger application (section 4 records
  why neither machine applies here).

## 0. Inheritance, board, concurrent work

- Closest proved mechanisms: the sheet involution `x -> -x` conjugating `3n+1` to `3n-1`
  ([THM-4474 (E)](../../01-canon/theorems/THM-4474-strategy-cube-provability-by-cycle-densities.md), THM-4610's
  `F_p(x) = -C_p(-x)`, klein's [Antihydra A1](fixed_twists_conductors_reverse_potential_20261025/antihydra/ANTIHYDRA.md));
  the root-27 packet's carry-cylinder theorem (C3), tree language (C1–C2), leaf-expansion law (L1–L2) and credit
  `q_J = 1 + 3J - K_J` ([anchor](root27_machine_bridges_20261026/anchor/README.md)); THM-2228 (Mahler carry tails);
  the aliquot note's D26 ([collatz_aliquot_lehmer_five_20260927](collatz_aliquot_lehmer_five_20260927.md): untouchables
  as leaves, Pollack–Pomerance's `0.17`).
- Canonical hostiles: the unit translation (it conjugates the affine steps but not the Syracuse maps); mod 15, where
  both prime factors alternate separately but with incompatible shifts; the halting sources `5, 9, 13, ...` whose
  carry prefixes are complete trees (positive controls for the rank law); the odd number 3, touchable only through
  the prime square 4.
- Corrected near misses: MISTAKE-529 (a repeated numeral read as a bridge) governs the `{3,5}/{2,7}` item; the
  root-27 audit's warning that only the canonical frontier guarantees the first mismatch is a leaf.
- Least-used sidecars: the order of `2` modulo `m` (it decides the alternation); the Haar measure of carry cylinders
  (it turns klein's golden model into a theorem for random sources); the semiprime-only sources of odd values.
- Board: sheet / modulus / halving count; carry word `kappa(z)` / first mismatch / credit; Goldbach partition /
  non-semiprime source / residue class.
- Concurrent work: the root-27 packet (2026-10-26) already contains (C3) and (L2); Theorem B restates them as a
  single rank law, adds the exhaustive all-shape controls, and adds the halting-share theorem. Nothing in Theorem B's
  first half is claimed as new.

## 1. Theorem A: the two sheets modulo an odd modulus

For odd `m`, `k >= 1`, and odd `x`, write the edge tables
`E_k^+(m) = {(x mod m, (3x+1)·2^{-k} mod m)}` and `E_k^-(m) = {(x mod m, (3x-1)·2^{-k} mod m)}`. The halving count `k`
is not a function of `x mod m`, so a modular table of a sheet is a family indexed by `k`; this is the honest form of
"the arithmetic table mod odd `n`".

**Theorem A.**
1. (Affine step only.) `x -> 3x-1` is the conjugate of `x -> 3x+1` by the unit translation `x -> x-1`, over `Z` and
   over every `Z/m`. The translation does not commute with halving, so the Syracuse maps `S^+`, `S^-` are not
   conjugate by it (counterexamples among odd `x < 2000`: check T1b).
2. (Negation.) For every odd `m` and every `k`, `E_k^-(m) = {(-x, -y) : (x,y) in E_k^+(m)}`. This is the sheet
   involution read modulo `m`.
3. (Alternation by a halving shift.) Suppose `-1` lies in the subgroup generated by `2` in `(Z/m)^*`, and let `t` be
   the half of the order of `2`, so `2^t = -1 (mod m)`. Then for every `k`,
   `E_{k+t}^-(m) = {(-x, y) : (x,y) in E_k^+(m)}`:
   the `3n-1` table with `t` more halvings is the `3n+1` table read from the negated source, with the same target
   residue. Modulo 3 (`t = 1`) the tables alternate with the parity of `k`: the plus sheet lands in residue 2 for odd
   `k` and 1 for even `k`, the minus sheet the reverse. Modulo 5, `t = 2`; modulo 9, `t = 3`; modulo 27, `t = 9`.
4. (Where it fails.) If `-1` is not a power of `2` modulo `m`, no halving shift relates the two tables in this way:
   `m = 7` (`2` has order 3), `m = 15` (the shifts `t = 1` mod 3 and `t = 2` mod 5 are incompatible), `21, 23, 31, 45, 63`.

*Proof.* (1) `(3(x-1)+1)+1 = 3x-1`. (2) `(3(-x)-1)·2^{-k} = -(3x+1)·2^{-k}`. (3) With `2^{-t} = -1`,
`(3(-x)-1)·2^{-(k+t)} = (-(3x+1))·2^{-k}·(-1) = (3x+1)·2^{-k}`. (4) Equality of the tables for all odd `x` forces,
comparing two sources, `2^{-k'} = -2^{-k}` modulo `m`, i.e. `-1` a power of `2`. The finite checks (T1–T4, odd `m <= 101`,
`k <= 6`) are controls, not the proof. `QED`

**The `{3,5}` versus `{2,7}` split, typed.** The prime 2 is the halving prime: the shift `t` lives in the group
generated by `2`. Modulo 3 and modulo 5 the sheets alternate (`t = 1, 2`); modulo 7 they do not, because `2` has odd
order 3; modulo 15 they do not, because the two alternations are out of phase. This is an exact statement about
which odd moduli see the two sheets as one table, and it separates `{3,5}` from `7` for a reason (the order of `2`).
It is not a statement about 14 or 15 as numbers, and no identification with LRC(14), the fifteens of session 65, or
the C15 symmetry of the double-triple packet is claimed; those are different objects (MISTAKE-529 applies).

## 2. Theorem B: the shape-free rank on the root-27 ray

Setting (root-27 anchor, inherited). `T(z) = ceil(3z/2) = (3z+c)/2` with carry bit `c = z mod 2`; credit `q -> q+2`
on `c = 1` and `q -> q-1` on `c = 0`, from `q = 1` at the blank source `z = 7`. A halt certificate is a full
ternary-tree word `w in {0,1}^n` (`Tree = 0 | 1 Tree Tree Tree`) with `v_2(Q_z(w)) >= n`, where
`Q_z(w) = 3^n z + C_w`, `C_w = sum_{i<n} c_i 2^i 3^{n-1-i}`. The open target is `v_2(Q_7(w)) < |w|` for every tree word.

Let `kappa(z) in {0,1}^N` be the carry word of `z`: the successive bits `c` along `z, T(z), T^2(z), ...`. For a binary
word `w` let `LCP(w, kappa(z))` be the length of its longest common prefix with `kappa(z)`.

**Theorem B.**
1. (Carry word is length-free.) For every `n` and `j <= n`, `3^n z + sum_{i<j} kappa_i 2^i 3^{n-1-i} = 0 (mod 2^j)`.
   Thus one infinite word serves every length `n`; the congruences of (C3) are homogeneous under multiplication by 3.
2. (Rank law.) For every binary word `w` of length `n`, tree word or not,
   `v_2(Q_z(w)) = LCP(w, kappa(z))` if this is `< n`, and `v_2(Q_z(w)) >= n` iff `w` is the length-`n` prefix of
   `kappa(z)`. Hence the **rank** `rho_z(w) = |w| - v_2(Q_z(w))` depends on the shape of `w` only through the position
   of its first disagreement with the source's carry word. Replacing a leaf at position `p` by `1000` changes the rank
   by exactly `+3 + (h - min(p,h))` when `p != h`, and by `3 - delta` with `delta >= 1` when `p = h` (this is the
   packet's (L2)). No sequence of shape changes can lower the rank below the frontier's.
3. (Equivalent target.) The root-27 inequality holds iff no prefix of `kappa(7)` is a complete tree word, i.e. iff the
   walk with steps `+2` on a `1` and `-1` on a `0`, started at `1`, never returns to `0` along `kappa(7)`. Along the
   frontier the rank equals `q - LCP(0^q, kappa(7) after the prefix)`, where `q` is the current credit: the rank is
   the number of pending slots not yet confirmed empty.
4. (Halting share of random sources.) With respect to Haar measure on the odd 2-adic integers, the set of sources
   whose credit walk returns to `0` has measure exactly `phi^{-3} = sqrt(5) - 2 = 0.2360...`; conditioned on the
   forced prefix `1110` of the blank source, the measure is `phi^{-6} = 9 - 4·sqrt(5) = 0.0557...`, klein's model
   number. The finite census gives `4709` of the `20001` odd sources `z <= 40001` halting within 500 symbols, share
   `0.2354`; `465` of `2001` for `z <= 4001` within 20000 symbols. The halting lengths are `1 mod 3`, with `2500`
   sources halting at length 4 (exactly the residue class of `z` modulo 16 with carry prefix `1000`).

*Proof.* (1) If `3^n z + sum_{i<j} kappa_i 2^i 3^{n-1-i} = 0 (mod 2^j)` then multiplying by 3 gives the same
congruence for `n+1`; the base case `n = j` is the cylinder theorem (C3): the first `j` carry bits of `z` are the
unique bits with `Q_z = 0 (mod 2^j)`. (2) Reducing `Q_z(w)` modulo `2^j` involves only the first `j` bits of `w`,
with odd weights `3^{n-1-i}`; by (1) the congruence modulo `2^j` holds iff the first `j` bits of `w` equal those of
`kappa(z)`. So `v_2(Q_z(w)) >= j` iff `w` and `kappa(z)` agree on the first `j` positions. The leaf law is a
corollary: an expansion at `p` changes nothing before position `p`, writes a `1` at `p`, and shifts the rest. (3) The
frontier word after reading `t` symbols is the prefix of `kappa(7)` followed by `q` leaves; it is a complete tree
word iff `q = 0`. A certificate of any shape has, by (2), its prefix of length `v_2(Q)` equal to a prefix of
`kappa(7)`, and a certificate needs `v_2(Q) >= |w|`, so it is a prefix of `kappa(7)` that is itself a complete tree
word. (4) By (C3) the cylinder of a carry prefix of length `j` is one residue class modulo `2^j`, so under Haar
measure the carry bits are independent fair bits. The walk from `1` hits `0` with probability `p` solving
`p = 1/2 + p^3/2`, `p = 1/phi`; an odd source has first bit `1` and credit `3`, giving `p^3`; the blank source's forced
`1110` leaves credit `6`, giving `p^6`. `QED`

Controls (script, 30 checks): (1) for seven sources and all `n <= 40`; (2) on all `9525` tree words with at most
seven internal nodes for seven sources, including the halting sources `5` and `1`, and on `20000` random binary words;
(L1) and the rank law on every leaf expansion of words of length `<= 16`; the frontier table of `7` for sixty symbols;
200 halting sources whose carry prefixes are verified to be complete tree words violating the inequality; the
blank sources `7` and `27` reach credit above `10^5` after `2·10^5` symbols with minimum credit `3`.

**What this does and does not do.** It shows that "arbitrarily changing tree shapes" cannot create an arithmetic
repair: the rank is a function of the first mismatch with the fixed source's carry word, so the universal
finite-word inequality is literally the orbit statement, and the only free parameter is the source. It does not
bound `kappa(7)`. The share `sqrt(5) - 2` is a statement about random sources; `7` is one source, and its first
`2·10^5` symbols keep the credit far from `0`.

## 3. Theorem C and the census: untouchable numbers

`s(n) = sigma(n) - n`; `m >= 2` is untouchable if `m != s(n)` for every `n >= 2`.

**Theorem C (classical, assembled with its full case list).**
1. `s(n)` is odd iff `n` is an odd non-square or `n` is an even square or twice a square. The odd values from
   semiprimes are `s(pq) = 1 + p + q` (`p < q` odd primes); `s(p^2) = p + 1` is even for odd `p`; `s(4) = 3`,
   `s(8) = 7`, `s(p) = 1`.
2. If `m >= 9` is odd and `m - 1 = p + q` with `p != q` primes, then `m = s(pq)` is touchable. Hence a second odd
   untouchable number `m` requires `m - 1` to be an even number with no partition into two distinct primes.
3. Every even `e` in `[8, 10^6]` has such a partition (check U5), and for `8 <= e <= 4·10^18` the verified minimal
   Goldbach partition uses the prime `p(e) <= 9781 < e/2` once `e > 19562`, so the partition is distinct (CITED,
   Oliveira e Silva–Herzog–Pardi; record `9781`). Therefore **any second odd untouchable number exceeds `4·10^18`**,
   and exhibiting one is equivalent to exhibiting a distinct-Goldbach exception that is also missed by every
   non-semiprime source of (1).

*Proof.* (1) `sigma(n)` is odd iff `n` is a square or twice a square; `s(n) = sigma(n) - n` is odd iff exactly one
of `sigma(n)`, `n` is odd. (2) `pq` has proper divisors `1, p, q`. (3) is the cited computation plus U5. `QED`

**Census to `10^6` (FINITE-EXACT, proved complete).** Write `n = p·k` with `p = spf(n)`, `k = n/p`. For composite `n`,
`1` and `k` are proper divisors, so `s(n) >= 1 + k` and `s(n) <= N` forces `k <= N - 1`. The script examines every
such `n`: `p <= 10` (direct `sigma` sieve to `10^7`), `k = 1` (primes, `s = 1`), `k` prime (semiprimes and prime
squares, through a Goldbach convolution of the prime indicator and `m = p + 1`), and `k` composite with `spf(k) >= p > 10`
(`2,509,482` explicit pairs). Results:

| quantity, `m <= 10^6` | value |
|---|---:|
| untouchable numbers | `150,232` (A057978; density `0.15023`) |
| odd untouchable numbers | `{5}` |
| odd `m in [3, 10^6]` that are touchable | `499,998` |
| of these, with a source that is not a semiprime or prime square | `443,731` |
| **Goldbach-dependent odd `m`** (every source is `pq` or `p^2`) | `56,267` (`11.25%` of odd touchables) |
| even untouchables `0 / 2 mod 4` | `83,603 / 66,628` |
| untouchables `0 / 1 / 2 mod 3` | `65,349 / 41,873 / 43,010` |
| untouchables `0,1,2,3,4 mod 5` | `27,675 / 35,525 / 31,933 / 33,931 / 21,168` |
| independent brute force to `3000` (all `n <= 2999^2`) | agrees, `318` untouchables |

The first Goldbach-dependent odd numbers are `3, 9, 11, 17, 19, 23, 25, 27, 29, 35, ...` (3 only through the prime
square 4); the largest below `10^6` is `999,999`. So about one odd number in nine below `10^6` owes its
touchability entirely to Goldbach; the other eight in nine are also hit by `n` with three or more prime factors or by
`2k^2`-type sources. The residue tables are shares in a finite range, typed FINITE-EXACT; no asymptotic residue law
is claimed. Two features worth a later hostile probe: untouchables are over-represented at `0 mod 3` and
under-represented at `4 mod 5`.

**The exact share as `N -> infinity` (CITED status).** The lower density is positive (Erdős 1973), at least `0.0324`
(te Riele 1976), at least `1/48` (Banks–Luca), at least `0.06` (Chen–Zhao 2011). The density is not known to exist;
Pollack–Pomerance (2016) give a heuristic value near `0.17`, and the counts of A057978 reach `171,128,671,374` at
`10^12`, density `0.1711`. Our `10^6` value `0.1502` is the A057978 term. Nothing here improves these.

**The Chowla map.** For `s'(n) = sigma(n) - n - 1` (the owner's `F = U + S` carrier), `m` is `s'`-touchable iff `m+1`
is `s`-touchable; the `s'`-untouchables are `1, 4, 51, 87, 95, 119, ...`, the `s`-untouchables shifted down by one. The
Goldbach burden moves from the odd values to the even ones (`s'(pq) = p + q`), and the odd `s'`-untouchables inherit the
positive density of the even `s`-untouchables. This is definitional, recorded so that the two threads are not treated
as different problems.

**Creative tricks, with verdicts.**
- *Search for a second odd untouchable.* Pointless below `4·10^18` by Theorem C.3; above it, a candidate would first
  be a distinct-Goldbach exception, and Goldbach exceptions are expected not to exist. Verdict: the question is
  Goldbach-complete in practice; the only exact content is the reduction and the `11.25%` dependence statistic.
- *Rule it out.* Would need every even `e >= 8` to be a sum of two distinct primes, or coverage of the exceptions by
  non-semiprime sources. The second route cannot be uniform: the Goldbach-dependent odd numbers have positive
  share in the range computed. Verdict: out of reach; recorded as a stopping reason.
- *Exact share.* Contrast with Theorem B.4: the Antihydra halting share is an exact golden number because carry
  cylinders are residue classes (a measure-preserving coding), while the aliquot preimage structure has no such
  coding, which is exactly why Pollack–Pomerance's `0.17` is a heuristic and the density's existence is open.

## 4. Frobenius defect and clone mergers: why they do not enter

- The session-65 Frobenius defect concerns a vertex permutation's unordered-pair orbits and their endpoint-return
  signs; here there is no permutation acting on pairs. The carry word `kappa(z)` is a 2-adic coding, not an orbit
  structure. Recorded as not applicable.
- The clone-merger test (deficiency-three ownership) identifies isomorphic presentations of a tournament class by
  merging light gaps. Tree words of equal length and shape class cannot be merged: the packet's rotation law (C6)
  shows that even a rotation changes the source congruence, and Theorem B.2 shows the rank is pinned to the literal
  first mismatch. Recorded as not applicable, with (C6) as the hostile.

## 5. Open and next

- The two over/under-represented residue classes of untouchables (`0 mod 3`, `4 mod 5`) deserve a Pollack–Pomerance
  style heuristic computed per class and a `10^7` census (the script runs in 27 s at `10^6` with `L = 10^7`; a `10^7`
  run needs a `10^8` sieve or a segmented one).
- Theorem A's shift `t` is a function of `ord_m(2)`; the moduli where `-1` is a power of `2` are exactly those with
  no prime factor `p` having `ord_p(2)` odd and with compatible 2-parts of the orders. A clean classification by the
  2-adic valuations of `ord_p(2)` over `p | m` is a one-paragraph exercise not written here.
- Nothing in this note moves the root-27 inequality; the owner's "source-preserving rank" is now an exact object
  (`rho_z(w)`, Theorem B.2) with a one-line law, and the only lever left is the source.

## 6. Reproduction

```bash
python3 04-computation/experiments/sheet_tables_mod_odd_kps_20261025.py
python3 04-computation/experiments/antihydra_tree_rank_kps_20261025.py
python3 04-computation/experiments/untouchable_census_goldbach_dependence_kps_20261025.py 1000000
```

The census needs numpy and about 30 s and 1 GB; the other two are standard library and run in seconds.
