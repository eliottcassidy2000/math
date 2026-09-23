# Loops through 1 of every length, the record structure of log_2 3, and the escapes from 1, 1/2 and the dyadic threads

**Status. PROVED (hand proofs below, each also checked by code where marked):**
- the loop equation and the sharp lower bound. Every E-loop through 1 with `s>=2` multiplications has `2^(K+1)>3^(s+1)`, so the 1-escape at depth `k` needs exactly `K=K0(k-1)=floor(k log_2 3)`.
- optimality of the 1-escape. Within `k` moves it is the only way `1+3^k u` can descend.
- the cost formula, the splicing lemma, and the reduction of HYP-9122 to one base loop and one cycle per record denominator of `log_2 3`.
- the height lower bound. Hence splicing a finite set of loops and cycles gives ratio `<3` for only finitely many lengths. Along any arithmetic progression of lengths `H(s)/s` is unbounded, and along the convergent records `H(q-1) >= q^2/(3 ln 2)*(1-o(1))`.
- the 1/2-transfer lemma and its sharpness (price exactly `2c(k-1)/3`, which lies in `(1,2)`), and the two-move exits.
- `1/2, 43/32, 59/64, 145/128, 209/256, 371/256, 499/512` all lie in `Bad_inf`.
- the explicit escape routes of `43/32` and `59/64` with their descent inequalities.

**FINITE-EXACT:**
- `minK(s)=K0(s)` for every `2<=s<=6000` (DP, values `<=3*10^7`), with the exact minimal heights `H(s)`.
- explicit loops for every `2<=s<=2500`, each checked by an independent exact verifier. They are generated from 13 base loops and 11 hub cycles.
- hub cycles realising all record pairs `(q, ceil(q log_2 3))` with `q<=306`.
- the level-12 certificate landscape (see-saw numbers).

**OPEN:**
- HYP-9122 for all `s`. No explicit infinite family was found. The precise failure is in section 6.
- the 1/2 escape beyond one certificate.
- closure of the see-saw.

Collatz OPEN. No independent-agent audit yet.

Session `collatz-procgen-20260922`, research sub-lane of lane one ([choice_ladder](collatz_procgen_20260922_choice_ladder.md), sections 1, 4, 5). Scripts: `04-computation/experiments/collatz_procgen_20260922_loops_*`. Output: [collatz_procgen_20260922_loops.out](collatz_procgen_20260922_loops.out). Reproduction is in section 10.

## 0. Setting

A **reverse move** is `x -> (2^k x - 1)/3` with `k>=0`. It is legal iff the result is a positive integer (for the 3-adic points: a 3-adic unit) not divisible by 3.

A **loop of length `s`** is a legal reverse path `1=x_0 -> x_1 -> ... -> x_s=1`. Read forward, it is an E-cycle through 1 with `s` arrows `n->3n+1` and `K=sum k_i` halvings. The *multiplication points* are `x_1,...,x_s` (with `x_s=1`).

Notation used throughout:
- `L = log_2 3`.
- `K0(s) = floor((s+1)L)`, the largest `K` with `2^K<3^(s+1)`.
- `c(s) = 2^{K0(s)}/3^s`.
- `eps(n) = ceil(nL) - nL`, which lies in `(0,1)`.

Then `c(s) = (3/2)*2^{eps(s+1)}`, which lies in `(3/2,3)`.

A **record** is an `n` with `eps(n) < eps(n')` for all `n'<n`. These are the upper best-approximation denominators of `L`:

`1, 3, 5, 17, 29, 41, 94, 147, 200, 253, 306, 971, 1636, 2301, ..., 5626, 6291, ..., 14936, 15601, ...`

The records from 306 to 14936 are `306 + 665j`, and `15601` is the next convergent.

## 1. The loop equation and the sharp lower bound (PROVED)

**Proposition 1.1 (loop equation).** Let `k_1,...,k_s>=0` and `K_i = k_1+...+k_i`. These are the moves of a legal loop through 1 iff

`2^K = 3^s + sum_{i=1}^{s} 3^{i-1} 2^{K-K_i}`.

Equivalently, `2^K = sum_{i=0}^{s} 3^i 2^{e_i}` with `e_i = K - K_{i+1}` nonincreasing and `e_{s-1}=e_s=0`. In particular legality is automatic: any word in `{n->3n+1, n->n/2}` whose affine map fixes 1 is a legal E-cycle.

*Proof.* The recursion gives `3^i x_i = 2^{K_i} - B_i` with `B_i = sum_{j<=i} 3^{j-1}2^{K_i-K_j}`. So every `x_i` lies in `Z[1/3]`, and `x_s=1` iff the identity holds.

Conversely, assume the identity. Running the forward maps `y -> (3y+1)/2^{k_i}` from `x_s=1` shows that every `x_i` also lies in `Z[1/2]`. Hence every `x_i` is an integer. It is positive (the forward maps preserve positivity) and satisfies `x_{i-1} = 2^{-k_i}(3x_i+1)`, which is `≢ 0 (mod 3)`. The intermediate halvings act on `2^{k_i-j}x_{i-1}`, which is an integer, hence even whenever `j<k_i`. ∎

**Proposition 1.2 (lower bound).** Every loop with `s>=2` has `2^{K+1} > 3^{s+1}`, that is `K>=K0(s)`.

Hence the 1-escape lemma at depth `k` (which needs `2^K<3^k`) is available **iff a loop of length `k-1` with exactly `K0(k-1)` halvings exists**, and HYP-9122 is equivalent to `minK(s)=K0(s)` for all `s>=2`. The minimal ratio is then exactly `c(s)`.

*Proof.* We have `2^K = sum 3^i 2^{e_i} >= sum_{i<=s} 3^i = (3^{s+1}-1)/2`, with equality iff every `e_i=0`, that is iff `3^{s+1}-1=2^{K+1}`.

Now `3^n-1` is a power of 2 only for `n=1,2`:
- if `n` is odd, `3^n-1 ≡ 2 (mod 4)`;
- if `n=2n'`, then `3^{n'}-1` and `3^{n'}+1` are powers of 2 differing by 2, so `n'=1`.

So for `s>=2` we get `2^{K+1} > 3^{s+1}-1`, hence `2^{K+1} > 3^{s+1}` (the two sides have different parity). ∎

This sharpens lane one's bound `13/9` to `c>3/2`. For `s>=2` the only ratio below 3 that a loop can have is `c(s)`, the unique value of `2^K/3^s` in `(3/2,3)`.

**Proposition 1.3 (cost formula).** `2^K/3^s = prod_{v}(1+1/(3v))` over the `s` multiplication points.

In particular, the initial climb `1->4->13->...->y_n = (3^{n+1}-1)/2` costs exactly `(3/2)(1-3^{-(n+1)})`. *Proof:* `x_{i-1} = x_i*(3/2^{k_i})(1+1/(3x_i))`; multiply around the loop. ∎

## 2. The 1-escape is optimal; no other exit shapes within `k` moves (PROVED)

**Proposition 2.1.** Let `m=1+3^k u` with `k>=2`, `u>=1`, `3∤u`.
- (i) Every legal reverse path from `m` of length `i<=k-1` is the lift of a legal path from 1: `x_i(m) = x_i(1) + 2^{K_i}3^{k-i}u`, with `x_i(1)>=1`. All its values exceed `m`.
- (ii) A legal path of length `k` ends below `m` **iff** it is a loop of length `k-1` with `2^K<3^k`, followed by the move 0. It then ends at `2^K u`.

*Proof.* `x_j(1) = x_j(m) - 2^{K_j}3^{k-j}u` is an integer. For `j<=k-1` it is `≡ x_j(m) ≢ 0 (mod 3)`. By induction it is `>=1`, because `x_j(1)=0` would need `x_{j-1}(1)=1` and `k_j=0`, giving `0 ≡ x_j(m) (mod 3)`.

Hence `2^{K_i} = 3^i x_i(1) + B_i > 3^i`, and `x_i(m) > 1 + 3^k u = m`.

At `j=k` put `y = x_k(1) >= 0`:
- If `y>=1`, the same bound gives `x_k(m) > m`.
- If `y=0`, then `x_{k-1}(1)=1` and `k_k=0`, so the first `k-1` moves form a loop. Its endpoint is `x_k(m) = 2^K u`, and `2^K u < 1+3^k u` iff `2^K<3^k`. ∎

**Consequences.**
- "Robust alternatives" of other shapes do not exist within `k` moves. For example, a loop of length `k-2` followed by a two-step exit either appends the trivial move `k=2` (so it *is* a loop of length `k-1`) or ends above `m`.
- Beyond `k` moves, every exit depends on the residue of `u`. As an example, suppose the minimal loop were missing and only a loop with `K0+1` halvings existed (ratio `2c`). Then `y=2^{K0+1}u` descends with one more move in two cases: `y ≡ 4,7 (mod 9)`, which always works; and `y ≡ 2,8 (mod 9)`, which works iff `c<9/4` (for large `u`). So two of the six unit residue classes mod 9 are covered always, and four of the six when `c<9/4`.
- All remaining exits belong to the see-saw (section 8).

## 3. The loop table (FINITE-EXACT)

**DP.** `loops_height` is a layered DP over values `<=3*10^7` with lexicographic objective (min `K`, then min height). It is valid because a `K0`-loop is `K`-minimal at every intermediate state (Proposition 1.2).

It gives **`minK(s)=K0(s)` for every `2<=s<=6000`**, together with the exact minimal height `H(s)` (the least possible maximum multiplication point). Cross-checks:
- the lane-one program, run with `s<=80` and values `<=2*10^7`, gives the same `minK`;
- the reconstruction DP `loops_dp` agrees for `s<=80`.

Peak memory is 361 MB and the run takes 8 minutes.

**Explicit loops.** For every `s<=80` a minimal-height loop is printed in the `.out`, as a forward word (`M^a` = `a` consecutive `3n+1`, `H^b` = `b` halvings) and as its reverse moves. Examples:
- `s=4`: `M3 H3 M H4` (`1->4->13->40->20->10->5->16->...->1`)
- `s=16`: `M3 H M2 H M2 H4 M2 H2 M2 H M H4 M2 H3 M H2 M H8`, ending `...->28->85->256->...->1`
- `s=40` (height 3949): climbs through the multiplication points `1,4,...,1093` to 3280, wanders with multiplication points in `[227,3949]`, and ends `341=(4^5-1)/3 -> 1024 -> 1`

**Table for `s<=80`.** Columns: `K0(s)`, `c(s)`, minimal height `H(s)`, and how the loop is built in the splicing family of section 4. `Cq+L(t)` means "splice the cycle `C_q` into loop `t`"; `base` means a primitive base loop. For `s=1` the minimal loop is the trivial one, with `K=2<K0(1)`.

| s | K0 | c(s) | H(s) | family | s | K0 | c(s) | H(s) | family |
|---|---|---|---|---|---|---|---|---|---|
| 1 | 2 | 1.3333 | 1 | trivial | 41 | 66 | 2.0231 | 209 | C5+L(36) |
| 2 | 4 | 1.7778 | 1 | base | 42 | 68 | 2.6974 | 61 | C5+L(37) |
| 3 | 6 | 2.3704 | 1 | C1+L(2) | 43 | 69 | 1.7983 | 317 | C5+L(38) |
| 4 | 7 | 1.5802 | 13 | base | 44 | 71 | 2.3977 | 61 | C5+L(39) |
| 5 | 9 | 2.1070 | 13 | C1+L(4) | 45 | 72 | 1.5985 | 821 | C5+L(40) |
| 6 | 11 | 2.8093 | 13 | C1+L(5) | 46 | 74 | 2.1313 | 209 | C5+L(41) |
| 7 | 12 | 1.8729 | 13 | C3+L(4) | 47 | 76 | 2.8417 | 61 | C5+L(42) |
| 8 | 14 | 2.4972 | 13 | C3+L(5) | 48 | 77 | 1.8945 | 317 | C5+L(43) |
| 9 | 15 | 1.6648 | 61 | C5+L(4) | 49 | 79 | 2.5260 | 61 | C5+L(44) |
| 10 | 17 | 2.2197 | 13 | C5+L(5) | 50 | 80 | 1.6840 | 533 | C5+L(45) |
| 11 | 19 | 2.9596 | 13 | C5+L(6) | 51 | 82 | 2.2453 | 209 | C5+L(46) |
| 12 | 20 | 1.9731 | 29 | C5+L(7) | 52 | 84 | 2.9937 | 61 | C5+L(47) |
| 13 | 22 | 2.6308 | 13 | C5+L(8) | 53 | 85 | 1.9958 | 317 | C5+L(48) |
| 14 | 23 | 1.7538 | 61 | C5+L(9) | 54 | 87 | 2.6611 | 61 | C5+L(49) |
| 15 | 25 | 2.3385 | 29 | C5+L(10) | 55 | 88 | 1.7741 | 497 | C5+L(50) |
| 16 | 26 | 1.5590 | 533 | base | 56 | 90 | 2.3654 | 209 | C5+L(51) |
| 17 | 28 | 2.0786 | 61 | C5+L(12) | 57 | 91 | 1.5770 | 1397 | C17+L(40) |
| 18 | 30 | 2.7715 | 29 | C5+L(13) | 58 | 93 | 2.1026 | 317 | C5+L(53) |
| 19 | 31 | 1.8477 | 61 | C5+L(14) | 59 | 95 | 2.8035 | 61 | C5+L(54) |
| 20 | 33 | 2.4636 | 29 | C5+L(15) | 60 | 96 | 1.8690 | 497 | C5+L(55) |
| 21 | 34 | 1.6424 | 209 | C5+L(16) | 61 | 98 | 2.4920 | 209 | C5+L(56) |
| 22 | 36 | 2.1898 | 61 | C5+L(17) | 62 | 99 | 1.6613 | 821 | C17+L(45) |
| 23 | 38 | 2.9198 | 29 | C5+L(18) | 63 | 101 | 2.2151 | 317 | C5+L(58) |
| 24 | 39 | 1.9465 | 61 | C5+L(19) | 64 | 103 | 2.9534 | 61 | C5+L(59) |
| 25 | 41 | 2.5954 | 61 | C5+L(20) | 65 | 104 | 1.9690 | 317 | C5+L(60) |
| 26 | 42 | 1.7302 | 209 | C5+L(21) | 66 | 106 | 2.6253 | 209 | C5+L(61) |
| 27 | 44 | 2.3070 | 61 | C5+L(22) | 67 | 107 | 1.7502 | 533 | C17+L(50) |
| 28 | 45 | 1.5380 | 949 | base | 68 | 109 | 2.3336 | 317 | C5+L(63) |
| 29 | 47 | 2.0507 | 61 | C5+L(24) | 69 | 110 | 1.5557 | 2081 | C29+L(40) |
| 30 | 49 | 2.7342 | 61 | C5+L(25) | 70 | 112 | 2.0743 | 317 | C5+L(65) |
| 31 | 50 | 1.8228 | 209 | C5+L(26) | 71 | 114 | 2.7657 | 209 | C5+L(66) |
| 32 | 52 | 2.4304 | 61 | C5+L(27) | 72 | 115 | 1.8438 | 497 | C17+L(55) |
| 33 | 53 | 1.6203 | 533 | C5+L(28) | 73 | 117 | 2.4584 | 317 | C5+L(68) |
| 34 | 55 | 2.1604 | 61 | C5+L(29) | 74 | 118 | 1.6390 | 949 | C29+L(45) |
| 35 | 57 | 2.8805 | 61 | C5+L(30) | 75 | 120 | 2.1853 | 317 | C5+L(70) |
| 36 | 58 | 1.9203 | 209 | C5+L(31) | 76 | 122 | 2.9137 | 209 | C5+L(71) |
| 37 | 60 | 2.5604 | 61 | C5+L(32) | 77 | 123 | 1.9425 | 497 | C17+L(60) |
| 38 | 61 | 1.7070 | 497 | C5+L(33) | 78 | 125 | 2.5900 | 317 | C5+L(73) |
| 39 | 63 | 2.2759 | 61 | C5+L(34) | 79 | 126 | 1.7266 | 821 | C29+L(50) |
| 40 | 64 | 1.5173 | 3949 | base | 80 | 128 | 2.3022 | 317 | C5+L(75) |

**Hard cases.** These are the lengths just below a record, `s=q-1`. The column "proved LB" is the lower bound on `H(s)` from Proposition 5.1. All rows are FINITE-EXACT.

| q (record) | s=q-1 | K0(s) | c(s) | eps(q) | H(s) | proved LB |
|---|---|---|---|---|---|---|
| 41 | 40 | 64 | 1.517293 | 0.016537 | 3,949 | 896 |
| 94 | 93 | 148 | 1.514128 | 0.013525 | 11,501 | 2,946 |
| 147 | 146 | 232 | 1.510970 | 0.010512 | 24,001 | 6,184 |
| 200 | 199 | 316 | 1.507818 | 0.007500 | 55,109 | 12,018 |
| 253 | 252 | 400 | 1.504673 | 0.004487 | 107,885 | 25,701 |
| 306 | 305 | 484 | 1.501534 | 0.001475 | 376,657 | 95,112 |
| 971 | 970 | 1538 | 1.501469 | 0.001412 | 1,257,905 | 325,498 |
| 1636 | 1635 | 2592 | 1.501403 | 0.001349 | 2,309,861 | 577,594 |
| 2301 | 2300 | 3646 | 1.501338 | 0.001286 | 3,384,637 | 854,461 |
| 5626 | 5625 | 8916 | 1.501010 | 0.000971 | 10,735,213 | 2,777,900 |

The ratio `H(s)*eps(s+1)*3 ln2/s` (which is `>=1` asymptotically by Proposition 5.1) equals 2.5 to 4.3 at the records. It is 3.40 at `s=40`, 3.79 at `s=305`, and between 3.85 and 3.99 at every record `s>=1635`. The median height grows linearly: about 2,000 for `s<=1000` and about 21,000 for `5000<s<=6000`.

The earlier "thin margin" worry is resolved: `c(s)` comes arbitrarily close to 3 (for example `s=52`, `c=2.99374`), but `K0(s)` is always the right `K` and is realised. The hard cases are the other end, `c(s)` close to `3/2`.

## 4. Structure: splicing, records, and the continued fraction of log_2 3

**Lemma 4.1 (splicing, PROVED).** Suppose a loop of length `s` visits `h`, and `C` is a cycle through `h` with `q` multiplications and `p` halvings. Splicing `C` in at `h` gives a loop of length `s+q` with `K+p` halvings, and its ratio is multiplied by `2^p/3^q`.

Now let the loop have `K0(s)` halvings and let `p=ceil(qL)`, so that the ratio of `C` is `2^{eps(q)}`. Then the spliced loop realises `K0(s+q)` iff `eps(s+1)+eps(q)<1`, and in that case `eps(s+q+1)=eps(s+1)+eps(q)`.

*Proof.* `K0(s)+p = (s+q+1)L + eps(s+1) + eps(q) - 1`. This equals `floor((s+q+1)L)` iff the fractional correction lies in `(-1,0)`. ∎

**Hub cycles (FINITE-EXACT, verified exactly).** Every loop through 1 visits 4. Every loop with `s>=2` and ratio `<=128/81` climbs at least two steps and so visits 13. Proof: a climb of length 1 forces either the multiplication points `1,2,7`, with cost `>= 4/3*7/6*22/21 = 1.63`, or a return to 1, with cost `>= (4/3)^2`, since every loop costs at least `4/3`. The cycles found:
- `C1`: trivial, at 1, `(1,2)`.
- `C3`: at 4, `4->13->40->20->10->5->16->8->4`, `(3,5)`.
- `C5`: at 13, `13->40->20->61->184->92->277->832->...->13`, `(5,8)`.

Run the minimal-cycle search `loops_hubs` through the climb points (`q<=120` or `306`, values `<=3-4*10^6`). Through `13, 40, 121, 1093`, as `q` increases, each new minimum of `eps` after the shortest available cycle is a record pair `(q, ceil(qL))` of `log_2 3`, provided the hub is high enough for that record's ratio `2^{eps(q)}`. The shortest available cycle is `q=3` at 13 and `q=10` at the higher hubs. Hub 4 behaves differently:
- Through 4, the successive minima are `(3,5), (15,24), (27,43), (39,62), (92,146)`. These are exactly `(q-2, ceil(qL)-3)` for the records `q = 5, 17, 29, 41, 94`, with ratio `(9/8)*2^{eps(q)}`. So they are the record base loops `B_q` viewed as closed walks at 4 (see section 6). The one exception is `(3,5)`, which is both `B_5` and `C3`.
- Through 13: `(3,5)`, `(5,8)`, `(17,27)`.
- Through 40: `(17,27)`, `(29,46)`. The record 41 is unavailable there, and `(82,130) = 2*(41,65)` appears instead.
- Through 121: `(17,27)`, `(29,46)`, `(41,65)`, `(94,149)`.
- Through 1093: every record `q=17,29,41,94,147,200,253,306` is realised, with ratio `2^{eps(q)}` down to `1.001023`. Through 364 the record `q=306` needs `p=486`, so that hub is too low.

**Theorem 4.2 (reduction to records, PROVED).** Suppose that:
- for every record `q>=3` there is a loop `B_q` of length `q-1` with `K0(q-1)` halvings;
- for every record `q` there is a cycle `C_q` with `(q, ceil(qL))` through a point `h_q`;
- each `B_q` visits `h_{q'}` for every record `q'<=q`.

Then every `s>=2` has a loop with `K0(s)` halvings. In other words, HYP-9122 reduces to the records, which are a sparse set.

*Proof.* Let `q*` be the largest record `<=s+1`. If `s+1=q*`, take `B_{q*}`.

Otherwise put `n=s+1-q*`. Since `s+1` is not a record and `q*` realises `min_{n'<s+1} eps(n')`, we have `eps(s+1)>eps(q*)`. Hence `eps(n)=eps(s+1)-eps(q*)` (additivity mod 1).

Greedy decomposition: if `n` is not a record, subtract the largest record `q'<n`; again `eps(n-q')=eps(n)-eps(q')`. This writes `n=q_2+...+q_t` with records `q_j<=q*` and `sum eps(q_j)=eps(n)`, hence `sum ceil(q_j L)=ceil(nL)`.

Splice `C_{q_2},...,C_{q_t}` into `B_{q*}` at their hubs. The length is `s`. The halvings number `floor(q*L)+ceil(nL) = floor((s+1)L)`. ∎

**FINITE-EXACT instance.** With 13 base loops at `s = 2, 4, 16, 28, 40, 93, 146, 199, 252, 305, 970, 1635, 2300` (that is, `s=q-1` for the records `q<=2301`) and the 11 cycles `C1, C3, C5, C17, ..., C306`, the splicing closure produces loops for **every `2<=s<=2500`** (`loops_family.py`, dynamic hub checks). All 2,499 are verified independently (`loops_verify.py`: reverse legality, forward E-simulation, the identity of Proposition 1.1, and `K=K0(s)`).

With only `C1, C3, C5` the base set is exactly `{2} ∪ {s : c(s)<=128/81}`: 1 + 189 = 190 lengths up to 2500, all realised and verified.

**Answer to "recursive pattern?"** Yes. The lengths `s` for which new loops are needed are exactly `s=q-1` for the record (upper best-approximation) denominators of `log_2 3`. Every other length arises by splicing the record cycles, in the Ostrowski-type greedy order of Theorem 4.2. The increments `K0(s+1)-K0(s)`, which form the Beatty word:
- `=2` is splicing `C1`;
- `=1` requires the other record cycles.

## 5. Heights, and why no cheap infinite family exists (PROVED)

**Proposition 5.1.** A loop of length `s>=2` whose multiplication points are all `<=H` has

`H >= [(s - 1 - log_3(2H+1))/(eps(s+1) ln 2) - 1]/3`.

*Proof.* The climb of length `n` costs `(3/2)(1-3^{-(n+1)})`. The first multiplication point after the climb is `<=y_n/2`, and its factor is at least `1+(4/3)3^{-(n+1)}`, which compensates the deficit. The remaining `s-n-1` points each cost `>=1+1/(3H)`.

So `2^{eps} >= (1+1/(3H))^{s-n-1}`, where `eps=eps(s+1)`. Use `ln(1+x) >= x/(1+x)`, and `y_{n-1}<=H`, which gives `n<=log_3(2H+1)`. ∎

**Consequences.**
- `H(s) >= (s-O(log s))/(3 ln 2)` for every `s`. So a fixed finite set of loops and cycles yields ratio `<3` for only finitely many lengths, since each splice multiplies the ratio by at least the least cycle ratio.
- Along an arithmetic progression `s=a+bj`, the values `eps(s+1)` are equidistributed in `(0,1)` (Weyl), so they come arbitrarily close to 0 and **`H(s)/s` is unbounded** along any AP family. Such a family must contain loops with ratio arbitrarily close to 3/2, which are the hard cases of the all-`s` statement.
- At the upper convergents `q` of `L` we have `eps(q) < 1/q_next <= 1/q`, so `H(q-1) >= q^2/(3 ln 2)*(1-o(1))` along the records.
- The data show the bound is sharp up to a factor of about 4 (section 3).

## 6. The explicit infinite family: not found (OPEN), and where it fails

Theorem 4.2 reduces HYP-9122 to one base loop `B_q` and one cycle `C_q` per record `q`.

**The record base loops cannot come from anything smaller by splicing (PROVED).** Suppose a loop of length `q-1` with `K0(q-1)` halvings were a loop `L'` of length `s'>=2` with a positive cycle `(q',p')` spliced in. Proposition 1.2 gives `ratio(L') >= c(s')`, so `c(q-1) >= c(s')*2^{p'}/3^{q'}` and hence `eps(q) >= eps(s'+1) + (p' - q'L) > eps(q-q')`. This contradicts the record property.

Only the host of length 1 remains, and it always works in a trivial way. Every loop of length `s>=2` starts `1->4` and ends `...->4->2->1` (the last multiplication lands on `4^j`, `j>=1`). So the loop is the trivial loop with a closed walk at 4 of length `s-1` spliced in, and its ratio is `(4/3)` times that of the walk. Hence `B_q` is the same object as a closed walk at 4 with `(q-2, ceil(qL)-3)`, whose ratio is `(9/8)*2^{eps(q)}`. For example `B_5` is the trivial loop with `C3` spliced in.

So the record problem is equivalent to finding closed walks at 4 whose ratio exceeds `9/8` by the factor `2^{eps(q)}`. By Proposition 5.1 these must climb to height at least about `q/(3 eps(q) ln 2)`.

What a proof would need is an operation that *raises* a loop by an order of magnitude while changing `(s,K)` by a lower-convergent pair, such as `(665,1054)`, whose ratio is below 1. Positive cycles cannot do this.

Attempted mechanisms, all of which failed:
- **Pumping a cycle (PROVED to fail).** Write `f_Q(v) = q* + rho(v-q*)` with `rho = 3^{s_Q}/2^{K_Q}`, which is not 1. If `P Q^r R` fixes 1 for two values of `r`, then `rho^{r_1}(u-q*) = rho^{r_2}(u-q*)` with `u = f_P(1)`, so `u = q*`. Then `Q` is a cycle at `u` and the family is plain insertion, whose ratio tends to infinity.
- **The 2-adic highways.** These are the lifts of the cycles at `-1` and `+1`: `2^a*alpha-1 -> 3^a*alpha-1` (ratio `(3/2)^a`) and `4^a*beta+1 -> 3^a*beta+1` (ratio `(3/4)^a`). Chaining them to 1 requires exact coincidences between the climb `(3^{n+1}-1)/2` and numbers like `R_j=(4^j-1)/3`. These are S-unit equations. For a bounded number of blocks the loop equation has only finitely many non-degenerate solutions (the S-unit theorem, CITED from general knowledge, not re-derived here), and no degenerate infinite family was found.
- **Descending from `R_j` along the 3-adic digits of `j`.** This exhausts the digits after about `log j` steps, while about `1.26j` steps are needed.

The hard base loops found (for example `s=40`, `93`, `305`, `970`) share one shape:
- a climb of 7 to 11 steps;
- a wander at heights about `s/eps`;
- a landing on a large `R_j` (`341 = R_5`, `5461 = R_7`), then a fall down `4^j`.

## 7. Escapes (part b)

### 7.1 The point 1 (recap and sharpening)

`m=1+3^k u` descends in exactly `k` moves to `2^{K0(k-1)} u`, with factor `c(k-1)/3`. This holds iff the loop of length `k-1` with `K0(k-1)` exists (Proposition 2.1). The loop exists for all `k<=6001` (DP) and explicitly for `k<=2501`.

At level 12 the certified factor of the classes `1+3^k u` is **exactly** `c(k-1)/3` for every `k<=12` (`.out` section 6).

### 7.2 The point 1/2 (PROVED lemmas, FINITE-EXACT statistics)

The reverse paths from `1/2` are the paths from 1 with the first move increased by one: `(2^{k+1}/2-1)/3 = (2^k-1)/3`. The move `k=1` gives 0 and is illegal.

**Lemma 7.1 (1/2-transfer).** Let `k>=2`, `m=(1+3^k w)/2`, `3∤w`. Let `(k_1,...,k_s)` be a loop of length `s=k-1` with `K` halvings. Then the reverse path from `m` with moves `(k_1+1, k_2, ..., k_s, 0)` is legal and ends at `y=2^K w`.

*Proof.* `x_i(m) = x_i(1/2) + 2^{K_i}3^{k-i}w`, where `x_i(1/2)` is the loop value for `i>=1`. So the path is congruent to the loop mod 3 up to step `s`, and `x_s(m)=1+3*2^K w`. The move 0 then gives `2^K w`, which is prime to 3. Checked exactly on 1,552 random instances.

**Lemma 7.2 (sharpness).** Let `k>=3`. Every legal reverse path from `m` of length `<=k` stays above `m`. The least endpoint at length `k` is `2^{minK(k-1)} w`, which lies in `(m, 2m)` when `minK=K0`.

So the **escape price of the 1/2-thread is exactly `2c(k-1)/3 ∈ (1,2)`**, since `y/m = 2^{K+1}w/(1+3^k w)` tends to that value from below as `w->infinity`. A direct descent would need `c<3/2`, which is impossible for `k>=3`.

*Proof.* As in Proposition 2.1 with `2^{K_i-1} = 3^i x_i(1/2) + B_i`. A virtual endpoint `>=1` gives a value `>3^k w`. A virtual endpoint 0 forces the loop, and then `2^{K+1}>3^k` (Proposition 1.2) gives `2^K w >= m`, with equality excluded. ∎

For `k=2` the trivial loop gives a direct descent to `4w`, factor `8/9`.

The route "1/2 -> 1 by the move 3, then a 1-loop of length `k-2`, then the move 0" is the special case of Lemma 7.1 whose loop starts with the trivial loop. Its price is `(2/3)(4/3)c(k-2) = 8c(k-2)/9`. This is never below `2c(k-1)/3`, and equals it iff `K0(k-1) = K0(k-2)+2`. So routing through 1 first never helps.

**Lemma 7.3 (two-move exits; required loop ratios).** With `y=2^K w`, `K=K0(k-1)` and `c=c(k-1)`:
- `y ≡ 4,7 (mod 9)`: the move 0 gives `(y-1)/3 < m` for every `c<9/2`. So this always works.
- `y ≡ 2,8 (mod 9)`: the move 1 gives `(2y-1)/3`, which is `<m` iff `2^{K+2}<3^{k+1}`, i.e. **`c(k-1)<9/4`**. This holds for a fraction `log_2(3/2) = 0.585` of the `k`.
- `y ≡ 1,5 (mod 9)`: this needs longer continuations.

In general, `m` descends as soon as the class of `y` has a certificate of factor `f<3/(2c)`, since the endpoint is `< f*y < f*(2c/3)m`.

FINITE-EXACT, level 12: `y` runs over all units as `w` does. Among the 1,062,882 unit classes mod `3^13`, the number whose best level-12 certificate fails this threshold (including the 30 exceptional classes) is between **56 and 285 for every `3<=k<=40`**, a fraction at most `2.7*10^{-4}`. In this range the counts repeat with period 12 in `k`. So one certificate after the transfer finishes the escape except on landing classes within `3^{-13}` of the hostile threads and their pre-images (for example `y ≡ 2 (mod 3^13)`, where `2 -> 1` lands on the thread of 1).

### 7.3 Dyadic hostile points (PROVED)

**Criterion.** For a reverse path from `q` ending at `x_s`, the multiplier is `(x_s + B/3^s)/q` with `B/3^s >= (1-3^{-s})/2`.

`loops_escapes.py` enumerates the finite non-integer part of the reverse tree of `q`, with exact multipliers, and checks three things:
- every non-integer node has multiplier `>1`;
- every integer node `>=2` gives `>= (2+1/3)/q > 1`;
- the value 1 appears only at depths `s` with `1+(1-3^{-s})/2 > q`, which follows from where the nodes of the form `4/2^k` sit.

Result: **`1/2, 43/32, 59/64, 145/128, 209/256, 371/256, 499/512` are in `Bad_inf`**. For `371/256>13/9`, the case `s*=3` is needed. The non-integer parts of the reverse trees:
- `43/32`: `{13/4}`;
- `59/64`: `{17/8, 19/2, 5/2, 1/2}`;
- `209/256` and `499/512` also pass through `1/2`.

All of them are also exceptional classes at level 12, and so are lane one's other reconstructions `661/512, 715/512, 1241/1024, 1753/2048, 4235/4096`. So the dyadic threads are pre-images of the threads of 1 and 1/2, at multiplier `>1`.

**Routes and descent inequalities (PROVED, checked on integers).**
- **`43/32`.** Take `32m=43+3^k W` with `3∤W`. The moves `3,2,0` give `1+3^{k-3}W` (the thread of 1 at precision `k-3`, route multiplier `32/27`). The 1-escape then ends at `2^{K0(k-4)}W`. This is `<m` iff `2^{K0(k-4)+5}<3^k`, i.e. **`c(k-4)<81/32`**, which holds for a fraction `log_2(27/16) = 0.755` of the `k` (for example it fails at `k=10`, where `c(6)=2.809`).
- **`59/64`.** Take `64m=59+3^k W`. The moves `3,2,0` give `(1+3^{k-3}W)/2`, the thread of 1/2. So `59/64` is **second order**: price `(32/27)(2c(k-4)/3) = 64c(k-4)/81 ∈ (1.19,2.37)`, then a continuation as in 7.2.

## 8. See-saw assessment (part c; FINITE-EXACT numbers, OPEN conclusion)

**Exact certificates at level 12** (`loops_seesaw.c`: classes mod `3^13`, best multiplier `2^K/3^s` over paths of length `<=12` and their prefixes, stored exactly):
- There are **30 exceptional classes**. They include the classes of `1, 1/2, 43/32, 59/64, 145/128, 209/256, 371/256`. All 29 classes other than 1 are `≡ 5 ≡ 1/2 (mod 9)`: every hostile point except 1 sits in the mod-9 class of 1/2.
- **Worst certified factor overall: `2^19/3^12 = 0.98654`.** It is attained on the thread of 1 at `k=12` (the 1-escape `c(11)/3`). By Proposition 2.1, the class `1+3^r u` has best factor exactly `c(r-1)/3` at level `r`, so the worst certified factor at level `r` is at least `max_{k<=r} c(k-1)/3`, whose `limsup` is 1.
- Only 255 certified classes (`2.4*10^{-4}`) have factor `>=1/2`.
- By distance to the exceptional set (`d` = the largest `j` for which the residue mod `3^j` is shared with an exceptional class), the worst factor is:
  - `d<=1`: `2/3`. These are the classes `2,4,7,8 (mod 9)`; each has a one-move descent at rate at least `ln(3/2) = 0.405` nats per consumed digit.
  - `d<=4`: `8/9`, the 1/2-thread at `k=2`.
  - `d<=7`: `2^11/3^7 = 0.936`.
  - `d>=8`: `0.98654`.
- The best achievable descent rate per digit falls to `0.0589`, `0.0094` and `0.0011` nats in those bands.

**Prices of the hostile threads.**
- The thread of 1 needs no see-saw: it is a direct descent with factor `c(k-1)/3<1`.
- The 1/2-type threads cost more than 1: `2c(k-1)/3` for 1/2, and `64c(k-4)/81` for 59/64. The largest price per consumed digit is `ln(2c(3)/3)/4 = 0.114` nats, attained by 1/2 at `k=4`.

**Conclusion (OPEN).**
- A one-round Applegate–Lagarias see-saw (escape, then one certificate) does not close at level 12, or at any fixed level.
- The reason is that the landing class `2^K w` is arbitrary, so the relevant landing factor is the worst over all certified classes, 0.98654. Then `(2c(k-1)/3)*0.98654 > 1` except when `c(k-1) < 1.5205`, i.e. at the lengths `k-1 = 40, 93, 146, ...` just below records. At a higher level `r` the worst factor is at least `max_{k<=r} c(k-1)/3 >= c(2)/3`, while the price approaches 2.
- It does close except on a fraction `<=2.7*10^{-4}` of landing classes, and those are again hostile-thread neighbourhoods.
- So the induction reduces to chains of hostile transfers. An adversarial `m` can arrange such chains digit by digit, since each 1/2-transfer consumes `k>=3` digits of `w`. With the prices of the known threads, a chain can raise the value by at most about `m^{0.104}` (`0.114` nats per digit over `log_3 m` digits) before the digits of `m` run out. The infinitely many other dyadic threads have not been priced.
- What remains is a point `2^K w_final` with `w_final` small. Its descent is governed by the 3-adic digits of `2^K`: a transversality statement of the Erdős ternary-digits type, not a finite certificate problem. **Refined (wave 3, 2026-09-23):** the chains read the **low** (3-adic) digits of `2^K w`, controlled by the kappa formula `v_3(2^K - r) = 1 + v_3(K - kappa(r))`, not the full or top expansion that Erdős's problem concerns. Every 3-adic-window form of Erdős is FALSE ([transversality foundry](collatz_procgen_20260922_transversality_foundry.md) §3; [Q1 mirror](collatz_procgen_20260922_q1_mirror.md) §5).
- So the see-saw, as such, does not look closable. A proof needs a digit-weighted potential that beats 0.114 nats per digit on hostile stages, plus control of the ternary digits of `2^K` at the end of hostile chains.

## 9. Hypotheses

- **HYP-9122 (sharpened, OPEN).** For every `s>=2` there is a loop with exactly `K0(s)` halvings. This is equivalent (Proposition 1.2) to the original form. By Theorem 4.2 it follows from the **record form**: for every record `q` there are a base loop `B_q` and a cycle `C_q`, compatible at hubs.
  - Evidence: all records `q<=2301` realise both objects (the cycles for `q<=306`, which is all that `s<=2500` needs), and the base loops exist for all records `q<=5626`.
- **Height law (conjecture).** `H(q-1) ~ 1.3 q/(ln 2 * eps(q))` at records. The observed constant is within 3% of 1.3 for every record `1636<=q<=5626`, and ranges over 0.85-1.44 for the small records.
- The 1/2-thread cannot be escaped at a price below 1 (proved). Its escape is inherently a two-stage process.

## 10. Reproduction

```
bash 04-computation/experiments/collatz_procgen_20260922_loops_run.sh > 05-knowledge/results/collatz_procgen_20260922_loops.out
```

Measured run time 10 min 14 s on one core; peak memory about 700 MB (the height DP 361 MB; the checkpointed base loop `s=2300` about 670 MB). Programs:

| program | what it does |
|---|---|
| `loops_dp.c` | layered min-`K` DP with reconstruction, 16-bit `K` |
| `loops_height.c` | lexicographic min-`K`/min-height DP |
| `loops_recon1.c` | checkpointed single-loop or hub-cycle reconstruction, memory `~2 sqrt(s)` layers |
| `loops_hubs.c` | cycles through a hub |
| `loops_family.py` | the splicing closure |
| `loops_verify.py` | independent exact verifier |
| `loops_table.py` | the loop table |
| `loops_escapes.py` | `Bad_inf` proofs, integer checks of the transfer lemmas, price tables |
| `loops_seesaw.c` | level-`r` certificates |
