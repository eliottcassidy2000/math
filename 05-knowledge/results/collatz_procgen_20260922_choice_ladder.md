# Choice collapses the Collatz exceptional set: the E-graph descent game, its hostile rationals, and the `6 mod 8` excursion

**Status: PROVED for the scoped lemmas stated as such (hand proofs below);
FINITE-EXACT for every table (two independent code paths agree where
marked); CITED for Applegate--Lagarias (primary source read, sections 2 and
4). OPEN: the E-SCC conjecture (Q1 and Q2), HYP-9120--9122. Collatz OPEN.
No independent-agent audit yet.**

Session `collatz-procgen-20260922` (machine `mac-mini`), lane one of the
procedural-approach session. Scripts: `04-computation/experiments/collatz_procgen_20260922_*`;
reproduce everything with
`bash 04-computation/experiments/collatz_procgen_20260922_choice_ladder_run.sh`
(output `collatz_procgen_20260922_choice_ladder.out`, about four minutes).
The runner stops the DFS at `m=32`. The counts quoted at `m=35,36,40` come
from `e_forward_dfs 26 40 q1bad_dump.txt`, which takes about 70 billion
DFS nodes (about 20 minutes), and the `3^17` count from
`e_reverse_dump 16`.

## 0. Inheritance

* Closest proved mechanism: **Applegate--Lagarias, "The 3x+1 semigroup"
  (J. Number Theory 117 (2006); arXiv math/0411140), CITED.** Their proof of
  the weak 3x+1 conjecture is an induction in which every residue class
  mod `2^j` descends by bounded lookahead *after multiplying by wild
  integers*, except the single class `-1 mod 2^j`, which "resisted
  elimination for `12<=j<=30`" and provably can never be eliminated because
  "all iterates of `-1`, times multipliers, remain negative" (their section
  2); a separate escape lemma (their Lemma 2.2, multiplier `(2^j+1)/3`,
  `T^j(mx)=x+(x+1)/2^j`) handles it.
* Inherited object: the graph `E` (arrows `n->n/2`, `n even`; `n->3n+1`,
  all `n`) and the conjecture E-SCC = Q1 and Q2 of
  [extended_collatz_scc](collatz_mod6_20260917_extended_collatz_scc.md)
  (Q1: every `n` reaches `1`; Q2: `1` reaches every `m` with `3 does not divide m`).
  The earlier session tested Q2 only with the *greedy* 3-adic map `G`
  ([three_adic_g_map](collatz_mod6_20260917_three_adic_g_map.md)).
* **Prior work on `E` itself (found by the barrier-atlas lane).** Le &
  Smith, "Observations on cycles in a variant of the Collatz Graph"
  ([arXiv 2109.01180](https://arxiv.org/abs/2109.01180), unrefereed; the
  abstract was confirmed this session), define the *Loosened Collatz
  Function*: `3x+1` may be applied to even and odd `x`. That is exactly
  `E`. They study its cycle tuples and conjecture that every `n` with `3`
  not dividing `n` lies on a cycle (Conj. 1, per the lane's reading of the
  paper). They do not state Q1, Q2 or strong connectivity. Q1 and Q2
  together give a cycle `1->n->1`, which implies their Conj. 1. **With Q1
  known below `2^71` (every Collatz orbit is an `E`-path; Barina 2025,
  CITED) and Q2 verified here below `2.02*10^13`, Le--Smith's Conj. 1 holds
  for every `n<2.02*10^13` with `3` not dividing `n` (FINITE-EXACT plus
  CITED).**
* Canonical hostile: the rising family `n=2^k u-1` (forward) and
  `m=3^j u+1` (backward, `244->...->256`).
* Least-used sidecar: the choice itself. Every earlier certificate in the
  repository was deterministic.

## 1. The descent game and its exceptional set

For a (possibly nondeterministic) system of affine moves on integers,
a **certificate** for a residue class is a legal path whose legality is
decided by the class and whose total multiplier is `<1`. The
**exceptional set** at level `m` is the set of classes with no certificate
of precision `m`; its intersection over `m` is a closed subset `Bad_inf`
of `Z_2` (forward games) or `Z_3^x` (backward games).

* **Forward E (Q1).** From odd `x`: forced `3x+1`. From even `x`: halve
  (cost `-log 2`, one bit consumed) or excursion `x->3x+1->9x+4`
  (cost `2 log 3`, no bit consumed). Value function over classes mod
  `2^m`: `W_m(x)=min over moves of cost+min(0,W(next))`, `W_0=+inf`;
  exceptional iff `W_m>=0`.
* **Backward E (Q2).** From a 3-adic unit `x`, move `k>=0` with
  `2^k x=1 mod 3`, to `y=(2^k x-1)/3`, legal iff `3 does not divide y`,
  cost `k log 2-log 3`. Classes mod `3^(r+1)`.
* **Partial choice `E_S`.** The excursion is allowed only at even `x` whose
  residue mod `2^J` lies in `S`. `S` empty is Collatz, `S` = all evens is `E`.

Multiplicative descent implies actual descent for every positive member of
the class on the backward side, and for all members above an explicit
threshold on the forward side; conversely actual descent of a positive
integer implies multiplicative descent. Hence **Q1 holds iff no integer
`n>1` lies in `Bad_inf(E)`, and Q2 iff no integer `m>1` lies in
`Bad_inf` of the backward game** (strong induction; small cases by the
earlier `10^6`/`10^7` verifications).

## 2. The ladder (FINITE-EXACT)

Number of exceptional classes:

| system | level | exceptional classes | growth |
|---|---|---|---|
| Collatz `T` (forward, no choice) | mod `2^26` | `1,037,374` | `log2(count)/m=0.769`, rising |
| `E_S`, `S={0 mod 16}` | mod `2^22` | `92,684` | as Collatz |
| `E_S`, `S={0 mod 8}` | mod `2^22` | `90,844` | as Collatz |
| `E_S`, `S={2 mod 8}` | mod `2^22` | `56,550` | slight gain |
| `E_S`, `S={0 mod 4}` | mod `2^22` | `15,796` | moderate |
| `E_S`, `S={6 mod 8}` | mod `2^22` | `3,238` | **most of the gain** |
| `E_S`, `S={2 mod 4}` | mod `2^22` | `1,452` | nearly full |
| `E` (all evens) | mod `2^22` / `2^26` / `2^36` | `782` / `777` / `908` | polynomial-looking |
| greedy `G` (backward, no choice) | mod `3^16` | `12,404` | `log3(count)/J=0.572`, rising |
| backward `E` (Q2) | mod `3^15` / `3^17` | `36` / `52` | tiny |

At the levels where the forward-`E` count drops (the drops follow the
Beatty pattern of `log_2 3`) the counts are
`124, 391, 369, 255, 561, 454` at `m=16,21,24,27,32,35` (to `m=40`:
`1030`). Collatz's count doubles almost every level. Controls: the forward
DP and an independent DFS (bit-precision tracking, branch-and-bound) give
**identical** exceptional sets at `m=27` (`255` classes); the backward DP
and an independent Python DFS give identical lists at `r=6,8,10`
(`14,18,32` classes).

**Exact dimension of the no-choice exceptional set (PROVED).** For the
shortcut Collatz map on `Z_2`, `Bad_inf(T)` is the set of `x` whose parity
sequence has at least `s log_3 2` odd steps in every prefix of length `s`.
The Bernstein--Lagarias parity-vector map is a bijection mod every `2^s`,
hence a 2-adic isometry, so
`dim_H Bad_inf(T) = h(log_3 2) = 0.94995...` bits, where `h` is binary
entropy: the upper bound covers by cylinders with at least `ps` ones; the
lower bound takes Bernoulli(`p'`) paths with `p'>p=log_3 2`, which stay
above the line with positive probability after a prefix `1^N`, and lets
`p'` decrease to `p` (Besicovitch--Eggleston). The finite counts above
approach this exponent slowly. No priority is claimed for this standard
consequence. On the backward side the inherited sharp rate
`exp(-I(c))=0.758751` gives greedy-`G` exponent `log_3(3*0.758751)=0.748`.

So the ladder is: no choice, dimension `0.95` (forward) and `0.75`
(backward); `E`-choice, counts that grow slowly; Applegate--Lagarias
multipliers, the single point `-1`. **The `E` exponent is not
determined.** Local minima of the forward counts (`124, 255, 454` at
`m=16, 27, 35`) grow by about `0.1` bit per level. The backward counts
(`14` at `r=6` to `52` at `r=16`) grow by about `0.12` trit per level.
This fits a small positive dimension, near `0.1`, as well as
polynomial growth (about `m^2`). `dim Bad_inf(E)=0` is therefore
OPEN, not supported. Deciding it needs `m` near `80`, beyond the
present DFS.

## 3. Where the gain comes from: the rising-run excursion (FINITE-EXACT)

Allowing the excursion only at evens `=6 mod 8` removes `96.5%` of the
Collatz exceptional classes at `2^22`; `2 mod 4` removes `98.4%`, close to
all of `E`. Every even `=6 mod 8` that follows an odd step is `3x+1` with
`x=7 mod 8`, the start of a run of at least three rises; halving it
continues the run (`(3x+1)/2=3 mod 4`). The collapse is therefore driven by
**the freedom to interrupt a rising run**, i.e. to leave the 2-adic
neighbourhood of `-1`. Excursions at `0 mod 8` do nothing measurable.
HYP-9121 records the conjecture that `S={6 mod 8}` already has
subexponential exceptional growth.

**Greedy fingerprint (FINITE-EXACT, `2^20`,
`collatz_procgen_20260922_greedy_choice_fingerprint.py`).** The script
adds even classes mod `64` to `S` greedily, each time taking the class
that most reduces the exceptional count. The order is `54, 62, 60, 28, 22,
14, 50, 30`, with counts `14509, 7643, 4503, 2369, 1631, 1282, 1046, 893`,
against `27328` for Collatz and `664` for all evens. Eight of the 32 even
classes mod 64 recover `97%` of the full collapse. The first choice is
`54=-10 mod 64`, the image `3x+1` of odd `x=39 mod 64`, the entry of a
three-rise run; it alone removes `47%`, more than the deep `-1`
neighbourhood `62 mod 64` removes first. Five of the eight classes are
`6 mod 8` (`54, 62, 22, 14, 30`), two are `4 mod 8` (`60, 28`) and one is
`2 mod 8` (`50`).

**The additive-choice zoo (FINITE-EXACT, `2^22`,
`collatz_procgen_20260922_additive_choice_zoo.c`).** Suppose each odd step
may use `3x+b` for any `b` in a set `B`. Then every `B` with two distinct
odd elements has **zero** exceptional classes: `{1,-1}`, `{1,3}`, `{1,5}`,
`{1,9}`, `{1,17}`, `{1,33}`, `{1,-3}`, `{1,7}`, `{1,-7}`, `{1,5,9,13}`.
Singletons `{1}` and `{-1}` give `93,222`. All odd `b` are 2-adically
conjugate (`x->lambda x`), so switching between two copies moves the orbit
between disjoint hostile sets, and nothing stays hostile. For `{1,-1}` the
reason is elementary: exactly one of `3x+-1` is divisible by `4`, the case
the barrier-atlas lane calls the trivial one-player form of Althöfer's
`3n+-1` game. The ladder of freedoms is therefore: none (dimension
`0.95`), branch choice `E` (thin), additive choice (empty). The pasted
snippet's `3n+sgn(n)` is a fixed rule, not a choice, so it gains nothing
here.

## 4. The hostile points are rationals over the other prime

**Forward `E` (2-adic), PROVED membership.**

* `-1` is in `Bad_inf(E)`. Every `E`-path from `-1` stays in the negative
  integers. A prefix with `a>=1` multiplications and `b` halvings ends at
  `(-3^a+B)/2^b<=-1`, with carry `B=sum 3^(a-1-i)2^(b_i)>=3^(a-1)`, so
  `3^a>=2^b+B` and the multiplier satisfies `r=3^a/2^b>=1+r/3`, i.e.
  `r>=3/2`. This is Applegate--Lagarias's argument, and it transfers
  verbatim to `E`.
* `-13/9` is in `Bad_inf(E)`. The first move is forced (`-13/9` is 2-adically
  odd) to `-10/3`, and all later values are negative. If the second move
  halves, the third is forced and `B/2^b>=(5/9)r`, so `r(13/9)>=1+(5/9)r`,
  giving `r>=9/8`. If the second move is an excursion, `B/2^b>=(4/9)r`,
  giving `r>=1`, and `r=1` is impossible. Earlier prefixes have `r=3, 3/2, 9`.

**FINITE-EXACT structure.** The `908` exceptional classes mod `2^36` were
reconstructed as rationals of small height. Every reliable reconstruction
is negative with a power-of-3 denominator, in about `[-1.6,-1]`:
`-1, -13/9, -35/27, -97/81, -113/81, -275/243, -307/243, -355/243, -371/243,
-793/729, ...`, accumulating at `-1` from below. The family is not bounded by `-3/2`:
`-371/243=-1.527` and `-43/27=-1.593` also show no descent. Exact-rational DFS
from each listed point finds no descent before a 2-million-node limit, whereas nearby
rationals `-5/3, -11/9, -37/27, -41/27, -7/5, -9/7` descend within 2 to 5 halvings
(`collatz_procgen_20260922_hostile_rational_check.py`). An exact membership
criterion for the family is OPEN (HYP-9120 item 1). The positive
reconstructions have height products near `5*10^7`, where about ten chance
matches are expected among `908` classes; they are not claimed. The only
exceptional class mod `2^36` containing an integer of absolute value below
`1.5*10^8` is `-1`. Hence **every positive integer below `1.5*10^8` lies in
a class with a multiplicative certificate of precision 36**. On the
forward side, a multiplicative certificate becomes actual descent only
above the path's threshold `B/(2^b-3^a)`. Q1 itself is implied
pointwise by Collatz, since Collatz orbits are `E`-paths. It therefore
holds for every `n<2^71` by the published Collatz verification (CITED,
Barina 2025, per the barrier-atlas lane). Q2 is not implied by
Collatz (the arrows point the other way); section 5 gives its own
verification (`2.02*10^13`).

**Backward `E` (3-adic), PROVED membership.**

* `1` is exceptional. Reverse paths from `1` stay at positive integers,
  since the move `k=0` from `1` gives `0` and is illegal. So
  `2^K=3^s x_s+B>=3^s+B>3^s`.
* `1/2` is exceptional. The move `k=1` gives `0` and is illegal, and every
  legal first move lands on an integer `>=1`. So `2^K>=2(3^s+B)` and
  `r>2`.

**FINITE-EXACT structure.** The reconstructed exceptional threads are `1`,
`1/2`, and positive dyadic rationals in about `(0.8, 1.5)`: `43/32, 59/64,
145/128, 209/256, 371/256, 499/512, 661/512, 715/512, 1241/1024, 1753/2048,
4235/4096`. These are provisional at modulus `3^17`. Several classes agree
with `1/2` only to lower precision.

**Duality.** On the forward, 2-adic side the hostile rationals have
3-power denominators and are negative. On the backward, 3-adic side they
have 2-power denominators and are positive. Negation maps the plus sheet
to the minus sheet (`T_+(-n)=-T_-(n)`), so for the minus-sheet `E_-` the
forward hostile rationals lie on the positive side: `+1` and non-integers.
This matches the inherited observation that `Q1_-` holds to `10^6`. The
plus sheet's negative Collatz cycles are exactly where the relaxed
plus-sheet problem stalls.

## 5. Escapes

**Lemma (1-escape, PROVED).** Let `m>=2` with `v_3(m-1)=k>=2`. Suppose
there is a forward `E`-cycle through `1` with `s=k-1` multiplications and
`K` halvings, `2^K<3^k`. Then `m` has a legal reverse path to
`2^K (m-1)/3^k<m`, which is not divisible by `3`.

*Proof.* Write `m=1+3^k u` with `3` not dividing `u`, and run the loop's
reverse moves from `m`. The residues agree with those from `1` modulo
`3^(k-i)` at step `i`, so each of the `s=k-1` moves is legal. The endpoint
is `1+3*2^K u`, which is `4` or `7 mod 9`. Hence the move `k=0` is legal
and gives `2^K u`, which is not divisible by `3`. Finally
`2^K u<1+3^k u` iff `2^K<=3^k`.

**FINITE-EXACT.** A loop of every length `s<=40` exists with minimal ratio
`2^K/3^s` in `[1.517,2.96]`, always `<3`, with values bounded by
`2*10^7`. So every `m` with `2<=v_3(m-1)<=41` descends, which covers the
hostile point `1` for all `m<3^42+1`. The margin is thin: `s=11` gives
`2.9596`, and no loop with ratio `<3/2` appears for `s>=2`. The lemma may
therefore fail at a rare `s`; a robust version must allow other exit
words (HYP-9122). Every loop through `1` has ratio at least
`13/9` for `s>=2`, from the last two carry terms.

**Uniform descent off two classes (PROVED; elementary; Lean-checked).**
Let `m>1` with `3` not dividing `m` and `m` not `1` or `14 mod 27`. Then at
most two reverse `E`-moves reach an integer `y<m` with `3` not dividing
`y`. The certificate depends only on `m mod 27`:

| `m` | route (forward reading `y -> ... -> m`) | `y` | factor |
|---|---|---|---|
| `4, 7 mod 9` | `y -> 3y+1 = m` | `(m-1)/3` | `1/3` |
| `2, 8 mod 9` | `y -> 3y+1 = 2m -> m` | `(2m-1)/3` | `2/3` |
| `10, 19 mod 27` | `y -> 3y+1 -> 4m -> 2m -> m` | `(4m-4)/9` | `4/9` |
| `5, 23 mod 27` | `y -> 3y+1 -> 8m -> 4m -> 2m -> m` | `(8m-4)/9` | `8/9` |

Kernel-checked in core Lean 4.30 (no Mathlib) as `q2_descent_off_1_and_14`
in `04-computation/lean/standalone/collatz_procgen_20260922_q2_mod27_descent.lean`
(`lean <file>`; axioms `propext`, `Quot.sound` only). The DP and the
exact-rational search (`collatz_procgen_20260922_q2_uniform_descent_check.py`)
confirm it independently. With lookahead 6 the worst class still has
factor `8/9` (`..._q2_worst_factor.c`), so deeper search does not beat this
table. **Mirror theorem for the `3n-1` sheet** (same Lean file, also
kernel-checked with axioms `propext`, `Quot.sound` only):
`q2_minus_descent_off_26_and_13` covers graph `E_-` (`n->3n-1`, `n->n/2`)
for `m>4`, `3` not dividing `m`, `m` not `26` or `13 mod 27`. These are
the negations of `1` and `14`. The hypothesis `m>4` is necessary. At
`m=4` the route `(8m+4)/9` returns `4` itself, along the `E_-` cycle
`4->11->32->16->8->4`. Its negation is the plus-sheet cycle
`-4->-11->-32->-16->-8->-4`, the inherited `G`-cycle `{-4,-11}`
([g_negatives_joint_carry](collatz_mod6_20260921_g_negatives_joint_carry.md)).
So the two mirror lemmas differ exactly at a known signed cycle.
Direct BFS confirms that `1` reaches every `m<=5000` prime to `3` in `E_-`
(`collatz_procgen_20260922_q2_minus_small.py`).

In the two excluded classes the two-move routes fail. From
`14 mod 27`, `(8m-1)/3=1 mod 9` and the next `k=0` move gives a multiple
of `3`. From `1 mod 27`, `(4m-1)/3=1 mod 9` fails in the same way.
**Consequence: Q2 reduces exactly to the two hostile neighbourhoods
`m=1 mod 27` and `m=14=1/2 mod 27`.** The thread computation to depth 27
agrees. The `134` exceptional classes mod `3^28` are the class of `1` and
`133` classes `=1/2 mod 27`; the latter include the dyadic points
`1/2+3^j/2^e`, such as `43/32=1/2+27/32`, `209/256=1/2+81/256` and
`4235/4096=1/2+2187/4096`. The counts at depths `12,19,24,27` are
`30,67,94,134`.

**Finite verification of Q2 to `7.1*10^12` (FINITE-EXACT).** The
thread DFS (`collatz_procgen_20260922_e_reverse_dfs.c`) agrees with the DP
at every depth `r<=16`. Pushed to depth `30`, the exceptional counts are
`67, 94, 134, 190` at `r=19, 24, 27, 30`. The smallest positive
representative of an exceptional class other than `1` grows like
`3^(r-4)`: it is `1.86*10^6` at `r=16`, `1.96*10^8` at `r=19`,
`7.17*10^12` at `r=30`, and `2.02*10^13` at `r=31` (about `2*10^10` DFS nodes). On the backward side a multiplicative certificate
always gives actual descent (the carry only subtracts), and values stay
positive integers prime to `3`. Strong induction therefore shows that
**`1` reaches every `m<2.02*10^13` with `3` not dividing `m` in `E`**. The
inherited verification was `10^7`, by greedy `G`.

**The `1/2` neighbourhood costs at least `8/3` (PROVED).** Every legal
reverse path from the 3-adic point `1/2` first lands on an integer
`y>=1`. The move `k=1` gives `0` and is illegal. So
`y=(2^K/3^s)(1/2)-B/3^s`, which gives `2^K/3^s=2(y+B/3^s)>=2(1+1/3)=8/3`.
Hence an integer `m` with `v_3(2m-1)=j` that follows a path of `1/2` for
`s<j` steps multiplies its offset `m-1/2` by at least `8/3`. Routing through
`1` with the 1-escape lemma never descends directly: the offset multiplier
is at least `(8/3)(13/9)/3>1`. Closing Q2 near `1/2` therefore needs an
amortized, Applegate--Lagarias-type argument. It must charge each escape
(bounded cost) against the 3-adic digits it consumes, and use the
uniform `8/9` lemma on fresh classes. This is the precise remaining gap
for Q2 (HYP-9120).

**The exits of `1/2` are the Collatz trunk (PROVED, and checked on 2000
random cases).** Let `2m-1=3^j w` with `j>=2`. The legal first reverse
moves from `m` are exactly `k=2i+1` with `i` not `0 mod 3`, and they land at
`N_i+4^i 3^(j-1)w`, where `N_i=(4^i-1)/3=1,5,85,341,5461,...` is the
inherited trunk `R^(i-1)(1)`
([row_braid_typing](collatz_mod6_20260917_row_braid_typing.md)). The
illegal ones, `N_3=21, N_6=1365, ...`, are multiples of `3`. So the hostile
point `1/2` is the common 3-adic limit of the trunk's reverse images, and
every escape passes near a trunk integer. The cheapest route goes through
`N_1=1`. Its offset cost has floor `(8/3)(4/3)/3=32/27` at `j=3` and
`(8/3)(13/9)/3=104/81` for `j>=4`. The post-escape value `X=2^(K+2)w`
depends only on the loop's total halving count `K`, which a ratio `<3`
allows in at most two values. An adversarial `w` can put one choice near
`1/2` again and the other near `1`. So Applegate--Lagarias's near-free
escape (their cost is `1+2^(-j)`) has no analogue here. The amortization
must use the descent near `1`, or steer by exits through the higher trunk
points `N_i`.

**The forward half does not reduce to one class (FINITE-EXACT,
`collatz_procgen_20260922_q1_worst_factor.c`).** Outside `x=-1 mod 4` every
forward class descends at once (factor at most `3/4`). But at `2^22` there
remain `6, 18, 48, 120` exceptional classes outside `-1 mod 8, 16, 32, 64`,
and the worst certified factor outside `-1 mod 2^J` is `1.5` for `J>=5`.
The exceptional classes mod `2^36` have `v_2(c+1)` spread from `2` to `36`
(bulk `4..9`). The forward hostile family reads `-1-2^i/3^j` (`-13/9=-1-4/9`,
`-35/27=-1-8/27`, `-97/81=-1-16/81`, `-113/81=-1-32/81`, `-275/243=-1-32/243`,
`-307/243=-1-64/243`), accumulating 2-adically at `-1` from depth `i=2`.
The backward family `1/2+3^j/2^e` starts at depth `3`, and the backward
move `k=0` contracts by `1/3`, which no forward move matches. This is why
Q2 is the cleaner half.

**Chained hostile landings are common (FINITE-EXACT,
`collatz_procgen_20260922_half_chain.c`).** Take every
`m=14 mod 27` up to `3*10^7` (`1,111,111` values) and search for a shortest
reverse path below `m` (iterative deepening; excursions capped at
`1000m`; depth at most `16`):
* depths are mostly `4..7` (histogram `4:493827, 5:164609, 6:274348,
  7:118887, ...`), with maximum `16` (at `m=2520518`);
* `16` values need more, all deep in the `1/2` neighbourhood, e.g.
  `m=(3^11+1)/2=88574` needs `14`;
* the largest excursion on a shortest path is about `959m`;
* **`31.4%` of shortest paths pass through another hostile class**
  (`1` or `14 mod 27` at depth at least `3`).

So chains of `1/2`-type landings are typical, not exceptional, for
integers. Any amortized proof of Q2 must pay for them. This is the empirical
side of the gap analysed in the synthesis (section 5, item 1).

**The price of escaping `-1` on the forward side (FINITE-EXACT upper bounds).**
This is the least multiplier over `E`-paths from `-1` that consume exactly
`b` bits, counting the forced multiplication after the last halving. It
stays in `[4.05,34.2]` for `b<=70`, oscillates with `frac(b/log_2 3)`, and
shows no growth. Cheap long negative `E`-cycles, such as the `{-5,-7}` ride
at `9/8` per three bits and the seven-cycle at `2187/2048` per eleven,
keep it bounded. Plain Collatz pays `(3/2)^b`.

## 5b. Which half of Collatz the class-level method can see (PROVED, elementary)

For the shortcut map, a positive integer `n` whose orbit enters a positive
`3n+1` cycle descends multiplicatively. A positive cycle satisfies
`x(2^K-3^L)=B>0`, so its loop multiplier `3^L/2^K` is `<1`, and repeated
loops drive the path multiplier to `0`. Hence **the plus-sheet exceptional
set `Bad_inf(T)` contains no positive integer that is eventually periodic.
Its positive integer points, if any, have divergent orbits.** Nontrivial
plus-sheet cycles are invisible to multiplicative descent. They sit
exactly at the actual-descent thresholds `B/(2^K-3^L)`, which is the
inherited cycle gate.

On the minus sheet the loop multiplier is `>1`: `3/2, 9/8, 2187/2048` for
`{1}, {5,7}` and the seven-cycle. The cycle minima `1, 5, 17` are
multiplicatively hostile (checked over 60 odd steps). Non-minimal cycle
points such as `7` descend to the minimum.

So the class-level instruments of this lane address only the divergence
half of Collatz. The sign enters as follows: positive plus-sheet cycles
contract multiplicatively, positive minus-sheet cycles expand. This is a
concrete *order* statement, the sign of `2^K-3^L`, of the kind the
inherited word-function theorem requires.

## 6. A negative control: the undirected Collatz game

Allowing backward moves along `T` itself gives a game that is
**equivalent to Collatz (PROVED)**. Every `m>1` has an undirected
`T`-path to a smaller integer iff every component contains `1`, by taking
the minimum of a component. It does **not** collapse the exceptional set
(argued, and supported FINITE-EXACT: `collatz_procgen_20260922_undirected_game.c`,
output `..._undirected_game.out`). A backward `T`-move to an odd predecessor
`(2^k y-1)/3` requires `k>=1` and yields nothing at multiples of `3`.
From a class `=1 mod 3` the cheapest predecessor grows by `4/3`, so the
3-adic hostile set of the backward `T`-game has positive measure, for
example `3Z_3`. The one move `T` lacks is `E`'s `k=0` move, the even
predecessor `(y-1)/3`. This locates the relaxation's power exactly.
A BFS on the undirected shortcut graph over all `3<=m<=10^6` looks for the
first value `<m`, with depth at most `60` and values at most `64m`. The
mean depth is `2.21`: all evens and all odd `m=2 mod 3` take one move, via
`m/2` or `(2m-1)/3`. The tail does not collapse. The maximum depth grows
(`40, 43, 59` at `N=10^3, 10^4, 10^5`), and the unresolved fraction stays
near `0.2%` (`2077` at `10^6`), much like forward stopping times.

## 7. What this says about the anchor

* **Mechanism.** Choice turns the `0.95`-dimensional exceptional set into a
  far thinner one: finite-level counts fall by three orders of
  magnitude. Its rational points lie over the other prime, and its
  exact dimension is open (section 2). The Collatz obstruction is
  therefore not a density phenomenon but a *no-choice* phenomenon. A proof
  must show that positive integers avoid a `0.95`-dimensional 2-adic
  Cantor set. This is a transversality problem of the same type as
  Erdős's ternary digits of `2^n` and Mahler's `Z`-numbers.
* **Where the missing freedom sits.** The relaxation's gain is concentrated
  in interrupting rising runs (`6 mod 8`). Collatz offers no such move, so a
  proof must control long rises with order or archimedean information.
  This is consistent with the inherited word-function theorem.
* **Choice is generic (correction).** The same choice collapses the `5n+1`
  relaxation: its mass falls to about `4e-5` at `2^64`, and choice games
  tip only near `q=10`
  ([sibling ladder](collatz_procgen_20260922_sibling_dimension_ladder.md) §4;
  MISTAKES 2026-09-22). So the ladder locates the difficulty but does not
  carry the drift or the sign, the two features that make Collatz
  special. A proof strategy suggested by the relaxation must add both.
* **E-SCC proof program (HYP-9120).**
  1. Finite certificate tables at a fixed level.
  2. Escape lemmas for the hostile families `-1` (bounded price, section 5),
     the `-p/3^j` family, `1`, `1/2` and the dyadic family.
  3. An Applegate--Lagarias see-saw induction closing when the
     bounded-lookahead descent factor times the escape price is `<1`.

  Open are the exact description of `Bad_inf`, uniform escape lemmas for
  its infinitely many rational points, and the bookkeeping of repeated
  escapes. For Q2 the first step is done: the Lean-checked mod-27 lemma
  leaves only `1` and `14 mod 27`. The last step is the real core: near
  `1/2` escapes cost at least `32/27`, not `1+2^(-j)`, and chains of
  hostile landings are typical (section 5).

## 8. Research cards used / candidate

Used: "Separate unbounded local support from a height-bounded modular cover",
"Test structured adversaries, not only random samples", "Search the
statement before the method" (it found Applegate--Lagarias's `-1`
obstruction), and "Use redundant code paths as detectors" (DP versus DFS).

**Candidate card: "Measure the exceptional set of the relaxation ladder."**
When an all-orbits statement resists proof, compute the exceptional-set
growth for the no-choice system, a partial-choice family and a
full-choice relaxation. A large drop in exceptional growth pins the
missing freedom (here the `6 mod 8` excursion). The hostile points of the
relaxed system are its escape obligations.
Counterindication: equivalent reformulations with choice, such as the
undirected game, can keep the full dimension. Evidence: this lane, and
Applegate--Lagarias's single-class obstruction (a distinct thread). Not
promoted.
