# Clone clocks, 3H+2 and Collatz on the negative integers

klein, 2026-11-13 (session `chanpak-881-negcollatz`). **PROVED IN NOTES / FINITE-EXACT; Chan–Pak remains OPEN.**
Owner directive: merge the `462·2^k−5 → 1386·2^k−13 → 4158·2^k−37` chain and the relation of `3H+2` to
Collatz on the negative integers; escape the fixed-block obstruction at 881 with changing cores and selective
ports so that `405·2^k+71` becomes a repeatable ready return; connect 3±1, 4±1, 8±1, the non-multiples of
three, and the Jacobian / Hadwiger–Nelson / Hodge lanes.

| Result | Status | Where |
|---|---|---|
| `3H+2 = T(2H+1) = 2T(H)+1` on odd `H`. In `r=(H+1)/2` the three C2 digits are `3r−1, 3r, 3r+1`: the `3x−1` (negative-Collatz) step, pure tripling, and the `3x+1` step. The non-multiples of three are exactly the `r` reached by a Collatz-type step. `r(7)=3·1+1`, `r(21)=3·4−1`. | PROVED | [COLLATZ §1](COLLATZ.md) |
| Dyadic rays shadow Collatz cycles with Størmer exchange rates. In particular `T^{3t}(462·8^t−5)=462·9^t−5` on the Catalan cycle `{−5,−7,−10}`, and the chain's anchors `−5→−13→−37` end on the unbalanced `−17` cycle. | PROVED / VERIFIED | [COLLATZ §2](COLLATZ.md) |
| `4x±1` and `8x±1` specify two and three odd-branch applications of `3x∓1`. These are strict ascents for positive sources, except the minus-map fixed point1. | PROVED (classical; corrected indexing and boundary) | [COLLATZ §3](COLLATZ.md) |
| Every transitive-clone ternary clock `H_m=A3^m−σ` has `A` even and **σ odd**, and is constant mod 4. Clone clocks never realize `3H` or `3H+4`. | PROVED | [CLOCKS §1](CLOCKS.md) |
| Complete census of ternary clocks of every digit through order 9: `0,4,30,206,1309`. Digit-1 amplitudes are `60,100,156,180`, extending Session 197's minimum. Near-transitive order 10 adds `260,300,468,540`. None of 416 digit-1 clocks is ready at its first member. | FINITE-EXACT | [CLOCKS §1](CLOCKS.md) |
| **Two-clock skeleton Π (9 vertices):** `H(Π[u←M,v←N]) = 60·P_M(3)·P_N(2) + H(M)·H(N)` for all tournaments `M, N`, i.e. `1+60·2^k·3^j`: exact repeatable returns `H↦3H−2` and `H↦2H−1` by one-vertex growth. | PROVED (complete table) | [CLOCKS §3](CLOCKS.md) |
| **One arc reversal realizes the negative-integer Collatz step:** in Π, reversing the arc from the last `v`-clone to `c` gives `H↦(3H−1)/2`, i.e. `T(−H)=−H'`. | PROVED | [CLOCKS §4](CLOCKS.md) |
| **Digit-3 two-clock Ψ (10 vertices):** `H = 160·P_M(3)P_N(2) + 40·P_M(3)H(N) − 3·H(M)H(N)`, i.e. `3^j(160·2^k+40)−3`: an exact `H↦3H+6` return for every `k`, with a genuine dyadic parameter. | PROVED (complete table) | [CLOCKS §3b](CLOCKS.md) |
| Two-clock census, 21,047 lines in five universes, complete through order 9: the C2 shape `(σ,σ_2)=(1,0)` never occurs. Pure `v`-sectors occur with digits `−11,−9,−5,−3,−1,3`. The stronger claim pushed earlier, constant term never `−1`, is **refuted**: products such as `(60·3^j−1)(2^k+1)` and one genuine family. | FINITE-EXACT; hypothesis OPEN; stronger form REFUTED | [CLOCKS §5](CLOCKS.md) |
| **Infinitely many seeded unit lines:** 217 certified unit two-clocks `H=1+(α3^j−σ_2)2^k` give, for every `j`, an explicit ray of the line `1+c·2^s`. That seeds infinitely many odd `c` for HYP-9658 (A), e.g. `c=30·3^j−1`, whose slopes follow the C2 map `c↦3c+2`. All 868 two-clock tuples are certified by complete tables. | PROVED (certificates) | [CLOCKS §6](CLOCKS.md) |
| A ready repeatable 405 return cannot be a double-cloning family (from the complete 881 fibre). Uniform selective ports reduce to modules. **Concurrent result:** the [complete intrinsic inverse](../target_directed_cores_20261112/TARGET_INVERSE.md) shows `405·2^17+71` has no in-degree-one unit-ready presentation at any order, so no construction makes the whole ray unit-ready. Non-unit interfaces, direct target surgery, or a change of arithmetic (Π, Ψ) remain. | PROVED corollary + route map; ray closure CITED (concurrent) | [RAY405](RAY405.md) |
| 2–3 axis, planar Jacobian shears, Hodge bidegrees, Hadwiger–Nelson shape, and the old rational/irrational, odd/even, positive/negative, addition/multiplication dualities | typed bridges; analogies flagged | [BRIDGES](BRIDGES.md) |

**What this changes.** The question "can `3H+2` be repeated exactly by a physical move?" has a sharp answer in
the clone language. One clock can do it: Session 197's `60·3^m−1`, which is never ready in the screened range.
Two clocks can do it exactly for other digits. On the +1 anchor (Π) the moves are `x↦3x+2`, `x↦2x+1` and the
Collatz step **on the negative integers** `x=−H`, the last realized by a single arc. On the −3 anchor (Ψ) the
move is `H↦3H+6` with a free dyadic parameter. The digit-1 shape `H+1=3^jF(k)`, which is the C2 line, never
appears. The obstruction, if real, is specific to digit 1 together with a pure `v`-sector: the constant `−1` alone
does occur, for instance in products like `(60·3^j−1)(2^k+1)`. With the inherited 881
fibre, a ready repeatable `405·2^k+71` return must leave the clone/port language entirely. Concretely: cores
that change with `k` in a non-uniform or numeration-driven way, or an arithmetic other than exact `3H+2`.

The same census pays a dividend to the existing 2-adic-line program. Each unit-offset two-clock seeds a whole
sequence of lines `1+c·2^s` with `c=α3^j−σ_2`, so infinitely many slopes `c` beyond every finite census get
explicit rays. In the family `c=30·3^j−1` the C2 map acts on the slope instead of the count.

**Not claimed.** No new Chan–Pak range, no all-order existence of the missing C2-line two-clock or its
nonexistence, no Collatz, Jacobian, Hodge or Hadwiger–Nelson statement. Literature priority was searched only
inside this repository. The `3x−1` ⇔ negative-integer equivalence, the parity-vector lemma and the ascent classes
are classical.

## Inheritance and credit

Session 197 (minimum ternary clock 60, run expansion, matching bound), Session 299 §4, the Session 295
excursion-factorial filtration, the path-bin calculus (P, E, S), the six-state funnel compiler
(berggren_target_return §3; reimplemented and checked against brute force here, `scripts/funnel.py`), the 881 fibre
and finite-atlas theorems (Session 296, mrb_seven_cutoffs), the multicore F3 chain and inflation exclusions,
MODULE_RETURN, the same-night H17 exact-repair packet, HYP-9339 (2–3 axis), and the dyadic-pair template census
(kind-pasteur S21 `modfam.c`, mac-mini dyadic_cones). That census covers modules visited at most twice; the two-clocks
here are its ternary-by-dyadic continuation. Concurrent changing-core, growing-t, ternary-grammar and surplus work is
mapped in [INCOMING](INCOMING.md).

Concept board at close: exact clone clocks / sector digits `(σ,σ_2)` / Collatz exchange rates / tail-law
offsets / readiness. The decisive next test is a single two-module table with `(σ,σ_2)=(1,0)`. It would refute
the C2-line exclusion and give an exactly repeatable `3H+2` with a free dyadic parameter. The alternative is a proof
that a pure `v`-sector forbids digit 1; the tail-law screen of CLOCKS §7 (pure responses always carry offset
digit 3) shows where such a proof should look. META-PATTERNS card used: "type every analogy and implication". A
candidate lesson for the cards, from the refuted constant-term hypothesis: run the closure operations a
coefficient respects (products, converse) before proposing an obstruction for it. It is not promoted from one
session.

## Reproduce

    python3 -B 05-knowledge/results/clone_clocks_negative_collatz_20261113/verify.py        # ~30 s, needs cc
    python3 -O -B 05-knowledge/results/clone_clocks_negative_collatz_20261113/verify.py
    python3 -B 05-knowledge/results/clone_clocks_negative_collatz_20261113/verify.py --full # + order-10 readiness

Census regeneration (C, `gentourng` from nauty): `c/ternscan.c` and `c/ternscan_all.c` read `gentourng -q n` on
stdin; `c/backscan.c n maxb part nparts`; `c/extscan.c e`; `c/extone.c part nparts` reads skeleton lines
(`n u v rows…`) and extends each by one vertex; `c/certify.c` computes complete tables. Outputs are in `data/`
(gzip). Tail-law screens: `python3 -B scripts/tail_c2_multi.py 3` and `scripts/tail_constants.py` (about 45 s each). Exploratory scripts
that found nothing are kept with their stated scope: `scripts/ternary_funnel_search.py`,
`scripts/twoclock_search.py`, `scripts/eig3.py`, `scripts/atlas.py`.
