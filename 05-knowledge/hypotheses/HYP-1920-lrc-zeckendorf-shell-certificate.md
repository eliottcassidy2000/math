---
id: HYP-1920
status: OPEN
source: codex-2026-05-31-S460
related:
  - THM-357
  - THM-360
  - THM-365
  - THM-366
  - HYP-1853
  - HYP-1881
  - HYP-1890
  - HYP-1900
  - HYP-1902
  - HYP-1910
---

# HYP-1920: The n=14 LRC proof has a Zeckendorf-shell certificate

## Statement

Write every speed column by its Zeckendorf support, i.e. as an independent
set in the Fibonacci path.  In this coordinate, the `n=14` Lonely Runner gate
branch is not just a composite-denominator anomaly; it is a finite
Fibonacci-cube row-cover problem.

The mandatory `14`-gate is

```text
14 = F6 + F1 = 13 + 1.
```

Its local endpoint invoice factors exactly:

```text
all minimum lower covers
  = {1,3,5,7,9,11,13} union {one even bridge in 2,4,6,8,10,12}.
```

Thus the private endpoint rows force the whole odd lower fan, and the only
local freedom is a single even Zeckendorf-cube bridge.  The proposed `n=14`
proof should branch on the forced odd fan, quotient the six even bridges as a
local fiber, and then prove coarse denominator rows plus exported owner-debt
rows give positive dual weight on every bridge fiber.

## Evidence

`lrc_zeckendorf_bridge_s460.py` connects three older threads:

```text
Zeckendorf:  independent sets in the Fibonacci path, x=1.
Tournament:  H(T)=I(Omega(T),2), x=2.
LRC:         endpoint rows covered by speed columns.
```

The script re-audits the mandatory gate rows for `n=14,15,16` using lower
columns `1..n-1`:

```text
n=14 exact size 8
  forced = (1,3,5,7,9,11,13)
  minimum covers = 6
  free parts = (2), (4), (6), (8), (10), (12)

n=15 exact size 10
  forced = (1,2,4,7,8,11,13,14)
  minimum covers = 8

n=16 exact size 9
  forced = (1,3,5,7,8,9,11,13,15)
  minimum covers = 1
```

So `n=14` is special: its exact cover family is a product of a forced odd fan
and a one-dimensional even bridge choice.  The forced odd fan has Zeckendorf
digit load

```text
F1:2 F2:1 F3:2 F4:2 F5:2 F6:1.
```

The denominator shell also lines up with the current LRC frontier:

```text
14 = F6+F1
15 = F6+F2
16 = F6+F3
18 = F6+F4  (Lucas/min-gap pair)
21 = F7     (pure Fibonacci reset)
```

This suggests that the `14,15,16,18,21` progression is a sequence of
Zeckendorf-shell payload moves, not merely a list of adjacent runner counts.

Finally, the owner-debt rows from S440 become clearer in Zeckendorf support:
the S380 gate ladder exports its largest exposed labels to high-index owners

```text
154 = F11+F5+F2
168 = F11+F7+F3
182 = F11+F8+F3+F1.
```

Thus gate-heavy repairs increase a Zeckendorf height/payload potential while
they move endpoint debt deeper in the `2 x 7` product tree.

## Interpretation

The old Zeckendorf/tournament bridge says

```text
I(P_m,1) = Fibonacci count
I(P_m,2) = Jacobsthal/tournament path-conflict count.
```

The LRC endpoint program adds a third object: a row-cover incidence matrix
whose columns are natural numbers.  Zeckendorf supports give a canonical
coordinate on those columns.  In the natural-operation language, no-carry
addition is support union in the Fibonacci path, while multiplication is a
carry-heavy sparse shadow.  The `n=14` gate fan is the first place where this
coordinate gives an exact cover factorization.

## Predictions

1. The six `n=14` even bridge choices should remain locally equivalent for
   the gate endpoints but become inequivalent after adding coarse denominator
   rows, primitivity rows, or exported owner-debt rows.
2. A dual certificate for `n=14` should assign the same base weight to the
   forced odd fan and then a bridge-dependent correction that is always
   positive globally.
3. `n=15` should behave as a two-free-column Zeckendorf-shell transition,
   while `n=16` should behave as a rigid dyadic shell; this matches the
   existing Hall-deficit and dyadic-cover evidence.
4. Future LRC feature extractors should include Zeckendorf top index, digit
   load, support overlap, and no-carry/carry status alongside the existing
   `2^h*odd_core`, Bruhat-Tits depth, and product-sum features.

## Sources

- `04-computation/lrc_zeckendorf_bridge_s460.py`
- `05-knowledge/results/lrc_zeckendorf_bridge_s460.out`
- `07-reflections/lrc-zeckendorf-shell-bridge-s460.md`
- `04-computation/zeckendorf_tournament.py`
- `07-reflections/zeckendorf-non-consecutivity-pairing.md`
- `07-reflections/summand-graph-fermat-zeckendorf.md`
- HYP-1900
- HYP-1910
