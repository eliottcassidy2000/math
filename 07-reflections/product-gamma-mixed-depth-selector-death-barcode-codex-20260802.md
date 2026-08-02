# Product-Gamma mixed-depth selector: exact persistence death barcode

**Status:** `FINITE-EXACT REFLECTION / THEOREM-GRADE CANDIDATE / NO GMC OR NC2 CLOSURE`

**Date:** 2026-08-02

**Exact companion:**
[`04-computation/gmc_mixed_depth_selector_death_barcode_scout.py`](../04-computation/gmc_mixed_depth_selector_death_barcode_scout.py)

**Stored transcript:**
[`05-knowledge/results/gmc_mixed_depth_selector_death_barcode_scout.out`](../05-knowledge/results/gmc_mixed_depth_selector_death_barcode_scout.out)

## 1. Inheritance: why depth must become a random variable

THM-3137 classifies probability laws on the eight physical one-pole currents
at support `(1,3)`, bank `I2`.  Such a law repairs every degree 5 through 9,
but degree 10 has an exact two-upset Farkas wall.  A law on the 33 physical
two-pole currents is also killed at degree 10 by one strictly negative upset.

The tempting conclusion would be that randomization has failed.  The correct
conclusion is subtler: **fixed deletion depth** has failed.  The union of the
eight depth-one and 33 depth-two states has a new convex direction unavailable
inside either face separately.

The selector is now a probability law on a stopping state, not just on a pole
value.  This changes the object and must be typed explicitly.

## 2. A three-state law repairs degree ten

On the 41-state simplex, take

```text
lambda_(1)   =2/9,
lambda_(2)   =2/9,
lambda_(2,2) =5/9.                                           (1)
```

The state `(r)` means virtual subtraction of one physical pole letter; the
state `(r,s)` means virtual subtraction of the two-letter multiset.  The two
copies of pole 2 in `(2,2)` are physically available in the reduced pole
multiset.

The exact upset counts for degrees 5 through 10 are

```text
10, 27, 47, 168, 573, 3588,                                 (2)
```

with respectively

```text
8, 25, 45, 166, 571, 3586                                   (3)
```

nonzero 41-coordinate response rows.  Law `(1)` pairs **strictly positively**
with all 4,401 nonzero rows.  Its smallest numerator after multiplying by 9 is

```text
802795680, at N=5 and upset {(5)}.                            (4)
```

Exact max flow equals exact demand in every degree 5 through 10.  Thus mixing
prefix depths is load-bearing: neither pure-depth face meets the degree-10
Hasse cone, while their joint convex hull does.

## 3. The displayed law dies at degree eleven

At degree 11, law `(1)` has

```text
flow   =37492336246227928606600,
demand =37492554756606147580488.                             (5)
```

The deficit is

```text
218510378218973888,                                          (6)
```

at first unresolved type `(8,1,1,1)`.

This does not mean that the degree-11 selector slice is empty.  Conflating
failure of one persistent law with emptiness of a later fibre was the first
numerical near-miss in this scout.

## 4. Degree eleven has its own exact selector

There is an exact four-state probability law supported on

```text
(1), (1,2), (2,2), (7,8).                                   (7)
```

Its rational weights are printed in full in the stored transcript; their
decimal sizes are approximately

```text
0.29871105728189945,
0.15190393645883254,
0.5473981219351144,
0.00198688432415355.                                        (8)
```

All 19,540 nonzero degree-11 upset rows pair nonnegatively with this law;
exactly three are active equalities, and exact max flow saturates demand.

The same law has degree-pass profile

```text
N=5,6,7,8,9,10,11:
  pass, pass, pass, fail, fail, fail, pass.                   (9)
```

So the individual feasible slices are not nested in degree.  A selector may
disappear for three degrees and reappear at degree 11.

## 5. The persistent section dies exactly at degree eleven

Put `P_N` for the probability laws on the 41 states whose degree-`N` current
is Hasse-positive, and define the cumulative feasible set

```text
C_d = intersection_{N=5}^d P_N.                              (10)
```

Law `(1)` proves `C_10` is nonempty.  Four exact upset rows prove `C_11` is
empty.  They are:

```text
F8       =P_8 \ {(1^8)},
F10      =P_10\ {(1^10)},
F11three =P_11\ {(1^11),(2,1^9),(2,2,1^7)},
F11one   =P_11\ {(1^11)}.                                   (11)
```

Take their positive integer combination with coefficients

```text
(461295, 3948, 1, 22).                                       (12)
```

The resulting 41-coordinate response row is strictly negative in every
coordinate.  Its exact range is

```text
-213093302419995239808
 <= coordinate <=
-340278874544980224.                                        (13)
```

If one law belonged to every `P_N`, each of the four rows in `(11)` would
pair nonnegatively with it, so their positive combination would too.  But a
probability vector pairs strictly negatively with `(13)`.  Contradiction.

Therefore

```text
C_10 is nonempty,       C_11 is empty,       P_11 is nonempty. (14)
```

This is a genuine persistence/death-barcode statement.  What dies is not the
degree-11 fibre; it is the possibility of one common selector section across
all scales 5 through 11.

## 6. Holotopy reading

The stochastic selector simplex is a finite parameter space, and each degree
cuts it by an upset halfspace arrangement.  The family `{P_N}` is therefore a
sequence of convex chambers.  Individual chambers can move nonmonotonically,
as `(9)` shows.  Their cumulative intersections `{C_d}` are nested and carry
the honest persistence information.

In this language:

- the depth-one and depth-two faces are individually disjoint from `P_10`;
- their join meets `P_10` at the explicit law `(1)`;
- `P_11` has a separate four-state chamber `(7)--(8)`;
- the four-row dual certificate `(11)--(13)` obstructs any global section
  through degrees 5 to 11.

This is the useful holotopy invariant: not merely whether each scale has a
positive object, but whether the positive objects admit a common section under
the identity comparison of selector coordinates.

## 7. Exact typing boundary

What is preserved:

- the fixed support `(1,3)`, bank `I2`, and distinguished residual alphabet
  `Q`;
- physical multiplicities of all removed pole letters;
- exact zero-mass Young currents and every partition-upset observable;
- one fixed probability law when discussing a cumulative set `C_d`.

What is absent:

- a positive decomposition of the original scalar response or PSD operator;
- a sequential Markov rule whose stopping distribution is `(1)`;
- compatibility with the deterministic THM-3120 pole order;
- a canonical comparison between selector coordinates at different supports;
- arbitrary-degree or arbitrary-radial NC2.

The exact result is about convex averages of virtual-prefix currents.  Calling
it a stochastic stopping rule is justified at the finite current level, but
not yet at the original operator level.

## 8. Next decisive tests

1. **Depth-three restoration.**  Add the physically allowed three-pole
   submultisets and test whether `C_11` revives.  The strict Farkas vector
   identifies exactly which new state directions must cross the wall.

2. **Coupled changing law.**  Replace a common selector by a transport kernel
   from `P_N` to `P_(N+1)`.  The reappearance `(9)` suggests a zigzag or
   sheaf-style persistence object rather than an ordinary nested filtration.

3. **Original-response realization.**  Search the product-Gamma integral for
   a latent stopping variable whose one-/two-pole law is `(1)`.  Mean data are
   insufficient; the THM-3137 variance hostile remains load-bearing.

4. **Cross-support obstruction.**  Repeat the 41-state construction at the
   support `(1,2),I2` Farkas wall.  This is the cheapest test of whether random
   depth repairs the earlier portability failure or is special to `(1,3)`.
