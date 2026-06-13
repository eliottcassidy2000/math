---
source: codex-2026-06-03-S585
status: exact translated-AP audit + recursive proof synthesis
tags: [LRC, additive-circuits, extremal-family, index-stratification, fractal, HYP-2118]
---

# Index-stratified extremal families and the recursive summand fractal

The new thought is:

```text
the extremal family stratifies by index.
```

This is not just a poetic way to say "there are variants of the AP."  In the
additive branch it can be made exact.

Take the translated interval

```text
I(k,q) = {q, q+1, ..., q+k-1}.
```

For fixed `k`, translating by `q` preserves the pair-sum multiplicity profile
up to a shift.  Therefore it preserves the whole 4-term collision count.  But
the visible 3-term folds depend on whether the sum lands back inside the
original window.  In coordinates `q+i`, `q+j`, the sum is visible exactly when

```text
i + j <= k-1-q.
```

So the visible-fold count is a clipped prefix of a fixed triangular pair-sum
profile:

```text
F(k,q) = #{0 <= i < j < k : i+j <= k-1-q}.
```

That is the precise index stratification.

## What S585 saw

At `k=13`, corresponding to the `n=14` additive model, S585 found:

```text
q:        1  2  3  4  5  6  7  8  9 10 11 12 13
folds:   36 30 25 20 16 12  9  6  4  2  1  0  0
4terms:  125 throughout
M/delta: 1.000 -> 4.789
```

The formula has zero mismatches.  The AP at `q=1` is the floor face.  As `q`
increases, the same fold mass is pushed out into hidden shells above the
row.  Raw 4-term energy sees none of this movement.

The hidden node with maximal pressure keeps multiplicity `6`, but its shell
walks outward:

```text
q=1:  s=14, shell=1
q=6:  s=23, shell=5
q=13: s=37, shell=12
```

Adding that hidden sum lowers `M` for every translated AP tested.  The fold is
still there; it has simply moved one layer up in the summand completion.

## The fractal part

The recursively completed object is not the original speed set.  It is the
summand tower

```text
S
S+S
(S+S)+(S+S)
...
```

For interval supports under distinct-pair completion, the support widths obey

```text
w_{r+1} = 2*w_r - 3,
w_r = 2^r*(k-3)+3.
```

For `k=13`:

```text
13 -> 23 -> 43 -> 83 -> 163 -> 323.
```

Each layer repeats the same mechanism: a pair-sum tent profile is clipped by a
window.  The visible 3-fold layer is only one slice of this tower.  A 4-term
relation is the next slice; adding the hidden sum node turns it back into
visible 3-fold data.

This is the precise fractal nature I think we want: not a visual Cantor set,
but an iterated clipped-convolution object whose address is an index vector.

## How it connects to the older indices

The translated-AP index `q` is only one coordinate.

Opus HYP-2116 adds the vertical 2-adic self-affine index:

```text
v -> 2v,
t -> t/2,
v2(v) records the scale.
```

Opus HYP-2117 adds the binary IFS/doubling-map index for translated APs:

```text
D(x)=2x,
T(x)=2x+1,
x -> 2x mod n degenerates exactly on the even frontier.
```

S585 adds the additive summand index: translation clips the same pair-sum tent
into visible and hidden shells.  These are complementary coordinates of the
same extremal tower.

HYP-2091 adds a parity index:

```text
even LRC n -> clean polygon
odd LRC n  -> wall/tie mesh
```

HYP-2092 adds the `C=2n-1` arithmetic shell index.  At `n=14`, `C=27`, so the
antipodal summand shells split as

```text
gcd 1: 9 shells
gcd 3: 3 shells
gcd 9: 1 shell
```

HYP-2079 adds a dyadic depth index:

```text
gap halves,
boundary debt doubles,
gap*debt stays fixed.
```

HYP-2115 adds the virtual completion depth:

```text
4-term collision -> hidden sum node -> two visible 3-folds upstairs.
```

So the proof address of a row should be something like:

```text
(parity lane,
 C-gcd shell,
 visible-fold clip q,
 hidden-sum shell,
 virtual completion depth,
 dyadic debt depth,
 endpoint-owner labels).
```

That is a lot of coordinates, but the payoff is that each coordinate has a
known proof gate.  The problem becomes recursive routing rather than blind
speed enumeration.

## Proof instinct

The AP is not merely "the extremal."  It is the zero-index face where the
pair-sum tent is clipped as far left as possible and the low folds land inside
the original row.  Translating the interval keeps the 4-term shadow but moves
the actual fold payload out of the visible window.  This is why shifted APs
are safe despite maximal 4-energy.

A proof of Lemma A should use this: if the visible fold clip is empty, the
fold payload has been exported to hidden shells.  Unless a hidden shell is
promoted to a real speed, it should behave like a diagnostic rather than a
blocker.  A proof of Lemma B should use the complementary branch: once the
hidden node is visible, the row is back in literal fold/shield territory.

## Tournament Analysis note

S585 used index lenses as tournament vertices:

```text
Phi_gap_value
C_gcd_shell_index
visible_fold_index
hidden_sum_shell_index
dyadic_debt_depth
raw_4term_energy
```

The observable was

```text
(certificate_rank, recursion_rank, preserves_labels, -cost, maturity).
```

The resulting ranking was transitive:

```text
Phi > C-gcd shell > visible fold index > hidden sum shell
    > dyadic debt depth > raw 4-term energy.
```

This is a useful warning.  Raw 4-term energy is the weakest coordinate because
it cannot see the deformation index.  The hidden-shell coordinate is the first
place where the deformation stops hiding information.
