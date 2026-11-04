# An exact inverse for one variable binary run

**PROVED within the declared 34-core family — [independent audit](FAMILY_AUDIT.md) passed.** This decides every
one-variable-run ready-funnel template in an explicitly supplied finite core
bank. There is no prefix or suffix length cutoff. Arbitrary clock-dependent
words, several variable runs and other cores are outside this language.

## Inherited ready funnel and the smaller state

The [inherited ready-count formula](../../two_entry_funnels_20261103/TRANSPORT.md)
uses actual core data `(h,p2,e,b,F,G,C_joint)`. Write

```
c=G+p2,       K=C_joint+F+b.
H(word)=c L(word)+h B(word)-e D(word)-K.              (1)
```

The tail word is empty when the core entry is proper, or nonempty ending1.
Here L is the inherited tail L(2), and B,D are B(1),D(1). A final converse
and the indegree-one ready birth produce an actual start-ready auxiliary
whose sole out-neighbour is the retained root. Its start gap is the core
entry count b. The named construction is retained by the selector.

For a word prefix p of length n define five integers `(T,B,eta,L,Z)`, initially
`(1,1,1,1,1)`. Appending a bit a gives

```
T'   = 2T-a,
B'   = B                         if a=0,
       B+T-eta                   if a=1,
eta' = a,
L'   = L+aZ,
Z'   = 2Z.                                             (2)
```

Then Z=2^n, L=1+v(p), where v reads p in little-endian order, and
`T=2^n-v(reverse(p))>=1`. The recurrence for L is immediate from its
weighted-bit formula. For B, specialize the inherited six-state tail
automaton to m=1. There are T assignments with a permitted first after-root
vertex, ignoring the last before-root condition; one assignment has no
after-root vertex. Assigning a new0 to the before-root side is inadmissible,
while putting it after permits exactly the old B good assignments. For a
new1, all T before assignments are good, and the after assignments add B
minus the no-after assignment exactly when its old last bit was good.
That indicator is eta, with eta=1 for the empty prefix. This gives (2).
For nonempty words, D is1 precisely when the first bit is0; D(empty)=1.

This state computes the same scalar as the inherited full endpoint carrier.
It does not authorize deletion, identify physical names, or replace joint
root-response checks on a selected source.

## Complete slope factorization

Consider every template

```
p a^k q,       a in {0,1},       k>=1,                (3)
```

with arbitrary fixed finite p,q and a valid final tail bit. If a=0 then q
must be nonempty ending1. If a=1, q may also be empty. Put X=2^k. Iterating
(2) along the variable run gives exact affine expressions:

| Run | T | B | L | Z |
|---|---|---|---|---|
|0^k|X T_p|B_p|L_p|X Z_p|
|1^k|X(T_p-1)+1|X(T_p-1)+B_p-T_p+2-eta_p|X Z_p+L_p-Z_p|X Z_p|

The last-bit state after the run is the fixed bit a. Each subsequent
application of (2) is affine, so (1) is exactly `A 2^k+B0`, not a fit to
sampled clocks. Its leading coefficient has the factorization

```
A = v(q)     (c 2^|p| + h T_p)          for a=0,
A = (1+v(q)) (c 2^|p| + h(T_p-1))       for a=1.      (4)
```

For zero runs, each1 of q adds `2^i T_p` to the B slope and `2^i Z_p`
to the L slope. For one runs the run already contributes one unit of each
slope before those suffix contributions. This proves (4), including empty q.

For a prescribed positive A, both parenthesized factors are at least
`c 2^|p|`. Thus `c 2^|p|<=A` is a proved prefix-length bound. For each such
p and each bit a, test divisibility by that factor. Its quotient uniquely
determines v(q), hence the entire q (positive integers have a unique finite
little-endian word ending1). The zero-run quotient must be positive; the
one-run quotient is at least1. Then evaluate the exact constant term by
the affine recurrence and compare it to the requested B0.

This procedure is complete over all finite p,q for each actual bank entry.
It may list several presentations of the same family; no uniqueness or
isomorphism quotient is claimed. Its completeness does not extend to two
independently growing runs or to words selected afresh at each k.

## Exact requested families and a positive control

The standalone [family source](../../../../04-computation/nonzero_preparation_20261103/family.py)
uses the 34 actual entry choices listed in
[surgery_ready.py](../../../../04-computation/nonzero_preparation_20261103/surgery_ready.py).
It tests the current family `756*2^k-83` and all fourteen formal descendants
through three requested2/4 steps, with A multiplied by3 and B0 sent to
`3B0+d`. No family has a template of either ready gap in this bank and
one-variable-run language. All prefix and suffix lengths are covered by (4).

The inherited common family `90*2^k+5` is a positive control: the exact
inverse recovers both gaps. The small-state recurrence is independently
compared to the inherited six-state formula on every binary word through
length6. Every slope-divisibility candidate is also checked at k=1,2,3,5
using the original word recurrence. The all-k result comes from the affine
derivation, not those checks.

Normal replay passes604,826 always-active gates. Receipt SHA256:
`7bcff67f4d7c675c8b7f94b009892c45ec41bc7b6de0a524e31158b2c44798a4`.

```text
python3 04-computation/nonzero_preparation_20261103/family.py
python3 -O 04-computation/nonzero_preparation_20261103/family.py --output /tmp/nonzero_family_optimized.json
```

The positive finite preparations found at H1429 remain valid. This theorem
explains why those isolated selector hits do not automatically promote to
a fixed-core, one-variable-run preparation of the whole original clock
family. The next construction must change the core bank, use a genuinely
adaptive word, or retain more than one changing run.
