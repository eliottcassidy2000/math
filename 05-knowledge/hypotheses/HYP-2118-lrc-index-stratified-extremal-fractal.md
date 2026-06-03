---
id: HYP-2118
status: OPEN proof route; S585 gives exact translated-AP index evidence
source: codex-2026-06-03-S585
related: [HYP-2117, HYP-2116, HYP-2115, HYP-2114, HYP-2113, HYP-2112, HYP-2092, HYP-2091, HYP-2079]
---

# HYP-2118: LRC extremal families stratify by recursive summand index

## Claim

The additive extremal family should not be treated as a single AP point.  It is
an indexed recursive family whose floor face is the AP and whose translated
faces preserve the same 4-term shadow while moving fold information into hidden
summand shells.

For the interval family

```text
I(k,q) = {q, q+1, ..., q+k-1},
```

the pair-sum multiplicity profile is translation-invariant, so the 4-term
count is fixed.  The visible 3-fold count is the clipped prefix

```text
F(k,q) = #{0 <= i < j < k : i+j <= k-1-q}.
```

Increasing the index `q` does not destroy the fold pattern.  It exports it
from visible folds in `I(k,q)` to hidden pair-sum shells above the row.  Adding
the appropriate virtual sum node recovers the fold one layer higher, exactly
as in HYP-2115.

Thus the useful proof coordinate is a recursive index vector, not raw additive
energy:

```text
(parity lane, C-gcd shell, visible fold clip, hidden sum shell,
 dyadic debt depth, endpoint-owner labels).
```

## Evidence

S585 audits the translated AP row at `k=13`, the LRC `n=14` additive model.

The exact table is:

```text
q=1..13:
  visible fold counts = 36,30,25,20,16,12,9,6,4,2,1,0,0
  formula mismatches = []
  4-term counts = 125 throughout
  M/delta rises from 1.000 to 4.789
```

The strongest hidden virtual node always has multiplicity `6`, while its shell
index walks outward:

```text
q=1:  s=14, shell=1
q=6:  s=23, shell=5
q=13: s=37, shell=12
```

Adding that hidden node lowers the augmented maximin in every tested translated
row, showing real virtual-fold pressure:

```text
q=1:  M=1/14, M(+s)=1/15
q=8:  M=0.2857, M(+s)=0.2286
q=13: M=0.3421, M(+s)=0.2600
```

This sharpens HYP-2114/HYP-2115: raw 4-term energy cannot see `q`; visible-fold
and hidden-shell indices can.

S585 also records the distinct-pair summand-completion width recursion for an
interval support:

```text
w_{r+1} = 2*w_r - 3,
w_r = 2^r*(k-3)+3.
```

For `k=13` the widths are

```text
13 -> 23 -> 43 -> 83 -> 163 -> 323.
```

So the precise fractal object is an iterated clipped tent profile: each
summand-completion level repeats the same pair-sum clipping mechanism at a
larger dyadic scale.

## Proof Route

Use the index vector as a recursive proof address.

1. If the visible fold clip is nonzero, route to Lemma B / fold-shield /
   delta-clock machinery.
2. If visible folds vanish but hidden shells remain, keep hidden sum nodes as
   labels and route to Lemma A / gap-structural randomness.  Do not use raw
   4-term energy as a certificate.
3. If a hidden node becomes real under completion, hand it to the visible
   fold branch.
4. If the `C=2n-1` shell is nonunit, use HYP-2092's gcd-stratum descent and
   endpoint-owner labels.
5. If the dyadic depth increases, use the HYP-2079/HYP-2073 conserved product
   law: gap halves, debt doubles, and the packet address survives.

This complements both HYP-2116 and HYP-2117.  HYP-2116 supplies the `x2`
vertical 2-adic index; HYP-2117 supplies the binary IFS/doubling-map view of
translated APs; this hypothesis supplies the clipped pair-sum/summand-shell
index.  Together they suggest an induction on the lexicographic proof address
rather than on speed size alone.

## Assumption Challenge

Tournament vertices need not be runners.  For this branch the useful vertices
are index lenses: `Phi` gap values, visible fold clips, `C`-gcd shells, hidden
sum shells, dyadic debt depths, and raw 4-term energy as a negative control.

The quotient preserves: where fold mass lives in the recursive summand
completion and which proof gate can consume it.

The quotient destroys: individual runner labels and exact endpoint owner
provenance until they are reattached as fibre labels.

The challenged assumption is that extremality is one flat family.  The AP is
better read as the index-zero face of a self-similar summand completion tree;
translated APs are higher-index faces with the same 4-term shadow but safer
visible geometry.

## Files

- `04-computation/lrc_index_stratified_extremal_fractal_s585.py`
- `05-knowledge/results/lrc_index_stratified_extremal_fractal_s585.out`
- `07-reflections/lrc-index-stratified-extremal-fractal-s585.md`
