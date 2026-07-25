# LRC14: marginal spectra are addresses; simultaneous geometry is the proof object

Source: codex-2026-07-25, the affine-chord scout, THM-2258, and the
THM-2263/THM-2145/THM-2270 synthesis.

Status:

- **PROVED + VERIFIED-EXACT:** THM-2258 independently excludes the scalar
  depth-three profile `(3,4,5)`. Its literal five-comb residual is at least
  `317/1155`, while THM-2243 requires at most `8907541/62377224`.
- **FINITE-EXACT / SUPERSEDED IN CONSEQUENCE:** the affine-chord companion
  closes all `1680` mixed projective direction triples in its safe Bellman
  enlargement, leaving `48` triples in one projective direction class.
  THM-2250 later proves the stronger literal equality of normalized cores;
  THM-2257 and THM-2258 close the whole profile. The companion is retained
  for the root-incidence mechanism and its quotient audit, not as a current
  theorem dependency.
- **PROVED incoming connection:** THM-2263 applies the same compatibility
  lesson to the remaining depth-one owner ledger. Gapwise pair extrema must
  obey one multiplicative ladder, and the unique strict worst profile is
  `(1,3,5)`, not three independent copies of the global pair maximum.
- **PROVED invariant from canonical THM-2145, with a simultaneous
  strengthening in THM-2270, not a closure:** bounded relations cross every
  labelled `6+7` cut, and one larger relation can cross all `1716` cuts at
  once. This is more than relation rank, but it does not exclude a row.

Nothing in this reflection decrements the current ledger beyond the
independent THM-2257/THM-2258 closure. LRC(14) remains open with the 165
first-depth-one scalar rows and the routed non-scalar frontiers.

## 1. The affine chord is the one-step object

For one thirteen-unit core `u`, normalize a thirteen-root fibre by the odd
guard `H`. Conditional on the future danger bit, the current danger roots
are:

```text
future bit 0:  {beta,beta+a},
future bit 1:  one endpoint of that pair,

a=-H u^(-1) mod 13.                                  (1)
```

The direction `a` is fixed by the core. The affine anchor `beta` is not: it
depends on the parent phase and previous root history. This is the exact
coordinate that an arbitrary-coupling Bellman discards.

The safe enlargement used by the companion fixes the three directions
throughout the recursion but allows every anchor and endpoint orientation
independently at every history. It is therefore an upper bound on every
actual arithmetic root process. Exact thirteen-root counting then gives:

```text
oriented direction triples                 12^3=1728
mixed projective classes closed                  1680
one projective direction class left                48
worst mixed bound                  567556891/5676327384
margin below 961/6930              36232889557/936594018360. (2)
```

The projective quotient is forced, not cosmetic. A two-point chord with
direction `a` is the same family as one with direction `-a` after moving its
anchor. The hostile orientation probe `(2,2,11)` is the minimal warning:
`11=-2 mod 13`, so treating these as different chord types manufactures a
false distinction. The exact state is

```text
F_13^*/{+-1}, together with the affine anchor sidecar. (3)
```

### The failed simplification

Setting `beta=0` would turn the chord into a direction-only object and makes
several attractive incidence inequalities appear stronger. It is invalid:
fixing a direction does not fix which two roots are selected, and the
continuation may choose a different translate at every parent. The first
failed implication is

```text
fixed normalized core
  -/-> fixed affine root chord.
```

The strongest safe survivor is instead: fixed core gives a fixed
**projective direction**, while anchors remain adversarial.

## 2. THM-2258: near-maximal charges force a large collision tax

Let

```text
A_i=C_H intersect D_(q_i),          c_i=measure(A_i),
p=measure(C_H minus union_i D_(q_i)).                (4)
```

THM-2243 gives `p<=B=8907541/62377224` for a hypothetical `(3,4,5)`
survivor. The union bound therefore forces

```text
sum_i c_i>=5/7-B=5092517/8911032.                    (5)
```

This looks like a marginal optimization problem. It is not. THM-2080 first
turns (5) into a finite address:

```text
four global largest distinct-ratio charges:
  5/42, 9/77, 4/35, 4/35;

forced fifth charge:
  at least 122343163/1143582440;

tail bound:
  every reduced product ab<=50;

exact bank:
  13 ratios, 1287 five-subsets, 15 marginal survivors. (6)
```

Restoring the literal five-way union destroys all 15 rows. Define

```text
Omega=sum_i c_i-measure(union_i A_i)
     =integral_(C_H)(owner multiplicity-1)_+.         (7)
```

The exact atlas gives

```text
p>=317/1155,
Omega>=29/220,
p-B>=150561431/1143582440.                           (8)
```

At the worst row the marginal surplus over (5) is only
`183527/1143582440`, but the collision tax is `29/220`. The arithmetic
ratios that make every edge individually good force the union to be very
inefficient.

This is the structural meaning of THM-2258:

```text
near equality in a marginal relaxation
does not imply near equality in its simultaneous realization.           (9)
```

## 3. THM-2263: the pair extrema have to live on one ladder

THM-2255's global ramified pair maximum is `25/1183`, uniquely at reduced
ratio `1:169`. Applying that maximum independently to all three blocker
pairs forgets the cocycle

```text
c_3/c_1=(c_3/c_2)(c_2/c_1).                          (10)
```

The three edge ratios are not independent labels. THM-2263 retains their
valuation gaps and proves the complete parity-sensitive pair spectrum. The
compatible strict extremum occurs at

```text
(lambda_1,lambda_2,lambda_3)=(1,3,5),
(c_1,c_2,c_3)=g(1,169,28561).                        (11)
```

Two adjacent edges attain the gap-two maximum; the outer edge must have
gap four. The exact compatible ledger improves the strict exclusive-owner
mass to

```text
15041431/197927730,
```

and its labelled expiration image to

```text
15041431/70270200>1/7.                               (12)
```

The missing implication is now sharply located. THM-2261 proves that raw
owner expiration is surjective, so no one-core future carrier contains the
image. The next useful object must retain the named root sheets or the
handoff target supplied by the full scalar cover. A better mass comparison
without that target cannot close a row.

## 4. THM-2145/2270 beyond relation rank: no balanced cut annihilates the space

THM-2264 is an exact duplicate of THM-2144 and is superseded. THM-2265 is
a valid but weaker duplicate of canonical THM-2145 and must not be used as
the dependency. THM-2164 already proves unconditionally

```text
dim W_105>=2
```

for the distinct zero-safe rows. THM-2145 contributes a different
coordinate. Put `W=W_584` and choose a basis from bounded relation
generators

```text
m^(1),...,m^(r),             ||m^(j)||_infinity<=584.
```

Define the thirteen column charges

```text
b_i=v_i(m_i^(1),...,m_i^(r)) in Z^r.                 (13)
```

Every basis row is a relation, so

```text
sum_i b_i=0.                                         (14)
```

THM-2145 says exactly that

```text
sum_(i in F)b_i !=0             for every |F|=6.     (15)
```

Indeed, if the vector sum vanished, every element of `W` would have zero
weighted transfer across that cut, contradicting the height-`584` crossing
relation supplied for `F`. By (14), no seven-subsum vanishes either.
THM-2270 further proves that one height-`504576` relation crosses all
`1716` cuts simultaneously; its support is at least eight.

Rank alone does not imply (15): a high-rank relation space may lie inside
the kernel of one cut functional. The new invariant is a zero-sum vector
configuration with no balanced zero subsum. Its useful finite side is that
the coefficient directions

```text
(m_i^(1),...,m_i^(r))
```

come from a bounded integer alphabet, while the speeds provide radial
scales. The next probe should classify low-rank direction configurations
or prove that some six-subsum must vanish under the additional scalar
valuation structure.

This is a hypergraph/cut object, not naturally a tournament: the intrinsic
observable lives on `6+7` cuts, and forcing it into pairwise orientations
would discard the partial-sum predicate.

## 5. Reusable research move

The successful move across THM-2258 and THM-2263 is:

```text
1. derive a near-extremal sum of local capacities;
2. rank the exact local spectrum and use a tail bound;
3. reduce every local label to a finite high-capacity bank;
4. restore the simultaneous constraint:
     literal union, ratio cocycle, shared root anchor, or cut functional;
5. measure the compatibility debt instead of optimizing each edge again. (16)
```

Use it when a proof is within a small numerical margin and its local
inequalities have rigid equality cases. Do not use it when labels are
genuinely independent, the target factors through the marginal data, or no
tail bound makes the high-capacity spectrum finite.

The affine-chord scout supplies the caution paired with the move: retain
every translation/anchor coordinate needed by the next operation. A
direction spectrum without its affine anchor is another marginal quotient,
not the simultaneous object.

## Artifacts

- `01-canon/theorems/THM-2258-depth-three-uniform-five-charge-spectrum-overlap-exclusion.md`
- `04-computation/lrc14_depth_three_uniform_five_charge_overlap_thm2258.py`
- `05-knowledge/results/lrc14_depth_three_uniform_five_charge_overlap_thm2258.out`
- `04-computation/lrc14_depth_three_affine_chord_projective_companion.py`
- `05-knowledge/results/lrc14_depth_three_affine_chord_projective_companion.out`
- `01-canon/theorems/THM-2263-thirteen-adic-gap-pair-spectrum-and-profile-sharp-owner-floor.md`
- `01-canon/theorems/THM-2145-two-block-spectral-crossing-and-6-plus-7-carry.md`
- `01-canon/theorems/THM-2270-simultaneous-balanced-cut-relation-and-six-uniform-orientation.md`
