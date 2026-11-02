# Independent audit: ordinal profile dichotomy and restricted strong zero rigidity

**Verdict: ACCEPT, all orders in the stated restricted domains.** Reviewed
`zero_rigidity/PROFILE_DICHOTOMY.md` and its ordinal-zero dependency
`zero_rigidity/ORDINAL_ZERO_CLASSIFICATION.md`; exact versions and supporting
certificate inputs are pinned in `PROFILE_DICHOTOMY_SOURCE_PINS.json`.

The proof establishes profile comparability for completions of a deleted
vertex in `TT_k[C3,...,C3]`, and scalar zero rigidity after cyclic substitution
**only when the replacement vertex retains its original orientations to the
other two outer modules**. It does not establish the arbitrary strong-quotient
zero-rigidity conjecture. The actual Q4 hostile is a counterexample to extending
the *profile dichotomy*, not a counterexample to scalar zero rigidity.

## 1. Independent reconstruction of the all-order proof

### Correct carrier and positivity

If `f_W(j)` counts ordered covers by j nonempty directed paths, selecting
j occupied colors in increasing color order gives

```
P_W(m) = sum_j binom(m,j) f_W(j).
```

There is no missing j! factor here: ordering the paths is already part of
f. In the colored odd-cycle-family expansion, each selected directed odd
cycle is monochromatic, each unused vertex has an independent color, and
the cycle receives its inherited OCF weight 2. Thus a vertex-disjoint cycle
family contributes exactly

```
2^number_of_cycles * m^(n - sum_C (|C|-1)).
```

This proves the displayed nonnegative monomial expansion. Conversion by
`m^e=sum_j j! S(e,j) binom(m,j)` preserves nonnegativity. In particular,
retaining only cycle families confined to parts of any vertex partition
gives `f_W >= f_(W1=>...=>Ws)` entrywise, regardless of cross arcs.

Two distinct positivity statements are used correctly: partition inclusion
first gives a nonnegative *monomial* difference; the later finite-core
comparison only gives a nonnegative *ordered-profile* difference and is
propagated by the positive ordinal-join law. The proof does not silently
identify those two cones.

The ordinal-join coefficient

```
(p+q-j)! / [j! (p-j)! (q-j)!]
```

is correct. Forget local path orders, select j A-to-B component matches,
concatenate them, and order the resulting global paths. The p! and q!
normalizations cancel to give precisely that coefficient. Consequently,
joining with C3 multiplies a first nonzero coordinate-2 difference by
`2*f_C3(1)+f_C3(2)=2*3+6=12`, and does not create a coordinate-1 difference.

### Exhaustive low-state and positive-core split

For root-pair bits `(a,b)` with `u->w`, reversal changes the bits to
`(1-b,1-a)`, reverses the intact-block order, and changes each intact
indicator t to 3-t. Reversing every directed path preserves each ordered
cover count. This validates every use of duality in the case split.

- For the cyclic pair `(1,0)`, all original triangle families survive.
  Any wrong earlier or later incidence creates a new cross-block triangle.
  Its contribution makes H strictly larger unless the cut is canonical.
- For middle pair `(0,1)`, canonical outside incidences give a low state.
  Otherwise there is a wrong later incidence, after reversal if needed;
  retaining that intact triangle gives one of the three first six-vertex
  cores.
- For pair `(0,0)`, any selected earlier vertex gives one of the other
  three six-vertex cores. If none exists, failure of the declared monotone
  pattern is exactly the existence of later blocks `i<j` with `t_i>0`
  and `t_j<3`; retaining those two blocks gives one of nine nine-vertex
  cores. Pair `(1,1)` is the reversed case.

Thus the fifteen-core list is exhaustive, not inferred from testing a
maximum k. Deleting intervening ordinal blocks preserves every retained
cross orientation; rotations of each intact C3 act transitively on both
its one-element and two-element subsets. The indicators therefore retain
all information needed for these cores.

For the low states, no mixed intact block gives an ordinal collection of
k-1 triangles and three transitive singleton vertices. With
`p(m)=m^3+2m`, its exact difference is `-2m*p^(k-1)`.
One mixed intact block gives the unique strong four-vertex tournament,
with colored polynomial `m^4+4m^2`, and the surviving root pair is a separate
transitive summand. The exact difference is then `-4m^2*p^(k-2)`.
Both have nonpositive ordered-profile coefficients and strictly negative
H difference. The second case requires k>=2, as the manuscript states.

For a positive core, partition the actual completion into that core and
all remaining intact triangles. Partition inclusion followed by the
positive ordinal-join law gives entrywise domination over the original
ordinal tournament. A positive first coordinate propagates strict H;
the sole zero-first-coordinate six-vertex core instead propagates its
positive second coordinate by successive factors of 12. The intermediate
partition inequality can only add nonnegative coordinates. This proves
exactly the three alternatives of Theorem A, including the uniqueness of
full-profile equality.

### Finite cores are exact local certificates

The independent implementation reconstructed all fifteen core tournaments
from their cut descriptions, checked the listed literal masks, and counted
all coordinates by a separate algorithm. The six-vertex differences are
exactly

```
(6,60,144,96,0,0), (4,56,144,96,0,0), (0,24,72,48,0,0),
(6,60,144,96,0,0), (16,128,288,192,0,0), (24,168,360,240,0,0).
```

The nine nine-vertex differences have the claimed common lower envelope

```
(4,648,8088,35040,67920,60480,20160,0,0),
```

attained at intact indicators (1,2). All seventeen complete profiles in
`profile_leafs.txt`—the fifteen cores and two ordinal controls—match, not
just their signs or first coordinates. The existing literal-permutation
certificate is therefore independently confirmed rather than trusted as
an unexplained numerical table.

During audit, the producer corrected one harmless metadata label from
`parent_masks` to `core_masks`: that field contains the completed core U,
not the deleted card. All masks and numerical vectors are unchanged. The
final source pins and both independent replays use the corrected label.
The final manuscript additionally appends a producer whole-theorem replay
paragraph; its preceding 9,497 bytes, containing the theorem and proof,
remain exactly identical to the version originally reviewed.

## 2. Why cyclic strongification is strictly zero-rigid here

With A and B fixed and nonempty, maximal module-run contraction expresses
`H(C3[M,A,B])` as a nonnegative linear functional of the ordered profile of
M. The run word `M,A,B` gives a positive coefficient of `f_M(1)`;
`M,A,B,M` gives a positive coefficient of `f_M(2)`. Both use ordinary
Hamilton paths in A and B, which exist in every nonempty tournament.

A low completion is weakly smaller in every coordinate and strictly smaller
in coordinate 1, so its ambient H strictly decreases. A high completion is
weakly larger in every coordinate and strictly larger in coordinate 1 or 2,
so its ambient H strictly increases. Only the canonical completion gives
equality. This supplies Theorem B for arbitrary nonempty A,B, not just the
small test modules.

The hypothesis retaining every x-to-A and x-to-B orientation is indispensable
to this proof: without it there is no unchanged one-module substitution
functional. No conclusion for arbitrary outside cuts follows. The original
and edited ambient tournaments are strong because the outer quotient is C3,
but strongness alone does not restore the lost functional hypothesis.

## 3. Ordinal scalar-zero dependency

The companion theorem's explicit capacities have the correct endpoint
orientation and normalization. In the cut formula the target increment is
`3L=2h`, and all terms are nonnegative. Two out-of-order intact indicators
already exceed this budget: two mixed blocks cost at least `2L+5L/3`, while
an extreme indicator makes the cross term at least `5L`. Thus intact
indicators must be monotone.

The cyclic pair spends all 3L immediately and forces canonical outside
incidences. Pair00/11 cannot hit the target. The middle pair must make exactly
one wrong incidence in an adjacent triangle, costing 2L plus the triangle's
mixed internal cost L; nonadjacent wrong incidences cost at least 6L. Hence
the only noncanonical zeros are the stated previous/next adjacent defects.
The k=1 case is separate and correct. Counts of labelled restoring cuts are
1 for k=1, 4 at an endpoint root, and 7 at an interior root.

The scalar-zero theorem is a closed all-order theorem for a nonstrong ordinal
quotient. It is neither a strong-quotient counterexample nor a proof that all
strong quotients are zero-rigid.

## 4. Actual strong-Q4 mixed-profile hostile

The independent checker rebuilt pair-bit code4, its four C3 substitutions,
the deletion of vertex0, and all three literal cuts. It verified the quotient,
original tournament, deleted parent, and all completed children are strong.
Every whole ordered-profile difference agrees with the producer receipt.

For cut996 the original H is 241137 and the new H is 237411. The difference is

```
(-3726,146796,1919304,10247376,31474800,61699680,
 79591680,66286080,32659200,7257600,0,0).
```

It has both signs. Thus replacing the ordinal hypothesis in Theorem A by
“arbitrary strong quotient” is false. Yet its H difference is -3726, not zero,
so it does not refute the strong-quotient scalar-zero conjecture.

As an independent literal context check, putting this module into
`C3[M,K1,K1]` reverses the H-change sign:

```
3*(-3726)+146796 = 135618 > 0.
```

The full ambient H was separately counted by Hamilton-subset DP; this was not
merely evaluation of the displayed expression. Cuts867 and995 also reproduce
exactly: their H differences are +13002 and +23058 but several higher profile
coordinates are negative. All three are actual mixed-profile controls, not
formal capacity networks and not zero witnesses.

No all-parameter claim about the separate two-parameter strongification
search is certified here. Its finite-box sign observations are not needed
for either accepted all-order theorem.

## 5. Independent implementation and replay

`audit_profile_dichotomy.py` imports no producer code. It reads the compact
core/leaf/hostile receipts as asserted data, but reconstructs the graphs and
computes expected values independently:

1. Backward last-endpoint DP computes H for every induced subset.
2. An unordered path-partition recurrence forces the least vertex into its
   chosen part, eliminating part-order multiplicity. Weight each part by its
   independently counted H; multiply the r-part total by r! to obtain the
   ordered cover coordinate.

This is different from the producer's forward permutation/separator
polynomial recurrence. Six-vertex controls additionally agree with literal
word/separator enumeration. Every gate raises explicitly and remains active
under Python optimization.

Completed exact scope in each execution:

- 15 finite cores and 17 full leaf profiles;
- all 228 block-profile/root cases for ordinal k=1..3;
- all 836 literal labelled cuts for k=1..3, including the exact zero counts;
- 108 literal outer cyclic-context checks using singleton, transitive-pair,
  and directed-triangle side modules;
- 669924 purely combinatorial bounded-core case-reduction states for k=1..8,
  reaching all fifteen allowed core types;
- all three actual strong-Q4 hostile vectors and their literal singleton
  cyclic-context responses.

The 669924-state check tests the stated case split; it is **not** an expanded
strong-quotient census or a proof by finite extrapolation. The universal
quantifiers are supplied by the audited structural proof above.

Run from the independent directory, preserving sibling `zero_rigidity/` data:

```sh
python3 -B audit_profile_dichotomy.py > replay.normal.json
python3 -B -O audit_profile_dichotomy.py > replay.optimized.json
cmp replay.normal.json profile_dichotomy.normal.json
cmp replay.optimized.json profile_dichotomy.optimized.json
sha256sum -c PROFILE_DICHOTOMY_SHA256SUMS
```

Normal and optimized outputs are byte-identical. The exact receipt digest,
script digest, and reviewed-source versions are bound by the dedicated
manifest and source-pin file. Existing independent stability files and their
aggregate manifest were not changed. No Git mutation, fetch, commit or push
was performed.
