# Dichotomy Recursion Mode Atlas

This pass went back through the repo's recurring binary splits: odd/even,
positive/negative, addition/multiplication, `+1` versus `/2`, `*2` versus `+2`,
and the broader sum/product/fraction/recursion grammar.

The main consolidation is that these are not decorative analogies.  They are
quotient guardrails.  Each split says which label survives and which coordinate
gets destroyed if we turn it into a scalar too early.

After the concurrent HYP-3002/LTI-152 curried packet functional tower landed,
this atlas is best read as its mode-switch companion: HYP-3002 says a proof
quotient is a partial evaluation, while HYP-3004 says these dichotomies are the
recursion modes that make the partial evaluation theorem-safe or lossy.
After the S171 rebase, HYP-3003/LTI-153 is the narrower summand/multiplicand
Farey-basis merge; this atlas treats that add/multiply split as one subcase
inside the wider recursion-mode switchboard.

## What Changed

The strongest older sentence is from the orbit-functor work:

```text
odd/multiplicative = orbit propagation
even/doubling      = projection defect
additive/fold      = denominator pincer
CRT                = where the projection defect lives
```

The signed-LRC work supplies the complementary sign law:

```text
positive/negative is observer-blind but pair-visible.
```

Signs do not change `M(S)`, but they select which pair clocks become sums.  So
the sign dichotomy is a laboratory for the additive pair-sum face, not a witness
machine by itself.

The treebolic `+2`/`*2` work supplies the geometry:

```text
*2 = vertical descent in the 2-adic tree
+2 = horizontal line motion weaving through 2-adic depths
```

The Collatz work supplies the affine/valuation version: `+1` is the bounded
defect after the right logarithm, while division by powers of two is the
valuation drop that counts recursion depth.

The recent HYP-2998/HYP-2999/HYP-3000/HYP-3001/HYP-2984/LTI-151 stack adds the
current LRC interpretation.  Farey sum is theorem-safe because it is affine in
`M`; Farey product is a scheduler only after the unit-excess gate; smoothing is
chosen only after packet route.

## S170 Script Result

`04-computation/dichotomy_recursion_mode_atlas_codex_s170.py` compares twelve
mode carriers with Tournament Analysis over proof carriers rather than runners.

The transitive path is:

```text
additive_pair_sum_face
> sign_cut
> triune_fraction_recursion
> parity_fold
> zeckendorf_path_normal_form
> smoothing_switchboard
> dyadic_doubling_tree
> farey_product_scheduler
> multiplicative_unit_orbit
> farey_sum_affine_check
> collatz_affine_halving
> plus_two_line_motion
```

That ranking is not a universal hierarchy.  It reflects the current LRC14 proof
state: pair-sum owners are the most concrete retained coordinate; signs expose
them; triune recursion tells us how to repair scalar projections.

## New Working Rule

Every dichotomy should be recorded as:

```text
binary split + preserved predicate + destroyed coordinate + recursion law.
```

For LRC packet work, add fields:

```text
parity_block
sign_cut_status
additive_pair_sum_owner
multiplicative_unit_orbit
recursion_boundary_state
smoothing_route
```

This is the cleanest next test for HYP-2963: a hard packet should have one
primary mode and named side-channel debts.  Multiple primary modes are not a
failure; they are exactly where a zipper theorem should live.
