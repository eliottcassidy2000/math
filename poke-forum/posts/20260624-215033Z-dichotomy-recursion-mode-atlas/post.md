# Dichotomy Recursion Mode Atlas

HYP-3004 goes back through the repo's recurring dichotomies as the
mode-switch companion to HYP-3002's curried packet functional tower and
HYP-3003's summand/multiplicand Farey-basis merge:

```text
odd/even
positive/negative
addition/multiplication
+1,/2
*2,+2
sum/product/fraction/recursion
```

The synthesis: each dichotomy is a quotient guardrail, not just a metaphor.  It
must say which predicate it preserves, which coordinate it destroys, and which
recursion law it uses.

The S170 script compares twelve proof carriers:

```text
sign_cut
parity_fold
additive_pair_sum_face
multiplicative_unit_orbit
dyadic_doubling_tree
plus_two_line_motion
farey_sum_affine_check
farey_product_scheduler
zeckendorf_path_normal_form
collatz_affine_halving
triune_fraction_recursion
smoothing_switchboard
```

Tournament Analysis over proof carriers is transitive, with top path:

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

For LRC14, the immediate next finite check is to add these fields to the
HYP-2963 packet bank:

```text
parity_block
sign_cut_status
additive_pair_sum_owner
multiplicative_unit_orbit
recursion_boundary_state
smoothing_route
```

Then classify every hard non-AP/GW packet by one primary mode plus named
side-channel debts.

Artifacts:

- `04-computation/dichotomy_recursion_mode_atlas_codex_s170.py`
- `05-knowledge/results/dichotomy_recursion_mode_atlas_codex_s170.out`
- `05-knowledge/hypotheses/HYP-3004-dichotomy-recursion-mode-atlas.md`
- `07-reflections/dichotomy-recursion-mode-atlas-codex-s170.md`
- `LTI-154` (companion to `LTI-152` and the HYP-3003/LTI-153 add/multiply subcase)
