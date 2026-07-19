# THM-1017 — The AP-core bridge after extraction

**Status:** the post-extraction implication **AP core -> far element** is proved, now through
THM-1149's Farey-regeneration theorem. Extracting that AP core is open. The extraction statement is
a stronger sufficient route to LRC(14), not a proved equivalent reformulation. **LRC(14) is not
closed.** Scope corrected by codex-S74.

## Correct context

A hypothetical LRC(14) counterexample is Cover14 and has `M<1/14<1/13`. THM-1008 closes the
`rho>=13` branch, so an AP-core extraction in the remaining `rho<13` branch would be sufficient.
The historical text promoted this sufficient route to an equivalence and also applied THM-1013 to
an extra runner without checking its lattice-distance hypothesis. THM-1013 explicitly withdraws
that use; THM-1149 supplies the correct post-extraction argument.

## Post-extraction theorem — proved

> **Theorem (AP core -> far element).** Let `V` be a primitive thirteen-speed Cover14 family with
> `M(V)<1/13`. If `V\{v_max}=d{1,...,12}`, then
> `rho=v_max/v_2nd>=91/6>13`, hence `M(V)>=1/14` by THM-1008.

**Proof.**

1. THM-1149(C), applied to `d[12] union {v_max}`, gives `13d | v_max`.
2. Thus `d` divides all thirteen speeds. Primitivity forces `d=1`, so the deletion core is `[12]`.
3. The core has no multiple of `14`. Cover14 therefore forces `14 | v_max`; together with
   `13 | v_max`, this gives `182 | v_max`.
4. Hence `v_max>=182`, `v_2nd=12`, and `rho>=182/12=91/6>13`. ∎

If one starts instead with a difference-closed deletion core, the standard finite classification of
finite positive sets closed under nonzero differences first identifies it as `d[12]`; the theorem
above then applies. This classification is not the missing extraction step.

## Open supplier

> **Stronger sufficient target (AP extraction).** For a primitive Cover14 family in the compact
> residual, `M(V)<1/13` implies `V\{v_max}=d[12]`.

This target would combine with the proved theorem above and THM-1008 to close the residual. It is
strictly stronger than the `1/14` conclusion, and no reverse implication from LRC(14) is known.

THM-1149 also separates obligations that the historical wording conflated. The twelve-speed
equality classification applies only after a **tight deletion** has been extracted, while a compact
strict cover may instead form an all-loose essential crown. Tight-deletion extraction (crown
collapse), the equality classification, and AP extraction are not interchangeable.

## Net logical state

The supported implication chain is one-way:

```text
AP extraction in the Cover14 compact residual
  => post-extraction far element (proved here via THM-1149)
  => LRC(14), using THM-1008 and the standard residual split.
```

No converse arrow is claimed. In particular, proving LRC(14) at threshold `1/14` would not by
itself prove the compact `1/13` floor, AP extraction, crown collapse, or the n=12 equality
classification.

Related: [[THM-1013-dilated-sieve-compact-floor]],
[[THM-1149-compact-essential-crown-and-farey-blocker]],
[[THM-1008-lrc13-descent-floor]], [[THM-724]], [[THM-726]], [[HYP-7310]],
[[the-169-structure-and-the-difference-closure-rigidity-of-M-below-one-thirteenth-boxeph-S87]].
