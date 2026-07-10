# The difference-primitive case collapses the conclusion, not the wall

*kind-pasteur-2026-07-10-S127. Owner: "prove TightRigidity for the difference-primitive case." I proved
the piece that is provable — the conclusion collapses to `{1,…,13}` — and I am reporting honestly that the
restricted rigidity is still the conjecture. This note explains why the restriction moves nothing.*

---

## The hope, and why it is misplaced

`TightRigidity` says: `μ(S) = 0 ⟹ S is a dilated interval c·{1,…,13}`. It has a *dilation freedom* — the
scale `c` is a free parameter — and one might hope that pinning `c` (by restricting to difference-primitive
families, where `gcd` of the differences is `1`) removes enough freedom to make the statement provable.

It does not, and the reason is a clean separation of the two halves of the statement:

- **The dilation freedom lives entirely in the *conclusion*.** "Dilated" = "`= c·{1,…,13}` for *some* `c`".
- **The hard content lives entirely in the *hypothesis→shape* implication:** `μ = 0 ⟹ the absolute speeds
  form an arithmetic progression at all`.

Primitivity only touches the first. It fixes `c = 1`, collapsing "dilated" to the single set `{1,…,13}`.
The implication `μ = 0 ⟹ AP-structure` — the wall — is untouched, because AP-ness is scale-free: it is
exactly as hard to prove a primitive family is `{1,…,13}` as to prove a general family is `c·{1,…,13}`.

## What I proved (kernel-pure)

In `LRCTightRigidity.lean`, foundational axioms only:

```lean
def Primitive (v : Fin 13 → ℤ) : Prop := ∀ d : ℤ, 2 ≤ d → ¬ (∀ i, d ∣ v i)

theorem dilated_primitive_eq_range (hprim : Primitive v) (hdil : DilatedFamily v) :
    (univ.image fun i => |v i|) = Finset.Icc 1 13
```

The collapse. If `v` is primitive and dilated, the scale is forced to `c = 1`. Proof: every `|vᵢ| = c·kᵢ`
sits in the dilated image, so `c ∣ |vᵢ| = |vᵢ|`, hence `c ∣ vᵢ` for *every* `i`; primitivity forbids
`c ≥ 2`, so `c = 1` and the image is `(Icc 1 13).image (1·) = Icc 1 13`. This genuinely *sharpens*
`TightRigidity`'s conclusion on the primitive class: the only primitive tight family the rigidity permits
is the arithmetic progression itself, no dilate.

And the honest boundary marker:

```lean
def PrimitiveTightRigidity : Prop :=
  ∀ v, (∀ i, v i ≠ 0) → Primitive v → volume (safePeriod v) = 0 →
    (univ.image fun i => |v i|) = Finset.Icc 1 13   -- ← still the conjecture, NAMED not proved

theorem primitiveTightRigidity_of_tightRigidity : TightRigidity → PrimitiveTightRigidity
```

The last theorem is the point stated in reverse: the full rigidity *gives* the primitive rigidity via the
collapse, confirming that restricting to primitive families throws away no essential content. It is a
restatement of the wall, not a way over it.

## Why the restriction still implies LRC(14)

The residual class carries difference-primitivity as one of its nine defining clauses. So `μ = 0 ⟹
{1,…,13}` on the primitive class is enough to run the same measure-floor argument
(`safeMeasureFloor_of_tightRigidity`): a residual family is scale-gapped (ratio `> 13`), so it is *not*
`{1,…,13}` (whose ratio is exactly `13`), so its safe set cannot have measure zero. `PrimitiveTightRigidity`
is therefore `≥ LRC(14)` — precisely as hard as the unrestricted rigidity.

## The search evidence

Over all `13`-subsets of `[1,N]` with `Vmax = N` and `gcd = 1` (primitive), `N ≤ 20`: the only measure-zero
family is `{1,…,13}` itself. No primitive non-AP tight family exists in range — consistent with
`PrimitiveTightRigidity`, not a proof of it. (The dilates `c·{1,…,13}` with `c ≥ 2` are all *non*-primitive,
so they never counterexample the primitive statement — which is exactly why the collapse is clean.)

## The lesson

A restriction that removes a *parameter* from the conclusion (here, the dilation scale) can look like
progress on the *theorem*, but it only trims the conclusion. The difficulty of an implication lives in the
map from hypothesis to structure, and a scale-free structural claim (arithmetic progression) is invariant
under exactly the kind of normalization that primitivity provides. When asked to prove a restricted case,
the honest first question is: *does the restriction touch the hard implication, or only the shape of its
output?* Here it was the latter. I proved the shape-collapse and named the untouched wall.

*Files: `LRCTightRigidity.lean` (`Primitive`, `dilated_primitive_eq_range`, `PrimitiveTightRigidity`,
`primitiveTightRigidity_of_tightRigidity`), `lrc14_primitive_tight_search_kps_S127`. Continues
[[the-residual-is-one-measure-floor-kps-S127]]; the wall is still
[[describing-the-wall-the-scale-gap-is-the-separator-kps-S127]].*
