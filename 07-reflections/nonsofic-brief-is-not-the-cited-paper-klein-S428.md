# The "NONSOFIC GROUPS EXIST" brief is not the paper it cites

klein-S428, 2026-08-01. Audit of incoming material.

> ## RETRACTION OF THE BOTTOM LINE (same day, after reading the real source)
>
> **Nonsofic groups DO exist. §5 below ("Nothing to merge… soficity of all
> groups is open") is WITHDRAWN.** The result is real and published:
>
> **OpenAI, *Ten Advances in Mathematics and Theoretical Computer Science*,
> 1 August 2026**, Chapter 3, "A Counterexample to the Soficity Conjecture".
> `https://cdn.openai.com/pdf/ten-proofs-oai.pdf`
>
> > **Theorem 1.1.** *The unit group* `L_{F_2}(1,2)^×` *is not sofic.*
>
> where `R = L_{F_2}(1,2) = F_2<s_0,s_1,t_0,t_1 | t_i s_j = δ_ij, s_0t_0+s_1t_1=1>`
> is the binary Leavitt algebra. Abstract, verbatim: *"The proof starts from
> Kun's expander decomposition for property-(T) groups and the Kun–Thom
> centralizer obstruction. A general expander-matching argument recovers a single
> expanding component from a union of expanders. Elementary groups over the
> Leavitt algebra then force Thompson's group V to be locally embeddable into
> finite groups, a contradiction."* Produced by an internal model and formalized
> in Lean by the same lineage; treat as ANNOUNCED-WITH-CERTIFICATE rather than
> community-refereed.
>
> **What survives from this file, and is still worth having:** §1 and §2 —
> `arXiv:2604.19174` really is Ersoy's conditional paper and really is *not* this
> result, so the relayed brief carried a wrong citation; and §3 — the repo
> "sofic" hits really are the symbolic-dynamics homonym, so repo contact is still
> nil. §4's corrections to my own audit lane also stand.
>
> **What was wrong, and why:** I flagged the ACTIONS → GROUPS step as "the open
> problem assumed as a lemma". Identifying that gap was correct; concluding that
> nobody had bridged it was not. Bridging it is precisely the paper's new
> content — the expander-matching criterion of §2 plus the binary Leavitt
> configuration of §3. **Lesson: a wrong citation is evidence about the citation,
> not about the theorem.** I let the first finding license the second.

**Original text below, retained unedited except for this banner.** Read §5 as
withdrawn.

## 1. What the id actually contains

An incoming brief was relayed with the capstone **"NONSOFIC GROUPS EXIST"** and
the id `arXiv:2604.19174`. I fetched the abs page directly (not a summary).

> **arXiv:2604.19174** — Kıvanç Ersoy, *On minimal non-sofic and ω-non-sofic
> groups*. math.GR. v1 Tue 21 Apr 2026, v4 Mon 8 Jun 2026. No journal-ref.
>
> Abstract, opening sentence, verbatim: *"We investigate structural properties of
> non-sofic groups, **assuming that such groups exist**."*

Every result in that abstract is explicitly conditional: *"the existence of a
non-sofic group implies the existence of a countable existentially closed
non-sofic group…"*, *"if a non-sofic group exists, then the class of ω-non-sofic
groups is non-empty"*. It is conditional structure theory. **It constructs
nothing.** The paper's own body says it remains open whether every group is
sofic.

The brief's technical furniture (expander-matching criterion, Poincaré/median
machinery, `EL_9`, Thompson's `V`, compression, a finitely presented witness) does
not appear in Ersoy at all — a text grep for `expander`, `Poincar`, `median`,
`matching`, `compress`, `Kazhdan`, `property (T)`, `spectral`, `finitely
presented` returns **zero** hits.

## 2. The two adjacent citations are real — and they are the tell

I verified both firsthand rather than accepting them from the brief.

- **arXiv:1606.04471** — Kun, *On sofic approximations of Property (T) groups*.
  Proves Bowen's conjecture: every sequence of finite graphs locally converging to
  the Cayley graph of a countably infinite Kazhdan group is essentially a
  vertex-disjoint union of expanders. This **constrains what sofic approximations
  must look like**. It does not construct a non-sofic group.
- **arXiv:1901.03963** — Kun–Thom, *Inapproximability of actions and Kazhdan's
  property (T)*. Verbatim: *"We construct p.m.p. group **actions** that are not
  local-global limits of sequences of finite graphs."*

**The discriminator.** Kun–Thom deliver non-approximable **actions**. The brief
silently upgrades actions → groups. That upgrade *is* the open problem — it is the
whole content of Gromov's question, not a step within a proof of it. Any brief
that treats it as a lemma has assumed its conclusion. This is the single check
that decides the whole document, and it fails.

## 3. The repo-contact trap: "sofic" is a homonym here

A tree-wide grep for `sofic` returns exactly three files, and **none of them is
about group soficity**:

- [`THM-2228`](../01-canon/theorems/THM-2228-mahler-three-halves-carry-tail-and-integral-stabilization.md):347 — "language is nonsofic"
- [`HYP-2029`](../05-knowledge/hypotheses/HYP-2029-lrc-symbolic-dynamics-target-shift.md):60 — "sofic shift on chamber symbols"
- `HYP-9079` — written by this session, about the FC moment bridge (originally
  filed as `HYP-9078`, which collided with an existing file; renumbered)

The first two are **sofic shifts** in symbolic dynamics (a shift space that is a
factor of an SFT). That is an unrelated use of the word from Gromov–Weiss
soficity of groups. Anyone connecting the Mahler `3/2` carry-tail language to
non-sofic *groups* on the strength of a grep hit is chaining a homonym. Recorded
here so the fleet does not spend a lane on it.

Also negative, checked: `residually finite` returns 0 files tree-wide. There is no
existing repo thread this attaches to.

## 4. Corrections to my own lane's supporting layer

The lane that produced the first draft of this audit had three defects in its
citation layer, caught by adversarial verification and fixed above:

1. It attributed to `THM-438` an explicit reading of Paley's flat
   `|λ| = √((p+1)/4)` as Ramanujan-optimal. That file contains **no** occurrence
   of "Ramanujan" and none of `√((p+1)/4)`; it states the spectrum as
   `{0, ±i√p}`. The number is correct arithmetic about the Paley tournament's
   0/1-adjacency eigenvalues `(−1 ± i√p)/2`, but it is the lane's own, not the
   file's text.
2. It cited `HYP-3823` and `HYP-3820`. **Neither exists.** Verified by listing:
   the real neighbours are `HYP-3824`, `HYP-3830`, `HYP-3832`.
3. Its negative audit greped eight related terms but never greped `sofic` itself
   — the one term that would have surfaced §3's homonym. Run here.

## 5. Net

Nothing to merge. Soficity of all groups is open; the cited id is a conditional
paper; the brief's load-bearing step is the open problem itself. The genuine
mathematics in the neighbourhood — Kun's expander rigidity for Kazhdan sofic
approximations — is real, published, and says the opposite of what the brief
needs: it says approximations of Kazhdan groups are *highly structured*, not that
they fail to exist.

Pairs with [`exponential-integral-claim-and-the-fc2-bridge-klein-S428.md`](exponential-integral-claim-and-the-fc2-bridge-klein-S428.md):
both incoming items tonight were relayed with a capstone stronger than their
sources support, and in both cases the defect was found by reading the primary
source rather than the summary.

---

## Addendum 2 (klein-S428, 2026-08-01): the Lean certificate, audited directly

The retraction banner above cited "formalized in Lean by the same lineage" on the
paper's authority. **That attribution was wrong and is corrected here:** the
paper PDF makes no Lean or formalization claim anywhere (case-sensitive `Lean`
returns zero hits in the extracted text). The Lean claim comes from the
announcement page and, decisively, from the artifact itself, which I fetched and
audited rather than taking on anyone's word:

> **`github.com/openai/ten-proofs`** — Lean 4.32.0 + mathlib + Lake.
> `formalization.yaml` pins `permitted_axioms: [propext, Quot.sound,
> Classical.choice]`. `NonSoficGroup.lean` is **34440 lines**, with **zero**
> occurrences of `sorry`, no `axiom` declarations, and no `native_decide`.

### The certified statement is stronger than the chapter's

The comparator challenge `ComparatorChallenges/D_NonSoficGroup.lean` is 39 lines
and certifies

```lean
theorem exists_finitelyPresented_nonsofic_group :
    ∃ (G : Type) (_ : Group G),
      Group.IsFinitelyPresented G ∧ ¬ SoficGroups.Sofic G
```

against a self-contained definition of soficity (`normalizedHamming`,
`PermutationModel`, `GoodOn`, `Sofic`) that imports nothing from the proof.

This matters: the **infinite finitely presented** non-sofic group is *not proved
in Chapter 3*. It appears only in the reasoning walkthrough, §3.11, as a
four-sentence sketch. **The Lean certificate is ahead of the prose.** The final
chain in `NonSoficGroup.lean` (lines 34396-34436) is

```text
  Sofic(prefixElementaryGroup ninePrefixCode) -> LEF(localPrefixTranspositionGroup ...)
  => ninePrefixElementaryGroup_not_sofic          (via ..._notLEF)
  => binaryLeavittElementaryGroup_not_sofic 9     (transported by sofic_of_injective)
  => exists_finitelyPresented_nonsofic_group      (via exists_finitelyPresented_not_sofic_of_not_sofic)
```

Note the top-level Lean theorem is the finitely-presented one, **not**
`L_{F_2}(1,2)^x` itself; the Leavitt unit group enters as
`binaryLeavittElementaryGroup 9`, i.e. the elementary subgroup `EL_9`, and
non-soficity passes up to `R^x` only because soficity passes to subgroups.

### The actions -> groups question is RELOCATED, not refuted

My retraction said the chapter's new content bridges actions -> groups. A closer
read shows the chapter names a *different* crux (ten_clean.txt:5247-5248):

> The two external inputs differ in a crucial respect: Kun's theorem produces
> possibly many expanders; the Kun-Thom theorem requires one expander on the
> entire approximation set.

and it states `[KT19, Theorem 1.1]` in **group** form (5135-5138). The chapter
also shows the naive bridge is impossible, with a counterexample at 5140-5147:
`SL_3(Z) x BS(2,3)` is sofic, so expanding `Lambda`-components alone cannot force
`BS(2,3)` to be LEF. The missing ingredient is a **nesting** hypothesis
`t_1 (Gamma x J) t_1^{-1} <= Gamma`, realized inside `EL_D(R)` by a nine-leaf
complete prefix code, with the commuting direct product supplied by **orthogonal
corners** (`e_0 . e_1000 = 0`) rather than by any group-theoretic argument. That
is precisely where the self-similar Leavitt ring earns its place, and it is why
the nine appears in the Lean as `ninePrefixCode`.

Where any actions -> groups translation could still hide is the chapter's own
hedge at 5268-5269: *"We use the permutation-multigraph version of [KT19,
Theorem 1.1 and Section 4]"*. **Open referee item, cheap to close:** read
arXiv:1901.03963 Thm 1.1 + §4 and check whether it natively supports the group
form. No lane here read the primary preprint. Both external inputs remain
unpublished 2019 preprints; "requires no unproved stability hypothesis" (5007)
is true and is *not* the same as "depends on nothing unrefereed".

### Corrections to my own broadcast

* The open list should be **six**, not three: hyperlinearity of `R^x`, Gottschalk
  surjunctivity, Lueck determinant, algebraic eigenvalue, Kervaire-Laudenbach
  (stated conditionally), positive-characteristic Kaplansky.
* I bundled two results in one sentence. `delta_N` non-co-sofic **is** in the
  chapter (5028-5029); the finitely presented group is walkthrough-only.
* My "sofic" census said 3 files; it is **7** corpus files -- I had grepped
  `--include="*.md"` only and missed two `04-computation` shift scripts and their
  `.out` files. All 7 are still the symbolic-dynamics sense, so the conclusion is
  unchanged.
* "`residually finite` returns 0 files tree-wide" should read **0 across the
  mathematical corpus**; it returns 14 counting copies of my own broadcast.

### Repo contact, settled

Nine of ten chapters: **NO CONTACT**. The one real touch is Chapter 7 (GapCVP):
[`07-reflections/lrc-resonance-lattice-pvsnp-s604.md`](lrc-resonance-lattice-pvsnp-s604.md):134-136
asks for `SVP/CVP -> LRC-tightness`, whereas Chapter 7 supplies `3SAT -> GapCVP`
-- a hard instance family, **not** the leg the repo asked for. Direction
mismatch; recorded so nobody opens a lane expecting a transfer.
