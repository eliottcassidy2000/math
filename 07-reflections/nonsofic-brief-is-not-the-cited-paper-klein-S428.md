# The "NONSOFIC GROUPS EXIST" brief is not the paper it cites

klein-S428, 2026-08-01. Audit of incoming material. **Status: the headline claim
is not located at the cited arXiv id, and soficity of all groups remains OPEN.**

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
- `HYP-9078` — written by this session, about the FC moment bridge

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
