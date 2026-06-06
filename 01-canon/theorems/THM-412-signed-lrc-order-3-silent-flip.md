# THM-412 — Signed-LRC order-3 silent flip: the AP_n sign-orbit degenerates through the point x=(2n−1)/3

**Status:** PROVED (degree-parity lemma + Euler-decomposition corollary; both verified exhaustively).
**Source:** monad-explorer-2026-06-06-S703. Sharpens and corrects S699's T4 ("sign-orbit
collisions are config automorphisms" — FALSE). Builds on the signed-LRC theory of S699
(T1 gauge-invariance, T2 sign=cut, T3 zero-clock⟺shell-partner) and the prime-3 thread
(THM-407, n=14 ⟺ C=27=3³).

---

## Setup (the signed-LRC cut-spectrum)

Config = the runners `V = {1,…,n−1}` of the arithmetic progression `AP_n` (observer = speed 0),
shell modulus `C = 2n−1` (THM-401). A sign vector `ε∈{±1}^{n−1}` is a 2-coloring (a **cut**) of the
runners; consider cuts up to global swap, so there are `2^{n−2}` of them. The **folded clock** of a
pair `{i,j}` (S699c) is
```
   clock(i,j) = |i−j|              if i,j monochromatic (same color)
              = ρ(i+j)            if i,j bichromatic   (a cut edge),   ρ(s) := min(s mod C, C−(s mod C)).
```
Each cut yields a **folded clock-multiset**; the **sign-orbit** is the number of distinct such
multisets. (Mono clocks are never folded: `|i−j| ≤ n−2 < C/2`.)

A cut `χ` and the cut `χ⊕flip(x)` obtained by flipping the color of one runner `x` are
**folded-equal** ("a silent flip of `x`") iff flipping `x` leaves the whole clock-multiset
unchanged. Flipping `x` toggles every incident edge `{x,y}` between its mono value `|x−y|` and its
cut value `ρ(x+y)`. Define the **value-multigraph** `G_x` on the value set with one edge
`{|x−y|, ρ(x+y)}` for each `y ∈ V\{x}`.

---

## Lemma A (silent flip ⟺ `G_x` Eulerian) — PROVED

> A silent flip of `x` (a coloring `χ` with `χ` and `χ⊕flip(x)` folded-equal) **exists iff every
> vertex of `G_x` has even degree** (i.e. `G_x` is Eulerian).

**Proof.** Flipping `x` replaces, for each `y`, the contribution of edge `{x,y}` by the *other*
member of the value-pair `{|x−y|, ρ(x+y)}`. The multiset is preserved iff we can choose, for each
incident edge, which member is "before" vs "after" so that the before-multiset equals the
after-multiset. That is exactly a 2-coloring of the edges of `G_x` in which every value-vertex has
equal degree in both color classes. Such a balanced edge 2-coloring exists iff every vertex of `G_x`
has even degree: necessity is immediate; sufficiency follows by decomposing the (all-even) `G_x`
into closed Eulerian walks and 2-coloring each walk alternately, which balances every vertex. The
non-incident edges are untouched, and the colors of `x`'s neighbours are free, so the required
coloring `χ` is realizable. ∎

## Lemma B (degree parity of `G_x`) — PROVED

> For every runner `x∈{1,…,n−1}`, the odd-degree vertices of `G_x` are **exactly** `{x, ρ(2x)}`.
> Hence `G_x` is Eulerian iff these two coincide, i.e. iff `x = ρ(2x)`.

**Proof.** Fix a value `v`. The degree of `v` in `G_x` counts the runners `y∈V\{x}` with `|x−y|=v`
(mono side) plus those with `ρ(x+y)=v` (cut side). In `ℤ/C`:
`|x−y|=v ⟺ y ≡ x∓v`, and `ρ(x+y)=v ⟺ x+y ≡ ±v ⟺ y ≡ −x±v`. The four residues are
`±(x−v)` and `±(x+v)`. The runner set `V={1,…,(C−1)/2}` is a **half-system** of `ℤ/C`: every
nonzero class `{r,−r}` has exactly one representative in `V`. So each of the two `±`-pairs
contributes exactly one in-range runner — giving `deg(v)=2` — *except* for two corrections:
(i) the pair `±(x−v)` degenerates to `{0}` (contributing 0, not 1) when `v ≡ ±x`, i.e. `v = ρ(x)=x`;
(ii) an in-range representative is the excluded runner `x` itself, which happens precisely when one
of `±(x−v), ±(x+v)` equals `x`, i.e. when `v = ρ(2x)`. Each correction subtracts exactly 1, flipping
the parity at that single value. Thus `deg(v)` is even for all `v∉{x,ρ(2x)}` and odd at
`v∈{x,ρ(2x)}` — *unless* `x=ρ(2x)`, in which case the two corrections land on the same vertex and
cancel, leaving all degrees even. (Consistent with the handshake lemma: a multigraph always has an
even number of odd-degree vertices.) ∎

**Verified:** Lemma B checked for all `(n,x)`, `n=4..40` (777 pairs, 0 mismatches),
`signed_lrc_thm411_proof_check_s703e.out` [filename predates the THM-411→412 renumber].

## Theorem (order-3 silent flip) — PROVED

> `x = ρ(2x)` ⟺ `2x ≡ −x (mod C)` ⟺ `3x ≡ 0 (mod C)` ⟺ `x = C/3` (since `1≤x≤(C−1)/2 < C`).
> Therefore: **`G_x` is Eulerian iff `x = (2n−1)/3`**, which exists as a runner iff `3 ∣ (2n−1)`.
> Consequently, when `3 ∣ (2n−1)` the order-3 torsion runner `x=(2n−1)/3` admits a silent single
> flip, so the `AP_n` folded sign-orbit is **strictly smaller than `2^{n−2}`**.

The witnessing runner is the **order-3 point** of `ℤ/C` (`3x≡0`): the signed-LRC face of the
prime-3 structure (THM-407). At `n=14`, `C=27=3³`, the silent runner is `x=9`; `27` being the
3-richest small modulus, its sign-orbit is the most degenerate (4096→4027, i.e. 69 colliding
groups; the 8 single-flip collisions all flip `x=9`).

---

## Corollaries / scope

- **Correction of S699 T4.** S699 reported the sign-orbit deviating from `2^{n−2}` and ascribed the
  collisions to "config automorphisms." This is false: the **unfolded** cut-spectrum (clocks taken
  with no mod-`C` fold) is **faithful for every config tested** (all `2^{n−2}` cuts distinct;
  verified exhaustively `n≤7` all speed-sets in `[1,14]`, plus adversarial symmetric/Sidon-violating
  sets). So no automorphism is involved — *every* collision is a **fold-only** artifact created by
  reducing clocks mod `C=2n−1`. THM-412 is the clean sufficient cause when `3∣C`.

- **Sufficient, not necessary, for orbit degeneration.** `3∣C` ⟹ orbit `< 2^{n−2}` via a *single*
  silent flip. Composite `C` with **no** factor 3 (e.g. `C=25` at `n=13`, `C=35` at `n=18`,
  `C=49`) still has a smaller orbit, but those collisions are all **multi-flip** (no single runner
  has `G_x` Eulerian — verified `eulx=∅` for `3∤C`). The full law is conjectural: see HYP-2267.

---

## Honest status

- **Proved:** Lemma A (silent flip ⟺ `G_x` Eulerian); Lemma B (odd-degree vertices `={x,ρ(2x)}`);
  the order-3 theorem `G_x` Eulerian ⟺ `x=(2n−1)/3`. Hence `3∣(2n−1) ⟹` orbit `<2^{n−2}`.
- **Verified:** Lemma B (n≤40); the `x=C/3` Eulerian property (n≤199, all `3∣C`); unfolded
  faithfulness (exhaustive n≤7 + adversarial).
- **Conjecture (HYP-2267):** `AP_n` folded sign-orbit `= 2^{n−2}` **iff `2n−1` is prime**
  (verified n=3..22, C=5..43). The `3∣C` direction is THM-412; other composite factors give
  multi-flip collisions whose count law is open.
- **Not claimed:** any statement about the observer gap `M` (which is sign-invariant by T1); this is
  purely about the finer signed (pair-clock) invariant.

**Artifacts:** `04-computation/signed_lrc_orbit_collisions_s703.py`,
`signed_lrc_faithfulness_s703b.py`, `signed_lrc_fast_confirm_s703c.py`,
`signed_lrc_collision_mechanism_s703d.py` (+`.out`s),
`05-knowledge/results/signed_lrc_thm411_proof_check_s703e.out`,
`07-reflections/signed-lrc-fold-collisions-and-the-order-3-silent-flip-s703.md`. HYP-2267, T757.
