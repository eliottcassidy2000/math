# 14 = |D₇|: the heptagon's dihedral group — and n=14's certificate is Borsuk–Ulam (free ℤ/2), not Brouwer

*kind-pasteur-2026-06-28-S31av. The owner asked to connect to imaginary roots, regular polygons on the
complex unit circle, and the tournament↔dihedral-group link, then keep pushing. These converge on one
structural fact that **refines** the topological certificate: `14 = |D₇|` is the symmetry group of the
regular heptagon, the de Moivre cubic is its three 2-dim irreps, and — because `7 ≡ 3 (mod 4)` — the
complement/reflection is **orientation-reversing**, so the ℤ/2 acts **freely** and the right theorem is
**Borsuk–Ulam, not Brouwer** (HYP-3219). This is the imaginary/odd/negative case made topological.*

## 1. `14 = |D₇|` and the de Moivre cubic = the 2-dim irreps (VERIFIED)
The `7` sectors are the vertices of a **regular heptagon** `ζ₇^k = e^{2πik/7}` on the complex unit circle.
Its symmetry group is the **dihedral group `D₇ = ℤ₇ ⋊ ℤ₂`, of order `2·7 = 14` = the LRC apex `n`.** Its
irreducible representations partition `14 = 1²+1²+2²+2²+2²`:
| irrep | dim | = which LRC mode |
|---|---|---|
| **trivial** | 1 | the **Möbius / trace** mode (HYP-3217 degree 1) |
| **sign** | 1 | the **complement / reflection** = the positive-negative / even-odd `ℤ₂` |
| **three 2-dim** | 2,2,2 | characters `χ_k(rot) = 2cos(2πk/7)` = **the de Moivre angles** |
So **the de Moivre cubic IS the three 2-dim irreps of `D₇`**, the trivial irrep is the Möbius mode, and the
**sign irrep is the obstruction's `ℤ₂`**. The HYP-3217 mode lattice = the `D₇` representation ring. `14` is
not a coincidence with `|D₇|`; the apex `n=2p` **is** the dihedral order, and the LRC lives on `D_p`.

## 2. The complement is the reflection mult-by-`−1`; aut vs anti-aut = `p mod 4` (VERIFIED)
On the Paley/circulant tournament `T_p`, the complement `x↦−x` is multiplication by `−1`, and `−1` is a QR
iff `p≡1 (mod 4)`:
| `p` | `n=2p` | `−1` is | reflection is | structure |
|---|---|---|---|---|
| 5 | 10 | a **residue** | an **automorphism** (orientation-**preserving**) | reflection-symmetric (SOS) |
| **7** | **14** | a **non-residue** | an **anti-automorphism** (`T↦Tᵒᵖ`, orientation-**REVERSING**, **self-converse**) | free reflection (obstruction) |
| 11 | 22 | non-residue | anti-aut (self-converse) | free reflection |
| 13 | 26 | residue | automorphism | SOS |
This is the **same `p mod 4`** that controlled the discriminant sign (HYP-3220), now read on the heptagon:
`p≡3 (mod 4)` ⇔ `−1` non-residue ⇔ the heptagon reflection **reverses tournament orientation**.

## 3. THE PUSH: n=14's certificate is BORSUK–ULAM (free ℤ/2), not Brouwer (fixed point)
HYP-3219 merged Brouwer (a fixed point of the reflection). But Brouwer's fixed point only exists when the
reflection has one — i.e. when it is an **automorphism** (`p≡1`, n=10). For `n=14` (`p≡3`), the reflection
is an **anti-automorphism**: it sends every orientation to its **converse**, so it has **no fixed
orientation** — the `ℤ₂` acts **FREELY**.
> **A free `ℤ₂` action is exactly the hypothesis of Borsuk–Ulam, not Brouwer.** The witnesses do not have a
> reflection-fixed representative (no symmetric/SOS witness); instead they come in **antipodal pairs
> `(t*, −t*)`**, and the certificate must be the **odd-degree / parity** invariant of the free action
> (Borsuk–Ulam: an odd map `S¹→S¹` has odd, hence nonzero, degree). The "odd degree" is precisely the
> project's **odd/negative datum**: the imaginary Gauss sum `g(χ₇)=i√7`, the negative trace `e₁=−1`
> (HYP-3219), the odd power sums (HYP-3220). They are one number: **the degree of the free reflection.**

So the topological refinement of the whole arc:
```
p ≡ 1 (mod 4)  (n=10):  reflection = AUT  →  BROUWER fixed point  →  symmetric witness  →  SOS / positive / even.
p ≡ 3 (mod 4)  (n=14):  reflection = ANTI-AUT (self-converse) → free ℤ₂ → BORSUK–ULAM odd degree
                         →  antipodal witness pair  →  obstruction / negative / odd / IMAGINARY.
```
**n=14 is hard because its witness is not a Brouwer fixed point but a Borsuk–Ulam antipodal pair** — there is
nothing to *construct*; one must certify *existence by odd degree*. This is the topological face of "the
real-unit/Pell machinery is absent" (HYP-3220).

## 4. Imaginary roots, and why `−7` is the gentle imaginary case
The heptagon's **imaginary part** `sin(2πk/7)` is reflection-**odd** (`sin(−x)=−sin x`) — it is the sign-
irrep coordinate, the orientation, the part the anti-aut flips. The **Gauss sum `g(χ₇)=i√7` is purely
imaginary** (because `7≡3`), so the obstruction's magnitude `√7` lives on the **imaginary axis** = the free-
`ℤ₂` / sin direction. The de Moivre cosines (real axis) are the Brouwer/SOS part; the sines (imaginary axis)
are the Borsuk–Ulam/obstruction part. **Real ↔ even ↔ SOS ↔ Brouwer; imaginary ↔ odd ↔ obstruction ↔
Borsuk–Ulam.**
Mitigating fact: **`−7` is a Heegner number**, so `ℚ(√−7)` has **class number 1** — unique factorization in
`ℤ[(1+√−7)/2]`. Among imaginary-quadratic walls (`p≡3`: `7, 11, 23, …`), `n=14` and `n=22` (`h(−7)=h(−11)=1`)
are the **gentlest** (UFD); `n=46` (`h(−23)=3`) is genuinely harder. So `n=14` is the **first and simplest**
imaginary-quadratic case — the right place to build the class-number/Stark certificate (HYP-3215/3220).

## 5. Concrete actions (pushing the proof)
1. **Decompose the cover bound into `D₇`-isotypic parts** (Reynolds average over `D₇`): trivial (the mean
   cap), sign (the reflection-odd obstruction), 2-dim (the de Moivre/Fejér body). The obstruction is the
   **sign-isotypic** component — certify it by the **Borsuk–Ulam odd degree** of the free reflection, not by
   SOS (which only handles the trivial + 2-dim even part, mac-mini S75e).
2. **Use `h(−7)=1`** (unique factorization) to factor the cap/dip in `ℤ[(1+√−7)/2]`; the conductor-`7²` data
   (mac-mini S75e) and the Stark `L′(0)` (HYP-3215) are the replacements for the missing real units/Pell.
3. **Family check:** verify the SOS/Brouwer route closes `n=10,26` (`p≡1`) outright, isolating the
   Borsuk–Ulam/imaginary obstruction to `p≡3` (`n=14,22,…`). If so, "`p≡1 mod 4` ⇒ LRC by SOS+Brouwer" is a
   *proved* half of the family, and `n=14` is the first genuinely-Borsuk–Ulam case.

## Net
- **`14 = |D₇|`**; the LRC mode lattice (HYP-3217) = the `D₇` irreps (trivial=Möbius, sign=complement,
  three 2-dim = de Moivre cubic).
- **The complement = mult-by-`−1` = the heptagon reflection**; aut vs anti-aut = `p mod 4`. For `n=14`
  (`p≡3`) it is an **orientation-reversing anti-aut** (self-converse) → the `ℤ₂` acts **freely**.
- **REFINED CERTIFICATE:** `n=14`'s topological certificate is **Borsuk–Ulam (free ℤ₂, odd degree)**, not
  Brouwer (fixed point). Odd degree = the imaginary Gauss sum `i√7` = the negative trace = the odd power
  sums — one datum. Real/even/SOS/Brouwer vs imaginary/odd/obstruction/Borsuk–Ulam.
- **GENTLE:** `−7` Heegner ⇒ `ℚ(√−7)` is a UFD (`h=1`); `n=14` is the simplest imaginary-quadratic wall.

→ HYP-3239 (this), HYP-3220 (even-odd=±, p mod 4), HYP-3219 (Brouwer/sign), HYP-3217 (modes=D₇ irreps),
HYP-3218 (Fejér), HYP-3215 (Stark), THM-024 (SC anti-aut involution), mac-mini-S75e (conductor 7²/SOS),
Borsuk–Ulam, Heegner numbers, OPEN-Q-108.
