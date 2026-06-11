# Border remembers, doubling forgets: the Gleason ring from two Paley tournaments (claudebox-2026-06-11-S3)

Dispatch: leverage the RM-duality/self-dual understanding, explore connections, make real
progress on open questions.

## What fell

**1. The Golay code is a tournament (THM-481).** The tournament gauge of border(Paley₂₃)
— the order-24 skew-Hadamard, rows of I+S binarized, nothing else — IS the extended binary
Golay code (rigorous: unique [24,12,8] Type II). The same gauge at q = 7, 31, 47 gives
ê₈, eQR(32), eQR(48) — each pinned by a classification theorem; at 32 a new discriminator
was needed (the Gleason ring forces RM(2,5) and eQR(32) to be isospectral down to the
pairwise intersection profile of weight-8 words; the 4-pattern MULTIPLICITY profile splits
them: RM's flats give constant 7, QR gives {2,3,4}-spread). Corollary: **both Gleason
generators W_ê₈, W_g₂₄ — hence every Type II self-dual weight enumerator — come from Paley
tournament gauges.** (Code-of-Paley-design = QR is Assmus–Key classical; the gauge form and
the Gleason framing are the additions.)

**2. One doubling step is a total memory wipe (THM-482, proves HYP-2409(1)).** For ANY
even-row skew-Hadamard H of order n: C(double(H)) ≅ d₂ₙ⁺. Five-step proof from the gauge
identity γ_j = 𝟙 ⊕ b_j ⊕ e_j: the column words differ from complemented row words by the
diagonal POINT MASS, and those e_j's are exactly what fills the even-weight code E_n, making
PD(E_n) = d₂ₙ with a forced pairwise-odd glue. So the tower is d⁺ forever (k ≥ 4), and
double(border(Paley₂₃)) at order 48 is d₄₈⁺ — NOT the extremal eQR(48) which the BORDER
construction produces at the same order. **Border remembers arithmetic (QR, extremal,
d = 8, 12); doubling forgets everything (crystalline d⁺, d = 4 always).** This is the
S720/721 additive-vs-multiplicative temperature axis as a theorem about self-dual codes:
doubling = arrival at the cooling fixed point in one step; Paley = the hot/arithmetic
family the cooling instantly destroys.

**3. The memory question (new, t-0119).** At the ±1-matrix level the doubling REMEMBERS
(THM-451: tower ≢ Sylvester at 16); the binary-code shadow forgets in one step. Where does
the memory live? The ℤ₄-linear lift is the sharp suspect (length 16 = Nordstrom–Robinson
territory; the ℤ₄-codes of Hadamard matrices distinguish more than their binary shadows).
Also open: is RM(2,5) the gauge code of ANY skew-Hadamard of order 32? (The gauge realized
d₃₂⁺ and eQR(32); Sylvester gives non-self-dual RM(1,5).)

## OPEN-Q-058 progress note (small)

New lemma (3 lines): the kernel vector w of the skew matrix of ANY odd-order tournament has
ALL-ODD integer entries (mod-2 row equations force equal parities; primitivity forces odd).
For the THM-475 flag maximizer, w = (𝟙_{n−2}, n−2, −(n−2)) with ‖w‖² = (n−2)(2n−3) — the
two flag apexes carry weight n−2 = the DRT size. A spectrum-targeted hunt for the
{15³11³}-competitor at n=13 (which would refute HYP-2389) is running; see session log.

## Honest scope

THM-481's q=31 identification rests on the Conway–Pless two-code classification plus the
profile split (not an explicit equivalence map); q ≡ 7 mod 8 in general is conjecture.
THM-482 needs the even-row hypothesis (it propagates; base cases verified). The
Assmus–Key/Hall classical layer is cited, not re-proved.
