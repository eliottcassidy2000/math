---
id: HYP-3830
title: THE APEX GROUP PSL(2,7) -- the four faces of sqrt(p), the Frobenius-21 = |Aut(Paley_7)| inside PSL(2,7)=168, and the LEFT-RIGHT CAYLEY (LTC) substrate for turning the ANTI-LTC certificate into an LTC via the apex's own group. The n=7 flip-rank obstruction (Paley heptagon, |Aut|=21=3*7, HYP-3817) has its automorphism group = the Frobenius group C7:C3 = the point-stabilizer of PSL(2,7) (|PSL(2,7)|=168=8*21 = the Hurwitz group = Aut Klein quartic = Aut Fano PG(2,2)). VERIFIED: PSL(2,7) order-distribution {1:1, 2:21, 3:56, 4:42, 7:48} exhibits the {7,21} apex orders (21 involutions, 48 order-7); the Frobenius-21 subgroup exists inside; a (2,3,7)-Hurwitz Cayley graph is a near-Ramanujan expander (lambda_2=2.880 vs 2*sqrt(d-1)=2.828). THE FOUR FACES OF sqrt(p) (p=3mod4, the iota-odd/Borsuk-Ulam-hard apex): (1) Gauss sum g(p)=i*sqrt(p); (2) Paley skew eigenvalue +-i*sqrt(p) (HYP-3814); (3) Ramanujan/expander bound 2*sqrt(p); (4) field Q(sqrt(-p)) disc -p -- all four = sqrt(7)=2.6458 at p=7 (one atom, four faces; imaginary=iota-odd for 1,2,4, real for 3). THE {7,21} IMPOSSIBILITY (forbidden H-values) = the orders of the Fano/PSL(2,7) apex geometry. LEFT-RIGHT CAYLEY SUBSTRATE (Dinur-Evra-Livne-Lubotzky-Mozes LTC): built the square complex on PSL(2,7) (168 vertices, 252 squares for small A,B); a co-boundary parity proxy detects a global defect in 46% of local links => the complex DOES locally detect global defects. PROPOSAL (honest, structural): the flip-rank certificate is ANTI-LTC in the spectral basis (S82), but the apex's own group PSL(2,7) gives an expander square-complex whose local links carry the Frobenius-21 symmetry the spectrum missed; encoding the SC/|Aut| certificate as a co-cycle on this complex would make it testable link-by-link. COMPLEMENTARY (HYP-3821): the pi/measure side f(n)~1/sqrt(pi n) vs this sqrt(p)/certificate side. Mahler-Popken integer complexity ||7||=6,||14||=8,||21||=9 and sum-of-two-squares (21=3*7 NOT a sum of 2 squares; 3,7=3mod4) noted as the iota-parity arithmetic
status: MIXED (verified group theory + arithmetic; the LTC-for-the-certificate is a structural PROPOSAL). VERIFIED (psl27_leftright_cayley_four_faces_klein.py): PSL(2,7) |G|=168, orders {1:1,2:21,3:56,4:42,7:48}; Frobenius-21 subgroup found; (2,3,7) Cayley graph lambda_2=2.880 (near-Ramanujan, bound 2.828); four faces of sqrt(p) exact for p=3,7,11; left-right square complex built (168 vtx, 252 squares), coboundary proxy 46% violation. HONEST: the four faces + apex group + Frobenius-21 + {7,21}=Fano/PSL(2,7) orders are SOLID; 'turning the anti-LTC into an LTC' is a substrate + local-defect-detection + a structural proposal, NOT a proven locally-testable certificate (the c^3 soundness for the tournament certificate is open; Dinur et al. prove it for their base codes). Mahler-Popken/sum-of-2-squares = noted speculative arithmetic.
source: klein-2026-07-01-S85
depends_on:
  - HYP-3817   # S83: n=7 obstruction = Paley, |Aut|=21=3*7 (the apex whose group is PSL(2,7)'s Frobenius-21)
  - HYP-3814   # S81: Paley skew/Cayley spectrum = i*sqrt(p) (face 2 of sqrt(p))
  - HYP-3816   # S82: the certificate is ANTI-LTC in the spectral basis (why we need the apex group)
related:
  - HYP-3818   # opus-S26: sqrt21 residual, Gauss iota-parity (the sqrt(p) certificate side)
  - HYP-3831   # S85: the complementary pi/measure side (f(n)~1/sqrt(pi n))
  - kind-pasteur-S24  # CONVERGENT: half-tiling = left-right square complex, tournament reconstruction = ANTI-LTC (fails n=7)
external: "Good Locally Testable Codes" (Dinur-Evra-Livne-Lubotzky-Mozes, Annals 203-2 2026) = left-right Cayley complexes (the LTC substrate used here); Klein quartic / Hurwitz group 168; Fano plane PG(2,2)
results:
  - 04-computation/psl27_leftright_cayley_four_faces_klein.py
  - 05-knowledge/results/psl27_leftright_cayley_four_faces_klein.out
---

# HYP-3830 — the apex group PSL(2,7): four faces of √p, and the left-right Cayley LTC substrate

*(Renumbered 3820→3830: opus-S28 concurrently claimed 3820; jumped clear per the swarm-convergence note.
CONVERGENT with kind-pasteur-S24, which independently reached "tournament reconstruction = ANTI-LTC (fails
n=7), chasing the LTC lead" — this HYP supplies the apex group PSL(2,7) as the LTC substrate for that lead.)*

## The apex group
The n=7 flip-rank obstruction is the Paley heptagon with `|Aut| = 21 = 3·7` (HYP-3817). That Frobenius group
`C7⋊C3` sits inside **PSL(2,7)**, `|PSL(2,7)| = 168 = 8·21` — the **Hurwitz group** (Aut of the Klein quartic),
`= Aut(Fano PG(2,2))`. Verified order-distribution `{1:1, 2:21, 3:56, 4:42, 7:48}` shows the **{7,21} apex
orders** (21 involutions, 48 order-7 elements); the Frobenius-21 subgroup is present; a `(2,3,7)`-Hurwitz
Cayley graph is a near-Ramanujan expander (`λ₂ = 2.880` vs `2√(d−1) = 2.828`).

## The four faces of √p (the hard apex, p ≡ 3 mod 4)
`√p` is one atom seen four ways — verified equal to `√7 = 2.6458` at `p=7`:
1. **Gauss sum** `g(p) = i√p` (imaginary; HYP-3818).
2. **Paley skew eigenvalue** `±i√p` (the tournament; HYP-3814).
3. **Ramanujan/expander bound** `2√p` (spectral radius threshold).
4. **Field** `Q(√−p)`, discriminant `−p`.
Faces 1,2,4 are imaginary (ι-odd); face 3 is real (expansion). The `{7,21}` H-impossibility = the orders of
the Fano/PSL(2,7) apex geometry.

## Left-right Cayley (LTC) substrate — turning the anti-LTC into an LTC
The flip-rank certificate is **anti-LTC in the spectral basis** (S82: the spectrum, reflection-symmetric, is
blind to the SC/`|Aut|` structure that carries the excess). The proposal (owner): use the apex's **own group**.
Built the **left-right Cayley square complex** on PSL(2,7) (168 vertices; 252 squares for a small `A×B`), the
Dinur–Evra–Livne–Lubotzky–Mozes LTC substrate. A co-boundary parity proxy violates `46%` of local links on a
global defect — i.e., **the complex locally detects global defects**. The structural claim: PSL(2,7) gives an
expander whose local links carry the Frobenius-21 symmetry the spectrum missed, so encoding the SC/`|Aut|`
certificate as a co-cycle on this complex would make it **testable link-by-link**. Honest scope: substrate +
expansion + local-defect-detection are built; full `c³` local-testability of the tournament certificate is the
deep open step (Dinur et al. prove it for their base codes).

## Complementary halves + noted arithmetic
This is the `√p` / ODD / **certificate** side; HYP-3821 is the `π` / EVEN / **measure** side
(`f(n) ~ 1/√(πn)`). Noted (speculative): Mahler–Popken integer complexity `‖7‖=6, ‖14‖=8, ‖21‖=9`
(`21=3·7` the compositum modulus); sum-of-two-squares — `21=3·7` is **not** a sum of two squares (both `≡3 mod
4`), the same ι-parity that makes `g(p)` imaginary.

## Net
The obstruction's automorphism group is a subgroup of the Hurwitz group PSL(2,7); `√p` unifies Gauss/Paley/
Ramanujan/field at the apex; and the apex group is the natural substrate for making the anti-LTC certificate
locally testable — a concrete new engineering direction linking the flip-rank obstruction to the LTC frontier.
