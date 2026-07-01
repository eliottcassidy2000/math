---
id: HYP-3814
title: TOURNAMENTS LIVE IN THREE GEOMETRIES GLUED BY THE CAYLEY TRANSFORM; COMPLEMENT IS ONE REFLECTION IN ALL THREE; THE ODD/EVEN (Redei) DUALITY IS ONE FACT SEEN THREE WAYS. A tournament is (1) a TILING of the staircase delta_{n-2}, (2) a RUNNER-SET on the unit circle = the eigenvalues of the Cayley transform U=(I-S)(I+S)^{-1} of its skew +-1 matrix S=A-A^T (S skew => U orthogonal, eig(U) on |z|=1), and (3) the Cayley transform itself is the GLUE S|->U from the combinatorial/staircase side to the circle spectrum. COMPLEMENT = ONE INVOLUTION in three coordinates: staircase = grid reflection sigma:(x,y)->(n+1-y,n+1-x) (opus-S18); skew = negation S->-S (A->A^T); circle = conjugation theta->-theta (U(-S)=U(S)^{-1}=U(S)^T). Via Cayley these are the SAME reflection. ODD/EVEN = ONE FACT: a class is SELF-COMPLEMENTARY (SC) iff FIXED by the reflection; Redei (H odd) makes the fixed tiling-fiber ODD (SC merged node = H), the moved fiber EVEN (NS merged node = H(T)+H(T^op)=2H) -- the SC-odd/NS-even merged-fiber parity, read as (staircase) sigma-fixed-fiber-odd, (circle) palindromic-measure-odd-mass, (spectral) U perm-conjugate to U^{-1}. PLUS the skew/Cayley spectrum is ALWAYS complement-blind (spec(-S)=spec(S), skew spectra are +-symmetric) => the reflection acts TRIVIALLY on spectra, explaining the S72 spectral-weakness. PLUS n-parity: n odd <=> S singular (odd skew) <=> U has eigenvalue +1 (a runner PINNED at angle 0). Bonus (exact): the Paley-p (p=3 mod 4) Cayley eigenvalue sits at cos(theta) = -(p-1)/(p+1) (=-3/4 at p=7, -5/6 at p=11), encoding the GAUSS SUM +-i*sqrt(p) -- NOT the p-th roots of unity (those are the VERTEX loop, a DIFFERENT circle; two circles clarified). LRC relevance: the runner picture (circle) and the tournament picture (staircase) are one object; complement/converse is the antipode iota; the covering-min phase cloud (HYP-3813) lives on the same circle
status: VERIFIED (linear-algebra identities are exact + proved by hand; iso-class facts exact n=3,4,5). PROVED-BY-HAND: (a) S skew => U=(I-S)(I+S)^{-1} orthogonal (U^T U=I since (I+-S)^T=(I-+S) commute); (b) complement A->A^T => S->-S => U->U^{-1}=U^T => eig conjugated (reflection theta->-theta); (c) spec(-S)=spec(S) always (skew spectra +-symmetric) => Cayley spectrum complement-blind; (d) n odd => det S=0 (odd skew) => +1 in eig(U). VERIFIED (computed): n=3,4,5 all tournaments -- U orthogonal, |eig|=1, complement->U^{-1}, (n odd <=> +1 eig) ALL True; class counts |G_n|=2,4,12=A000568, #SC=2,2,8, NS-merged=0,1,2 (match canon), SC fiber H ODD + NS fiber 2H EVEN all True (Redei), spec(S)=spec(-S) all True; Paley p=7,11 cos(theta)=-(p-1)/(p+1) exact. HONEST: a SYNTHESIS/reframe unifying the tiling model, the OPUC/runner circle, and the skew-Cayley spectrum into three coordinate systems for one object, with complement=reflection and the Redei odd/even as the reflection's fixed-point parity. Not a new theorem; a coherent geometric dictionary. Corrects a RECOLLECTION (not the file): HYP-3802's 'roots of unity' is the vertex loop, not the Cayley spectrum.
source: klein-2026-07-01-S81
depends_on:
  - HYP-3802   # S70: runners on the loop / dihedral atoms / vertex-loop roots of unity (circle i)
  - HYP-1772   # Redei: H(T) odd (=> SC fiber odd) + |Aut| odd
related:
  - HYP-3813   # S80: the covering-min phase cloud lives on the same (circle) geometry
  - HYP-3804   # S72: skew-spectrum weakness -- EXPLAINED here (spectrum is complement-blind)
  - HYP-3811   # S78: complement = grid reflection sigma; quarter-fold Klein-four
  - HYP-3808   # S75: merged-metagraph fiber parity (SC odd / NS even) -- the reflection fixed/2-orbit
  - HYP-3801   # S69: Verblunsky/OPUC lonely measure on the circle
results:
  - 04-computation/three_geometries_cayley_complement_klein.py
  - 05-knowledge/results/three_geometries_cayley_complement_klein.out
---

# HYP-3814 — three geometries, one reflection, one odd/even fact

## The three geometries (one object)
A tournament `T` on `n` vertices is simultaneously:
1. **Staircase.** A tiling of the right-isosceles staircase `delta_{n-2}` (fix the base Hamiltonian path;
   tiles = non-consecutive arcs; the cube `Q_m`, `m = C(n-1,2)`).
2. **Circle.** A runner-set on the unit circle: the eigenvalues of the **Cayley transform**
   `U = (I - S)(I + S)^{-1}` of the skew `+-1` matrix `S = A - A^T`. `S` skew `=>` `U` orthogonal and
   `eig(U)` lie on `|z| = 1`. (Cayley maps the imaginary axis `spec(S) ⊂ iR` to the circle.)
3. **The Cayley transform is the GLUE** `S |-> U`: the combinatorial/staircase side maps to the circle
   spectrum. It is a bijection on the relevant data and intertwines the two complement actions.

## Complement is ONE reflection (three coordinates)
| geometry | complement action | why |
|---|---|---|
| staircase | grid reflection `sigma: (x,y) -> (n+1-y, n+1-x)` | opus-S18 (`sigma` = complement, linear involution) |
| skew matrix | `S -> -S` | `T^op: A -> A^T`, so `S = A - A^T -> A^T - A = -S` |
| circle | conjugation `theta -> -theta` | `U(-S) = (I+S)(I-S)^{-1} = U(S)^{-1} = U(S)^T`; `|lambda|=1 => 1/lambda = conj(lambda)` |

Via Cayley, `sigma <-> (S -> -S) <-> (theta -> -theta)` are the **same involution**, transported between
coordinate systems. Verified (n=3,4,5, all tournaments): `complement -> U^{-1}` holds exactly.

## The odd/even (Redei) duality is ONE fact
A class is **self-complementary (SC)** iff **fixed** by the reflection. Rédei (`H(T)` = Ham-path count is
**odd**) then gives the **merged-fiber parity**:
- an **SC merged node** carries a single class, tiling-fiber `H` — **ODD**;
- an **NS merged node** carries a complement pair, fiber `H(T) + H(T^op) = 2H` — **EVEN**.

This one fact reads three ways: **[staircase]** `sigma`-fixed fiber is odd; **[circle]** the reflection-fixed
(palindromic / self-conjugate) measure has odd Rédei mass; **[spectral]** SC `<=>` `U` permutation-conjugate
to `U^{-1}` (the reflection's fixed set). Verified (n=3,4,5): SC fiber `H` odd = True, NS fiber `2H` even =
True, class counts match `A000568` and the canon NS-merged `= 0,1,2`.

## Two corollaries (both verified)
- **Spectral blindness (explains S72/HYP-3804).** Skew spectra are `+-`-symmetric, so `spec(-S) = spec(S)`
  **always**: the reflection acts **trivially** on the spectrum. Hence the skew/Cayley spectrum **cannot
  detect complement** — the exact reason the skew-spectrum is a weak invariant.
- **`n`-parity is a pinned runner.** `n` odd `=>` `det S = 0` (odd-size skew) `=>` `+1 in eig(U)`: a runner
  **pinned at angle 0**. (n=3,4,5: `n` odd `<=> +1` eigenvalue, exactly; no `n=4` tournament has a singular `S`.)

## Bonus: the Cayley spectrum encodes the Gauss sum (two circles)
The **vertex loop** (circle i, HYP-3802) places a circulant tournament's vertices at roots of unity. The
**Cayley spectrum** (circle ii, this file) is different: the Paley-`p` (`p = 3 mod 4`) skew circulant has
`spec(S) = {0, +-i*sqrt(p)}` (the Gauss sum), so its Cayley eigenvalues sit at
> `cos(theta) = (1 - p)/(1 + p) = -(p-1)/(p+1)` (`-3/4` at `p=7`, `-5/6` at `p=11`) — verified exact,
an irrational multiple of `2pi` (**not** a root of unity; `U^p != I`). So "roots of unity" is a vertex-loop
statement, **not** a Cayley-spectrum statement — two circles, both with complement = reflection.

## LRC relevance
The runner picture and the tournament picture are **one object** in two coordinate systems glued by Cayley;
complement/converse is the antipode `iota`. The covering-min phase cloud (HYP-3813) lives on the same circle,
so the tournament reflection `sigma` and the runner antipode `iota` are the same map — the fold work (S78/S79)
and the runner work (S80) are two faces of one reflection.

## Net
Three geometries (staircase / circle / Cayley-glue), one reflection (complement), one odd/even fact (Rédei
fixed-point parity). The skew-spectrum's weakness and the `n`-parity's pinned runner both fall out; the Gauss
sum appears as the Cayley angle. A coherent geometric dictionary, not a new theorem.
