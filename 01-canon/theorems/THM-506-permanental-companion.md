---
id: THM-506
title: The permanental companion of the OCF; characteristic and permanental polynomials determine H exactly through order eight, with first failure at nine
status: >
  PROVED for the cycle-packing coefficient theorem and signed/unsigned split.
  COMPUTER-ASSISTED / FINITE-EXACT for recovery of H at n=8: two independent
  exhaustive computations of all 2^21 fixed-Hamiltonian-path presentations,
  with 3217 joint-polynomial fibers and no unequal-H fiber. The direct formula
  proves recovery for n<=7. An explicit independently checked n=9 collision,
  followed by source extension, proves nonrecovery for every n>=9.
  CORRECTED 2026-10-23: the earlier n=8 failure claim confused non-affine
  recovery with arbitrary functional nonrecovery. The n=8 rank gap survives
  only as an affine obstruction.
source: monad-explorer-2026-06-15-S6
depends_on:
  - THM-505   # OCF non-spectral defect: H = spectral skeleton + carriers; Sachs-basis form
  - THM-502   # closed-walk census; c6 non-spectral from onset; Witt transform
  - THM-499   # H = I(Omega,2); alpha_2 first non-spectral OCF ingredient at n=6
related:
  - HYP-2515
  - HYP-2514  # master cycle-packing polynomial Phi (spectrum & H as two specializations)
  - OPEN-Q-096
  - reflection: the-permanent-is-the-unsigned-twin-of-the-spectrum
  - reflection: the-master-cycle-packing-polynomial
---

# THM-506 — the permanental companion of the OCF

> **Correction, 2026-10-23.** The old first-failure threshold n=8 was false.
> The sharp threshold is recovery through n=8 and failure from n=9.
> Parts 1–2 survive unchanged. Part 4 below replaces the invalid rank inference.
> [Independent proofs and exact receipts](../../05-knowledge/results/session37_inheritance_audit_20261023/REPORT.md).

The master cycle-packing polynomial (HYP-2514, reflection `the-master-cycle-packing-polynomial`)
is `Phi(T;{y_k}) = sum_{L linear subdigraph} prod_{C in L} y_{|C|}`, with the spectrum the
**signed** vertex-graded face (`y_k=-x^{-k}`, Sachs) and `H=I(Omega,2)` the unsigned odd-only
fugacity-2 face. This theorem identifies the **unsigned, all-length, vertex-graded** face as a
second classical matrix function — the permanental polynomial — and uses it to bound exactly
how much of `H` is recoverable from the spectrum-plus-permanent.

Throughout, `A` is the 0/1 adjacency of a tournament `T` on `n` vertices; a *linear
subdigraph* `L` is a set of vertex-disjoint directed cycles (length `>= 3`, since tournaments
have no 2-cycles); `#cyc(L)` is the number of cycles; `|V(L)|` the number of covered vertices.

## Part 1 — the permanental coefficient theorem (PROVED)

> **per(xI + A) = sum_{m=0}^{n} e_m^{unsigned} x^{n-m},   e_m^{unsigned} = #{ L : |V(L)| = m }**,

the **unsigned, all-length** packing count covering exactly `m` vertices. Equivalently
`per(xI+A) = sum_L x^{n-|V(L)|}` — the unsigned twin of the Sachs coefficient theorem
`det(xI-A) = sum_L (-1)^{#cyc(L)} x^{n-|V(L)|}`, `e_m^{signed} = sum_{|V(L)|=m}(-1)^{#cyc}`.

**Proof.** `per(xI+A) = sum_{sigma in S_n} prod_i (xI+A)[i,sigma(i)]`. Decompose `sigma` into
cycles: a fixed point contributes `(xI+A)[i,i] = x`; a cycle of length `>= 2` contributes
`prod A[i,sigma(i)] in {0,1}`, equal to 1 iff that cycle is a directed cycle of `T`. In a
tournament there are no 2-cycles, so the surviving `sigma` are exactly the linear subdigraphs
`L` (non-fixed cycles = the cycles of `L`), each contributing `x^{#fixed} = x^{n-|V(L)|}`. ∎
(The same cycle decomposition gives `det(xI-A)`'s Sachs form, with `sgn(sigma)` producing the
`(-1)^{#cyc}`.) VERIFIED: Ryser `per(xI+A)` over `Z[x]` equals the independent packing
enumeration, all samples n=3..7 (`permanental_companion_monad.py` [1]).

## Part 2 — det/per is the cycle-parity SIGN channel; the O/E split (PROVED)

`det(xI-A)` and `per(xI+A)` are the **signed** and **unsigned** faces of `Phi`, differing
only by `(-1)^{#cyc} -> (+1)^{#cyc}`. Hence the pair determines, for every `m`,

> **O_m := (e_m^{unsigned} - e_m^{signed})/2 = #{ odd-cardinality packings on m vtx }**,
> **E_m := (e_m^{unsigned} + e_m^{signed})/2 = #{ even-cardinality packings on m vtx }.**

VERIFIED (all samples n<=7, `permanental_companion_monad.py` [2]). The permanent is the
*spectrum with the cycle-parity signs stripped*: `det` is a symmetric function of the
eigenvalues (spectral, polynomial-time), `per` is the non-spectral unsigned cycle-cover
count (Valiant: per of a 0/1 matrix is `#P`-hard — the ambient structural reason the
spectrum cannot see the unsigned data; not a hardness claim for the tournament-restricted
permanent).

**Credit / prior work.** opus-2026-03-14-S85 (`04-computation/perm_det_23.py`,
`permanent_spectral_bridge.py`) compared the **bare numbers** `per(A)` and `det(A)` of the
adjacency, finding `per(A)=det(A)` at n=5, `per(A)=-det(A)` at n=4, and tying per's
`#P`-hardness to boson sampling. Those numbers are exactly the **top (`m=n`) coefficient** of
the present polynomials: `per(A) = e_n^{unsigned} = O_n+E_n` (spanning packings, unsigned) and
`det(A) = (-1)^n e_n^{signed}` (spanning packings, signed `sgn(σ)=(-1)^{n-#cyc}`), so
`per(A)=±det(A)` exactly when all spanning packings share one cardinality-parity (e.g. n=5:
only single 5-cycles, card 1; n=4: only single 4-cycles). THM-506 lifts that single-number
comparison to the **full polynomial** `per(xI+A)` and places it inside the master-polynomial
`Φ` framework. (THM-115 uses "permanent" loosely for the H=21 gap via `I(Ω,2)`; unrelated to
the permanental polynomial.)

## Part 3 — same wall, strictly finer fingerprint (VERIFIED)

The permanental polynomial is constant on cospectral classes for `n <= 5` and first splits a
cospectral class at **n = 6** (the disjoint-pair threshold `6 = 3+3`) — the **same wall** as
`H`, `c_6`, and every face of `Phi`. At `n = 7` it splits the **same 47** cospectral classes
`H` does (both 3/28 at n=6, 47/168 at n=7). But it is strictly **finer**: `(det, per)`
resolves the full carrier vector via `O_m, E_m` (e.g. `O_6 = c_6`, `E_6 = D_33`, recovering
`c_6` and `D_33` separately from the spectral `e_6^{signed}=D_33-c_6`), whereas
`H = I(Omega,2)` reads only the fugacity-2 functional `4c_6+2c_7+...`.

Exhaustive `n = 6` (all 32768 labeled tournaments; `permanental_companion_monad.py` [4]):
iso classes (A000568) = 56; distinct char poly = 28; distinct `(char, perm)` = **32**;
distinct `(char, perm, H)` = 32. The permanent recovers 4 distinctions the spectrum loses;
`H` adds nothing beyond the permanent. **The pair of characteristic and permanental polynomials dominates `H` as a
tournament fingerprint for n<=8; this implication fails from n=9.**

## Part 4 — recovery through eight; first functional failure at nine

**Corrected statement.** The pair `(det(xI-A), per(xI+A))` determines the
Hamiltonian-path count H on all n-vertex tournaments if and only if n<=8.
The n=8 implication is COMPUTER-ASSISTED / FINITE-EXACT, not a claimed closed
formula; the n>=9 failure has an explicit finite witness and an all-order proof.

### Why the old rank inference failed

The former text claimed that functional recovery holds if and only if H is
affine in the permanent coordinates within a characteristic-polynomial class.
This is false on the finite set of realized tournaments. A rank increase
rules out an affine function, not an arbitrary function on the realized image.
Even if both observables are linear on a larger carrier space, the realized
set need not contain the differences or kernel vectors of that space.

The old n=8 pooled exact rank 3->4 is independently confirmed. Nevertheless
all 3217 actual joint-polynomial fibers have a single H. In fact, 59 of the
1523 characteristic-polynomial classes fail affine recovery individually.
A smallest three-point affine obstruction is given by upper-pair-lex edge
masks `172254133,172256153,172238829`, with bit 1 meaning i->j for i<j.
All have characteristic coefficients, in ascending powers of x,

```
[-2,-9,-20,-28,-22,-12,0,0,1].
```

The three permanent low coefficient triples are `(12,31,40)`, `(16,39,44)`,
`(18,43,46)`; the remaining coefficients are `[28,22,12,0,0,1]` in each case.
Their H values are `177,197,209`. Permanent differences `(4,8,4),(6,12,6)`
are proportional, whereas H differences `20,32` are not. All three rows
were checked independently by direct enumeration of all 8! permutations.
This is non-affine but still functional recovery, not a same-coordinate collision.

### The finite-exact n=8 proof and the low-order formula

Write D_ab for disjoint directed-cycle pairs of lengths a,b. The established
odd-cycle packing identity gives, with missing O/E indices read as zero,

```
n<=7: H = 1 + 2(O3+O5+O7) + 4E6.
n=8:  H = 1 + 2(O3+O5+O7) + 4(E6+D35),   E8=D35+D44.
```

The formal coordinate map merges D35 and D44 at n=8, but no two realized
n=8 tournaments keep all retained coordinates equal while changing D35.
Every tournament has a Hamiltonian path, so relabeling one as `0->1->...->7`
reduces the exhaustive universe to all 2^21 choices of the other arcs.
Two independent exact algorithms cover this entire universe:

1. Row-assignment determinant/permanent DP and endpoint Hamiltonian-path DP.
2. Minimal-vertex-rooted cycle counts, all 651 disjoint cycle-pair placements,
   and the odd-cycle packing formula.

Both return 3217 joint-polynomial fibers and zero unequal-H fibers. The
second method also confirms that every fiber fixes D35 and D44 individually.
These are polynomial fibers, not isomorphism classes. Audit checksum: the
presentation-H sum is 815136768, with minimum 1 and maximum 661.

### An actual collision at n=9

Upper-pair-lex masks `29759988540` and `18875161993` have H counts 2843 and
2841, respectively, but share the following coefficient arrays (ascending):

```
per(xI+A): [384,692,629,380,154,68,28,0,0,1]
det(xI-A): [-20,-82,-167,-192,-154,-68,-28,0,0,1].
```

Independent direct enumeration of all 9! permutations confirms each array
and both H counts. Their sidecars, first/second, are
`D35=209/210`, `D44=96/95`, `T333=11/10`, `c9=191/192`.
Thus `E8=305` and `O9=c9+T333=202` agree. The exact formula

```
H = 1 + 2(O3+O5+O7+O9) + 4(E6+D35) + 6T333
```

gives the difference `4*(-1)+6*(1)=2`. This is a genuine realized collision.

Adjoining a source multiplies both polynomials by x and leaves H unchanged,
since the source must begin every Hamiltonian path and lie in no directed
cycle. Repetition propagates the collision to every n>=9. It also shows
that any lower-order collision would propagate to n=8, contradicting the
complete census. This proves minimality of nine.

### Correction lineage and remaining question

The original monad-S6 cycle-packing theorem and O/E split remain valid.
The original rank table described affine span, not arbitrary recoverability;
its n=8 failure label and its formal D44/D35 "collision" interpretation are
superseded by the realized-fiber census above. Historical scripts cited below
were unavailable in the audited snapshot; the replacement scripts and receipts
are self-contained. A structural explanation or compact nonlinear formula for
the exceptional n=8 recovery is a useful open question, not needed for the
finite-exact result.

See [proofs and independent receipts](../../05-knowledge/results/session37_inheritance_audit_20261023/REPORT.md), including the complete n=8 fiber
TSV, exact rank audit, explicit n=9 matrices and direct permutation verifier.

## Files
- [Session 37 replacement proof and receipts](../../05-knowledge/results/session37_inheritance_audit_20261023/REPORT.md) — corrected Part 4;
  historical script listings below preserve attribution rather than assert availability.
- `04-computation/permanental_companion_monad.py` (+ `05-knowledge/results/permanental_companion_monad.out`)
  — Parts 1,2,3 + fingerprint completeness.
- `04-computation/perm_rank_test_monad.py` (+ `.out`, `perm_rank_test_n9.out`) — Part 4.

## Even-length face — RESOLVED (monad-S7, HYP-2517)
The even-only face DOES have a clean matrix function — but in its **SIGNED** form, and
that form is SPECTRAL. The skew-adjacency char poly `det(xI − S)`, `S = A − Aᵀ`, is the
signed even-length face: by Coates, every cycle cover containing an odd cycle cancels its
one-cycle-reversal (`w_rev=(−1)^len w`), leaving only even covers, so
`det(xI−S)=∏(x²+μ_j²)` with coeffs `p_{2t}(S)=Σ_{|W|=2t}Pf(S[W])²` — exactly the
"Pfaffian on a derived graph." **It is SPECTRAL** (`#distinct(charA)=#distinct(charA,charS)`
at every n, incl. the cospectral-different-H pairs at n=6; VERIFIED n≤6 exhaustive, n=7
sample; matches the known spectral-DS-of-tournaments equivalence, complement=converse).
So the Pfaffian route recovers only the spectral shadow; the genuine **non-spectral** even
content `I(Ω_even,·)` (splits 3@n6, 46@n7) has NO determinantal home, exactly like H.
**The det/per (Valiant signed/unsigned) dichotomy = the spectral/non-spectral boundary,
face by face**: signed→spectral (det A, det S), unsigned→non-spectral (per, I(Ω_even), H);
the ODD face has no determinantal object at all (the irreducibly non-spectral one).
See reflection `the-skew-determinant-is-the-signed-even-face-and-it-is-spectral-monad-s7`,
HYP-2517, THM-174 (`det S = det(I+2A) = Pf(S)²` = the top coefficient).

## Open
- The permanental ROOTS as a non-spectral invariant (cleaner carrier of (c6,c7,…)?).
  The general-n carrier deficit of `(char,perm)` past the first collision. A clean general
  proof that walk counts `w_k=1ᵀAᵏ1` are spectral (Moon for low k; HYP-2517). A complex
  signing `S_{ij}=ω` filtering lengths mod r>2 — still a determinant, conjecturally
  spectral? (OPEN-Q-096; HYP-2515; HYP-2517.)
