---
id: THM-3687
title: "Russell-cylinder exceptional-quartic actual-ring J0 lift"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  Over the
  irreducible quartic field cut out by THM-3683's sixth-debt polynomial, all
  four exceptional degree-eight folds have one uniform actual-target-ring
  lift through J_0=1.  A split finite-field fibre selects a 218-square only;
  the square is rebuilt over the quartic field, expanded to an 872-square
  over Q, solved exactly, and substituted into every one of the 219 full
  coefficient rows.  Split-prime cutoff-195 failure and cutoff-198
  feasibility are finite-field controls only.  No J_1, J_2, Keller-pair, or
  JC(2) claim is made.
source: jc_zero_debt_lift / exceptional-quartic actual lift, 2026-08-22
audit: >
  PASS -- root independently checked actual target-monomial typing, the
  G_1/F_1 column signs, the 218-field-square to 872-rational-square power-
  basis expansion, exact contraction, all 219 residual rows, and the strict
  finite-field-only boundary; a fresh normal replay matched the frozen
  hashes and transcript without correction.
depends_on:
  - THM-3683-russell-cylinder-sixth-debt-quartic-on-the-zero-fourth-parabola
related:
  - THM-3678-russell-cylinder-qdagger-actual-j0-lift
  - THM-3680-russell-cylinder-qdagger-coupled-stable-lift
script: 04-computation/jc2_russell_cylinder_exceptional_quartic_exact_j0_lift_thm3687.py
auxiliary_script: 04-computation/jc2_russell_cylinder_exceptional_quartic_modular_lift_thm3687.py
output: 05-knowledge/results/jc2_russell_cylinder_exceptional_quartic_exact_j0_lift_thm3687.out
script_sha256: ee89640ffa0d66a7ec044a9d4c870859e0fbf0c28df36d89cea4acb5d7433f39
auxiliary_script_sha256: 4630562017fb0378a24e65d5db4da10187172789d92202fe8512a3f6ae67c014
output_sha256: 16d3b7396213b2bedf4dd48d58ad363c417e3d9994f49889eb764691809cbcb2
hash_basis: raw LF bytes
---

# THM-3687 -- all four sixth-debt folds pass the actual `J_0` gate

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  The four
algebraic folds left by THM-3683 are not merely retained-jet survivors.  One
actual target-ring certificate over their common quartic field proves
`J_0=1` uniformly at all four complex embeddings.

All rings below have characteristic zero.

## 1. The four exceptional folds as one field-valued fold

Put

```text
F_6(T)=72783360T^4-77822208T^3-28419741T^2
                         +7849770T-1276420,             (1)

K=Q[alpha]/(F_6(alpha)).                                (2)
```

THM-3683 proves that `F_6` is irreducible and that its four complex roots are
exactly the degree-at-most-eight principal folds surviving the complete
retained order-six arbitrary-two-form screen.  In the notation of that
theorem, set

```text
P=x^2(x^2-1)^2,
R_1=P(1-x^2),                 R_2=P(4-9x),

p(alpha)=520alpha^2/9-1688alpha/81-5717/729,
Q_alpha=Q_6+p(alpha)R_1+alpha R_2.                     (3)
```

The exact companion independently gates

```text
Q_alpha(-1,0,1)=(-3,-3/4,-3),
Q_alpha'(-1,0,1)=(-9/2,1,9/2),                         (4)
```

as well as the zero-second, zero-fourth, and zero-sixth debt identities.  The
last is the literal field relation `F_6(alpha)=0`; none is inferred from a
floating-point root.

Use the Russell-cylinder target ring and compiler over `K`:

```text
R_K=K[b,c,e]/(c^2e-b(b+4)),

D=1+x^2q,
b=(D-1)(D+2)^2,       c=xD(D+2),       e=q(D+3),       (5)
```

and let `gamma` denote restriction through `q=Q_alpha(x)`.  The restriction
degrees are

```text
(deg gamma(b),deg gamma(c),deg gamma(e))=(30,21,18).   (6)
```

## 2. Exact actual-ring lift

Take zero-stable target coefficients

```text
U=c,                         V=e.                       (7)
```

There exist actual elements `F_1,G_1 in R_K` such that

```text
gamma(c)' gamma(G_1)-gamma(F_1) gamma(e)'=1
                                      in K[x].          (8)
```

Equivalently, for target functions beginning

```text
F^#=c+wF_1+O(w^2),             G^#=e+wG_1+O(w^2),      (9)
```

the constant stable coefficient of the pulled source Jacobian is the full
polynomial identity

```text
J_0=1.                                                     (10)
```

The words *actual elements* are load-bearing.  Use the raw target packet

```text
M_198={b^i c^j e^k:30i+21j+18k<=198}.                  (11)
```

It contains `187` monomials.  The two coefficient vectors in `(8)` are
indexed by these target monomials, not by arbitrary source polynomials.  The
frozen certificate has

```text
                  target support   restriction degree
F_1                    108                  198
G_1                    107                  195

delta(F_1) degree 190,              delta(G_1) degree 187. (12)
```

Their exact hashes are

```text
F_1 target: d3ac53e05a433a4d9d0a7d52b138d6c2e998341c8dc7477f20a79646221afc0c,
G_1 target: 0a37d8b98f192a78e5774c2d73209011615024e41654ad5a444db4373537f61d,

gamma(F_1):221deea140f9b010f6672eeb3bc4bbf6f126b6ac6f2e6b8263abf70b2476dcca,
gamma(G_1):f310fce20fecb88b151ad9dfda5300997a0287a84b0ed674d05327b71f562bf9. (13)
```

Each field coefficient is serialized in the ordered power basis
`(1,alpha,alpha^2,alpha^3)` before hashing.

## 3. Why the computation is characteristic-zero exact

The coefficient matrix for `(8)` has

```text
219 rows,                    374 columns.               (14)
```

The split fibre `(p,alpha)=(137,44)` selects `218` columns and `218` rows.
No coefficient is lifted from that finite field.  Every selected entry is
rebuilt in `K`, and multiplication by one element of `K` is expanded in the
power basis `(1,alpha,alpha^2,alpha^3)`.  Thus the selected field square
becomes an

```text
872 by 872 rational square.                             (15)
```

Exact rational solution of `(15)` proves that the selected square is
invertible over `K`.  Contracting the solution back to `K`, placing it on the
actual target monomials `(11)`, and substituting into **all** `219` rows gives
the zero residual in `(8)`.  Therefore the modular selector is only a sparse
index oracle; the proof itself is over `Q` and `K`.

Because `F_6` is irreducible, each of its four complex roots defines an
embedding `K -> C`.  Base-changing `(8)` along any of these embeddings gives
an actual-ring `J_0=1` lift for the corresponding complex fold.  This is the
precise all-four-root uniformity.

## 4. Hostile modular controls and their boundary

At each of the four roots in each completely split good prime `137` and
`163`, the raw finite-field packets give

```text
cutoff 195: operator rank 214, augmented rank 215,
cutoff 198: operator rank 218, augmented rank 218.       (16)
```

These eight independent residue fibres detect the same transition and guard
against a root-specific or bad-prime accident.  They are **FINITE-FIELD
CONTROLS ONLY**.  A modular rank failure does not prove characteristic-zero
infeasibility when a solution denominator may meet that prime, so `(16)` is
not a characteristic-zero minimality theorem.  The positive exact statement
is only feasibility at cutoff `198`, proved by Sections 2--3.

## 5. Strict stopping boundary

This theorem advances every THM-3683 exceptional fold through the first
actual target-ring gate.  It does not prove

- `J_1=0` or `J_2=0` over `K[x]`;
- a finite continuation through all stable orders;
- decomposability beyond the displayed first target coefficients;
- a noninjective polynomial Keller pair; or
- a counterexample to `JC(2)`.

The cheapest next test is the coupled `J_1,J_2` system over `K`.  Split-prime
experiments may select its squares, but only a characteristic-zero rebuild and
full residual substitution can promote that continuation.

## 6. Reproduction

```bash
python3 -B 04-computation/jc2_russell_cylinder_exceptional_quartic_exact_j0_lift_thm3687.py
python3 -O -B 04-computation/jc2_russell_cylinder_exceptional_quartic_exact_j0_lift_thm3687.py
```

Both modes must match the stored LF transcript byte-for-byte.  The companion
also runs the split-prime controls in `(16)` and contains no Python `assert`
statements.  **QED.**
