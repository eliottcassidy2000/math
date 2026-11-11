# Moving the Gaussian centre fills the mixed-derivative gap

Session269, 2026-11-11. **JC(2) remains OPEN.**

**New elementary result:** pointwise first-hit vanishing in both pure derivative
rays, throughout a nonempty open set of Gaussian centres, forces mixed finite
chaos. The derivative orders and the degrees of the membership witnesses need
not be bounded in advance. An explicit rational-certificate operation then
propagates the first hits to complete mixed quadrants.

**Conditional JC consequence:** taking the dimension comparison and finite-chaos
criterion of THM-5138 as working hypotheses, this proves detection with a moving
centre and collapses its universally quantified K-GMC/CV/SN ladder to JC(2).
THM-5138 is explicitly **PRE-AUDIT DRAFT** at the source pin below; its analytic
comparison is not re-proved here. The new propagation theorem does not depend
on that comparison, GMC, a classification of Keller maps, or D-module machinery.
A refinement handles one centre algebraically independent over the coefficient
field, but arbitrary prescribed-centre SN detection and JC(2) remain open.

## 1. Inheritance, niche, and concept board

Read-only source pin: `origin/main = 1ca9c5f83b55bcac7dacc0a6f5def38e02f245e7`.
The dirty shared checkout was not modified. Startup/navigation/protocol and
active guardrails were read from that pin; coordinator handles repository startup.

- Closest mechanism: THM-5138, `01-canon/theorems/THM-5138-stein-modules-of-keller-pairs-dimension-count-and-finite-chaos-criterion.md`, especially (A), (C), (D).
- Precise open target: `05-knowledge/hypotheses/HYP-9668-gaussian-ladder-for-keller-pairs-between-gmc2-and-jc2.md`.
- Underused sidecar: THM-4620, `01-canon/theorems/THM-4620-keller-dual-field-module-kernel-ladder-and-de-rham-cokernel.md`: the dual Weyl action and its inherited finite generation/cyclicity. Only the optional finite-test sharpening uses that cited fact.
- Recent operation: Session267's parameter-family charge rigidity, in `05-knowledge/results/session267_gaussian_planar_20261111/SYNTHESIS.md`. Here the family parameter is the actual Gaussian centre, not a charge-polynomial coefficient.
- Canonical hostile: the etale map on `C* x C`, `(u,v)=(x^2,y/x)`, where one pure direction is locally nilpotent but the other has nonvanishing residue moments. Source THM-5138(E).
- Near miss: at one centre, separate eventual row/column zeroes do not force a finite two-dimensional array. A diagonal array is the abstract hostile; it is NOT claimed realizable by a Keller pair.

Live board: (1) pure-ray Stein first hits; (2) mixed Hermite tails;
(3) centre differentiation; (4) constructible membership and countable covers;
(5) finite Weyl generators. The useful connection is 3+4: first hits at moving
centres become rational certificates, which can actually be differentiated.

## 2. Setup and main theorem

Let `B=C[x,y]` and `F=(f,g)` satisfy `{f,g}=1`. Set

    D = g_y partial_x - g_x partial_y,
    E = -f_y partial_x + f_x partial_y.

Thus `Df=Eg=1`, `Dg=Ef=0`, and `[D,E]=0`. For a centre `(a,b)` set

    U=f-a, V=g-b, T=D-V, S=E-U,
    M_(a,b)=T B + S B, Q_(a,b)=B/M_(a,b).

The quotient is a vector space, NOT a quotient algebra. We never multiply a
zero class by an arbitrary polynomial. Write `[h]_(a,b)` for its class.

**Theorem 1 (PROVED; moving-centre first-hit propagation).**
Suppose there is a nonempty Euclidean open set `Omega` in `C^2` such that for
every monomial `h=x^i y^j` and every `(a,b)` in `Omega` there are positive
integers `r,s` (depending on all this data) with

    [D^r h]_(a,b)=0,   [E^s h]_(a,b)=0.                 (FH)

Then for each monomial `h` there exist positive integers `R_h,S_h` and a
nonzero polynomial `P_h(a,b)` such that, at every centre with `P_h(a,b) != 0`,

    [D^m E^n h]_(a,b)=0
          whenever m >= R_h or n >= S_h.             (RECT)

In particular, outside a countable union of proper algebraic curves there
is finite chaos: for every `h in B` there is `N(h)` with

    [D^m E^n h]_(a,b)=0 whenever m+n >= N(h).

One can use `N(x^i y^j)=R_h+S_h-1`, and take a maximum over the finitely many
monomials of a general polynomial. The first-hit hypothesis is weaker than
assuming eventual pure-ray vanishing at each centre.

The same proof works if each required first-hit set is nonmeagre in `C^2`;
it suffices that it cannot lie in a countable union of proper algebraic subsets.
One must not replace this with bare Zariski density or with uncountability alone.

### 2.1 Finite membership is constructible

Fix a polynomial `k in B`. For integer `L>=0`, membership with degree-bounded
witnesses is the condition

    k = T p + S q,   deg_(x,y) p, deg_(x,y) q <= L.     (1)

In monomial bases this is a finite linear system

    A_L(a,b) c = v_k.

The entries of `A_L` are polynomials in `a,b`; include all target monomials
that occur in either side of (1). Solvability is the constructible condition
`rank A_L = rank [A_L | v_k]`, expressed by finitely many minors and rank strata.
Denote this constructible set by `Z_(k,L)`.

If `Z_(k,L)` is Zariski dense, then (1) has a solution over `K=C(a,b)`.
Indeed density forces equality of the two generic ranks: if the augmented
rank were larger generically, its maximal nonzero minor and one generic-rank
minor of `A_L` would give a nonempty open set of nonsolvability, contradicting
density of the constructible solvability set. Gaussian elimination over `K`
then supplies `p,q in K[x,y]`.

Fix `h` and the D direction. Hypothesis (FH) says

    Omega is contained in union_(r>=1,L>=0) Z_(D^r h,L).

If none of these constructible sets were Zariski dense, their proper
algebraic closures would be nowhere dense in the Euclidean topology. A
nonempty Euclidean open set, being a Baire space, cannot be covered by
countably many such sets. Consequently some `r=R_h,L` supplies a rational
certificate

    D^R_h h = T p_h + S q_h,   p_h,q_h in K[x,y].       (2)

Apply the same argument to the E direction to obtain

    E^S_h h = T p'_h + S q'_h.                         (3)

Crucially, no uniform exponent, witness degree, finite-dimensionality of
`Q_(a,b)`, flatness, or specialization theorem was used.

### 2.2 Flat centre connections and an explicit certificate compiler

On `K[x,y]` define

    nabla_a = partial_a + V,
    nabla_b = partial_b + U.

They commute with each other and with both `T,S`. For example,

    [partial_a,S]=1,   [V,E]=-1,

so `[nabla_a,S]=0`; the other identities are equally direct. Thus they
preserve `M_K=T K[x,y]+S K[x,y]` and induce a flat connection on its quotient.
For any parameter-independent `k in B`, integration by parts gives

    nabla_a[k]=[Vk]=[Dk],
    nabla_b[k]=[Uk]=[Ek].                             (4)

Here is a constructive version that does not rely on quotient notation.
If parameter-independent `k` has rational certificate `k=T p+S q`, then

    Dk = T(k+nabla_a p) + S(nabla_a q),                (5)
    Ek = T(nabla_b p) + S(k+nabla_b q).                (6)

For (5), apply `nabla_a` to the certificate and use
`Dk=Vk+Tk`; (6) is symmetric. Formulae (5)-(6) are the new local operations.
Every differentiated target remains parameter-independent. Repeatedly applying
them to (2) proves every mixed identity with `m>=R_h`; applying them to (3)
proves those with `n>=S_h`.

Choose `P_h` to clear all denominators of the four initial certificates.
Multiplication by `U,V` and differentiation of rational functions introduce
only higher powers of those same denominator factors. Therefore ALL the
compiled identities specialize at EVERY centre where `P_h != 0`. This proves
(RECT), not merely generic membership separately for each derivative.

There are countably many monomials, so the union of their denominator zero
sets is a countable union of proper algebraic subsets and cannot cover an open
set. Choose a centre outside it. Finite linear combinations establish the
last assertion of Theorem 1. QED.

## 3. A finite Weyl-generator strengthening

For a fixed centre, let `FC` be the vector subspace of all `h` with mixed
finite chaos. It is stable under `D,E`, and also under multiplication by `f,g`.
To see the last assertion without illegally multiplying quotient classes,
put `k_(m,n)=D^m E^n h`. Leibniz and integration by parts give

    [D^m E^n(fh)] = [k_(m,n+1)] + a[k_(m,n)]
                    + m[k_(m-1,n)],
    [D^m E^n(gh)] = [k_(m+1,n)] + b[k_(m,n)]
                    + n[k_(m,n-1)].                  (7)

Terms with a zero prefactor are omitted. If `h` has threshold `N`, then
`fh,gh` have threshold at most `N+1`. Hence `FC` is a submodule for the Weyl
action generated by `f,g,D,E`.

**Theorem 2 (PROVED conditional on the stated generating set).** If
`h_1,...,h_l` generate `B` for that Weyl action, it suffices to assume (FH)
only for these finitely many polynomials. Theorem 1's argument produces
rational first-hit certificates for each generator. On the nonempty Zariski
open set where their finitely many denominators do not vanish, every generator
has finite chaos; (7) makes the whole module have finite chaos there.

The inherited holonomicity/cyclicity assertion in THM-4620 supplies such a
finite set (indeed an existential cyclic generator), but this note does not
construct that generator from arbitrary forward Keller data. The theorem
is directly usable whenever a generating set is supplied. This is NOT a
claim that the algebra generators `x,y` are Weyl-module generators.

## 3a. A single centre that is generic over the coefficient field

There is a stronger precise single-centre version. Let `k` be the subfield of
`C` generated over `Q` by the finitely many coefficients of `f,g`. It is
countable. A pair `(a_0,b_0)` algebraically independent over `k` lies outside
every proper algebraic subset defined over `k`.

**Theorem 2a (PROVED; very-general-centre first-hit detection).** Suppose
`a_0,b_0` are algebraically independent over `k`. If for every monomial `h`
there are positive integers `r,s` such that

    D^r h, E^s h belong to M_(a_0,b_0),

then `F` has mixed finite chaos at that SAME centre.

Proof: fix a first hit and a finite witness degree `L`. Now regard the centre
variables as indeterminates `A,B`. The matrix `A_L(A,B)` and augmented vector
in (1) have coefficients in `k`. Evaluation at `(a_0,b_0)` embeds `k(A,B)`
into `C`; every nonzero polynomial minor stays nonzero. Consequently the
matrix and augmented matrix have exactly their generic ranks at this centre.
Solvability over `C` therefore implies solvability over `k(A,B)`. This gives
one rational certificate for each first hit, with coefficients in `k(A,B)`.
Apply the compiler (5)-(6). Every denominator is a nonzero polynomial in
`k[A,B]`, so none vanishes at `(a_0,b_0)`, even for the countably many monomial
certificates. The resulting rectangles establish finite chaos at the original
centre. QED.

This argument does not require the original membership witnesses to have
coefficients in `k(a_0,b_0)`: solvability of a finite linear system descends
from the extension field `C` by the rank criterion.

The word **very general** is essential: this excludes a countable union of
proper algebraic subsets depending on the coefficient field, not merely the
finite exceptional set/curves excluded by the inherited topological count.
Thus it is stronger than the open-centre theorem in its input requirement,
but it does not prove detection at an arbitrary prescribed centre. Under the
inherited THM-5138 finite-chaos criterion, it gives automorphism whenever this
same very-general centre is also eligible for that criterion. One can always
choose centres satisfying both conditions. This is a precise single-centre
version of HYP-9668(D), rather than a proof for every centre called generic
in a weaker geometric sense.

## 4. Consequences for the actual complex planar Jacobian question

**Conditional Corollary 3.** Assume THM-5138(B,C), with its stated comparison
hypotheses and eligible generic centres. Then any polynomial Keller pair
satisfying Theorem 1's moving-centre first-hit hypothesis is an automorphism.

Proof: choose a finite-chaos centre from Theorem 1 outside the algebraic
exceptional set excluded in THM-5138(B). Its comparison gives
`dim Q >= 2d-1`, while the elementary finite-chaos module argument gives
`dim Q <= d`, so `d=1`, and the inherited criterion concludes automorphism.
Equivalently, apply THM-5138(C) directly. No claim about maps on punctured
surfaces is substituted for this polynomial `C^2 -> C^2` statement.

**Conditional Corollary 4 (universal ladder collapse).** With the same
inherited input, the following universally quantified statements are equivalent:

1. JC(2).
2. K-GMC for every polynomial Keller pair and every Gaussian centre.
3. CV for every such pair and centre.
4. SN for every such pair and centre.
5. The two pure directions have a positive first hit for every monomial,
   for every such pair and centre.

The forward chain is inherited THM-5138(D), using the inherited GMC theorem
only for `1 => 2`; `4 => 5` is immediate and `5 => 1` is Corollary 3.
It is also enough to posit the universal statements at the single fixed
centre `(0,0)` FOR EVERY KELLER PAIR: apply that statement to the translated
pair `(f-a,g-b)` for each `(a,b)`. Their dual fields are unchanged and their
zero-centre kernel is exactly `M_(a,b)` for the original pair.

The last observation is a quantifier distinction. It does NOT say SN at
one centre for one fixed pair implies invertibility. That remains unproved
here for arbitrary prescribed centres; Theorem 2a covers a precisely specified
very-general centre. The new result removes the need for a separate one-centre detection
theorem in the UNIVERSAL K-GMC/CV/SN programme.

A useful contrapositive is explicit: subject to the inherited criterion,
a hypothetical noninvertible polynomial Keller pair has some monomial `h`
and one direction for which the first-hit locus is contained in a countable
union of proper algebraic subsets of centre space. Indeed, if every monomial
and both directions had a Zariski-dense finite-witness solvability stratum,
the certificate proof would give finite chaos and a contradiction. Thus some
pure ray is never killed at a very general centre. This is an existential
obstruction profile, not an effective algorithm selecting its monomial.

## 5. Controls, scope, and the lost-coordinate ledger

- Identity and every polynomial automorphism pass: dual fields are ordinary
  coordinate derivatives transported by the inverse, so every polynomial
  eventually has zero ACTUAL derivatives in both directions.
- A single pure direction is insufficient on general etale surfaces:
  `(u,v)=(x^2,y/x)` on `C* x C` has `E=x partial_y` locally nilpotent, while
  for `h=x^-2` the inherited residue functional gives
  `[D^m h]` detected by `(-1)^m m!/(2 a^(m+1)) != 0` whenever `a!=0`.
- Mere separate eventual vanishing of a numerical array misses its diagonal.
  The extra sidecar here is a certificate that varies rationally in the
  centre, not a scalar estimate at one centre.
- A countable Zariski-dense centre sample is insufficient for the Baire
  argument: it can itself be a countable union of closed points.
- A finite collection of tests cannot establish (FH) over an open set.
- Denominator powers may grow, but denominator LOCATIONS do not. This is why
  specialization commutes with the whole infinite compiled tail.
- No effective bound on `R_h,S_h,L` follows from the countable-cover argument.
- No external literature priority is claimed. The only new dependency used
  for the conditional JC consequence is explicitly inherited, not upgraded.

Connection passport: source = actual two-parameter Stein membership family;
target = mixed Hermite vanishing; map = rational-certificate extraction followed
by (5)-(6); preserved predicate = exact image membership; destroyed information
= effective bounds during the Baire extraction; restoration sidecar = explicit
initial certificates and their denominator polynomial; cheapest decisive test
= symbolic commutator/certificate checks, plus the punctured-surface residue
hostile. The source-flag and arbitrary high-dimensional special-image gaps in
Session267 are not bypassed or used.

## 6. Reproduction and next theorem

Run `python3 verify.py` in this directory; its output is `verification.txt`.
It checks the explicit compiler and commutators on a nonlinear shear-word
Keller automorphism using exact symbolic arithmetic, finite-chaos derivative
bounds in target coordinates, and the nonzero hostile residue formula. It
checks only the new algebraic mechanism, not the inherited comparison theorem.

Best next task: find forward-data rules that construct the first-hit rational
certificates for a supplied Weyl generator, or prove their existence throughout
a centre open set. A single certificate per direction is enough; one no longer
needs to guess an entire mixed Hermite tail. This is the remaining substantive
problem, and it is equivalent in the stated sense to the open JC programme.
