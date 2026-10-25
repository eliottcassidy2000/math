# Exact outer Picard rank 19 and hourglass Picard rank 58

**Status: PROVED over C; independently audited.** This is the current exact
conclusion. `LOWER_BOUND.md` remains a valid, independently closed lower-bound
snapshot; its two-value set is superseded here by the separately proved upper
bound. There is no remaining exact-rank blocker.

## Theorem

Let S_o be any of the six outer-deletion six-root K3 surfaces of the hourglass,
and let Z be the smooth K3 model of

    W²=E(E²−x²)(E²−z²)(E+x+z).

Let Y_7 be the minimal resolution of the seven-root hourglass surface. Then

    rho(Z)=rho(S_o)=19,
    rho(Y_7)=58,
    rank T(S_o)=3,
    rank T(Y_7)=20.                                    (1)

The seven six-root deletion K3 surfaces form **exactly two complex
isomorphism classes**: the central deletion of Picard rank 20, and the six
isomorphic outer deletions of Picard rank 19.

The elliptic fibration on Z with d=1−t² and equation

    Y²=(X−d)(X+d)(X−d(t+1))

has Mordell–Weil rank one over both C(t) and Q(t). Its section

    P=(dt,d²)

is nontorsion. This theorem does **not** assert that P is primitive or is a
Mordell–Weil generator; no height or complete integral-lattice computation is
being claimed.

## Proof, with independently checkable lower and upper halves

### A. Nineteen classes from the elliptic fibration

The affine substitution `t=x, X=−dz, Y=dw` identifies the generic fiber of
Z with the displayed cubic. Its minimal singular fibers are

    I2 at t=0,−2;  I0* at t=1,−1;  I2* at infinity.

Consequently its trivial lattice has rank
`2+1+1+4+4+6=18`. The section P specializes at the smooth fiber t=2,
after the integral shift u=X+3, to

    E_6: Y²=u³−36u,  Q=(−3,9),  2Q=(25/4,−35/8).

The duplication identity

    u(2R)=(u(R)²+36)²/[4u(R)(u(R)²−36)]

implies `v_2(u(2^n Q))=−2n` for every n>=1, so Q and therefore P are
nontorsion. Shioda–Tate gives rho(Z)>=19.

All minimality, infinity, specialization, and valuation details are proved
in [the lower-bound packet](LOWER_BOUND.md), with 87 exact checks in both
normal and optimized Python runs. In particular, the two lifts of z=−t
are P and −P, not two independent Mordell–Weil classes.

### B. Twenty classes would have incompatible reduction discriminants

The [independent upper-bound packet](../picard-upper/INDEPENDENT_UPPER_BOUND.md)
constructs a smooth proper K3 model over Z[1/6] by blowing up the branch
triples and then its ordinary crossings. It proves that all nine A1 and
eight D4 exceptional components are rational over the residue field. Together
with the pullback of a general line, they give a nondegenerate rank-18
rational lattice. Resolution adds precisely 17q to the double-plane point
count over F_q.

Two independently written exact counters agree on

| p | #Z(F_p) | #Z(F_(p²)) |
|---|---:|---:|
| 7 | 180 | 3384 |
| 13 | 408 | 32088 |

Removing the eighteen known Tate eigenvalues leaves a four-dimensional
orthogonal Frobenius problem. Its first two traces give nonzero middle
coefficients −42 and −234, respectively, which force the positive
reciprocity sign. The resulting full characteristic polynomials are

    (T−7)^20(T²+10T+49),
    (T−13)^20(T²+22T+169).                              (2)

The remaining normalized quadratic roots are not roots of unity: their
sums −10/7 and −22/13 are nonintegral rationals, whereas the sum of a root
of unity and its inverse is an algebraic integer. The proven Tate theorem
therefore gives geometric rank 20 for both reductions.

Artin–Tate and the square order of the finite Brauer group give reduction
Neron–Severi discriminant square classes −24 and −48, namely −6 and −3.
Indeed the normalized residual polynomials have values 24/7 and 48/13
at 1, and Artin–Tate multiplies by the corresponding prime. The negative
sign is the rank-20 intersection-form sign. The upper packet pins the
applicable Charles, Milne, and Liu–Lorenzini–Raynaud theorems, including
the latter's corrigendum.

If rho(Z_C)=20, specialization would embed its Neron–Severi lattice with
finite index in both rank-20 reduction lattices. Intersection preservation
would make both discriminants equal up to rational squares. Their ratio
2 is not a rational square. Hence rho(Z)<=19. With part A this proves
rho(Z)=19. This is a two-prime discriminant proof for the fixed surface,
not an assertion that the surface is generic in a family.

### C. Transfer back to the six-root and seven-root surfaces

The inherited even-sign quotient correspondence gives
`T(Z)_Q ~= T(S_o)_Q`; both surfaces have b2=22, so rho(S_o)=19 and
rank T(S_o)=3. This is a rational Hodge statement, not an integral-lattice
or Kummer-status identification.

The six outer deletions are isomorphic by the explicit hourglass symmetry.
The central deletion has Picard number 20. The established orthogonal
Hodge decomposition therefore gives

    rank T(Y_7)=2+6*3=20,
    rho(Y_7)=78−20=58.

The six outer surfaces form one complex isomorphism class by the explicit
D12 symmetry. Picard number is invariant under complex isomorphism, so none
of them is isomorphic to the central rank-20 deletion. Thus there are exactly
two classes among the seven deletion K3 surfaces. This strengthens the earlier
"at most two" conclusion without inferring nonisomorphism from different
nodal presentations.

Equivalently, the inherited exact identity is
`rho(Y_7)=6 rho(S_o)−56`. These inputs are in
`../deck-k3/HODGE_DECOMPOSITION.md`, Sections 2–4, and
`../hourglass-symmetry/RESULT.md`, Sections 3–4.

Finally, Shioda–Tate for the same fibration gives
`rank E(C(t))=19−18=1`. Since P is a nontorsion Q(t)-point, the subgroup
E(Q(t)) has rank at least one and at most that of E(C(t)), proving rank one
over Q(t) too.

## Exact bridge to Session57, with the center/root distinction retained

The smooth t=2 specialization is the congruent-number curve of area 6.
Its triangle map

    a=(u²−36)/Y, b=12u/Y, h=(u²+36)/Y

sends Q to (−3,−4,5), hence to the positive 3–4–5 triangle. Adding the
rational two-torsion point (0,0) gives `(12,36)`, precisely Session57's
chosen point on E_6; the double remains `(25/4,−35/8)`.

The central **square value** is c=h²=25, its **square root** is h=5, and
the fixed square gap is 4*6=24. Thus the associated square progression is

    (c−24,c,c+24)=(1,25,49)=(1²,5²,7²).

This explicit identification joins the K3 section to Session57's
`congruent-curve/PROOF.md`, Sections 2–3. It does not produce three
*different centers* in an arithmetic progression, a magic square of
squares, or a CM structure on the whole K3. A fiber with j=1728 is not a
proof that the surface is a CM or singular K3; here its Picard number is
nineteen, not twenty.

## Verification and boundary ledger

* **Lower half:** `check_lower_bound.py`, `LOWER_BOUND_CHECKS.json`, and
  `LOWER_BOUND_VERIFICATION_CERTIFICATE.json`; both Python modes pass.
* **Upper half:** `../picard-upper/check_counts.py` and its two matching
  result files; independent count agreement with the exploratory
  `probe_counts.py` in this folder.
* **Audit:** the [final independent rank audit](../audits/final/AUDIT.md)
  records PASS before this corollary was frozen. All good-reduction, component-splitting,
  reciprocity-sign, Tate, Artin–Tate, and discriminant-transfer checks passed.
* **No hidden sign guess:** exploratory primes with zero middle coefficient
  are explicitly left undetermined by two traces and do not enter the proof.
* **Not determined:** integral transcendental lattices, exact integral
  discriminants in characteristic zero, saturation/primitiveness of P,
  or arithmetic Picard rank over Q rather than geometric rank.
* **No new exploration:** this corollary only combines the closed lower
  proof, independently audited upper proof, and inherited Hodge identities.
