# Complementary regular norms determine every carried boundary jet

**Status: PROVED ANALYTICALLY + FINITE-EXACT + INDEPENDENTLY AUDITED.** A general leading-jet identity reduces the exact norm
multiplicity at x=-r to a regular complementary norm. The previously proved
regular certificates then settle seven boundary diagonals at all heights.
This is a polynomial-boundary theorem, not a new actual endpoint theorem or
a solution of general doubled-return noncancellation.

## 1. The incoming question and the recovered complement

The incoming `05-knowledge/results/third_20260906_trace.md` proves the exact
norm order h-1 at x=-h for all heights. Its overview explicitly asks for
complementary nonvanishing to settle other negative-integer boundaries.
The earlier all-height forced divisor gives only lower orders there.
Our `continuing4_20260906_regular_duality.md` identifies those complementary
rows after reflection and gives exact positive norm polynomials at heights
0 through6 (height0 is the empty determinant1). Combining the two retained
objects answers the incoming question on six additional unbounded diagonals.

The closest mechanism is the ramified small-root cluster. The canonical
hostile is old h=1,x=-1: the genuine carried quotient is -1/90720 at its
zero-root block, while naive raw specialization gives0. The corrected near
miss is identifying a stable trace perturbation with a stable determinant.
The underused sidecar is the norm of the complementary block, rather than
its dimension or individual first-row root geometry.

Source: the full carried multiplication operator. Target: its leading norm
jet at a negative parameter. Map: split the formal small block from the
coprime complementary factor, retain their two norms, and reflect only the
complement. Preserved: the exact coefficient, order and possible vanishing
of the determinant. Lost by either block alone: the other's zero or unit.
Needed sidecars: both first derivatives, the complement constant term,
all monomial powers and the factorial normalization. The finite tests below
compare the complete previously audited polynomial certificates, not just
observed signs.

## 2. Definitions and universal leading-jet identity

Fix h>=1 and1<=r<=h. Write H=h-r, z_0=2r+1 and delta=x+r. Use the full
carried rows with falling products

    p_(h,x)(t)=sum_(j=0)^h
      (2h+1)! (x+h)_(h-j) t^j /[(3h-3j)!(1+2j)!],
    q_(h,x)(t)=sum_(e=-1)^(2h)
      (2x+2h)_(2h-e) t^e /[(6h-3e)!(2+2e)!].

Let c_(h,h)(x) be the constant coefficient of det(wI-M), where M is
multiplication by the polynomial remainder of the full q modulo the monic
p. Polynomiality after cancellation is inherited from the proved carried
factor theorem. The inverse carry is part of M.

Let Phi_(H,z),Psi_(H,z) be the regular reflected rows from
`continuing4_20260906_regular_duality.md`, and let

    N_H(z)=det(-multiplication by Psi_(H,z) modulo Phi_(H,z)),
    N_0(z)=1.

Define the nonzero rational constants

    a0=(2h+1)! H! /[(3H)!(2r+1)!],
    A=(-1)^(r-1) (2h+1)! H! (r-1)!/(3h)!,
    B=2(2H)!(2r)!/(6h+3)!.

For every such h,r, with no restriction on H,

    c_(h,h)(-r+delta)
      = K_(h,r) delta^(r-1)+O(delta^r),

    K_(h,r)= B^r a0^(2r+1) N_H(2r+1)
                  /[A ((4h+2)!)^H].                  (1)

In particular the forced factor(x+r)^(r-1) has **exact** multiplicity r-1
if and only if N_(h-r)(2r+1) is nonzero. If the complementary norm vanishes,
the leading coefficient in(1) is zero and the order is higher; an identically
zero norm is interpreted as infinite order. The theorem does not assume
that this obstruction is absent at every complementary height.

## 3. Proof with both blocks retained

At x=-r, the first row is t^r Phi_(H,z_0)(t), whose complementary constant
is a0>0. The derivatives of its constant coefficient and of the inverse
carry are respectively A and B, directly from their consecutive falling
factors. Thus the two coprime factors t^r and Phi lift formally over
C[[delta]], giving a small block of rank r and a complementary block H.

Set delta=epsilon^r. The small roots have the form t_i=epsilon v_i(epsilon),
where the nonzero distinct leading values solve

    a0*v^r+A=0,       product_i v_i(0)=(-1)^r A/a0.

The inverse-carry coefficient is B delta+O(delta^2). Every other q
coefficient of index e<2r vanishes at delta=0; indices e>=2r are regular.
At a small root the inverse carry has epsilon order r-1. Nonnegative
indices e<2r have order at least r+e, and indices e>=2r have order at
least e>=2r. Hence the inverse carry uniquely supplies the leading term,
including r=1. Taking the characteristic norm of these r values gives

    det(-M_small)= (B^r a0/A) delta^(r-1)+O(delta^r).   (2)

Initially the error is a higher Puiseux order. The lifted small-block norm
belongs to C[[delta]], so its next possible exponent is the integer r.

At a nonzero complementary root rho the exact reflection identity is

    q_(h,-r)(rho)=rho^(2r) Psi_(H,z_0)(rho)/(4h+2)!.

The generic quotient and raw response agree on this complementary block;
this step does not specialize away the small inverse carry. Since Phi is
monic and product rho=(-1)^H a0, the complementary characteristic norm is

    det(-M_comp)|_(delta=0)
       = a0^(2r) N_H(z_0)/((4h+2)!)^H.               (3)

The product formula holds with multiplicities and needs no simplicity or
sign assumption on the complementary roots. Chinese remainders multiply
(2) and(3), proving(1). All its prefactors other than N_H are nonzero,
which proves the exact-multiplicity equivalence.

## 4. Seven uniformly exact boundary diagonals

The full regular characteristic certificate proves N_H(z)>0 for every
positive real z when1<=H<=6: its last coefficient is a product of positive
paired factors and a polynomial with strictly positive rational coefficients.
For H=0 the same statement is N_0=1. This particular positivity consequence
is purely algebraic and does not need the separately CITED real-root supplier
used by the actual Laurent first-return interpretation.

Therefore, for every h>=1 and r satisfying

    1<=r<=h,          0<=h-r<=6,

the norm has exact order r-1 at x=-r. Its leading coefficient has sign
(-1)^(r-1), since A has that sign and every other factor in(1) is positive.
For r=1 this means the norm value itself is nonzero, not a vanishing factor.
For H=0, a0=1 and(1) reduces to the inherited exact x=-h jet B^h/A.
The six diagonals H=1,...,6 are the additional conclusion. Each is unbounded
in h; this is not an extension of the finite actual endpoint table.

No primitive gcd hypothesis is required here. The reflected coefficient
identities and the norm certificates are valid without it. Coprimality enters
only the separate claim that a designated mass is an actual first support
return. Exact multiplicities for H>=7 still require complementary nonvanishing;
general response positivity and all-channel noncancellation remain OPEN.

## 5. Exact verification and stopping scope

The producer reads two pinned complete certificates: the old carried norm
polynomials in x-1 for h=6,...,10 and the regular residual polynomials for
H=1,...,6. It imports no producer implementation. For every1<=r<=h in those
heights with H=h-r<=6, it expands the complete old polynomial at x=-r,
checks all lower jet coefficients vanish, and compares the leading coefficient
with the full formula(1). All34 triples match exactly, including signs.
The symbolic all-height proof supplies the unbounded quantifiers.

The height-one carry control gives B=1/90720, so its small response is
-1/90720 and its characteristic norm is positive1/90720. The independent
referee reconstructs low-height norm polynomials directly from rational
Sylvester resultants, without either certificate as an oracle.

    python -B 04-computation/continuing5_20260906_complementary_norm_jets.py
    python -B -O 04-computation/continuing5_20260906_complementary_norm_jets.py

An optional --data argument selects the directory containing the two pinned
prior certificates. All276 always-active gates pass. Preserve those hashes
and the source/output pins when filing; no normalized output bytes are used
to claim raw-byte identity.

The proved connection is an exact reduction from one boundary problem to
the complementary response norm. A positive trace, an unmarked factor
degree, or a generic nonzero sample would not prove the needed unit. The
next unresolved input is N_H at the odd positive parameter with H>=7, or
another certificate that proves its nonvanishing there.

Independent [proof and exact referee](continuing5_20260906_complementary_norm_jets_audit.md) passes.

Filed checkpoint provenance: the [raw-byte manifest](continuing5_20260906_manifest.json)
pins the final report, source and output. Reviewed candidate-report hashes
above refer to the pre-promotion bytes. Source-location defaults and local
links were made portable where necessary; all emitted outputs were replayed
as raw bytes. The independent audit supplies the stated promotion basis.
