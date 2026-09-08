# Base-only Jacobians transport the quartic approximate root across the DG boundary

**Status: PROVED ANALYTICALLY + INDEPENDENTLY AUDITED; controls FINITE-EXACT.**
The canonical finite-pole mechanisms are inherited. The audited extension
checks their exact input and applies them at the previously unpaid boundary
of the specified DG surface. This does not close the quartic Keller stratum
or JC(2), and imposes no degree bound on the proposed second coordinate.

## 1. The actual inherited root and the proposed extension

Let P,Q belong to C[x,z], with deg_z P=4. The proposed algebraic statement is

    0!=J_(x,z)(P,Q) in C[x]
      implies that the canonical quadratic square approximate root
      of P is a polynomial in C[x,z].                         (1)

The Jacobian in (1) may vanish at finite base points and need not be a
constant. The conclusion is polynomiality, not monicity or a quadratic
member in the output pencil.

The inherited mechanism is
[THM-2158 / quartic-quadratic-deck-parity-and-exact-finite-pole-criterion](../../01-canon/theorems/THM-2158-quartic-quadratic-deck-parity-and-exact-finite-pole-criterion.md).
For

    P=V^2 z^4+beta z^3+gamma z^2+delta z+epsilon,

its uniquely specified approximate root with leading coefficient V is

    H=V z^2+[beta/(2V)]z+[4gamma V^2-beta^2]/(8V^3),
    deg_z(P-H^2)<=1.                                         (2)

Polynomiality is equivalent to the single congruence

    V^3 divides 4gamma V^2-beta^2.                            (3)

The explicit MISTAKE-248 correction is retained: (3) already forces V|beta,
so those are not two independent divisibility conditions. The later
[THM-2180 / quartic-first-pole-stratum-is-empty](../../01-canon/theorems/THM-2180-quartic-first-pole-stratum-is-empty.md)
and
[THM-2202 / uniform-all-degree-quartic-pole-closure](../../01-canon/theorems/THM-2202-uniform-all-degree-quartic-pole-closure.md)
close these finite poles for the reduced Keller branch. The all-degree
THM-2202 supersedes the earlier open-pole boundary wording in the older
notes. The remaining nonmonic terminal descent is not superseded.

The corrected near miss is to confuse ordinary polynomiality in x,z with
globality on a larger surface. The canonical hostile below separates them
exactly. The least-used sidecar is the first two Faber fluxes: their
constancy follows from fibre-degree zero of the Jacobian, even when the
base function on the right has zeros. The four live concepts are the
canonical root, the two fluxes, the second source chart, and the global
linear remainder. The map is an invertible affine change of the fibre
over C(r), followed by an actual finite-place pole test at r=0.
It preserves the root by uniqueness; it changes a nonzero constant
Jacobian to lambda r^2. The precise extension (1) pays that change.

## 2. Why the Faber pole inputs survive a base-only right side

First reduce the mate by constant polynomial target shears. If its current
positive degree is n and its leading coefficient is g_n(x), with A(x)
the leading quartic coefficient of P, the highest Jacobian coefficient is

    n A' g_n-4A g_n'=0.

Thus g_n^4/A^n is a nonzero constant. When 4|n, this implies
g_n=c A^(n/4) for c in C, so subtract c P^(n/4). This lowers the degree and
preserves the Jacobian. The process cannot terminate at a degree-zero
mate: then J(P,Q)=-P_z Q'(x), which cannot be nonzero and fibre-independent.
Hence its remaining degree n is positive and 4 does not divide n.
Taking multiplicities of the zeros of A in the leading-coefficient
identity now shows that A is a square. Choose V in C[x] with A=V^2.
There is no bound on n.

Work on the quadratic field or split deck U^2=V, put w=Uz, and depress
by a base-dependent translation Z=w+beta/(4U^3). The quartic becomes

    P0=Z^4+p Z^2+q Z+r.

The source Jacobian in these coordinates is kappa(x)/U, where
kappa(x)=J_(x,z)(P,Q). This is still independent of Z, even if it has
zeros or poles on the deck. For

    E_m=Pol_Z(P0^(m/4)),
    P0^(m/4)=E_m+c_(m,1)Z^(-1)+c_(m,2)Z^(-2)+...,
    Phi_m=4c_(m,1), Psi_m=4c_(m,2),
    Theta_m=4c_(m,3)+p c_(m,1),

[THM-2129 / quartic-faber-three-coefficient-boundary-classification](../../01-canon/theorems/THM-2129-quartic-faber-three-coefficient-boundary-classification.md)
gives the exact differential identity

    J_(x,Z)(P0,E_m)
       =(Z^2+p/4)Phi_m'+Z Psi_m'+Theta_m'.             (4)

The triangular subtraction also applies on this deck. A remaining leading
coefficient a_d(x) would contribute -4a_d' Z^(d+3). Since the right side
has fibre degree zero, a_d is constant. Subtract that constant times E_d
and continue, obtaining

    Q=H(P0)+sum_(4 not dividing m) c_m E_m,
    H in C[T], c_m in C.                              (5)

The constant field remains C in the algebraic function field of the deck.
In a split deck one works on a chosen component. Powers divisible by four
are precisely the global target shears already allowed.

Equation (4) and fibre-degree zero now force

    Phi_Q=sum c_m Phi_m in C,
    Psi_Q=sum c_m Psi_m in C.                         (6)

Only the final primitive equation changes:

    Theta_Q'=kappa(x)/U.                              (7)

No pole argument used below requires (7), or requires kappa to be a unit
at the place. This is the essential input audit; one cannot simply cite a
constant-Jacobian theorem after changing its hypothesis.

## 3. The complete finite-pole argument with these exact inputs

At every finite place dividing V, the original polynomial boundary

    Q(x,0)-H(P(x,0))

is regular. The constant Faber expansion (5), this regularity, and the two
constant fluxes (6) are all the inputs required by the pole mechanisms.
Here is their coverage, including the lowest degree and the odd cases.

1. **A pole of beta/V.** THM-2180, Sections 1--3, normalizes the depressed
   face to (T-1)^3(T+3). The top boundary term is uniquely deepest and has
   coefficient 4^n binom(n/4,n), which is nonzero whenever 4 does not divide
   n. Its contradiction uses the regular boundary, not the value of the
   Jacobian. Thus beta=V b with b polynomial.

2. **A remaining constant-coefficient pole.** Write D0=4gamma-b^2 and
   suppose at some place that nu(D0)<nu(V). Section 2 of
   [THM-2189 / nonsplit-quartic-deck-forces-the-remaining-pole-congruence](../../01-canon/theorems/THM-2189-nonsplit-quartic-deck-forces-the-remaining-pole-congruence.md)
   normalizes the negative scale to a nonconstant polynomial
   T0(X)=1+A X+B X^2+C X^3. The uniquely deepest boundary and successive
   fluxes require three consecutive coefficients n,n+1,n+2 of
   T0^(n/4) to vanish. If A!=0, THM-2129's recurrence makes this impossible
   for odd n and forces an exact quadratic square for n twice odd.
   If A=0 and C!=0, the same third-order recurrence descends to a zero
   constant coefficient, again impossible. If A=C=0,B!=0, the coefficient
   at the first even index in {n,n+1} is nonzero: it is a power of B times
   binom(n/4,ceil(n/2)), whose upper argument is not an integer. Thus a
   remaining bad place necessarily lies in the same exact square cage,
   and n=4r-2. This step uses regularity and (6); exact zero of one flux
   on a nonsplit deck is stronger than needed.

3. **Reduced degree n=2.** THM-2189, Section 5, applies before imposing
   deck parity. In its notation a,c,Lambda,Omega, a hypothetical bad place
   has nu(a)=-A<0, nu(c)=C with 0<C<2A, nu(Omega)=2C. For the unique
   possible lower odd seed coefficient k, the required quantities are

       (a^2 c+k a)/2,
       2a^3 Lambda+k a^2(c-1/2),
       a^4(2Omega-Lambda)+k a^3 Lambda.

   Regularity of the first forces C=A and k!=0; constancy of the second
   forces nu(Lambda)=A. In the third, -a^4 Lambda has order -3A while
   both possible competitors have order -2A. This is an uncancellable
   pole. The final primitive flux and the right-side value are absent.

4. **Every higher reduced twice-odd degree.** THM-2202, Sections 2--4,
   use only the same square-cage valuation data, constant coefficients
   in (5), the regular boundary, and (6). In the Lambda chamber the
   first flux forces an odd matching seed whose second flux is strictly
   shallower, leaving a unique negative second-flux term. In the other
   chamber the boundary forces an odd matching seed whose first flux
   is uniquely deepest. Its equal-order boundaries, possible vanishing
   Lambda, all lower even seeds, and all odd seeds are explicitly covered.
   No line uses Theta_Q' or the nonvanishing of kappa at the place.

This is an interface extension of the proved all-degree mechanisms, not
an inference from a bounded bank. The four cases exhaust every possible
finite pole. Consequently V divides D0, so (2) belongs to C[x,z]. This
proves (1); the [independent audit](planar_jc48_sep08_quartic_boundary_audit.md)
checks this explicit dependency transfer. Allowing a zero of kappa at the tested place
does not create a missing degree-two or higher-degree case.

A useful boundary control is

    P=x^2 z^4+z^2,  Q=xz^2,
    J(P,Q)=-2z^3,
    H=xz^2+1/(2x).

Here the Jacobian depends on the fibre, and the approximate root does
have a finite pole. Dropping fibre-degree zero from (1) is false.
Conversely the nonconstant-base examples in the verifier have kappa
vanishing where their leading coefficient vanishes and still have
polynomial canonical roots. They prevent a hidden unit assumption.

## 4. The root actually glues on the DG surface

Use the actual surface and two full charts of
[the DG filtration](planar_jc48_sep08_dg_quadratic.md):

    W=(P1_x x P1_z) minus {z=x^2},
    U0=A2_(x,t), Uinf=A2_(r,b),
    x=1/r, t=-r^2-r^4b, D={r=0},
    omega=dx wedge dt=r^2 dr wedge db.

Suppose F,G are global functions with dF wedge dG=lambda omega, lambda!=0,
and deg_t F=4. There is no degree assumption on G. On U0, the argument
above gives a polynomial canonical root H in C[x,t], with leading
coefficient V and

    F=H^2+L,       deg_t H=2, deg_t L<=1.              (8)

On the second chart, Fbar and Gbar are polynomials in r,b and satisfy

    J_(r,b)(Fbar,Gbar)=lambda r^2.                    (9)

Their degrees in b are their original degrees in t: the change is affine
with invertible slope -r^4 over C(r). The leading b^4 coefficient is

    r^16 V(1/r)^2=Vbar(r)^2,
    Vbar(r)=r^8 V(1/r).

Polynomiality of that coefficient forces deg V<=8, so Vbar is an actual
polynomial. Applying (1) to (9) produces a polynomial approximate root
in C[r,b], with the chosen leading coefficient Vbar.

The rational quadratic

    Htilde(r,b)=H(1/r,-r^2-r^4b)

has exactly that leading coefficient, and Fbar-Htilde^2 has degree at most
one in b. The canonical root is unique by matching degrees four, three,
and two. Thus Htilde is the newly obtained polynomial root. It follows
that H is regular on both full charts and hence is a global function on W.
Likewise L=F-H^2 is global. In the complete inherited filtration,

    H in L_2,      L in L_1.                           (10)

The proof checks regularity at the whole boundary chart, not only on its
function field or on a punctured divisor. All target shears used during
reduction are polynomials in the global F and preserve globality.

Also L is nonconstant. If it were constant then J(F,G)=2H J(H,G) could
not be a nonzero constant on the affine source plane. The global square
prefix (8)--(10) is the precise paid consequence. It does not say that H
or L is a linear output combination of F,G, that either has a constant
Jacobian mate, or that a root itself can be substituted as a target
coordinate. Thus the previously recovered exclusion of pairs with a
member in L_3 cannot be applied to H or L without another argument.

## 5. The proposed hostile is valid and outside the Keller hypotheses

The cheap example is

    F=t^4+2x^4t^2=t^4+2u^2, u=x^2t,
    H=t^2+x^4,          F-H^2=-x^8.

The first function is global, but in the full second chart

    Fbar=r^8(1+r^2b)^4+2(1+r^2b)^2,
    Hbar=r^4(1+r^2b)^2+r^(-4).

Hence its canonical root has genuine boundary pole order four. Its
boundary finite-pole congruence has numerator

    4gamma_bar Vbar^2-beta_bar^2=8r^20(1+r^8),
    Vbar=r^8,

which is not divisible by Vbar^3=r^24. The failure is exactly the one
that the added base-Jacobian interface excludes, not a deck issue.

This F cannot have any polynomial Keller mate already because both
F_x and F_t vanish on the whole line t=0. It is therefore an exact hostile
to root globality from F global alone, with no contradiction to (10).
The inherited monic quartic theorem also excludes a Keller mate, but is
unnecessary for this elementary control. The example does not refute
root globality under the full F,G hypotheses.

## 6. Exact controls and the remaining consumer

The [standalone verifier](../../04-computation/planar_jc48_sep08_quartic_boundary.py)
checks the canonical root and arbitrary affine fibre transport symbolically.
It reconstructs the Faber coefficient recurrence in degrees 1..12 with
independent p,q,r and p',q',r', and verifies the full identity (4),
including both flux derivatives. The unbounded property is inherited
from the exact identity and proofs above; these twelve degrees impose
no restriction on the mate.

Six polynomial pairs have nonconstant base-only Jacobian, with zeros of
that right side and of the leading coefficient allowed at the same place.
They include reduced degrees one and two and arbitrary target-shear
patterns increasing the raw mate degree. Other controls retain the
fibre-dependent-Jacobian counterexample, the literal global quartic with
nonglobal root, its actual critical line, the precise r^2 chart multiplier,
and global square prefixes with no mate asserted.

    python3 -B 04-computation/planar_jc48_sep08_quartic_boundary.py
    python3 -B -O 04-computation/planar_jc48_sep08_quartic_boundary.py

The remaining source consumer must use the simultaneous global pair
H in L_2 and L in L_1 in (8), together with the unrestricted mate G.
No full quartic Keller theorem or global pair exclusion follows merely
from membership of those two auxiliaries in already-excluded component
layers. The independent proof and source audit passes;
normal and optimized replays agree exactly on all 126 always-active gates.
The matching [planar_jc48_sep08_quartic_boundary.out](planar_jc48_sep08_quartic_boundary.out) contains 460 bytes.

- Source SHA256: `de2e1125bd6dc2446ac38b46b19e2f3e59ac71410b7965ff4b1716b49a86a439`.
- Output SHA256: `e3c2a76c515933f92a13b4305a05561910491cf6c1c5a887a1c21c6e7ac8a939`.
- Semantic SHA256: `6e989a1819d645fd97d7105a0ba44ee4c14935276ea15e102abda80b959fa80e`.

Source and output are frozen for the independent audit.
