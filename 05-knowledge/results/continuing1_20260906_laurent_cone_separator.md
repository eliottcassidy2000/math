# An actual carried Laurent response escapes every positive phase-weighted derivative window

**Status: PROVED relative to the stated inherited inputs; independently
audited exact certificates.** The [independent referee](continuing1_20260906_laurent_cone_separator_audit.md)
accepts the complete proof and reconstructs the literal coefficient and
spectral certificates with 646 always-active gates. Root independently read
and accepted the proof. The new statement is a limitation of a proposed sign mechanism on
an actual primitive trinomial support. It is not a counterexample to Laurent
noncancellation. General two-rung separation remains OPEN. The actual first
and doubled rows used below are both inherited objects, and their negative
doubled-row values at the first roots were already covered by earlier exact
work. No new endpoint-family closure or literature priority is claimed.

For the genuine support **(-21,1,12)**, the actual carried doubled response
cannot be expressed, even modulo its true first-row equation, as a
nonnegative combination of the coupled derivative windows. More strongly,
this remains impossible after multiplying each window by an arbitrary
nonnegative-coefficient Laurent polynomial in the positive phase s=-t.
All integer powers, including inverse powers, are excluded by one exact
linear separating functional. This is proved analytically; a bounded power
bank is only an independent arithmetic check.

**Concurrent recovery.** The incoming [endpoint33 proof](second_20260906_laurent.md)
independently contains the same endpoint15 constant-cone obstruction.
That quadratic obstruction here is a concurrent independent control. The
actual cubic all-integer Laurent-weight obstruction below is the stronger
new mechanism limitation, and the [positive amplitude repair](continuing1_20260906_laurent_real_power_boundary.md)
records its exact branch and coefficient-field boundary.

## 1. Inheritance and source-to-target contract

The nearest new input is the audited
[minimal coupled-window cone](creative_20260906_laurent_bridge.md).
It proves that, at a vanished interior coefficient of a real-rooted carrier,
all two-sided derivative windows reduce to min(k,n-k) nonnegative generating
rays. Its strict sign input is **THM-4440**,
`01-canon/theorems/THM-4440-signed-duplication-sos-and-real-rooted-laurent-return.md`.

The actual coefficient map comes from the proved
[A2B3 midpoint identities](overnight7_20260906_laurent_midpoint_transport.md)
and [alpha completion](overnight8_20260906_alpha_completion.md).
The original first-row real-root theorem is **THM-4436**,
`01-canon/theorems/THM-4436-complete-factorial-row-simple-negative-roots-and-trinomial-phase-wall.md`.
These are proved repository inputs. No new external theorem is invoked.

The canonical hostile and corrected near miss are the
[combined Euler/pencil model](combined_pencil_empty_core_morning_sep06.md),
which retains the carrier but permits a free relative scale, and the earlier
independent contiguous replacements, which lose its selected zero. Here the
entire genuine binomial coefficient law and the crossing factor 2t remain
fixed. The least-used sidecar is the action of the original phase variable
on the finite first-row quotient, rather than just the carrier's roots.

The live board is the actual carried midpoint row, the common zero
coefficient, the complete derivative-window cone, phase multiplication,
the quotient's three real embeddings, and mixed contiguous products.
The source is a literal trinomial constant-term row. The target is a cone
membership statement in its first-row quotient. The map below preserves
the actual coefficient normalization, both available carry endpoints and
every first root. Quotient reduction discards multiples of the actual first
row only; it does not discard its different positive roots or merge their
values. The sidecar needed for the dual certificate is the ordered three-root
embedding and its Lagrange weights. The cheapest decisive test is a small
quadratic separator; the cubic witness then defeats all positive Laurent
phase weights, not merely another finite derivative list.

## 2. The actual carrier and the original scale

Use the established A=2,B=3 progression with x=1, h>=1:

    q=h+1, k=2h+1, g=3h+2,
    support=(-6h-3,1,3h+3).

All three support exponents are congruent to one modulo g, so a support
return mass must be a multiple of g. The j=0 counts (1,3h,1) give an actual
return at mass g. Thus g is the first support-return mass. The complete raw
first and doubled rows are

    P_raw(t)=sum_(j=0)^h g! t^j /[(1+j)!(3h-3j)!(1+2j)!],
    Q_raw(t)=sum_(j=-1)^(2h) (2g)! t^j /[(2+j)!(6h-3j)!(2+2j)!].       (1)

The genuine lower carry j=-1 is retained. For example, setting the first
two Laurent coefficients to one and the third to a formal gamma gives
CT(f^g)=gamma P_raw(gamma^2) and
CT(f^(2g))=gamma^2 Q_raw(gamma^2). No coefficient has been independently
rescaled. The verifier reconstructs every return through mass 2g by literal
charge enumeration, independently of (1).

For s>0 set t=-s and define

    B(v)=sum_(i=0)^q (-1)^(q-i) binom(q+2i,q-i) v^i,
    C(v)=sum_(i=0)^(q-1) (-1)^(q-1-i) binom(q+1+2i,q-1-i) v^i,
    D(v)=sum_(i=0)^(q-1) (-1)^(q-1-i) binom(q+2i,q-1-i) v^i,
    K_B=s^q B(u^2/s), K_C=s^(q-1) C(u^2/s), K_D=s^(q-1) D(u^2/s),
    H=(1+u)^g K_B, H_C=(1+u)^g K_C, H_D=(1+u)^g K_D.

These are the original complete path kernels. In particular,

    [u^k]H=(-s)P_raw(-s),
    s Q_raw(-s)
      =s^-1 [u^(2k)]H^2 -2 [u^(2k-2)]H_C H_D.                       (2)

Equivalently, before dividing by the positive common power, the actual
response is s^-2{[u^(2k)]H^2-2s[u^(2k-2)]H_C H_D}. The factor s on the
mixed term is indispensable. This is the factor whose relative normalization
the combined-pencil model was allowed to change; it is fixed here.

The inherited unsigned full beta path polynomial has negative roots.
The reversal and alternating signs in the displayed B make its roots
positive. Thus its u^2/s pullback K_B, and consequently H, have only real
roots for s>0. Their endpoints are nonzero. The degree is
n=g+2q, with n-k=3q>k. At a root of P_raw(-s), the selected H coefficient
vanishes at its genuine interior index k.

Define the polynomial responses

    J_r(s)=s^-1 [u^(2k-2r)](D_u^r H)^2,  0<=r<k,
    T(s)=s Q_raw(-s).                                                (3)

Every displayed quotient by s is exact as a polynomial. The inherited
coupled-window theorem gives J_r(s)<0 at every first root. Its two-sided
cone reduction says that every allowed right/left window, modulo the
selected-zero equation, is a nonnegative combination of these k rays.

## 3. The first small obstruction and its positive repair

For h=1, support (-9,1,6), a constant coefficient suffices:

    T == (6305/1714)J_0  modulo P_raw(-s).

For h=2, support (-15,1,9), the monic first equation is

    p_2(s)=s^2-10s+1.

On a remainder a+bs define ell_2(a+bs)=-1376a-139b. Direct exact reduction
gives, for r=0,...,4,

    ell_2(J_r)=(872320,19900736,277245248,1956810240,4509872640),
    ell_2(T)=-126496.                                                (4)

Thus no nonnegative constant combination of all coupled windows represents
the actual T, even after imposing its true first equation. This is stronger
than finding negative coordinates in one chosen, unreduced polynomial basis.
It excludes every representation in the entire quotient cone.

There is nevertheless an exact positive phase repair at this support:

    T == (277182676568+75869733936 s)/191390856673 * J_0 mod p_2.       (5)

Both coefficients are positive. This reproduces the known negative actual
response through a richer observer certificate; it is not a new closure of
the already treated endpoint-15 family. It motivates testing arbitrary
positive phase weights rather than stopping at (4).

## 4. The actual cubic quotient and the infinite obstruction

Now fix h=3, so the primitive support is (-21,1,12), g=11, q=4, k=7,
n=19. Its actual kernels are

    B=v^4-10v^3+28v^2-20v+1,
    C=v^3-9v^2+21v-10,
    D=v^3-8v^2+15v-4.

The first row is

    P_raw(t)=110+4620t+9240t^2+330t^3,
    P_raw(-s)=-330p(s), p(s)=s^3-28s^2+14s-1/3.                       (6)

Work in the three-dimensional real algebra A=R[s,s^-1]/(p). The inverse
variable exists since p(0) is nonzero. Let C_pm consist of all finite sums

    sum_(r=0)^6 sum_(d in Z) c_(r,d) [s^d J_r],  c_(r,d)>=0.          (7)

Equivalently each ray may receive its own arbitrary nonnegative-coefficient
Laurent polynomial in the positive phase s. This includes every such
weighting of all 7*12 two-sided windows. It does not mean arbitrary
functions merely positive at the three first roots; that is a larger class.

**Theorem.** Neither [T] nor [T-J_0], the actual beta-skip payment, belongs
to C_pm or its closure.

Define ell(a+bs+cs^2)=y_0 a+y_1 b+y_2 c, where

    y=(-32559418467680575845,
       -813366686033280291,
       -19065454274144095).                                      (8)

Exact reduction gives

    ell(T)=-8569669792006914054574 <0,
    ell(J_0)=0, ell(s J_6)=0.                                    (9)

It remains to prove ell(s^d J_r)>=0 for every integer d, not just for the
finitely tested powers. The following argument supplies that step.

Let 0<alpha<beta<gamma be the roots of p. Exact disjoint rational intervals,
each with opposite endpoint signs of p, are:

| Root | Left endpoint | Right endpoint |
|---|---:|---:|
| alpha | 63077879/2516582400 | 504623033/20132659200 |
| beta | 380435/786432 | 31703/65536 |
| gamma | 241815037630759/8796093022208 | 30226879703845/1099511627776 |

These three intervals exhaust the degree and give three distinct simple
positive roots; in particular 1/3<beta<1. Rational interval Horner evaluation
on each interval gives J_r<0 for all seven r, T<0, and N<0, where

    N(z)=-32559418467680575845 z^2
           +910850350409022843369 z-433076656792870357777.           (10)

The signs J_r<0 also follow from the inherited real-rooted window theorem.
The direct enclosures separately verify them at this actual coefficient law.
The signs T<0 explicitly retain the distinction between cone failure and
actual noncancellation: the actual doubled row is negative at all three
first roots.

Lagrange interpolation gives ell(F)=sum_i lambda_i F(r_i), with

    lambda_i=N(r_i)/p'(r_i).

Indeed the Lagrange numerator is
s^2-(28-r_i)s+(14-28r_i+r_i^2), and applying (8) gives (10). Since p'
has signs (+,-,+) at its ordered roots, lambda has signs (-,+,-).

Fix r and put e_d=ell(s^d J_r), f_d=e_d/beta^d. Therefore

    f_d=A_r(alpha/beta)^d-B_r+C_r(gamma/beta)^d,
    A_r,B_r,C_r>0.                                               (11)

This identity holds for every integer d, including negative integers. Its
second discrete difference is strictly positive. The exact four-level
certificate in the JSON and transcript verifies, for every r=0,...,6,

    e_0>=0, e_1>=0, e_(-1)>3e_0, e_2>e_1.                        (12)

Together with 1/3<beta<1 this gives f_(-1)>f_0 and f_2>f_1. Strict
discrete convexity means f decreases up to the possible central minimum
and increases after it: for d<=-1 its value is at least f_(-1)>f_0>=0,
and for d>=2 its value is at least f_2>f_1>=0. The two central values
are themselves nonnegative. Hence e_d>=0 for every d in Z.

This proves separation from (7). Continuity of ell in the finite-dimensional
quotient also excludes its closure. Equation (9) gives the same negative
separator value for the actual skip payment T-J_0. No countable power bank
or numerical root plot substitutes for (11)--(12).

## 5. What can cross the wall, and what is still missing

The obstruction rules out adding further derivative orders, right windows,
positive powers, inverse powers, or nonnegative Laurent combinations of them
on this same H. It does not rule out a mixed-carrier inequality or a
coefficient-sign-changing multiplier certified positive on the actual root
set. For example the exact representative

    R(s)=142115562175391338833022911/115962939903341750549938130
      +(137850584919079100136401223/115962939903341750549938130)s
      -(825111792094668079242879/23192587980668350109987626)s^2

satisfies T==R J_0 modulo p. The same three rational intervals certify R>0
at every first root. Its negative quadratic coefficient is material: no
nonnegative Laurent-polynomial combination of the full original window
family can replace it. This is an exact repair for the named already-known
row, not an all-parameter positivity theorem for R or for interpolated ratios.

The recovered zero-preserving lowering in
[the Euler transport note, Section 4](laurent_transport_empty_core_next_sep06.md)
is also retained literally:

    A=(1+u)^(g-1)K_B,
    B_lower=(2/3)u(1+u)^(g-1)K_C,
    L=A+B_lower,
    k[u^k]H=g[u^(k-1)]L.

Normalize its degree 2k-2 quadratic coefficients by s^-1 as before. The
separator values of A^2, 2A B_lower, B_lower^2, and L^2 respectively are

    -2262024758344833167156,
    -10441745105422569984592/9,
     102255208380723434281160/27,
     9855304589145228814172/27.                              (13)

Thus isolated mixed pieces can cross the dual wall, while the complete
lowered square does not do so at phase degree zero. The mixed pieces are
not automatically sign certificates. In particular [u^(k-1)]A is nonzero
in the true quotient; removing B_lower loses the selected-zero predicate.
The cheap hostile A=(1-2u)^2, H=(1+u)A has H_2=0 but [u^2]A^2=24>0.
Real-rootedness alone cannot fix that deletion. The previously stated
interlacing hypothesis for real-rooted L is not silently promoted here.

The precise next obligation is an actual composition-recurrence or
divided-difference inequality controlling a mixed term such as (13), or a
uniform proof of rootwise positivity for an exact signed multiplier. The
separator determines which proposed extra generator must change sign under
ell before it can repair the current cone. This retains a useful positive
test instead of requesting more degree or phase-power enumeration.

## 6. Reproduction, finite universe, and frozen files

The standalone producer imports no repository mathematics, SymPy, numerical
root finder, or LP solver. Its three actual controls are exactly x=1,h=1,2,3.
It reconstructs literal charge fibres through mass 2g, all original kernel
coefficients, the first-zero identity, complete derivative squares, the
factor-s midpoint formula, and all allowed two-sided cone coordinates.
For the cubic separator it verifies the rational root enclosures and all
four-level inequalities. A separate forward/backward recurrence

    e_(d+3)=28e_(d+2)-14e_(d+1)+e_d/3

is compared with direct quotient multiplication for -16<=d<=16, on every
ray. This is an independent arithmetic path, not the all-integer proof.

[Source](../../04-computation/continuing1_20260906_laurent_cone_separator.py),
[output](continuing1_20260906_laurent_cone_separator.out),
and [full rational certificate](../../04-computation/continuing1_20260906_laurent_cone_separator_certificate.json):

```
python -B 04-computation/continuing1_20260906_laurent_cone_separator.py
python -B -O 04-computation/continuing1_20260906_laurent_cone_separator.py
python -B 04-computation/continuing1_20260906_laurent_cone_separator.py --certificate 04-computation/continuing1_20260906_laurent_cone_separator_certificate.json
```

Normal and optimized output is byte-identical LF, with 1,799 always-active
gates. Frozen SHA256 values are:

    source 5bd98bc1fe072fcfc88b34af46b2310e1a4e184dcefdff6c5a508a820d82d00b
    output 207b86469f1b5e8b5bb18a5cd1cd2a931c4c0c596a5fe66ba354923fe3069fee
    JSON   4b1ee5770b484e4164e692fbf2934f4099800b0b85d379ce01d8afef71040cc0

Exploratory files with `probe` in their names are not proof inputs. One
exploratory LP response failed an exact equality check and was rejected;
the final separation uses no LP verdict. The producer and proof remain
outside the repository pending parent integration.
