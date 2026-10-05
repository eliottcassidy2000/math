# Prime-ary addresses, a signed eta reader, and the missing height coordinate

2026-10-05. **PROVED:** the controller/shadow bijection, its two cost defects,
the signed enumeration, the source-weighted generating function, the
explicit golden-plane reduction and its lifted matrix orders, the
fixed-parent height completion theorem, and the elementary tree moment
identity below. **CITED:** the level11 eigenform and its Hecke recurrence.
Their combination identifies the signed enumeration with its coefficients
at powers of2. **FINITE-EXACT:** the accompanying controls. **OPEN:** a
fixed-seed Fourier bound, universal entry into grounded families, and Collatz.
These are scoped deductions; no priority claim is made.

Artifacts: [program](../../04-computation/experiments/prime_ary_controller_eta_20261005.py)
and [exact JSON output](prime_ary_controller_eta_20261005.json).

## 1. Inheritance, correction, and the working board

The closest proved mechanisms are the native five-letter language and
translation decoder in [carry interfaces](collatz_carry_interfaces_20261004.md),
sections2–3, and the branch isometry of the
[rooted refuel tree](finite_seed_refuel_tree_20261005.md), section5.
The canonical hostiles are the empty native interface LB and a residue
representative which is not the supplied source. The corrected near miss is
the Fourier package's nonprimitive-frequency tail artifact. Its current
[sections4l–4m](collatz_fourier_experiments_20261004.md) take precedence over
the original heavy-tail reading. The underused sidecar is the ordinary height
bound on a branch integer, alongside its successive ternary digits.

The anchor is source-specific root certification; the niche is signed
controller enumeration and the eta11 local factor; the wildcard is completion
of a prime-ary address tree. Our six concepts are **address / binary cost /
native guard / signed sum / ordinary height / grounded terminal**. The
META-PATTERNS cards used are “Controlled forgetting and unlabeled quotients
require a sidecar” and “Canonicalize mechanics; retain the predicate.”

The incoming packages through `2ffc27a9a2` were inherited, including the
rooted families, centroid stopping theorem, level11 torsion bridge and
two-clock Fourier energy. Fresh fetches during this session initially found
no additional commits. In particular, the uniform-start experiments give
finite evidence for different typical q=3 and q=5 rates, not a proved
asymptotic Lyapunov exponent or a fixed-integer estimate. No claim that
arithmetic deep digits cause the observed excess is inherited.

## 2. A five-letter native controller is a four-letter free shadow

Read words chronologically. The inherited library is:

| letter | action | native source | grade j | binary cost A | credit |
|---|---|---|---:|---:|---:|
| H | (729x+669)/1024 | 155 mod2048 | 3 | 10 | 2 |
| G | (9x+5)/8 | 11 mod16 | 1 | 3 | -1 |
| A | (81x+85)/128 | 187 mod256 | 2 | 7 | 1 |
| B | (81x+73)/128 | 7 mod256 | 2 | 7 | 1 |
| L | (9x-3)/16 | 219 mod256 | 1 | 4 | 4 |

Here j is half the ternary valuation of the slope numerator. H and L are
proved common-future dependencies, while G,A,B are actual odd-step blocks.
They must retain those types. Write N for the finite words with no LB,
including the empty word. The inherited native-language theorem says
exactly that these words have some legal positive source.

**Shadow theorem (PROVED).** Replacing every H by LB gives a bijection

    E : N -> {G,L,A,B}*,     E(H)=LB.

Its inverse C replaces every LB by H. Occurrences of LB cannot overlap.
An expanded H starts with L and ends with B, so it creates no extra LB
across its boundary. Every LB in E(w) therefore came from one H. Conversely
compressing all occurrences leaves no LB, and expanding recovers the word.
This proves both inverses, including arbitrary lengths and the empty word.

The map preserves j because j(H)=j(L)+j(B)=3. If h(w) counts H, it also gives

    A(w)=A(E(w))-h(w),
    credit(w)=credit(E(w))-3h(w).                  (1)

Thus it preserves the complete word as combinatorial information, but
does not preserve the binary clock or the funding semantics. In particular
an H is not an arithmetic identity LB. In fact

    F_LB(x)=(729x+925)/2048,
    H(x)=2 F_LB(x)-1/4.                           (2)

The native LB domain is empty; H has infinitely many native sources.

Transporting free concatenation defines an associative product on N:

    u star v = C(E(u)E(v)).

It is ordinary concatenation except at a boundary ending L / starting B,
where that one pair fuses to H. Associativity follows from the bijection,
not from an unproved interchange law. Equation(2) and the lost guard show
why this product is not composition of the arithmetic maps. Its single
boundary correction is a useful finite normal form, not a Collatz reroute.

### The two clocks and the exceptional fusion

There are two monomer letters G,L and two dimer letters A,B in the shadow.
Consequently the number c_j of native words of grade j has generating series

    sum c_j z^j = 1/(1-2z-2z^2),
    c_0=1, c_1=2, c_j=2c_(j-1)+2c_(j-2),
    c_j=sum_(r=0..floor(j/2)) binom(j-r,r) 2^(j-r). (3)

The first counts are 1,2,6,16,44,120,328,896. This is a literal bijection
with two-colour monomer/dimer tilings, not just agreement of recurrences.

The inherited weighted guard matrix is

    M = [[h+g+a+b,l],[h+g+a,l]],
    e_ordinary^T (I-M)^(-1) (1,1)^T
        =1/(1-h-g-a-b-l+bl).                      (4)

The shadow preserves a multiplicative weight exactly when h=bl. For
ternary grade alone h=z^3 and bl=z^3, giving (3). Keeping binary cost t too
gives instead

    1 / (1-z*t^3-z*t^4-2*z^2*t^7-z^3*t^10+z^3*t^11). (5)

The two cubic terms differ by precisely the fusion's one binary digit.
Adding a credit variable changes their exponents by three as in (1).
The grade cancellation therefore detects exactly which observer has
forgotten the arithmetic cost. It cannot be used after silently restoring
that cost. The grade-three witness is already H versus the shadow LB.

## 3. The eta11 coefficients at powers of2 are signed native counts

Let

    f(q)=q product_(n>=1)(1-q^n)^2(1-q^(11n))^2
        =sum_(n>=1) a_n q^n.

Its initial terms are

    q-2q^2-q^3+2q^4+q^5+2q^6-2q^7+0q^8-2q^9-2q^10+q^11+...

**CITED inputs.** The level11 normalized eigenform and curve identification
are recorded in [Elkies, conductor11](https://people.math.harvard.edu/~elkies/nature.html).
For a normalized weight2 eigenform with trivial character, at a prime p
not dividing the level the Hecke relation is

    a_(p^j)=a_p a_(p^(j-1))-p a_(p^(j-2)),  j>=2.

See [Ribet–Stein, Lemma16.1.7 and (16.1.4)](https://wstein.org/books/ribet-stein/main.pdf),
printed page152. Here p=2 is good and a2=-2, so

    sum_(j>=0) a_(2^j) z^j = 1/(1+2z+2z^2).       (6)

**Signed reader theorem (PROVED from those inputs).** Assign sign +1 to H
and -1 to each of G,A,B,L. For every j>=0,

    a_(2^j) = sum_(w in N, grade(w)=j) (-1)^(|w|-h(w)). (7)

Proof: under E the sign is (-1)^(shadow length), since H expands to two
negative letters. Each monomer and dimer has signed total -2, giving (6).
Equivalently substitute h=z^3, g=l=-z, a=b=-z^2 in (4); the two cubic
terms cancel. The empty word has sign1, matching a1, and grade1 has two
negative words, matching a2. This proves all grades, not only the checked
coefficients through q^1024.

The values are 1,-2,2,0,-4,8,-8,0,16,-32,32,0,... . The sign is an
enumeration observable; it is not primality, route credit, or a Fourier
phase. The whole eta product is not identified with a controller series:
(7) extracts only its coefficients at indices 2^j. Other eigenforms with
the same local2 factor have this reader too, so it does not characterize
level11 by itself.

For every good prime p the same Hecke identity has the explicit signed
tiling expansion

    a_(p^j)=sum_(r=0..floor(j/2)) binom(j-r,r) a_p^(j-2r) (-p)^r. (7a)

This can be viewed as a recursive tree with |a_p| monomer choices of the
sign of a_p and p negative dimer choices; the grade decreases by1 or2 at
every call. It is a finite signed counting proof, not a tree of positive
Collatz certificates. At2 it is exactly our two-plus-two shadow. At3,
a3=-1 gives one negative monomer and three negative dimers, hence the local
denominator1+z+3z^2; that is the polynomial already connected to the
[centroid return's lattice cube root](centroid_membership_20261005.md).
At the bad prime11 the Hecke factor is instead1/(1-a11*z)=1/(1-z): the
second-order term disappears. The three prime roles are therefore different
and specified by maps. The controller bijection is specifically the local2
case of this prime-indexed recurrence, not a claimed uniform identification
of all prime trees with one Collatz language.

### This closes an explicit loop to the golden clock

The signed recurrence uses

    K=[[-2,-2],[1,0]],       K^4=-4I.

Modulo3, K=J M J with

    J=[[0,1],[1,0]],   M=[[0,1],[1,1]];

the polynomial is X^2-X-1 and K has order8. This is exactly the matrix
already realized by Frobenius-at2 on the level11 curve's3-torsion in
[the level11 package, section4](level11_dessin_golden_torsion_20261005.md).
Thus (7) adds a specific signed controller reader to that earlier bridge.
The field F9 is the torsion operator field; it is not Z/9Z or the field
of coordinates of all torsion points.

Three recursions on related carriers now have distinct exact meanings:

| observable | characteristic polynomial | retained meaning |
|---|---|---|
| length counts on G/B/L | X^2-3X+1 | native golden sublanguage, eigenvalues phi^(+/-2) |
| unsigned ternary-grade counts | X^2-2X-2 | two-colour tilings, eigenvalues 1+/-sqrt3 |
| signed ternary-grade counts | X^2+2X+2 | local2 eta coefficients, eigenvalues -1+/-i |

The normalized trace law is also precise. For any2 by2 matrix K,
tr(K^(2r))=tr(K^r)^2-2 det(K)^r. The native golden matrix has determinant1;
the eta matrix has determinant2. Dividing its r-th trace by 2^(r/2)
recovers the same J->J^2-2 law. Hecke coefficients have initial value1,
whereas power traces have initial value2: at j=2 they are2 and0 here.
The shared polynomial alone does not identify these scalar sequences.

The four-vertex tournament's presentation/path counts are useful inspiration
for retaining a marking before counting. No bijection from those tournament
presentations to the present shadow letters or signed words is claimed.
The explicit preserved predicate here is native-word nonemptiness.

## 4. Restoring the source changes the eta cancellation

Let D_w be the full native source cylinder. For a nonempty native word
write F_w(x)=(P x+B_w)/Q, P=3^(2j), Q=2^A. The inherited closed-form guard
from [translation decoding](translation_phase_decoder_20261005.md) is

    P n+B_w = t Q mod s Q,
    (s,t)=(16,11) if w ends L, otherwise (2,1).     (8)

Because P is odd, this is one odd residue modulo sQ. Under normalized Haar
measure on odd2-adic integers its mass is

    mu(D_w)=2^-A * (1/8 if w ends L, else 1).      (9)

This is also its natural density within positive odd integers. Define

    Z_j(n)=sum_(w in N, grade(w)=j) sign(w) 1_(n in D_w).

Then the mean is exactly

    sum_(j>=0) E_odd[Z_j] z^j
       = (1+7z/128)/(1+3z/16+z^2/64-z^3/2048).    (10)

Proof: use letter weights h=z^3/1024, g=-z/8, l=-z/16,
a=b=-z^2/128 in (4), and terminal vector (1,1/8) instead of (1,1).
The numerator becomes 1-7l/8. Notice h is twice bl: retaining binary
cost has restored the cubic term. Equation(9) supplies the terminal factor.
For each coefficient the sum is finite, so no exchange of infinite sums
or limits is involved.

The first mean values are

    1, -17/128, 19/2048, 27/32768, -191/524288, ... .

At grade3, (7) is zero. In contrast Z_3(155)=1, with sole legal word H,
and Z_3(219)=-1, with sole legal word LGG. Over all32768 odd residues
modulo65536 the literal signed sum is27, independently confirming (10).
These are overlapping-cylinder signed sums, not union densities or a
claim that any counted word has a grounded terminal proof.

There is an exact operator reason the scalar cancellation cannot be lifted
unchanged. Put

    (U_i h)(n)=1_(n in D_i) h(F_i(n)).

Then U_L U_B=0, whereas U_H 1 is1 at155. A representation of the shadow
fusion as these native arithmetic operators would require U_H=U_L U_B,
which is false. The correction is to keep the operators themselves:

    Z_j = -(U_G+U_L) Z_(j-1) -(U_A+U_B) Z_(j-2) + U_H Z_(j-3),
    Z_0=1, Z_j=0 for j<0.                         (11)

Partition words by their first letter to prove (11). Invalid transitions
vanish by their guard indicators. This gives a faithful recurrence on the
source functions. It is generally an infinite-state operator problem;
the scalar eta clock is a useful reader, not its replacement.

### A golden plane survives the source weighting modulo3

There is a positive repair even after the scalar eta identity fails.
Put r_j=8*16^j E_odd[Z_j]. Equation(10) gives integers with series

    sum r_j z^j = (8+7z)/(1+3z+4z^2-2z^3).

Their recurrence matrix and initial vector are

    C=[[-3,-4,2],[1,0,0],[0,1,0]],
    (r_2,r_1,r_0)=(19,-17,8).

**PROVED.** Modulo3 its characteristic polynomial splits as

    X^3+X+1 = (X-1)(X^2+X+2).

More explicitly, with S=[[1,1,0],[1,1,2],[1,2,2]], det(S)=1 modulo3,

    C S = S diag(1,-M) mod3.                      (11a)

So the source-averaged clock has one fixed line and a golden F9 plane.
The initial vector reduces to the second column of S, so its fixed-line
component is zero. Equivalently, the numerator of the scalar series
cancels the factor1-z modulo3, giving

    r_j = 2*(-1)^j a_(2^j) mod3,
    E_odd[Z_j] = (-1)^j a_(2^j) mod3              (11b)

where the latter reduction is in Z[1/2]. This is an exact congruence at
every grade. It is stronger than noticing that two clocks have period8.
It still forgets higher ternary digits and the individual source. Already
r1=-17 differs from2*(-1)*a2=4 modulo9.

The complete source matrix has order8*3^(k-1) modulo3^k for every k>=1.
To prove it, (11a) gives exact order8 at k=1, and direct multiplication gives

    C^8 = [[505,-714,198],[99,802,-318],[-159,-378,166]].

The minimum3-adic valuation of an entry of C^8-I is1. If a matrix is
I+3^e D with e>=1 and D nonzero modulo3, cubing raises this minimum
valuation exactly by1, by the binomial expansion. Iteration proves the
order formula, including its lower bound. The extra fixed eigenvalue
modulo3 need not remain fixed at deeper precision: the characteristic
polynomial has its lifted root4 modulo9 instead of1. The one-plus-two
split is a concrete prime-local structure; no identification with an
arbitrary numerical triplet is needed.

This gives a new connection contract: the source-weighted cubic clock maps,
by reduction and the basis S, to a fixed line plus the negative golden
operator. It preserves this linear recursion modulo3 and destroys deeper
carry information. The residual ternary digits are the required sidecar.

## 5. Two prime trees and an ordinary-height gate

Each controller map is a similarity of Z3 onto a residue ball of measure
3^(-2j). Its translation tags modulo9 are, in H,G,A,B,L order, 3,4,2,5,6.
These balls are disjoint and avoid the seed0. This recovers the exact finite
translation history. The sum of their first-level masses is181/729.

At word length k the sum of formal cylinder masses is(181/729)^k, tending
to zero. The same bound holds after restricting to native words. An
infinite backward address therefore lies in a null subset of Z3, though
finite words remain lossless. This is a statement about translation space,
not the density of integer sources or a Collatz coverage bound.
For example

    F_(L^k)(0)=-(3/7)(1-(9/16)^k) -> -3/7 in Z3.

The limit has an infinite L address and is not a finite dyadic translation.
It parallels the denominator9 nonmember limit in
[centroid membership](centroid_membership_20261005.md); both require a
stopping condition even when all inverse letters are unambiguous.

The two clocks also give the elementary product identity

    |3^(2j)/2^A|_infinity * |3^(2j)/2^A|_2 * |3^(2j)/2^A|_3 = 1. (12)

Ternary contraction alone hides binary expansion and ordinary growth.
For G, G(n)+5=(9/8)(n+5): the anchor -5 restores the affine carry, and
the real multiplier9/8 still grows. Payment requires the actual source,
the carry and the funded potential, not just a prime-tree depth.

### Completing a fixed-parent decoder with a finite height bound

Fix a certified parent u=3 mod8 in the inherited refuel tree. Let
kappa(u) in1..729 solve 4^kappa(3u+1)=334 mod2187 and put

    C=1024*4^kappa(3u+1),  A0=4^729,
    psi_u(t)=(C A0^t-3031)/2187,  t>=0.            (13)

Every such child is certified. The inherited proof gives

    v3(psi_u(t)-psi_u(s))=v3(t-s),
    psi_u(t+e3^k)=psi_u(t)+e3^k mod3^(k+1).        (14)

The compatible residue bijections extend to an isometric homeomorphism
Z3->Z3: lift one ternary digit at a time for surjectivity, and use (14)
for uniqueness and continuity. For a fixed integer n, let t_k in[0,3^k)
be its unique depth-k branch address.

**Height completion theorem (PROVED).** The following are equivalent:

1. n=psi_u(t) for a nonnegative integer t;
2. t_k is eventually constant;
3. the sequence t_k is bounded in the ordinary absolute value.

The implication1->2 uses t_k=t mod3^k. The prefixes are nondecreasing
ordinary integers, so3->2. Under2, n-psi_u(t) is an integer divisible by
every3^k and is therefore0. For a nonmember the prefixes instead tend to
infinity, and psi_u(t_k) tends to positive infinity while converging to n
3-adically. This identifies exactly how residue completion can lose the
source integer.

There is a finite algorithm, not just a criterion about an infinite tail.
For supplied n find

    b=max{t>=0: psi_u(t)<=n},

using b=-1 if the set is empty. Strict ordinary growth of (13) permits
doubling and binary search with exact integer comparisons. If b=-1 reject.
Otherwise choose k with3^k>b, decode t_k, and accept exactly when
t_k<=b and psi_u(t_k)=n. Any solution lies in[0,b], where its residue
modulo3^k identifies it uniquely. This proves total acceptance/rejection
for this family. The required number of ternary digits is O(log log n)
for fixed u, because b grows logarithmically in n; this is not a claim
of logarithmic bit complexity for reading or comparing n.

For u=3, kappa=675 and the first child already has1353 binary digits.
The target155 has compatible branch prefixes

    2, 8, 8, 62, 143, 143, 143, 4517, 11078, 30761, ...

at precisions3,9,27,... . All finite congruences are realizable by certified
children. Nevertheless b=-1, so exact membership is rejected immediately.
The unbounded-prefix conclusion here is proved by the first-child bound,
not inferred from the displayed digits. This does not question convergence
of155: it is only a nonmember of this fixed-parent subfamily.

The new stopping rule is useful when switching constructors: return an
exact rejection of this candidate family and try another grounded rule,
rather than interpreting continuing digit refinement as progress toward
the same source. It does not prove that some family will eventually accept
every positive input. That universal entry statement is still OPEN.

## 6. A genuine prime-tree heavy tail, without misidentifying its source

There is an elementary model in which a small-prime effect can be proved.
For a prime p and0<rho<1, let H_j be normalized Haar measure on p^j Z_p,
and set nu=(1-rho) sum_(j>=0) rho^j H_j. A character of conductor p^ell
has Fourier coefficient rho^ell, because its integral over H_j is1
exactly when j>=ell. Thus this is an actual projective probability law.

For a uniform frequency U modulo p^n, and r>0, exact valuation counting gives

    E |nu_hat_n(U)|^(2r)
      = p^-n + (p-1)/p * sum_(j=0..n-1) p^-j rho^(2r(n-j)). (15)

The first term is the zero frequency. The other terms separate the
primitive conductor layers. Hence the moment decays like

    rho^(2rn) if p*rho^(2r)>1;
    n*p^-n if p*rho^(2r)=1;
    p^-n if p*rho^(2r)<1,

up to positive constants in the two strict cases. With
Y=rho^-n |nu_hat_n(U)|, for0<=j<=n,

    Pr[Y>=rho^-j]=p^-j.                            (16)

This is a truncated power-law tail with exponent log(p)/log(1/rho).
Smaller p gives the heavier tail at fixed rho. Conditioning on primitive
U removes it completely in this radial model: Y=1.

Equations(15)–(16) make a precise version of the prime-ary intuition and
explain why mixing conductor depths is dangerous. They do not identify
nu with the Syracuse law. The E5q correction found just such nonprimitive
frequency contamination in an earlier arithmetic ensemble. The surviving
E5r finite difference between uniform-start q=3 and q=5 towers remains a
different question about interlevel coupling. A tail model alone provides
neither that rate nor the required estimate at the fixed seed1.

## 7. Connection map and the next proof object

| source -> target | exact map and preserved predicate | loss / required sidecar |
|---|---|---|
| native controller -> free shadow | H expands to LB; bijective finite history, grade preserved | binary cost and credit corrections(1), native action(2) |
| signed shadow -> eta local2 | signed grade enumeration(7) | source, carry, funding and terminal proof |
| eta companion -> golden F9 clock | reduce K modulo3, conjugate by J | characteristic, torsion basis and scalar-output initial data |
| native word -> source function | full guard indicator in(11) | averaging loses the particular integer |
| branch address -> rooted child | isometry(14), actual constructor(13) | finite integer address needs the ordinary-height gate |
| p-ary frequency tree -> moments | conductor stratification(15) | aggregated moments forget which primitive layer contributed |

The proposed proof object is a finite **source-labelled certificate mesh**:
a shadow word with its fusion positions, its native source/target progression,
the actual source parameter, both clocks and affine carry, a height or paid
rank witness for each recursive call, and a checked terminal proof. Each
object maps to one integer; an integer may have many such objects. The
finite translation decoder can recover its controller history, and the
height-completed branch decoder can recognize a selected family exactly.
This combines existing checked constructors without assuming coverage.

The decisive targets now have separate types: prove a source-resolved
bound for operators such as(11), or prove every supplied integer enters a
grounded, well-founded certificate family. A small signed recurrence can
suggest cancellation, but the witness155 prevents us from substituting
its eta value for the required source-resolved statement. Complete residue
trees suggest addresses, but (13)–(14) plus ordinary height decide whether
they name the source. This is the reusable move: retain the coordinate in
which the obstruction is measurable, and make its rejection finite.

## 8. Reproduction and validation scope

Run from the repository root:

    python3 -B 04-computation/experiments/prime_ary_controller_eta_20261005.py
    python3 -B -O 04-computation/experiments/prime_ary_controller_eta_20261005.py

The experiment enumerates all28821 native words and all28821 free shadows
of grade at most10; checks boundary fusion and associativity; compares
independent recurrences and tiling sums through grade120; expands the eta
product through q^1024 by both literal factors and logarithmic derivative;
checks the explicit cubic/golden intertwiner and matrix orders modulo3^k
for k=1 through6;
checks every grade-three word at all32768 odd residues modulo65536;
constructs every branch residue through ternary depth6; checks 24-digit
addresses of0,1,155,223,233; and tests the height decoder on actual branches
and same-guard nonmembers. The radial-law census uses p=2,3,5,7,11,
depths1–5 and moment powers2,4. These are finite controls, not all-height
proofs by enumeration; the proofs above provide the quantified statements.
All 761,476 explicit checks passed. No Python assertion is load-bearing;
normal and optimized executions produced byte-identical JSON. The hashes are:

    script SHA256 e5ab42335154a68e4731a4bacc51a566a865f801599b879e00373290caac5e12
    JSON   SHA256 ab8821878f0adffc44b463776935f6079d13cbc07336934c9643d5e3667f0fd5

No Lean formalization or universal source coverage is claimed.
