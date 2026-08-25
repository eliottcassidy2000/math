---
id: THM-4037
title: "Centered binomial parity and singular fibres"
status: >
  PROVED + VERIFIED-EXACT. Centered reflection gives the exact sign of every
  binomial-basis atom. For every all-even rank packet, complete exact-period
  representation counts are divisible by the full sign-orbit size at even
  moduli and have one explicitly classified odd target at odd moduli. For
  ranks 2,4,6,8 that target is -3453/32768. The odd centered derivatives give
  a degree-24 critical-value eliminant; at eight explicitly listed primes the
  target fibre has exactly the displayed singular points and none lifts to
  level two. This is a local multiplicity and singular-fibre theorem, not a
  least-hole theorem or a classification of integer counterexamples.
source: codex sun-minimal-sturmian / odd-function bridge, 2026-08-24
depends_on:
  - LEM-004-tournaments-are-odd-functions
  - THM-4027-sun-two-four-six-eight-universal-modular-solubility
related:
  - THM-4026-sun-two-four-six-eight-binomial-counterexample
  - THM-4028-sun-two-four-six-eight-average-order-criticality
  - MISTAKE-225
parity_script: 04-computation/sun_2468_centered_parity_thm4037.py
parity_output: 05-knowledge/results/sun_2468_centered_parity_thm4037.out
singular_fibre_script: 04-computation/sun_2468_critical_singular_primes_thm4037.py
singular_fibre_output: 05-knowledge/results/sun_2468_critical_singular_primes_thm4037.out
---

# THM-4037 -- centered binomial parity and singular fibres

**PROVED + VERIFIED-EXACT.** There are two results here.  The first is a
general parity theorem for complete period boxes of even-rank binomial atoms.
The second applies the same centered involution to the critical fibres of
Sun's `2-4-6-8` polynomial.  The second result gives eight exact singular
primes for the counterexample target; it does not assert that those are all
singular primes.

## 1. Inheritance and the typed odd-function connection

For an involution `iota` without fixed points, LEM-004 identifies an
orientation with a sign function `sigma(iota d)=-sigma(d)`.  Its load-bearing
mechanism is the fixed-point obstruction: an odd `+-1` function cannot exist
at a fixed point.  Here the involution acts on a binomial top index, not on a
tournament difference set.  The binomial value is in the trivial
representation, while its first difference is in the sign representation.

THM-4027 supplies the exact period of `n -> C(n,k) mod m`.  That period is
essential below: using only `m` would be wrong at primes at most `k`.
MISTAKE-225 is the corrected near miss.  A centered index reflection,
tournament converse, torus inversion, and permutation sign are different
involutions unless an equivariant map preserves the target predicate.

The connection contract is therefore

```text
source:      odd/even functions on an involutive set (LEM-004)
target:      exact-period index sets for C(n,2r), and their shell increments
map:         n -> 2r-1-n;  Delta C(n,2r)=C(n,2r-1)
preserved:   trivial/sign C2 type, orbit pairing, and fixed-point parity
lost:        positive chamber, magnitude, bounded height, and common lift
sidecar:     exact period, role label, legal lower bound, and height box
hostile:     THM-4026 is locally represented at every modulus but has no
             bounded integral representation.
```

Thus the result is not the slogan "Sun's atoms are tournament odd
functions."  The lawful transfer is the fixed-point and sign-isotypic
calculus.

## 2. The centered binomial parity ladder

For `k>=1`, put

```text
B_k(n)=C(n,k),             iota_k(n)=k-1-n,            (1)
```

where `C(n,k)` is the generalized binomial polynomial at every integer `n`.
Then

\[
B_k(\iota_k(n))=(-1)^kB_k(n).                          \tag{2}
\]

Indeed,

\[
\prod_{j=0}^{k-1}(k-1-n-j)
=(-1)^k\prod_{j=0}^{k-1}(n-j),
\]

and division by `k!` proves `(2)`.  Pascal subtraction gives the second exact
identity

\[
\Delta B_k(n):=B_k(n+1)-B_k(n)=B_{k-1}(n).             \tag{3}
\]

Consequently an even-rank atom is centered-even, whereas its shell increment
is centered-odd about the half-step-shifted center:

\[
B_{2r-1}(2r-2-n)=-B_{2r-1}(n).                         \tag{4}
\]

In the centered coordinate

```text
u=2n-(2r-1),
```

the even atom factors as

\[
 B_{2r}(n)=Q_r(u^2),\qquad
 Q_r(s)={\prod_{j=1,3,\ldots,2r-1}(s-j^2)
          \over 2^{2r}(2r)!}.                          \tag{5}
\]

For Sun's four roles this reads

\[
\begin{aligned}
\binom w2&={u_2^2-1\over8},\\
\binom x4&={(u_4^2-1)(u_4^2-9)\over384},\\
\binom y6&={(u_6^2-1)(u_6^2-9)(u_6^2-25)\over46080},\\
\binom z8&={(u_8^2-1)(u_8^2-9)(u_8^2-25)(u_8^2-49)
              \over10321920}.
                                                               \tag{6}
\end{aligned}
\]

The legal chamber is

```text
u_2>=3,      u_4>=3,      u_6>=5,      u_8>=7,
```

with every `u` odd.  It chooses a positive section of the sign quotient and
is not itself invariant under `u -> -u`.  Equivalently, each even atom is the
positive-half prefix mass of its odd increment:

\[
\binom n{2r}=\sum_{j=2r-1}^{n-1}\binom j{2r-1}.         \tag{7}
\]

This boundary loss is why the centered symmetry alone cannot classify exact
integer holes.

## 3. Exact-period fixed points and the all-modulus parity law

Let `m>=2`, and write `p^a || m`.  By the exact prime-power period theorem in
THM-4027, the least period of `B_k(n) mod p^a` is

\[
P_{p^a,k}=p^{a+\lfloor\log_p k\rfloor}.                \tag{8}
\]

The factors for distinct primes are coprime.  Their product

\[
P_{m,k}=\prod_{p^a\parallel m}
 p^{a+\lfloor\log_p k\rfloor}                         \tag{9}
\]

is a period modulo `m` by CRT.  Conversely, any period modulo `m` is a period
modulo every `p^a`, hence is divisible by every factor in `(8)` and therefore
by `(9)`.  Thus `(9)` is the exact least period.  In particular,

```text
P_(m,k) is even  iff  m is even.                       (10)
```

Fix even ranks `k_i=2r_i`, `1<=i<=s`, and define the complete-period fibre

\[
R_m(t)=\#\left\{(n_i)\in\prod_i\mathbf Z/P_{m,k_i}:
       \sum_i\binom{n_i}{k_i}\equiv t\pmod m\right\}. \tag{11}
\]

On the `i`-th period apply

\[
\iota_i(n_i)=k_i-1-n_i\pmod {P_{m,k_i}}.              \tag{12}
\]

Periodicity and `(2)` show that `(12)` preserves the corresponding atom and
hence every fibre in `(11)`.

### Even modulus

If `m` is even, `P_{m,k_i}` is even, while `k_i-1` is odd.  The fixed-point
equation

\[
2n_i\equiv k_i-1\pmod {P_{m,k_i}}                     \tag{13}
\]

has no solution.  The independent coordinate reflections therefore define a
free action of `(Z/2)^s` on each fibre in `(11)`: a nonidentity group element
could fix a tuple only if every coordinate it flips had a fixed point.
Consequently

\[
2^s\mid R_m(t)\qquad\hbox{for every }t                 \tag{14}
\]

whenever `m` is even.

### Odd modulus

If `m` is odd, each period is odd and `(13)` has one solution `c_i`.  All
other indices occur in two-element orbits with equal atom value.  Therefore
the one-role histogram is, modulo two, a point mass at

```text
a_(2r)(m)=B_(2r)(c) mod m.                             (15)
```

The fixed value has the closed form

\[
a_{2r}(m)=(-1)^r\binom{2r}{r}16^{-r}\pmod m.           \tag{16}
\]

Here is the denominator-safe justification of `(16)`.  For each `p^a || m`,
the fixed class `c mod P_(p^a,2r)` agrees in `Z_p` with the half-integer
`r-1/2`.  Choose nonnegative integers in this period class converging
`p`-adically to `r-1/2`.  Exact periodicity makes all their binomial values
congruent to `B_(2r)(c) mod p^a`; polynomial continuity gives the limit

\[
\binom{r-1/2}{2r}
=(-1)^r{(2r-1)!!^2\over2^{2r}(2r)!}
=(-1)^r{\binom{2r}{r}\over16^r}.                      \tag{17}
\]

The denominator in `(17)` is a power of two and is a unit at every odd
prime.  Passing to the closed congruence class modulo `p^a`, then using CRT,
proves `(16)` even at the denominator primes `p<=2r`.

Convolution of the one-role histograms now gives the exact parity
classification

\[
R_m(t)\equiv
\mathbf1\!\left[t\equiv
 \sum_i(-1)^{r_i}\binom{2r_i}{r_i}16^{-r_i}\pmod m
\right]\pmod2                                          \tag{18}
\]

for every odd `m`.  The all-even-rank hypothesis is the failure boundary:
an odd-rank atom is anti-invariant under `(12)`, so its value fibre is not
preserved.

## 4. The Sun parity fingerprint

For ranks `2,4,6,8`, the distinguished target in `(18)` is

\[
-{1\over8}+{3\over128}-{5\over1024}+{35\over32768}
=-{3453\over32768}.                                    \tag{19}
\]

Hence, for every odd `m`, the complete-period representation count at an
integer `M` is odd exactly when

\[
m\mid 32768M+3453.                                    \tag{20}
\]

At the THM-4026 target

```text
N=896315812331399,
```

the exact factorization is

\[
32768N+3453
=5\cdot7\cdot42667\cdot199447\cdot98610539.           \tag{21}
\]

This fingerprint is not a counterexample test.  The target has even
complete-period count modulo `3`, `11`, and `33`, despite belonging to the
known sparse class modulo `33`; its low-pair character primes `17` and `23`
also do not divide `(21)`.  THM-4027 proves all those local counts are
positive.  Parity, local density, and bounded exact coverage are distinct
coordinates.

## 5. Odd derivatives and the critical-value eliminant

Work now over a field of characteristic `p>8`, so every factorial in `(6)` is
a unit.  With `s=u^2`, differentiation of `(5)` in the original top index
gives

\[
{d\over dn}Q_r(u^2)=4uQ_r'(u^2).                      \tag{22}
\]

The derivative is centered-odd.  Its critical squared coordinates are the
roots of

\[
sQ_r'(s)=0.                                           \tag{23}
\]

Direct differentiation and elimination give the following monic
critical-value polynomials:

\[
\begin{aligned}
V_1(Y)&=Y+{1\over8},\\
V_2(Y)&=(Y+{1\over24})(Y-{3\over128}),\\
V_3(Y)&=(Y+{5\over1024})
 \left(Y^2+{4Y\over243}-{1\over6075}\right),\\
V_4(Y)&=(Y-{35\over32768})
 \left(Y^3+{9Y^2\over640}-{3Y\over89600}
             -{1\over14049280}\right).
                                                               \tag{24}
\end{aligned}
\]

An exact certificate for `(24)` is

\[
V_r(Q_r(s))\equiv0\pmod {sQ_r'(s)}                    \tag{25}
\]

in `Q[s]`.  Both sides expose exactly `r` critical values with multiplicity;
the singular-fibre companion verifies all four polynomial divisions using
exact rational arithmetic.

Let `Z(V_r)` denote the multiset of complex roots of `V_r`, and define

\[
D(T)=\prod_{\alpha_r\in Z(V_r),\ 1\le r\le4}
 \left(T-\alpha_1-\alpha_2-\alpha_3-\alpha_4\right).
                                                               \tag{26}
\]

This is a monic degree-`24` polynomial in `Q[T]`.  Equivalently it is the
iterated additive resultant of the four polynomials in `(24)`.  The rational
linear/quadratic/cubic factor split in `(24)` gives eight additive-resultant
branches of degrees `1,3,2,6,1,3,2,6`; their degrees sum to `24`.

If the fibre

\[
F(w,x,y,z)=\binom w2+\binom x4+\binom y6+\binom z8=N  \tag{27}
\]

has a singular point over `F_p`, then every coordinate value is a root of the
corresponding `V_r`, so

\[
D(N)\equiv0\pmod p.                                   \tag{28}
\]

The exact branch evaluation verifies `D(N)!=0`, so `(28)` gives a finite,
effective prime candidate set.  It is only necessary over `F_p`: an
algebraic critical squared coordinate must also lie in `F_p` and be a square
there before it yields an actual centered coordinate `u`.  This field and
square sidecar cannot be dropped.

The all-centered branch of `(26)` is exactly

\[
T+{3453\over32768}.                                   \tag{29}
\]

Thus `(21)` is simultaneously the complete-period parity fingerprint and the
canonical centered-singularity fingerprint.

## 6. Eight exact singular fibres and their first-lift failure

The singular-fibre companion constructs all eight resultant branches by
Sylvester determinants and partially factors their values at `N`.  Among the
twenty exposed prime factors above `8`, it tests all roots of `(23)` in
`F_p`, retains the square-coordinate sidecar, and obtains the following exact
table.  Here

```text
s=(u_2^2,u_4^2,u_6^2,u_8^2),
c_p = number of singular points over F_p,
rho_p = (F(n_bar)-N)/p mod p
```

for the representatives `0<=n_bar_i<p`.  Every point at a given prime has
the displayed same nonzero residual.

| `p` | critical squared-coordinate row `s` | `c_p` | `rho_p` |
|---:|---|---:|---:|
| `31` | `(0,0,28,8)` | `4` | `7` |
| `499` | `(0,5,0,0)` | `2` | `356` |
| `35,837` | `(0,0,0,34,969)` | `2` | `27,183` |
| `42,667` | `(0,0,0,0)` | `1` | `25,517` |
| `199,447` | `(0,0,0,0)` | `1` | `69,875` |
| `508,471` | `(0,0,0,291,471)` | `2` | `151,075` |
| `98,610,539` | `(0,0,0,0)` | `1` | `70,179,963` |
| `931,932,289` | `(0,5,0,0)` | `2` | `529,232,077` |

The point counts follow from the exact finite-field roots: a nonzero square
coordinate has two signs, while a zero coordinate has one.  The program uses
`gcd(C(s),s^p-s)` to obtain every rational critical root, deterministic
finite-field splitting, and Tonelli--Shanks to enforce the square sidecar.
It then checks the four derivatives and the target equality directly.

To prove the nonlift assertion, let `n_i` represent a singular point and let
`t_i` be arbitrary.  Since `p>8`, Taylor expansion of the binomial polynomials
is denominator-safe and gives

\[
\binom{n_i+pt_i}{k_i}
\equiv\binom{n_i}{k_i}+pt_iB'_{k_i}(n_i)
\equiv\binom{n_i}{k_i}\pmod {p^2}.                    \tag{30}
\]

The sum at every lift is therefore congruent to the chosen representative.
Every `rho_p` in the table is nonzero, so none of the displayed singular
points lifts to a solution modulo `p^2`.

For the three all-centered primes above `8`, there is also a uniform form:

\[
\text{the centered residue tuple lifts modulo }p^2
\quad\Longleftrightarrow\quad
p^2\mid32768N+3453.                                   \tag{31}
\]

The factorization `(21)` is squarefree, so `(31)` independently forces the
three centered nonlifts.

If `S_p` is the number of solutions to `(27)` modulo `p`, each nonsingular
solution has exactly `p^3` lifts to the next prescribed target digit, while
the `c_p` singular solutions in the table have none.  In THM-4027's relative
local-density normalization this gives

\[
\sigma_{p^2}(N)=\sigma_p(N)-{c_p\over p^3},\qquad
\sigma_{p^a}(N)=\sigma_{p^2}(N)\quad(a\ge2).           \tag{32}
\]

The first two rows recover THM-4027's `p=31,499` hostile examples; the other
six are additional exact local sidecars.

## 7. Resultant and scope firewalls

The partial branch factorization exposes the twenty primes

```text
13,31,73,109,307,401,479,499,607,1753,3221,21617,
34667,35837,42667,199447,508471,5503241,98610539,
931932289.
```

The square/field audit rejects

```text
13,73,109,307,401,479,607,1753,3221,21617,34667,5503241
```

and retains exactly the eight rows in the table.  This is an exact
classification **within the exposed factors only**.  Several large branch
cofactors remain unfactored.  They may contain further primes, and a prime
dividing a resultant branch may still fail the field or square sidecar.
Accordingly:

1. the eight displayed primes are proved singular/nonlifting primes, not all
   singular primes of the target;
2. divisibility of `D(N)` is necessary, not sufficient, for an `F_p`
   singular point;
3. local singularity or local multiplicity parity is neither necessary nor
   sufficient for an exact integer hole;
4. `(18)` counts a complete modular period box, not the bounded positive box
   in THM-4026; and
5. no least-counterexample, exception-density, or all-counterexample
   classification follows.

These firewalls are load-bearing.  THM-4026 remains the decisive hostile to
every attempt to replace height by complete local data.

## 8. Exact companions

The parity companion checks `(2)--(4)` at `2,196` integer lanes and checks
every modulus `2<=m<=100`, all four exact role periods, all reflected lanes,
and the complete circular convolution (`185,308` modular lanes):

```text
python -B 04-computation/sun_2468_centered_parity_thm4037.py
```

The dependency-free singular-fibre companion uses only the Python standard
library.  It verifies the four exact critical-value divisions, constructs all
eight additive resultants by rational Sylvester determinants, checks the
partial prime valuations, obtains complete finite-field critical roots,
enforces the square sidecar, and checks every residual in the table.  Normal
and optimized runs are byte-identical:

```text
python -B 04-computation/sun_2468_critical_singular_primes_thm4037.py
python -B -O 04-computation/sun_2468_critical_singular_primes_thm4037.py
```

The computation certifies only the finite rows explicitly labelled as exact;
the parity and eliminant statements have the proofs above.  **QED.**
