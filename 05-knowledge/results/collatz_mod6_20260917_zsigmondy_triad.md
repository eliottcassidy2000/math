# Zsigmondy triad: Bang's 63, the critical orbit of x^2-7/4, and the 3-cycle parabola

**Status: PROVED (scoped) + FINITE-EXACT + CITED; two conjectures REFUTED with minimal witnesses; OPEN residue named.**
PROVED: Bang's `n=6` exception is exactly Catalan's `3^2-2^3=1`, equivalently no prime has `ord_p(2)=6`;
for `n<=6` the critical-orbit numerator `N_n` of `x^2+a/b` factors as `prod_(d|n) H_d(a,b)` with pairwise coprime
factors, so "no primitive prime at `n`" is the Thue equation `H_n(a,b)=+-1`; complete rational preperiodic graphs
and complete "which parameters make `x` preperiodic" lists by a valuation/escape argument; the parabola identities
`c=-7/4-s^2`, `lambda(s)`, `c(t)=-7/4-s(t)^2`. CITED: Bang 1886 / Zsigmondy 1892; PARI `thue` (Bilu-Hanrot,
unconditional flag) closing the Thue equations for `n=3,4,5` over all of `Q`; PARI `nfdisc`; Morton 1998,
Flynn-Poonen-Schaefer 1997, Stoll 2008; Krieger via the geometry note. FINITE-EXACT: Bang table `n<=40`, the
1601-parameter census to `n=8`, unit and Thue searches, the `H_4,H_5,H_6=+-1` boxes. REFUTED: "`f^n(0)=c/2^(2^n-2)`
at every real parabolic parameter" (witness `n=4`) and "`63 <-> -7/4`" as a structural map. OPEN: `H_6(a,b)=+-1`
beyond the finite box, everything at `n>=7`, and a uniform resultant theorem. **No novelty claim** for Bang/Catalan,
the Gleason product, or the Chebyshev conjugacy; the `n=3` Thue reduction, the parabola, the parabolic cycle and the
`Q(2cos 2pi/7)` identification are **inherited** from the geometry note and only re-verified here. No Collatz
statement is claimed. Session `collatz-mod6-20260917`, lane `zsigmondy_triad`, recovered and finalized 2026-09-21.

## Inheritance and concept board

Read and not re-derived: [the Collatz braid note](arithmetic_braids_20260917_collatz.md) (three-row typing,
row power-of-two exponents `1,5,3 mod 6` from `ord_9(2)=6`, inverse fibre `R(n)=4n+1`),
[the summand note](arithmetic_braids_20260917_summand.md), [the divisor note](arithmetic_braids_20260917_divisors.md),
and above all [the geometry note](arithmetic_braids_20260917_geometry.md): its (4) is the `eta`-parabola
`c=-(eta^2+7)/4` with `disc P=(eta^2+eta+7)^2` (here `eta=2s`), its (9) is the parabolic cycle
`P(x)=x^3+x^2/2-9x/4-1/8` with `Phi_3=P^2` and multiplier `1` at `-7/4`, its section 2 proves that the reversed
conductor-7 cycle of `x^2-2` is carried by `L(x)=-x-1/2` onto this parabolic cycle, its (11) proves
`f^3(0)=aQ(a,b)/b^4` with `Q=+-1` iff no new prime (census `b<=2000`), and its section 3 cites Krieger's
*Primitive prime divisors in the critical orbit of z^d+c* (Theorem 1.1: at most 23 exceptional indices for infinite
critical orbits; Theorem 6.1: bound 3 with the third-index exception exactly at `-7/4` on its explicit parameter set;
read there on 2026-09-17). Also [THM-4139](../../01-canon/theorems/THM-4139-rational-three-cycle-order-six-lift-and-horizontal-carrier.md)
(`PrePer(x^2-29/16,Q)` (2)-(3), marked chart (13), the multiplier `35/8` in (32), `Psi_6` of squaring (36)-(38)),
[THM-4146](../../01-canon/theorems/THM-4146-rational-three-cycle-order-six-lift-horizontal-divisor-fibre-firewall.md)
(`sigma`-parametrization (10)-(12), section 7 on why `63` is a lift degree and not a scalar period), and for Gaussian
squaring [THM-3341](../../01-canon/theorems/THM-3341-u-spine-square-hypotenuse-transplant-and-triangular-plane-torsors.md)
section 4 and [THM-3333](../../01-canon/theorems/THM-3333-gaussian-square-farey-pythagorean-triangular-light-cone.md)
(read for the `x=2A/C -> x^2-2` conjugacy only; **SCOPE:** no map found from their Pell-hypotenuse selectors to the
critical-orbit exception).

The closest proved mechanism is the geometry note's `n=3` numerator test (11); this note generalizes it by exact
pairwise resultants to `n<=6`. The canonical hostile is the `n=2` family `c=-1+-1/b`: infinitely many exceptions,
so no Zsigmondy-type finiteness holds on the quadratic side without excluding low indices. The corrected near
misses are (a) the recovered draft's claims that the `n=3` Thue reduction, the parabola and the `-7/4` parabolic
cycle were new (they are geometry (11), (4), (9)), (b) its incidence list for `-7/4` that stopped at `-77/16` and
missed `-93/16`, and (c) its ledger count "14" for the `n=1` census failures (true value 13). The least-used sidecar
is the resultant table `Res(H_d,H_n)=+-1`, the quadratic analogue of "non-primitive primes of `Phi_n(2)` divide `n`".

Notation: `f_c(x)=x^2+c`, `G_n(c)=f_c^n(0)` (Gleason polynomials), `N_n` the numerator of `f_c^n(0)` in lowest
terms; a *primitive prime divisor at `n`* is a prime of `N_n` dividing no `N_m`, `m<n`.

| Lane | Object / representation | Predicate tested | Lost coordinate / hostile |
|---|---|---|---|
| Multiplicative | `2^n-1 = prod_(d\|n) Phi_d(2)` | primitive prime at `n` iff `ord_p(2)=n` | `n=6`: `Phi_6(2)=3` divides `6` |
| Critical orbit | `N_n = prod_(d\|n) H_d(a,b)` (`n<=6`) | primitive prime at `n` iff `0` has exact period `n` mod `p` | `n=2`: infinitely many exceptions |
| Parabola | `c=-7/4-s^2`, `s=sigma+1/2` | rational unmarked 3-cycle iff `-7/4-c` is a square | rational points need Morton's `t` |
| Preperiodic graphs | `PrePer(f_c,Q)` by valuation + escape | which special `c` are points of which others | tail vs period must be kept |
| Chebyshev | `x=y+1/y`, `y -> y^2` | conductors `7,9` of the two 3-cycles of `x^2-2` | `y=2` (Fermat) is not a unit-circle point |

## 1. Bang's exception is Catalan, and it is why the rows live modulo 9

**FINITE-EXACT (`n<=40`).** `2^n-1` has a primitive prime divisor for every `1<=n<=40` except `n=1` (vacuous,
`2^1-1=1`) and `n=6` (`63=3^2*7`). Every non-primitive prime of `Phi_n(2)`, `n>2`, is the largest prime of `n` and
divides `Phi_n(2)` exactly once; in the table this happens at `n=6,18,20,21` with primes `3,3,5,7`. For `2^n+1`,
`n<=40`, the unique index without a primitive prime is `n=3` (`9=3^2`). **CITED:** Bang 1886 / Zsigmondy 1892:
`a^n-b^n` has a primitive prime divisor except `n=1` with `a-b=1`, `n=2` with `a+b` a power of two, and
`(a,b,n)=(2,1,6)`.

**PROVED.** `2^6-1=(2^3-1)(2^3+1)=7*9` and `Phi_1(2)Phi_2(2)Phi_3(2)Phi_6(2)=1*3*7*3`. A primitive prime of
`2^6-1` cannot divide `2^3-1`, so it divides `(2^3+1)/gcd(2^3+1,2^2-1)=9/3=3`, but `3 | 2^2-1`: impossible. The
exception *is* the identity `3^2-2^3=1` (Mihailescu 2004 is CITED for uniqueness of Catalan but is not needed for
this direction). Equivalently **no prime `p` has `ord_p(2)=6`**, because a prime with `ord_p(2)=n` is exactly a
primitive prime divisor of `2^n-1` (hostile sweep `p<10^5` in the output). Hence `9` is the least modulus in which
`2` has order `6` (`2^e = 2,5,8 mod 9` iff `e = 1,5,3 mod 6`), which is precisely why the inherited three-row
exponent law is a mod-`9`, period-`6` statement with no prime modulus available. THM-4139 (36)-(38) and THM-4146
section 7 say the same from the squaring-map side (conductors `9,21,63`, `ord_9(2)=6`). **SCOPE:** this is the
only contact between this lane and the Collatz rows; nothing dynamical about the rows is claimed.

## 2. The critical orbit of x^2-7/4 and the bounded census

`f^n(0) = -7/4, 21/16, -7/256, -114639/65536, ...` with numerators `-7`, `3*7`, `-7`, `-3*7*53*103`,
`7*419*563*3407`, `-3*7*47237*636069519847`, ...; the primitive parts (gcd-stripped, no factoring needed) are
`7, 3, 1, 5459, 803701079, 30046015909012739, ...` (the `n=7` numerator is `-7` times a 38-digit prime; the `n=8`
primitive part has 73 digits). **FINITE-EXACT (`n<=8`):** `N_3=-7` has no primitive prime
divisor; `n=1,2,4,...,8` all do. Also `f^3(0)/f(0)=1/64=2^-6`, the "63" the user sees.

Census (FINITE-EXACT): `c=a/b`, `b in {1,2,4,8,16,32,64}`, `|a|<=200`, `gcd(a,b)=1` (1601 parameters), orbit to
`n=8`. Zero is preperiodic only for `c in {0,-1,-2}` (tail/period `2/1`, `0/2`, `0/1`), reported separately.
Non-preperiodic parameters whose `N_n` has no primitive prime:

| `n` | count | parameters |
|---|---|---|
| 1 | 13 | `c=+-1/b`: `1, +-1/2, +-1/4, +-1/8, +-1/16, +-1/32, +-1/64` |
| 2 | 12 | `c=-1+-1/b`: `-3/2, -1/2, -5/4, -3/4, -9/8, -7/8, -17/16, -15/16, -33/32, -31/32, -65/64, -63/64` |
| 3 | 1 | `-7/4` |
| 4..8 | 0 | none |

**PROVED (all `c=a/b` in lowest terms; `n=3` inherited from geometry (11)).** `G_1=c`, `G_2=c(c+1)`,
`G_3=c*H_3` with `H_3=c^3+2c^2+c+1`. So `N_1=a`: no primitive prime iff `|a|=1`, i.e. `c=+-1/b` (this includes the
cusp `1/4`). `N_2=a(a+b)`, and primes of `a+b` are coprime to `a`: no primitive prime iff `|a+b|=1`, i.e.
`c=-1+-1/b` (this includes the period-doubling parameter `-3/4`), **infinitely many**. `N_3=a*F(a,b)` with
`F=a^3+2a^2b+ab^2+b^3`, `F = b^3 mod a` and `F = b^3 mod (a+b)`, so every prime of `F` is primitive and failure at
`n=3` iff `F(a,b)=+-1`. The three parabolic parameters `1/4`, `-3/4`, `-7/4` fail at their own periods for three
different reasons (`|a|=1`, `|a+b|=1`, a cubic Thue unit). The draft's expectation that `-7/4` is the only
non-preperiodic failure in the box is true for `n=3` only.

## 3. The n<=6 reduction: Gleason product, pairwise resultants, Thue closures

**PROVED (exact division, `n<=6`).** Define `H_n` recursively by `G_n=prod_(d|n) H_d`; the divisions are exact,
every `H_n` is monic, squarefree and irreducible over `Q`, of degree `1,1,3,6,15,27` for `n=1..6`, with
`H_n(0)=1` for `n>=2`, `H_n(-1)=1` for `n>=3`, and `H_n(-2) = -2, -1, -1, 1, -1, -1` for `n=1..6`. The exact
resultants are

| | `H_1` | `H_2` | `H_3` | `H_4` | `H_5` |
|---|---|---|---|---|---|
| `H_2` | `1` | | | | |
| `H_3` | `-1` | `-1` | | | |
| `H_4` | `1` | `1` | `-1` | | |
| `H_5` | `-1` | `-1` | `1` | `1` | |
| `H_6` | `-1` | `-1` | `1` | `-1` | `1` |

**Theorem (PROVED for `n<=6`).** Write `H_n(a,b)=b^(deg) H_n(a/b)`. For `c=a/b` in lowest terms,
`N_n=prod_(d|n) H_d(a,b)` exactly (each `H_d(a,b) = a^deg mod b` is coprime to `b`, and the degrees add to
`2^(n-1)`). If a prime `p` divided `H_d(a,b)` and `H_n(a,b)` with `d<n`, then `p` does not divide `b`, and `a/b mod p`
would be a common root of `H_d` and `H_n` modulo `p`, so `p | Res(H_d,H_n)=+-1`: impossible. Hence the factors are
pairwise coprime, every prime of `H_n(a,b)` is primitive at `n`, and every other prime of `N_n` divides some `N_d`
with `d|n`, `d<n`. Therefore, for `n<=6`,

```text
N_n has a primitive prime divisor   <=>   |H_n(a,b)| >= 2,
non-Zsigmondy at n                  <=>   H_n(a,b) = +-1   (a Thue equation of degree 1,1,3,6,15,27).
```

**`n=3`, the plastic unit (PROVED + CITED).** `F(a,b)=N(a+b*rho^2)` in `Q(rho)`, `rho^3=rho+1` the plastic number
(discriminant `-23`, squarefree, so `Z[rho]` is the full ring of integers): the minimal polynomial of `-rho^2` is
`t^3+2t^2+t+1=H_3`, so the airplane centre is `c=-rho^2=-1.75487766625` and `-7/4` lies within `0.0049` of it.
Units `+-rho^k=A+C*rho^2` with vanishing `rho`-coefficient, `|k|<=300`: `k=-14` (`4rho^2-7`), `-5` (`2-rho^2`),
`-1` (`rho^2-1`), `0`, `2`, giving `c=A/C in {-7/4,-2,-1,0}` and `b=0` (`c=infinity`, not a parameter). A direct
search `1<=b<=10^6` gives `[(-2,1),(-1,1),(0,1),(-7,4)]`. **CITED (PARI/GP `thueinit(F,1)`/`thue`, Bilu-Hanrot,
flag 1 = no GRH):** `thue(F,+1)=[[-7,4],[-1,1],[0,1],[1,0],[2,-1]]` and `thue(F,-1)` is its negative. **Hence over
all of `Q` the `n=3` non-Zsigmondy parameters are exactly `c in {0,-1,-2}` (preperiodic) and `c=-7/4`, and `-7/4`
is the coefficient ratio in `rho^-14=4rho^2-7`.** (The geometry note had this FINITE-EXACT to `b<=2000` and
explicitly left the all-height Thue classification open; it is now closed by the certified solver.)

**`n=4,5` (PROVED + CITED).** PARI `thueinit(H_4,1)`: `thue(+1)` gives the eight sign-variants of
`(1,0),(0,1),(-1,1),(-2,1)` and `thue(-1)=[]`; `thueinit(H_5,1)`: `thue(+1)=[[-1,1],[0,1],[1,0],[2,-1]]`,
`thue(-1)` its negative. All solutions are `b=0` or `c in {0,-1,-2}`. **Therefore for every rational `c` with
infinite critical orbit, `N_4` and `N_5` have primitive prime divisors, and the only non-preperiodic rational
parameter whose critical orbit misses a primitive prime at some `3<=n<=5` is `c=-7/4`, at `n=3`.** This closes
the draft's `n=4` OPEN item and is consistent with Krieger's Theorem 6.1 as cited in the geometry note.

**`n=6` (FINITE-EXACT, OPEN beyond).** `H_6(a,b)=+-1` was searched over all `|a|<=3b`, `b<=100`, and over `a`
within `2` of `b` times each real root (`-1.9964, -1.9668, -1.9073, -1.7729, -1.476`) for `b<=5000`: only
`(-2,1),(-1,1),(0,1)`. The same boxes for `H_4` (real roots `-1.9408,-1.3107`) and `H_5` (`-1.9854,-1.8608,-1.6254`)
agree with PARI. A side run of `thueinit` on the degree-27 `H_6` did not finish within this lane's time budget, so
`H_6(a,b)=+-1` is OPEN beyond the box. Nothing is proved for `n>=7` (no resultants computed), and whether
`Res(H_d,H_n)=+-1` for all `d<n` is OPEN (UNCITED-RECOLLECTION: Gleason's simple-root theorem gives it only
`2`-adically).

## 4. Parabolic mechanism at -7/4, the parabola, and the cycle fields

**PROVED (dynatomic; parabolic cycle inherited from geometry (9)).** `disc_x Phi_3(x,c)=-(4c+7)^3(16c^2+4c+7)^2`
and `Res_x(Phi_3,(f^3)'-1)=((4c+7)(16c^2+4c+7))^3`, so `Phi_3` has a repeated root iff
`Delta_3(c)=(4c+7)(16c^2+4c+7)=64c^3+128c^2+56c+49=0`, real branch `4c+7=0`. At `c=-7/4`,
`Phi_3=(8x^3+4x^2-18x-1)^2/64` and the collided cycle has multiplier `8*prod x_i=8*(1/8)=1` (Vieta): a parabolic
saddle-node 3-cycle, numerically `(-1.7469796, -0.054958132, 1.3019377)`, shadowing `0 -> -7/4 -> 21/16 -> -7/256`.

**PROVED (inherited: geometry (4) in `eta=2s`).** With `s=sigma+1/2` (`sigma` the cycle sum, THM-4146 (10)):
`c=-7/4-s^2` and `lambda(s)=1-14s-4s^2-8s^3`. The unmarked 3-cycle curve is a parabola double-covering the `c`-line,
branched exactly at the parabolic parameter (`s=0`, `lambda=1`); a rational unmarked 3-cycle exists iff `-7/4-c` is
a rational square. `-29/16=-7/4-(1/4)^2`: `s=-1/4` is the rational AP cycle (`lambda=35/8`, the value displayed in
THM-4139 (32)), `s=+1/4` the cycle `64X^3+16X^2-164X+23` (`lambda=-23/8`). `-2=-7/4-(1/2)^2`: `s=-1/2,+1/2` give
`sigma=-1,0`. `disc_X P_s=(4s^2+2s+7)^2`, so every rational-`s` cycle field is `Q` or a cyclic cubic. Exact
embeddings (polynomial-remainder certificates) and PARI `nfdisc` cross-checks `[49, 961, 81, 49, 961]`:

| `s` | `c` | cycle cubic | field (certificate) | `4s^2+2s+7` |
|---|---|---|---|---|
| `0` | `-7/4` | `8X^3+4X^2-18X-1` | `Q(2cos 2pi/7)`, root `3/2-y^2`, `y^3+y^2-2y-1=0`; nfdisc 49 | `7` |
| `-1/4` | `-29/16` | `64X^3+48X^2-132X-35` | splits over `Q` | `27/4` |
| `+1/4` | `-29/16` | `64X^3+16X^2-164X+23` | conductor 31, root `-1/4+y/2`, `y^3-y^2-10y+8=0`; nfdisc 961 | `31/4` |
| `-1/2` | `-2` | `X^3+X^2-2X-1` | `Q(2cos 2pi/7)`; nfdisc 49 | `7` |
| `+1/2` | `-2` | `X^3-3X+1` | `Q(2cos 2pi/9)`; nfdisc 81 | `9` |

The field equality between the parabolic cycle at `-7/4` and the conductor-7 cycle at `-2` is the geometry note's
section 2 (an explicit affine map), re-verified here by embedding; **no novelty claim**.

The identity `64G_3(c)-c=c(4c+7)(16c^2+4c+9)` explains `f^3(0)=c/64` at `-7/4`; the complex parabolic factor
`16c^2+4c+7` does not divide it. **REFUTED:** the guess "`f^n(0)=c/2^(2^n-2)` at every real period-`n` parabolic
parameter" holds for `n=1,2,3` (`4G_2-c=c(4c+3)`, and `c=1/4` trivially) but
`gcd(2^14 G_4-c, Delta_4)=1` exactly, where `Res_x(Phi_4,(f^4)'-1)=((4c+5)(16c^2-8c+5)(64c^3+144c^2+108c+135))^4`
and the real period-4 saddle node is `c_4 ~ -1.94055078898`. Minimal witness `n=4`. Moreover `2^n-2=2n` only for
`n=3`, which is the sole reason the exponent `6` in `2^-6` equals Bang's index `6`.

## 5. Families and Morton's chart

**PROVED.** Fixed points are rational iff `c=(1-r^2)/4` (`r=1,3,5/2,9/2 -> 0,-2,-21/16,-77/16`); 2-cycles iff
`c=-(3+s^2)/4` (`s=0,1,2,3,5/2,9/2 -> -3/4` (degenerate, multiplier `-1`), `-1,-7/4,-3,-37/16,-93/16`). The chart
THM-4139 (13) satisfies `c(t)=-7/4-s(t)^2` with `s(t)=(t^3+t^2-2t-1)/(2t(t+1))` and
`sigma(t)=(t^3-3t-1)/(2t(t+1))`, invariant under the relabeling `rho(t)=-1/(t+1)`, and `t=(p1-p2)/(p0-p1)`
(`t=1` iff AP iff `-29/16`); the numerator of `s(t)` is the geometry note's `eta(t)` numerator. `c(t)=-2` iff
`t^3+2t^2-t-1=0` (`t=1/(2cos 2pi k/7)`, the reciprocal of `q7`) or `t^3-3t-1=0` (`t=-2cos 2pi k/9`, the negation of
`q9`), and `Phi_3(x,-2)=(x^3-3x+1)(x^3+x^2-2x-1)` (geometry (6)). Rational periods `4,5,6` do not occur (Morton
1998; Flynn-Poonen-Schaefer 1997; Stoll 2008 under BSD: CITED).

## 6. Complete rational preperiodic graphs and the microcosm test

**PROVED algorithm** (odd `p`: `v_p(f(x))=2v_p(x)` when `2v_p(x) != v_p(c)`, so the denominator of `x` is `D` with
`den(c)=D^2`; real escape `|x|>(1+sqrt(1-4c))/2`; finite exact enumeration with revisit detection), the
THM-4146 section 1 argument for general rational `c`. Results (edges in the output):

| `c` | `PrePer(f_c,Q)` | structure |
|---|---|---|
| `0` | `{-1,0,1}` | `0` fixed; `1` fixed with `-1->1` (user's "x^2+0") |
| `-1` | `{-1,0,1}` | `0<->-1`, `1->0` (user's "x^2+1" is `x^2-1`) |
| `-2` | `{-2,...,2}` | `2,-1` fixed; `1->-1`; `0->-2->2` (user's "x^2+2" is `x^2-2`; geometry sec. 5) |
| `-3/4` | `{+-1/2,+-3/2}` | fixed `-1/2, 3/2`, each with one extra preimage |
| `-7/4` | `{+-1/2,+-3/2}` | 2-cycle `1/2<->-3/2`, `-1/2->-3/2`, `3/2->1/2` (same set as `-3/4`, non-isomorphic graph) |
| `-21/16` | `m/4, |m|<=7` | fixed `7/4,-3/4`; 2-cycle `{1/4,-5/4}`; `-7/4->7/4` |
| `-29/16` | `m/4, |m|<=7` | THM-4139 (3) |
| `-77/16` | `{+-7/4,+-11/4}` | fixed `-7/4, 11/4` |

Isomorphism classes among the twelve listed parameters: `{3/16,-3/4,-77/16}` and `{-13/16,-37/16}` coincide, all
others are distinct. `0,-1,-2` are not preperiodic for `-7/4`, nor is `-7/4` for itself (`v_2=-2 != -1`).

**Complete parameter lists for a given point (PROVED bounds `den(c)=den(x)^2`, `c<=min(1/4,|x|-x^2)`,
`c>=-(1+x^2)-sqrt(1+x^2)`, the last from `x^2+c>=-beta(c)`):**

| point `x` | parameters `c` with `x` preperiodic (tail, period) |
|---|---|
| `1/4` | `-29/16 (1,3)`, `-21/16 (0,2)`, `-13/16 (1,2)`, `-5/16 (1,1)`, `3/16 (0,1)` |
| `-3/4` | `-45/16 (2,1)`, `-37/16 (1,2)`, `-29/16 (2,3)`, `-21/16 (0,1)`, `-13/16 (0,2)`, `-5/16 (2,1)`, `3/16 (1,1)` |
| `-7/4` | `-93/16 (1,2)`, `-77/16 (0,1)`, `-37/16 (0,2)`, `-29/16 (0,3)`, `-21/16 (1,1)` |
| `-29/16` | `-1561/256 (1,2)`, `-1305/256 (0,1)`, `-633/256 (0,2)`, `-377/256 (1,1)` |

So the user's microcosm holds in exactly this form: the period-`n` parabolic parameter (`n=1,2,3`) is a periodic
point of exact period `1` and `2` of later parameters (`1/4`: `3/16,-21/16`; `-3/4`: `-21/16,-13/16`), and `-7/4`
additionally of exact period `3` (`-77/16,-37/16,-29/16`, each unique by the rational roots of `Phi_n(-7/4,c)`,
`n<=4`, with no rational period `4`); `-29/16` continues only to periods `1,2` (`-1305/256,-633/256`), never `3`.
The `-93/16` entry (tail `1` into the 2-cycle `{-11/4,7/4}`) is what a hand list stopping at `-77/16` misses.

## 7. Chebyshev, Gaussian squaring, and the typed analogies

`x=y+1/y` conjugates `y->y^2` to `x->x^2-2`; `y=2 -> 5/2`, whose orbit numerators are the Fermat numbers
`5, 17, 257, 65537, 4294967297` (pairwise coprime: every term has a primitive prime). Period-`n` points of `x^2-2`
are `2cos(2pi k/N)` with `N | 2^n-+1`; the exact-period counts are `2,2,6,12,30,54` for `n=1..6`, and for `n=3` the
two conductors are `7=2^3-1` and `9=2^3+1`, so `63=7*9` is the product of the conductors of the two 3-cycles of
`x^2-2`. Bang's `n=6` exception iff `2^3+1=3^2` has no new prime iff the `2^n+1` side of period 3 has no prime
conductor iff no prime of order `6` (section 1).

**Typed analogy A (Bang vs critical orbit).** Source: `2^n-1` at `n=6`, primitive part `Phi_6(2)=3`. Target: `N_n`
at `n=3`, primitive part `F(-7,4)=1`. Map on objects: **none** (multiplicative orbit vs polynomial critical orbit).
Shared predicate: "no prime sees the base point with exact period `n` mod `p`" (both PROVED: `p` primitive for
`2^n-1` iff `ord_p(2)=n`; `p` primitive for `N_n` iff `0` has exact period `n` under `x^2+c` mod `p`). Preserved:
divisor-indexed factorization `2^n-1=prod Phi_d(2)` vs `N_n=prod H_d(a,b)`, primitive part = the `d=n` factor.
Lost: Bang's uniform lemma "non-primitive primes of `Phi_n` divide `n`" (on the quadratic side the factors are
pairwise coprime for `n<=6`, but there is no uniform theorem here, and `n=2` has infinitely many exceptions).
Sidecar: the unit equation in `Z[rho]` at `n=3`. Decisive test: `3` versus `1`. **Verdict on "`63 <-> -7/4`":
REFUTED as a structural identification.** The only exact residues are `f^3(0)/f(0)=1/(63+1)` with exponent
`2^3-2` (`=2n` only for `n=3`), `63=7*9` attached to `c=-2` (the `s=1/2` point of the parabola centred at `-7/4`),
and the shared field `Q(2cos 2pi/7)`.

**Typed analogy B (Gaussian squaring, inherited).** Source: `z=(A+iB)/C -> z^2` on primitive Pythagorean points
(geometry (13)-(16), THM-3341 section 4, THM-3333). Target: `x=2A/C -> x^2-2`. Map: `x=z+z^-1` (the same Chebyshev
conjugacy, now on the unit circle). Preserved: denominators `C^(2^n)` and the pairwise coprime odd legs, so every
step of such an `x^2-2` orbit has a primitive numerator prime, with no exception at any depth (geometry (14)).
Lost: for `c != -2` there is no power-map conjugacy, which is exactly why exceptions like `-7/4` at `n=3` can
exist. Test: the Fermat orbit of `5/2` is the `|y|=2` instance, also exception-free. **SCOPE:** no map found from
THM-3341's Pell-hypotenuse selector to the critical-orbit exception.

## Reproduction

```text
python3 04-computation/experiments/collatz_mod6_20260917_zsigmondy_triad.py > 05-knowledge/results/collatz_mod6_20260917_zsigmondy_triad.out
python3 -O 04-computation/experiments/collatz_mod6_20260917_zsigmondy_triad.py   # identical stream except the timing line
```

Runtime about 6 s (the output's last line reports it); every check uses explicit `raise`. sha256 script
`e608d2de7d032b48fa8bc6372b4ec0d4ef006b853303989b94d6245d9219a3f5`, output
`b063b975b58ea05aaf69689297a9dd8714383bc0769ab6a2328c76425adcee4b` (the timing line is the only run-to-run
variation). PARI/GP (`gp` on this machine) supplies the
CITED completeness of the `F`, `H_4`, `H_5` Thue lists and the `nfdisc` cross-check; if `gp` is absent the script
prints that those items are OPEN and the remaining checks still run. Recovery note: the script recovered from the
transcript predated the agent's final fixes (float-vs-`Rational` equality in the family checks under SymPy 1.14, a
`\q` that broke the PARI call, a four-element `-7/4` incidence list missing `-93/16`, and a ledger count `14`); each
was repaired to the independently recomputed truth, never weakened.

## Stopping boundary / next question

Closed over `Q`: the critical-orbit primitive-divisor question at `n=3,4,5` (Thue, certified), and every incidence
list above. The next question with real content is the degree-27 Thue equation `H_6(a,b)=+-1` (a longer PARI run,
or a Baker-type bound with the explicit real roots), and whether `Res(H_d,H_n)=+-1` holds for all `d<n`; a proof
would give "non-Zsigmondy at `n` iff `H_n(a,b)=+-1`" uniformly and reduce Krieger's bound to a Thue-equation
family. The conductor law along the parabola (is the cycle-field conductor the `3`-adjusted squarefree kernel of
`4u^2+2uv+7v^2` for `s=u/v`? examples `7,9,31` fit, `s=-1/4` gives `27/4` and a split cubic) is untested beyond
these five points. No transfer to the Collatz rows beyond the mod-`9`/Catalan link of section 1 is claimed; the
row lane and the quadratic lane share the integer `9` and nothing else proved here.
