---
id: THM-2005
title: SUPPORT-DIRICHLET PROFILE EXTENSIONS AND THE AUTOMATIC/TOURNAMENT RECIPROCAL ATLAS
status: PROVED.  The full support Dirichlet profile and its abscissa, integer Abel--Dini lift, Egyptian conservation hyperplane, ordinary and centered polygonal digamma clocks, maximum-cyclic-triangle mass and condensation-hazard profile, Forcade reciprocal mass, certified tournament-census tail, fibbinary/Moser Mahler block laws, primitive-residue profile, and Sylvester remainder are proved.  The single-file referee passes normally and under Python optimization.  Numerical comparisons use mpmath; unnamed constants receive no unproved arithmetic classification
source: codex-2026-07-21 reciprocal-sequence continuation and union audit
depends_on: [THM-2000, THM-785, THM-819, THM-1360, THM-2016, HYP-3724, HYP-3008, HYP-3063]
related: [THM-488, THM-555, THM-874, THM-900, THM-1127, THM-1985, THM-1990, THM-2016, MISTAKE-209, MISTAKE-210, MISTAKE-217, MISTAKE-219, MISTAKE-220]
external:
  - "Lawrence Downey, Boon W. Ong, and James A. Sellers, Beyond the Basel Problem: Sums of Reciprocals of Figurate Numbers, College Mathematics Journal 39 (2008), 391--394, JSTOR stable 27646686, https://www.jstor.org/stable/27646686"
  - "Archived preprint: https://web.archive.org/web/20130529032918/http://www.math.psu.edu/sellersj/downey_ong_sellers_cmj_preprint.pdf"
script: 04-computation/support_dirichlet_automatic_tournament_atlas_thm2005.py
output: 05-knowledge/results/support_dirichlet_automatic_tournament_atlas_thm2005.out
formalization: 04-computation/lean/TournamentH7/TournamentH7/SupportHarmonicFigurate.lean
script_sha256: 9fa4f4d63e5b31c403cbf6cf9d3576d9b86d42a4468cd7a07966a42dcac80dee
output_sha256: ff95605d3ee48c8fdbc6a8948a033b6cb334c35a60a6db5422b8c494936e2472
formalization_sha256: 7a2d36331e924465887f5a19bebb2561fb3d5cb9d7b6d3bc69eb68166950b416
---

# THM-2005 -- support-Dirichlet profile extensions and exact sequence atlas

THM-2000 establishes the support/multiplicity correction, Abel--Stieltjes
criterion, master figurate surface, Faulhaber rows, Gauss theta product, and
reciprocal tournament-series reversal.  This companion identifies the larger
analytic object behind the scalar mass and closes several sequence families
which THM-2000 leaves at atlas level.

## 1. The underlying object is the support Dirichlet profile

Let `a=(a_n)` be a positive integer sequence, let
`nu(m)=#{n:a_n=m}`, and let `A={m:nu(m)>0}`.  THM-2000's collision law is
best written in extended nonnegative reals as

```text
sum_n 1/a_n
 =sum_(m in A)1/m+sum_m (nu(m)-1)_+/m.                (1)
```

Writing it as an equality, rather than a subtraction, avoids the undefined
form `infinity-infinity`.

The information-preserving transforms are

```text
D_A(z)=sum_(m in A)m^(-z),
G_A(w)=sum_(m in A)w^m.                              (2)
```

For real `z>0`, with equality allowed in `[0,infinity]`, Tonelli and the
gamma integral give

```text
D_A(z)=1/Gamma(z) integral_0^infinity
        u^(z-1)G_A(exp(-u)) du,
D_A(1)=integral_0^1 G_A(w)/w dw.                     (3)
```

The coefficients of either transform recover the support.  Coefficientwise
(and hence on every common absolute-convergence domain), Boolean support
operations obey the exact valuation

```text
D_(A union B)+D_(A intersect B)=D_A+D_B.             (4)
```

Thus `D_A(1)` is one lossy evaluation of a complete profile, not a sequence
fingerprint.

Let `q_1<q_2<...` enumerate an infinite support, and let
`N_A(x)=#{m in A:m<=x}`.  The number

```text
delta(A)=limsup_n log n/log q_n
        =limsup_x log N_A(x)/log x                   (5)
```

is the abscissa of convergence of the positive Dirichlet series.  If
`z>delta`, choose `alpha` strictly between them; eventually
`N_A(x)<=x^alpha`, so Abel--Stieltjes makes `D_A(z)` converge.  Conversely,
if `D_A(z)` converges, then

```text
n/q_n^z <=sum_(j<=n)q_j^(-z)=O(1),
```

so `delta<=z`.  Hence `D_A(z)` diverges for `z<delta`; either behavior can
occur on the boundary.  In particular `delta=1` is a region, not a verdict:
for `n>=2`, the support `ceil(n log n)` diverges at one, while
`ceil(n(log n)^2)` converges (finite initial repair is immaterial).

## 2. Abel--Dini is a recursive integer-support operation

Let `q_k` enumerate a divergent support and set

```text
H_k=sum_(j<=k)1/q_j.
```

Once `H_k>=1`, for any real `p>=0` define

```text
b_k=ceil(q_k H_k^p).                                 (6)
```

Then `b_k` is strictly increasing: the unrounded arguments increase by at
least `H_k^p>=1`.  Moreover `x<=ceil x<=2x` for `x>=1`, so the reciprocal
series is comparable to

```text
sum_k (1/q_k)/H_k^p.
```

The Abel--Dini theorem from THM-2000 therefore yields the exact support law

```text
sum_k1/b_k diverges for 0<=p<=1,
sum_k1/b_k converges for p>1.                        (7)
```

At `p=1`,

```text
sum_(k<=K)1/b_k=log H_K+O(1).                         (8)
```

Indeed the difference between
`(1/q_k)/H_k` and `log(H_k/H_(k-1))` is summable, as is the rounding error,
because `q_k>=k`.  Iterating the critical lift from `q_k=k` constructs
honest integer supports at the full Bertrand scales

```text
k log k log_2 k ... log_r k.                         (9)
```

So the boundary has no final member: it is a recursive operation on supports.

## 3. Egyptian splitting makes z=1 a conservation hyperplane

For `n>=2`, a collision-free support refinement in which both children are
fresh replaces one denominator by

```text
n -> {n+1,n(n+1)}.
```

Its Dirichlet-profile change is

```text
Delta_n(z)=n^(-z)[(n/(n+1))^z+(1/(n+1))^z-1].       (10)
```

The two bases in brackets sum to one.  Strict concavity/convexity of `x^z`
therefore gives

```text
Delta_n(z)>0  for 0<z<1,
Delta_n(1)=0,
Delta_n(z)<0  for z>1.                               (11)
```

Thus harmonic mass is not accidentally preserved: `z=1` is exactly the
conservation hyperplane of the Egyptian refinement flow.

Start with `{1,2,4,8,...}` and independently split any subcollection of
`2^j`, `j>=1`.  The new odd values `2^j+1` are distinct; the even values
`2^j(2^j+1)` have distinct 2-adic valuations; neither collides with an
unsplit power of two.  This gives continuum many supports with

```text
D_A(1)=2,                    N_A(x)=Theta(log x),
delta(A)=0.                                         (12)
```

The infinite triangular support and finite `{1,2,3,6}` also have mass two.
Even the pair `(D_A(1),delta(A))` is therefore non-rigid; the full profile
separates the supports.

## 4. Two polygonal digamma clocks and two confluences

### 4.1 Ordinary polygonal numbers

For side count `h>=3`,

```text
P_h(n)=[(h-2)n^2-(h-4)n]/2.
```

For `h!=4`, the finite harmonic-clock identity is

```text
1/P_h(n)=2/(h-4)
 [1/(n-1+2/(h-2))-1/n],                             (13)

sum_(n=1)^N1/P_h(n)
 =2/(h-4)[sum_(r=0)^(N-1)1/(r+2/(h-2))-H_N].        (14)
```

Consequently

```text
sum_(n>=1)1/P_h(n)
 =2/(h-4)[psi(1)-psi(2/(h-2))].                     (15)
```

This extends the even `h>=6` rows of Downey--Ong--Sellers to every odd
polygon without a separate argument.  At `h=4` the two clocks coalesce and
the divided difference becomes

```text
psi_1(1)=zeta(2)=pi^2/6.                             (16)
```

For comparison, if `u_n=n(alpha n+beta)/c`, where `alpha>0`, `c>0`, and
`alpha+beta>0`, then

```text
sum_(n>=1)1/u_n
 =c/beta[psi(1+beta/alpha)+gamma]   if beta!=0,
 =c/alpha zeta(2)                   if beta=0.       (17)
```

Distinct linear factors give digamma differences; repeated factors give
polygamma values.

### 4.2 Centered polygonal numbers

The centered `k`-gonal row from THM-1360 is

```text
C_k(n)=1+k n(n-1)/2,                 k>=3,n>=1.      (18)
```

Put

```text
Delta=sqrt(1-8/k),            r_+-r_-=Delta,
r_+/-=(1+/-Delta)/2.                                  (19)
```

Partial fractions and the same digamma clock give, for `k!=8`,

```text
sum_(n>=1)1/C_k(n)
 =2/(k Delta)[psi(1-r_-)-psi(1-r_+)].                (20)
```

For `3<=k<8`, the roots are conjugate and (20), on the principal branch, is
real.  At `k=8` they coalesce and

```text
C_8(n)=(2n-1)^2,
sum_(n>=1)1/C_8(n)=pi^2/8.                           (21)
```

At the first split real-root row,

```text
sum_(n>=1)1/C_9(n)=2pi/(3sqrt 3).                    (22)
```

The `n=1` summand stays one; every summand with `n>=2` decreases strictly with
`k`, so the masses decrease strictly to one.  The ordinary square at `h=4`
and centered octagonal row at `k=8` are parallel
digamma-to-trigamma confluences.

## 5. Exact tournament reciprocal rows

### 5.1 Maximum cyclic triples

For a tournament score sequence `d_1,...,d_n`,

```text
c_3(T)=C(n,3)-sum_i C(d_i,2).                        (23)
```

Convexity at fixed `sum_i d_i=C(n,2)` minimizes the last sum at regular
scores when `n` is odd and almost-regular scores when `n` is even.  These
score rows exist (for even `n`, delete one vertex from a regular tournament
of order `n+1`).  Hence

```text
M_3(n)=max_T c_3(T)
 =n(n^2-1)/24,                 n odd,
 =n(n^2-4)/24,                 n even.              (24)
```

The former claim `M_3(n)=C(n,3)` confused all triple slots with cyclic
triples.  Splitting parity and telescoping gives

```text
sum_(m>=1)1/M_3(2m+1)=18-24 log 2,
sum_(m>=2)1/M_3(2m)=3/4,                             (25)

sum_(n>=3)1/M_3(n)=75/4-24 log 2
                  =2.1144676665613125... .           (26)
```

The odd denominator is `sum_(j<=m)j^2`, so (25) is also the second
Faulhaber reciprocal mass in THM-2000.

The concurrent reducibility ceiling of THM-2016 turns the same row into a
temperature product.  With

```text
R_n=max{c_3(T):T reducible on n vertices}=M_3(n-1),
tau_c(n)=R_n/M_3(n),                                  n>=4,
```

one has the exact telescoping identities

```text
product_(j=4)^N tau_c(j)=1/M_3(N),
sum_(N>=3) product_(j=4)^N tau_c(j)
 =sum_(n>=4)1/R_n
 =75/4-24 log 2.                                    (26a)
```

The `N=3` product is empty.  The ratio itself exposes a second, exactly
harmonic support.  Formula (24) simplifies it to

```text
tau_c(n)=(n-1)/(n+2),     n even,
tau_c(n)=(n-3)/n,         n odd,

q_n:=3/(1-tau_c(n))
   =n+2,                  n even,
   =n,                    n odd.                    (26b)
```

Thus `(q_n)_(n>=4)=(6,5,8,7,10,9,...)`, the adjacent-pair reversal of the
cofinite harmonic support `{5,6,7,...}`.  If
`H_4^(z):=sum_(m=1)^4 m^(-z)`, then for `Re z>1`,

```text
sum_(n>=4)((1-tau_c(n))/3)^z
 =sum_(m>=5)m^(-z)=zeta(z)-H_4^(z).                 (26c)
```

This profile has abscissa exactly one.  Its boundary divergence is also
termwise explicit:

```text
(1/3)sum_(n=4)^N(1-tau_c(n))
 =H_(N+1)-H_4,                         N odd,
 =H_N-H_4+1/(N+2),                     N even,

(1/3)sum_(n=4)^N(1-tau_c(n))-log N -> gamma-H_4.    (26d)
```

There is an order-sensitive refinement.  Put `a_m=1-3/m`.  Sorting the same
hazard support gives

```text
S_4=1,
S_M=product_(m=5)^M a_m=24/[(M-2)(M-1)M],
sum_(M>=4)S_M=2.
```

The tournament order instead has prefix products
`P_3=1`, `P_N=product_(n=4)^N a_(q_n)=1/M_3(N)`.  Swapping every adjacent
pair changes only the intermediate prefix, so

```text
sum_(N>=3)P_N-sum_(M>=4)S_M
 =67/4-24 log 2
 =sum_(k>=2)72/[(2k-2)(2k-1)(2k)(2k+1)(2k+2)]>0.   (26e)
```

Thus the reciprocal maximum-`c3` mass is also the partition sum of cumulative
condensation ratios and the reciprocal mass of the reducibility ceilings.
The additive hazard profile depends only on its denominator support, while
the prefix-product partition function remembers the parity-shuffled word.

### 5.2 Forcade orders and the Mersenne clock

Let `E=sum_(p>=1)1/(2^p-1)` be the Erdos--Borwein constant.  The arc count at
the Forcade order `2^p` is

```text
d_p=C(2^p,2)=2^(p-1)(2^p-1).
```

Termwise,

```text
1/d_p=2[1/(2^p-1)-1/2^p],
```

and therefore the exact tournament-observable mass is

```text
sum_(p>=1)1/C(2^p,2)=2E-2
                       =1.2133903048305835... .       (27)
```

No claim about the subrow indexed by Mersenne primes is made; even its
infinitude is open.

The entire Forcade support profile has the Lambert expansion

```text
sum_(p>=1)d_p^(-z)
 =2^z sum_(r>=0) (z)_r/[r!(2^(2z+r)-1)],     Re z>0.       (27a)
```

Indeed, expand `(1-2^(-p))^(-z)` and interchange two absolutely convergent
series.  Formula (27) is its `z=1` slice.

### 5.3 A certified unlabeled-tournament census tail

Let `t_n` be the number of tournament isomorphism classes on `n` vertices.
Every orbit has size at most `n!`, so

```text
t_n>=2^C(n,2)/n!,
1/t_n<=a_n:=n!/2^C(n,2).                            (28)
```

Since

```text
a_(n+1)/a_n=(n+1)/2^n
```

decreases and is at most `1/2` for `n>=3`, every `N>=2` has the tail
certificate

```text
sum_(n>N)1/t_n
 <= a_(N+1)/[1-(N+2)/2^(N+1)]
 <=2(N+1)!/2^C(N+1,2).                              (29)
```

The same upper bound applies after support deduplication.  Thus any correctly
indexed A000568 prefix can be turned into a rigorous interval; its arithmetic
nature remains unnamed.

More generally, if `sigma=Re z>0`, then

```text
sum_(n>N)|t_n^(-z)|
 <= [(N+1)!/2^C(N+1,2)]^sigma
    /[1-((N+2)/2^(N+1))^sigma],             N>=2.   (29a)
```

For the canonical indexing `n=0,...,20`, Burnside's all-odd-cycle formula
gives

```text
(t_0,...,t_20)=
(1,1,1,2,4,12,56,456,6880,191536,9733056,903753248,
 154108311168,48542114686912,28401423719122304,
 31021002160355166848,63530415842308265100288,
 244912778438520759443245824,1783398846284777975419600287232,
 24605641171260376770598003978281472,
 645022068557873570931850526424042500096).
```

After deduplicating the repeated value one, the exact prefix is

```text
P_20=
21099328871173442978479709817904691186237927028914957892229218402318110941302184270517423289978542744036679862436756474093271800482440798353111262783654453239105958659179983400460013348652791077311743013830527
/11383296645907658352445029044902856765681508434252564284448144572805470654577744374492610675249889401912534749131567170501822811250421465546423774298097496376055720127690081166416741747792702633421074615848960,

0 < D_support(1)-P_20
  <=5568470782875
    /179343882456254548076397338228109815750132052193353662464.
```

Consequently

```text
1.85353413229011631733371526582380005497181657702917
 < D_support(1)
 < 1.85353413229011631733371526582380005497181660807831.  (29b)
```

The certified interval has width below `3.11*10^(-44)`.

The indexed sequence sum has the same interval shifted upward by two, exactly
the collision tax from the three initial ones.  To see that there are no later
collisions, add a unique dominating vertex to inject order-`n` classes into
order-`n+1` classes; for every order at least three a tournament with no
dominating vertex lies outside the image.  Thus `t_(n+1)>t_n` for `n>=2`.
The census support has abscissa zero, although its harmonic mass is not
presently recognized as a classical constant.

## 6. Automatic supports have exact block tails

Let `F` be the nonnegative fibbinary integers (binary words with no adjacent
ones) and let `M` be the Moser--de Bruijn integers (base-four digits zero or
one).  Their generating functions, including the zero term, satisfy

```text
G_F(w)=G_F(w^2)+wG_F(w^4),
G_M(w)=(1+w)G_M(w^4)=product_(j>=0)(1+w^(4^j)).      (30)
```

The first identity separates words ending in `0` from those ending in `01`;
the second chooses the lowest base-four digit.  In particular `M subset F`.

With `F_0=0,F_1=1` and `phi=(1+sqrt 5)/2`,

```text
#(F intersect [0,2^K))=F_(K+2),
#(F intersect [2^k,2^(k+1)))=F_(k+1),

#(M intersect [0,4^R))=2^R,
#(M intersect [4^r,4^(r+1)))=2^r.                  (31)
```

Therefore

```text
delta(F\{0})=log_2(phi),              delta(M\{0})=1/2.             (32)
```

Both harmonic masses converge, with exact block-tail certificates

```text
F_(K+3)/2^K
 <=sum_(m in F,m>=2^K)1/m
 <=F_(K+3)/2^(K-1),                                  (33)

3/2^(R+1)
 < sum_(m in M,m>=4^R)1/m
 <=2^(1-R).                                          (34)
```

For (33), sum the Fibonacci block counts with
`sum_(j>=0)F_(j+1)x^j=1/(1-x-x^2)`.  For (34), the `r`th block has `2^r`
values between `4^r` and `(4^(r+1)-1)/3`.

Writing `F_+=F\{0}`, `M_+=M\{0}`, their full profiles also obey

```text
(1-2^(-z))D_(F_+)(z)=sum_(m in F)(4m+1)^(-z),
(1-4^(-z))D_(M_+)(z)=sum_(m in M)(4m+1)^(-z).       (35)
```

The first analytic identity is asserted for `Re z>log_2(phi)` and the second
for `Re z>1/2`; coefficientwise they are formal identities, and on the
positive real axis below those thresholds both sides are understood only in
the extended nonnegative sense.  These equations are the faithful shadows of
the Fibonacci-cube and even-coordinate Boolean-cube carriers in
HYP-3008/HYP-3063.

## 7. Primitive residues and the Sylvester exact tower

For an LRC modulus `q>=2`, let the finite residue support

```text
A_q={1<=u<q:gcd(u,q)=1},
H_m^(z)=sum_(j=1)^m j^(-z).
```

Möbius inversion of the coprimality indicator gives the exact finite profile

```text
D_(A_q)(z)=sum_(d|q)mu(d)d^(-z)H_(q/d-1)^(z).        (36)
```

At `z=1`, THM-819's interval-core good measure is

```text
2/[q(q+1)] D_(A_q)(1).                              (37)
```

For prime `q`, this becomes `2H_(q-1)/[q(q+1)]`.  Primitive residues and
their inverse pairs, rather than runners, are therefore exact alternate
vertices for this LRC harmonic law.

At the opposite, doubly-exponential scale, let

```text
s_0=2,                 s_(j+1)=s_j^2-s_j+1.
```

Then

```text
2,3,7,43,1807,...,
sum_(j=0)^N1/s_j=1-1/(s_(N+1)-1),
sum_(j>=0)1/s_j=1.                                  (38)
```

The exact rational mass shows again that rapid growth predicts convergence,
not arithmetic complexity.

## 8. Tournament Analysis and challenged vertices

On a finite family of distinct finite-mass supports, ordering first by
`D_A(1)` and then lexicographically by the least element of a symmetric
difference creates a transitive tournament.  Its score multiset is
`{0,1,...,r-1}`, every SCC is a singleton, it has no directed cycle, and it
has one Hamiltonian path.  This is only mass-ranking telemetry.

We challenged sequence indices, values, multiplicity collisions, dyadic
blocks, automaton states, primitive residues, inverse pairs, Egyptian moves,
polygonal root clocks, and proof obligations as vertices.  Bounded-ratio
occupancy blocks with the Dirichlet profile as sidecar preserve the Abel
predicate.  Primitive-residue/inverse-pair vertices additionally preserve the
THM-819 LRC measure.  Scalar mass destroys both distinctions.

## 9. Verification and formalization boundary

The single-file referee contains no Python `assert` nodes and passes
byte-identically under ordinary and optimized Python.  It checks the finite
Abel and valuation laws, ordinary and centered polygonal splits, master
figurate factorizations, integer Abel--Dini lifts, automatic-language
recurrences and block counts, Egyptian profile orientation, Sylvester
remainders, primitive-residue Möbius profiles, Forcade identities,
tournament-tail ratios, the finite discrete-convexity/superadditivity core of
the repaired THM-2016 ceiling, its condensation product, harmonic-hazard
permutation and finite profiles, parity-shuffle tax, the Gauss product, and
both maximum-`c_3` parity sums.

The sorry-free Lean module `TournamentH7/SupportHarmonicFigurate.lean`
formalizes the master and ordinary-polygonal factorizations, their finite
reciprocal decompositions, the odd/even maximum-`c3` denominator and partial
fraction algebra, both condensation ratios and harmonic defects, and the
block sandwich.  Its isolated build passes, and the printed axiom audits
contain no `sorryAx`.  Infinite Abel/Mellin integrals,
digamma evaluation, automatic-language enumeration, the score-theoretic
maximum-`c3` argument, and the infinite parity sums remain paper proofs plus
independent exact referees; none is postulated as a Lean axiom.

This theorem deliberately makes no transcendence claim for census or
automatic masses and no completeness claim for the global Hamiltonian-path
value set.  ∎
