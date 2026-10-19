# Smooth arithmetic operations at the quadratic cusp

19 October 2026. **PROVED from the cited cusp-phase and sieve inputs.**
This is an elementary transfer theorem, not a prime-distribution theorem.
An [independent audit](review/independent_cusp_operations_review.md)
returned **PASS**. No literature-priority claim is made.

## 1. Result and source inputs

Fix a real escape radius R>1/2, and use the strict escape convention

\[
x_0=0,\qquad x_{k+1}=x_k^2+\tfrac14+\epsilon,\qquad
N(\epsilon)=\min\{k\geq0:x_k>R\}.
\tag{1}
\]

There are explicit real-analytic maps of the small positive parameter
epsilon that implement integer translation of N, and a real-analytic
involution that implements a prescribed complementary sum. The translation
maps are universal, independent of R; the reflection uses one convergent
endpoint constant C_R. They do not query an escape-band table or invert
the discontinuous integer-valued function N. The identities
for N hold outside a specified thin guard set. The maps' composition laws
hold exactly throughout their common domains, without any guard.

On a dyadic phase window of length M, the discarded proportion is
O_R(log M/M^2). Every bounded two-integer statistic then has an unguarded
integral-to-sum error O_R(||F||_infinity log M/M). A guarded, equally
weighted version is an **exact** sum identity, including exact positivity
equivalence. These assertions are uniform in the statistic, with the
displayed dependence on its supremum norm.

Combining this transfer with the prior every-target four-form rough-sieve
theorem gives, for every sufficiently large integer m, a set of parameters
of ordinary epsilon-measure at least c_R/(m^2(log m)^4) on which

\[
N(\epsilon)+N(\mathcal R_m(\epsilon))=m
\]

and each of the four integers

\[
6N(\epsilon)\pm1,\qquad
6N(\mathcal R_m(\epsilon))\pm1
\]

has at most 37 prime factors counted with multiplicity. The prime versions
remain unestimated arithmetic integrals; no twin-prime or Goldbach
conjecture is proved.

The two imported inputs are:

1. **Cusp phase.** The companion
   `parabolic_escape_precision_20261019.md`, Sections 1 and 4, defines
   C_R by convergent unperturbed endpoint sums and proves the strict
   escape-band convention. Its extension
   [next-order cusp phase](../parabolic_escape_next_order_20261019/parabolic_escape_next_order_20261019.md),
   Sections 1 and 6, proves the quantitative phase used below. The former packet is at
   [parabolic escape precision](../parabolic_escape_precision_20261019/parabolic_escape_precision_20261019.md).
2. **Rough twin-center sums.**
   [Every-target rough twin-center sums at power sieve depth](../prime_pairs_local_symmetry_20261019/power-depth-rough-pairs.md),
   Theorem 1 with sieve parameter 38. Its classical combinatorial-sieve
   input and target-uniform hypotheses are recorded and proved there.

All subsequent transfer arguments are included here. The convergence of
C_R, the next-order phase estimate, and the lower sieve are inherited
inputs rather than claims that the following elementary coordinate change
proves them anew.

## 2. The large-branch phase coordinate

Let epsilon_n be the unique positive solution of x_n(epsilon_n)=R and put

\[
t_n=\pi/\sqrt{\epsilon_n},\qquad A=\pi^2/3,
\qquad \tau(t)=t+C_R+A/t.
\]

The imported phase theorem states, with constants depending only on R,

\[
s_n:=\tau(t_n)=n+O_R(\log n/n^2).
\tag{2}
\]

Because escape is strict, the actual staircase is

\[
N(\pi^2/t^2)=n\quad\Longleftrightarrow\quad
t_{n-1}\leq t<t_n\qquad(n\geq2).
\tag{3}
\]

In particular its value at t_n is n+1, not n.

Choose once and for all T>sqrt(A), depending only on R, with
u_0:=tau(T)>1. Restrict to t>T, or equivalently to

\[
\mathcal D=(0,\pi^2/T^2).
\]

Here tau is strictly increasing, since tau'(t)=1-A/t^2>0. Its inverse is
the **large** quadratic branch

\[
t(u)=\frac{u-C_R+\sqrt{(u-C_R)^2-4A}}2,\qquad u>u_0.
\tag{4}
\]

Specifying the branch is essential: the quadratic equation generally has
two positive roots, and the smaller one does not lie in t>T. Define

\[
\Xi(u)=\frac{\pi^2}{t(u)^2},\qquad
U(\epsilon)=\tau(\pi/\sqrt\epsilon).
\tag{5}
\]

Then Xi:(u_0,infinity)->D is a decreasing real-analytic diffeomorphism,
and U is its inverse. Once n is sufficiently large, (3) becomes

\[
N(\Xi(u))=n\quad\Longleftrightarrow\quad
s_{n-1}\leq u<s_n.
\tag{6}
\]

### Lemma 1. One guard works uniformly above a scale

There exist K_R>0 and L_R>0 such that, for every real L>=L_R, with

\[
\eta_L=K_R\frac{\log L}{L^2}<\frac14,
\tag{7}
\]

every real u>=L/4 satisfying dist(u,Z)>eta_L lies in the domain of Xi and
obeys

\[
N(\Xi(u))=\lceil u\rceil.
\tag{8}
\]

The constants K_R and L_R are existential here. This packet neither
numerically calibrates them nor certifies a claimed finite onset.

**Proof.** By (2), choose B_R and n_0 such that
|s_j-j|<=B_R log j/j^2 for j>=n_0. Increase n_0 so that t_j>T there.
The function log x/x^2 decreases for x>=2. If n=ceil u and L is sufficiently
large, then n-1>=L/8>=n_0. Hence both boundary errors in (6) are at most

\[
B_R\frac{\log(L/8)}{(L/8)^2}
\leq64B_R\frac{\log L}{L^2}.
\]

Take, for example, K_R=128B_R+1 and then increase L_R so that eta_L<1/4
and L/4>u_0. The guard excludes integers. It gives

\[
s_{n-1}\leq n-1+\eta_L<u<n-\eta_L\leq s_n.
\]

Now (6) proves (8). No differentiability of the remainder in (2) has been
used. QED.

## 3. Analytic translations and complementary-sum involutions

For an integer q define

\[
\mathcal F_q(\epsilon)=\Xi(U(\epsilon)+q)
\tag{9}
\]

on the precise domain

\[
\mathcal D_q=\{\epsilon\in\mathcal D:U(\epsilon)+q>u_0\}.
\]

It is an increasing real-analytic diffeomorphism from D_q onto D_{-q},
with inverse F_{-q}. These maps form a **partial** integer action on D:

\[
\mathcal F_0=\mathrm{id},\qquad
\mathcal F_p\circ\mathcal F_q=\mathcal F_{p+q}
\tag{10}
\]

whenever u, u+q, and u+p+q all exceed u_0. A full action on the whole cusp
interval is not claimed, because negative translations can leave it.

In fact, the constant C_R **cancels from the translation map itself**.
Writing H(t)=t+A/t and taking its large-branch inverse gives the equivalent
formula

\[
\mathcal F_q(\epsilon)
=\frac{\pi^2}
 {\left[H_+^{-1}\!\left(H(\pi/\sqrt\epsilon)+q\right)\right]^2}.
\tag{10a}
\]

Thus these maps are independent of the escape radius on common cusp
domains. Extending q to a real parameter s, the same formula is the local
flow, for as long as it remains in D, of the explicit parameter equation

\[
\frac{d\epsilon}{ds}
=-\frac{2\epsilon^{3/2}}{\pi(1-\epsilon/3)}.
\tag{10b}
\]

This follows by differentiating U(epsilon(s))=U(epsilon(0))+s. The endpoint
constant locates the integer guard cells; it is not needed to evaluate
F_q. This is a flow of **parameter transformations**, not an analytic
embedding of the original recurrence x->x^2+1/4+epsilon.

For an integer target m with m-1>2u_0, define

\[
\mathcal R_m(\epsilon)=\Xi(m-1-U(\epsilon))
\tag{11}
\]

on

\[
\mathcal H_m=\{\epsilon\in\mathcal D:
u_0<U(\epsilon)<m-1-u_0\}.
\]

This is a decreasing real-analytic diffeomorphism of H_m onto itself.
The following are exact coordinate identities:

\[
\mathcal R_m^2=\mathrm{id},\qquad
\mathcal R_m\mathcal F_q=\mathcal F_{-q}\mathcal R_m,
\qquad
\mathcal R_m\mathcal R_\ell=\mathcal F_{m-\ell},
\tag{12}
\]

each on the common domain where every indicated intermediate point lies
in D and every indicated reflection is defined. They follow by applying
Xi to u+q and m-1-u. Thus the composition laws do not depend on a rounding
estimate, and do not assert that the original quadratic recurrence has a
global analytic flow.

In H coordinates the reflection reads H(t')=m-1-2C_R-H(t), so C_R does not
cancel from this operation.

### Theorem 2. Exact integer operations away from the guard

There is M_R such that the following hold.

**Translations.** Let M>=M_R be an integer, q an integer with q>=-M/2,
and u in [M,2M] with dist(u,Z)>eta_M. Then

\[
N(\mathcal F_q(\Xi(u)))=N(\Xi(u))+q
=\lceil u\rceil+q.
\tag{13}
\]

The constants do not depend on q; in particular they cover every fixed
integer q once M is sufficiently large. There is no upper bound on q.

**Reflection.** For an integer m>=M_R set

\[
a_m=\lfloor m/3\rfloor,\quad b_m=\lfloor2m/3\rfloor,
\quad J_m=[a_m,b_m],\quad
I_m=\{a\in\mathbb Z:m/3<a\leq2m/3\}.
\tag{14}
\]

If u in J_m and dist(u,Z)>eta_m, then Xi(u) lies in H_m and

\[
N(\Xi(u))=\lceil u\rceil\in I_m,\qquad
N(\mathcal R_m(\Xi(u)))=m-\lceil u\rceil.
\tag{15}
\]

In particular the two actual escape integers sum to m.

**Proof.** For translations, u>=M and u+q>=M/2; integer q preserves the
distance to Z. Lemma 1 therefore applies at both inputs, and
ceil(u+q)=ceil u+q.

For reflection, both u and v=m-1-u are at least m/3-1>=m/4 for all large
m, and dist(v,Z)=dist(u,Z). Write u=n-1+theta with 0<theta<1. Then

\[
m-1-u=m-n-\theta\in(m-n-1,m-n),
\]

so ceil(v)=m-n. Lemma 1 proves (15). The domain requirements follow by
increasing M_R. The full unit cells (n-1,n), n in I_m, are exactly the
noninteger parts of J_m. This accounts for all floor and endpoint
conventions in (14). QED.

At an integer u=k, the ceiling identity would instead give
ceil u+ceil(m-1-u)=m-1. Those points are deliberately excluded. Replacing
m-1-u with m-u would yield the wrong sum m+1 at every noninteger phase.

## 4. Integral identities, counting, and positivity

Write B=||F||_infinity for a bounded function F:Z^2->C. The same arguments
apply to real-valued statistics. Let

\[
G_M=\{u\in[M,2M]:\operatorname{dist}(u,\mathbb Z)>\eta_M\}.
\]

The complement has Lebesgue measure exactly 2M eta_M. On each cell
(n-1,n), n=M+1,...,2M, its good part is
(n-1+eta_M,n-eta_M), of length 1-2eta_M.

### Theorem 3. Translation statistics

For all M and q allowed in Theorem 2,

\[
\int_{G_M}F\bigl(N(\Xi(u)),N(\mathcal F_q(\Xi(u)))\bigr)\,du
=(1-2\eta_M)\sum_{n=M+1}^{2M}F(n,n+q).
\tag{16}
\]

Moreover,

\[
\left|\int_M^{2M}F\bigl(N(\Xi(u)),N(\mathcal F_q(\Xi(u)))\bigr)\,du
-\sum_{n=M+1}^{2M}F(n,n+q)\right|
\leq4BM\eta_M
=O_R\!\left(B\frac{\log M}{M}\right).
\tag{17}
\]

The bounds hold for arbitrary bounded F, including F that depends on M or
q. The implied constant depends only on R, not on F, M, or q.

**Proof.** Theorem 2 identifies the pair throughout each retained cell,
which proves (16). Outside G_M compare the actual integrand with
F(ceil u,ceil u+q); their pointwise difference has absolute value at most
2B. The latter integrand has full integral equal to the displayed sum;
integer endpoints are null sets. Multiplication by |[M,2M]\G_M|=2M eta_M
proves (17). QED.

For reflection put

\[
G_m^{\mathrm{ref}}=\{u\in J_m:
\operatorname{dist}(u,\mathbb Z)>\eta_m\},\qquad
L_m=b_m-a_m=|I_m|.
\]

The same proof gives

\[
\int_{G_m^{\mathrm{ref}}}
F\bigl(N(\Xi(u)),N(\mathcal R_m(\Xi(u)))\bigr)\,du
=(1-2\eta_m)\sum_{a\in I_m}F(a,m-a),
\tag{18}
\]

and

\[
\left|\int_{J_m}
F\bigl(N(\Xi(u)),N(\mathcal R_m(\Xi(u)))\bigr)\,du
-\sum_{a\in I_m}F(a,m-a)\right|
\leq4B L_m\eta_m
=O_R\!\left(B\frac{\log m}{m}\right).
\tag{19}
\]

In (18), F may depend on m. Thus the normalized guarded measure
du/(1-2eta) counts each retained arithmetic cell **exactly once**.

If F is nonnegative, positivity of (16) or (18) is equivalent to a
positive corresponding arithmetic sum. For a predicate, this is exactly
the existence of a satisfying pair within the stated integer range.
The unguarded estimates (17) and (19) alone do **not** give this positivity
equivalence. A tiny integral supported entirely on the bad set could be
a false positive; an unguarded integral greater than the stated error
bound would be sufficient.

Here is a concrete reason not to omit that qualification. Consider the
model staircase with thresholds s_n=n+(-1)^n rho, 0<rho<1/2, and value n
on [s_{n-1},s_n). For every even n and
u in (n-rho,n+rho), its values at u and u+1 differ by 2, although the
arithmetic translation to be modeled has difference 1. The indicator of
this wrong difference has a positive unguarded integral and a zero
arithmetic sum. This model is a rounding control, not a claim about the
specific cusp thresholds. An arbitrarily small threshold error still
allows the issue.

## 5. Ordinary parameter measure and the exact calibrated measure

The pushforward from phase to epsilon has positive density

\[
w(u):=-\Xi'(u)
=\frac{2\pi^2}{t(u)^3(1-A/t(u)^2)}
=\frac{2\pi^2}{u^3}\left(1+O_R(u^{-1})\right).
\tag{20}
\]

The derivative formula follows directly from (4)-(5). The asymptotic
uses t(u)=u-C_R+O_R(1/u); no threshold-remainder derivative is involved.
For each fixed 0<c<C, there are positive constants d_{R,c,C},D_{R,c,C}
such that for all sufficiently large M and u in [cM,CM],

\[
d_{R,c,C}M^{-3}\leq w(u)\leq D_{R,c,C}M^{-3}.
\tag{21}
\]

Consequently, in the epsilon window Xi([M,2M]):

\[
|\Xi([M,2M])|\asymp_R M^{-2},\qquad
|\Xi([M,2M]\setminus G_M)|
=O_R(\log M/M^4).
\tag{22}
\]

The discarded relative epsilon-measure is O_R(log M/M^2), the same order
as its relative phase measure. Analogous statements hold for J_m and
G_m^ref, since all their phases are comparable to m. Image intervals in
(22) are measured by length, irrespective of Xi's reversed orientation.

For any subset S of {M+1,...,2M}, the good parameter cells indexed by S
satisfy

\[
|\{\Xi(u):u\in G_M,\ \lceil u\rceil\in S\}|
\asymp_R M^{-3}(1-2\eta_M)|S|.
\tag{23}
\]

The constants are uniform in S, including the empty set in the sense of
the two non-strict inequalities. The same statement with m in place of
M holds for subsets of I_m in the reflection window.

If an exactly calibrated measure is desired, write

\[
d\nu_R(\epsilon)=|U'(\epsilon)|\,d\epsilon
=\frac\pi{2\epsilon^{3/2}}
\left(1-\frac\epsilon3\right)d\epsilon
\qquad(\epsilon\in\mathcal D).
\tag{24}
\]

The factor is positive because T>sqrt(A) implies epsilon<3. By change of
variables, dnu_R corresponds exactly to du. Thus (16) and (18) can be
read as epsilon-integral identities with this explicit density. Although
the location of each phase window depends on C_R, the calibration density
does not. Ordinary epsilon measure is not uniform counting measure:
it weights the nth good cell by an amount comparable to n^-3.

## 6. An unconditional every-target almost-prime application

For a positive integer m put z_m=(6m)^(1/38), and let S_m be the set of
a in I_m for which all four integers

\[
6a-1,\quad6a+1,\quad6(m-a)-1,\quad6(m-a)+1
\tag{25}
\]

have least prime factor at least z_m. The imported rough twin-center
theorem provides absolute c>0 and m_0 such that

\[
|S_m|\geq c\,m/(\log m)^4\qquad(m\geq m_0),
\tag{26}
\]

for **every** integer target m, with no exceptional target set.

### Theorem 4. Positive-measure rough complementary escape pairs

For each fixed R>1/2 there exist c_R>0 and m_R such that, for every integer
m>=m_R, the set

\[
\mathcal E_m=\{\Xi(u):u\in G_m^{\mathrm{ref}},\
\lceil u\rceil\in S_m\}
\tag{27}
\]

is a finite union of open intervals, contained in Xi(J_m), and

\[
|\mathcal E_m|\geq\frac{c_R}{m^2(\log m)^4}.
\tag{28}
\]

Every epsilon in E_m satisfies the exact complementary identity
N(epsilon)+N(R_m(epsilon))=m. Each of
6N(epsilon)±1 and 6N(R_m(epsilon))±1 has least prime factor at least z_m
and at most 37 prime factors counted with multiplicity.

**Proof.** Each a in S_m contributes the entire good phase cell
(a-1+eta_m,a-eta_m), of length 1-2eta_m>1/2. On it (15) identifies the
escape pair with (a,m-a), so both the additive identity and the roughness
condition hold pointwise. Formula (21) bounds its ordinary parameter
length below by a constant times m^-3. The cells are disjoint; summing
and using (26) proves (28). Their images are open intervals because Xi is
a diffeomorphism.

For all sufficiently large m every form in (25) is greater than 1 and
strictly less than 6m. If one had at least 38 prime factors, with
multiplicity, and each factor were at least (6m)^(1/38), it would be at
least 6m, a contradiction. Thus Omega<=37 for each individual form.
QED.

Equivalently, the calibrated measure of E_m is exactly
(1-2eta_m)|S_m|, and its ordinary measure is comparable to m^-3|S_m|.
The constants and the onset in this theorem are not numerically
calibrated. The number 37 is inherited from the cited conservative sieve;
it is not presented as an optimized almost-prime bound or a new sieve.

## 7. What this establishes, and what it leaves untouched

Taking F(a,b)=1_{a prime}1_{b prime}, formula (16) with q=2 gives an exact
guarded integral for the number of twin-prime pairs with first member in
{M+1,...,2M}. The twin-prime conjecture is therefore equivalent to
positivity of these guarded integrals for arbitrarily large dyadic
M=2^j. This is an equivalence, with no lower bound for the prime integral.

For the same prime predicate, (18) counts prime Goldbach representations
whose first prime lies in I_m. This is a **middle-third, balanced** count;
it must not be equated without qualification to unrestricted Goldbach
representations. Taking instead the predicate that all four forms in
(25) are prime counts balanced twin-center complementary pairs. Theorem 4
proves a rough almost-prime version, not positivity of that four-prime
integral. Analytic smoothness and exact composition laws do not estimate
any of these prime predicates.

The additional content beyond relabeling a staircase by its exact
thresholds is specific and limited: (4), (9), and (11) are closed
real-analytic formulas, with the translations independent of C_R, exact
algebraic relations, and a proved shrinking relative exceptional measure
for the actual dynamics.
They require no threshold lookup, no floor in the map formulas, and no
evaluation of N to construct the maps. The nonlinear dynamical work is
the inherited quantitative threshold estimate; the almost-prime work is
the inherited sieve. These operations are not an independent mechanism
for overcoming the prime-distribution obstruction.

Finally, C_R in the formulas means the exact endpoint constant. A fixed
numerical approximation with error at most kappa cannot silently replace
it while retaining a guard shrinking to zero. In coordinates built from
that approximation, the true phase differs by a constant of absolute
value at most kappa, so the rounding arguments require the guard width
eta_L+kappa instead. A practical certified finite-scale algorithm would
also require a numerical bound for K_R and a verified onset. Neither is
supplied by finite tests in this packet.

## 8. Reproducible controls and their limits

The companion script `verify_cusp_arithmetic_operations_20261019.py`
checks exact rational ceiling and guard identities, a hostile example
for unguarded positivity, numerical real-analytic composition identities,
the derivative in (20), and a finite set of actual escape pairs using
directed interval arithmetic with the previously certified enclosure

\[
-2.463836664<C_2<-2.463836643.
\]

The directed-interval checks use only this inherited enclosure and the
original recurrence, with no threshold table. They certify their listed
finite samples, not the existential asymptotic guard constant or onset.
The proof of the general theorem is Sections 2-6 above.
