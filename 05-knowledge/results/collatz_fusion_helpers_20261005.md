# Four receipt helpers: conformal repairs, a paid section return, and capacity

2026-10-05. **PROVED, scoped, author-audited:** finite-pool conformal repair;
cofinal neutral stock and its exact component classification; the canonical 3-parent section;
the J common-future controller, credit, guard language and decoder;
coupled prime/composite transport; a finite-seed depth/height bound.
**CITED:** Applegate--Lagarias's semigroup theorem and conformal Graver
decomposition. **FINITE-EXACT:** 631,890 checks in the accompanying program.
**OPEN:** universal grounded entry, an effective universally successful
edge-pool enrichment rule, and Collatz. These are scoped research results,
not an independent canon audit or a claim of literature priority.

Artifacts: [program](../../04-computation/experiments/collatz_fusion_helpers_20261005.py),
[exact JSON](collatz_fusion_helpers_20261005.json).
The four questions are inherited from
[universal weak receipts, section 8](collatz_universal_weak_receipts_20261005.md).

## 1. Inheritance, concept board, and incoming work

The anchor is the passage from universal weak receipts to grounded
certificates. The niche is conformal integer augmentation, recovered from
AMM12592. The wildcard is common catalytic context and finite-address capacity.

- Closest proved mechanism: W1, D1/D2, A1/A3 and G1/G2 in the
  [weak-receipt note](collatz_universal_weak_receipts_20261005.md): universal
  prime balance, labelled endpoint defects, source anchoring, and grounded
  deficit extraction.
- Canonical hostile: k copies of `3 -> 5 -> 1` certify `3^k` multiplicatively
  but have no outgoing stock at `3^k` for k>1. Neutral addition cannot fix
  that missing 3-source.
- Corrected near miss: a common-future signed packet has value n/d, not n;
  its open endpoint is not by itself a prime-balanced composite defect.
- Least-used relevant sidecar: retain the **sign orthant of every defect
  coordinate together with nonnegative edge stock**. This is exactly what
  the AMM slack-kernel argument needs before partial moves remain feasible.

The live board is **edge stock / endpoint defect / prime charge / 3-parent
section / paid guard / ordinary height**. Each experiment below changes at
least one of these coordinates without silently forgetting the others.

Concurrent commit `7de2181a21`,
[the four-question receipt note](collatz_receipts_four_questions_20261005.md),
adds the useful observation that the published semigroup construction inserts
multipliers at particular orbit labels. A local insertion has endpoint defect

\[
I_m(y)=[my]-[y]+\partial w_m,\qquad V(w_m)=1/m.
\]

This entire vector is prime-balanced. The moving label [my] is important,
but [y] and the fixed wild-certificate endpoints remain part of the vector;
one insertion is not automatically a single K-coordinate. This suggests
repairing insertion families one at a time. Section 4 supplies a paid
subfamily of the published multiplier-11 class.

## 2. Q1: a complete well-founded repair rule inside a finite edge pool

Let E be a finite set of actual positive odd Collatz edges. Write

\[
U(u)=(3u+1)/2^{v_2(3u+1)},\quad
Bc=\sum_{u\in E}c_u([u]-[U(u)]),\quad A=\Lambda B,
\]

where Lambda records the valuations at odd primes. Include n and 1 among
the vertex coordinates even if they are absent from the edge endpoints.
For a nonnegative weak receipt c for n, put b_n=[n]-[1] and D=Bc-b_n.
Then Ac=v(n) and D is prime-balanced with augmentation zero.

**Proposition R1 (conditional finite-pool completeness).** Form the integer
matrix

\[
M=\begin{pmatrix}A&0\\ B&-I\end{pmatrix}.
\]

If E contains a nonnegative receipt r with Br=b_n, then whenever D is
nonzero there is a Graver move (z,Bz) of M such that

\[
c+z\ge0,\quad Az=0,\quad Bz\text{ is conformal to }-D,\quad
\|D+Bz\|_1<\|D\|_1.                                      \tag{1}
\]

Consequently, repeatedly taking any move from the complete augmented
Graver bank satisfying (1) terminates at zero defect. If that bank has no
applicable strict move at a nonzero defect, no zero-defect receipt exists
in this pool.

**Proof.** The integer kernel vector `(r-c,-D)` is a conformal sum of
Graver atoms of M. Their second coordinates sum to -D, so at least one
has a nonzero second coordinate. Every partial sum of their edge
coordinates lies coordinatewise between c and r after addition to c;
it therefore preserves nonnegative stock. A nonzero second coordinate
conformal to -D reduces its integer l1 norm strictly. The same argument
applies at every new weak receipt because r remains feasible. The last
assertion is the contrapositive. The decomposition used here is standard:
[Onn, *Convex Discrete Optimization*, Lemma 4.2](https://arxiv.org/pdf/math/0703575).
The finite Graver bank follows from the orthant-wise Dickson lemma. □

The rank is **the norm of the full labelled endpoint vector**, not merely
the norm of its composite-coordinate coefficients. The theorem neither
supplies r nor bounds a useful edge pool uniformly in n. In this finite
deterministic graph, existence of r is also ordinary reachability of 1
from n, by G2; graph search is cheaper than computing a Graver basis.
The gain is a complete stock-preserving repair language, suitable for
seeking symbolic families. It is not a new reachability oracle.

### 2.1. A fully explicit repair, and its stopping boundary

Take E={3,5,7,9,11,13,17}. Starting with c_0=2e_3+2e_5 for 9, use

\[
z=-2e_3-e_5+e_7+e_9+e_{11}+e_{13}+e_{17}.
\]

Its multiplicative value is 1, and
`Bz=[9]+[1]-2[3]=R(3,3)`. The resulting stock is exactly the path
`9 -> 7 -> 11 -> 17 -> 13 -> 5 -> 1`; defect norm falls from 4 to 0.
Exhausting every conformal subvector finds no proper nonzero prime-balanced
subrelation. This verifies primitiveness in this pool, not an unrestricted
Graver-basis computation.

Prime balance in the same pool forces every receipt for `3^k` to have

\[
c_3=k-2t,\quad c_5=k-t,\quad
c_7=c_9=c_{11}=c_{13}=c_{17}=t,\quad0\le t\le\lfloor k/2\rfloor.
\]

Indeed balance at primes 7,11,17,13 makes the five t-coordinates equal;
balance at 5 and 3 then gives the first two formulas. Each available z
reduces the defect norm. For k>=3 it stops nonzero: the pool has no edge
out of `3^k`. For k=3 the norm sequence is 6,4; for k=4 it is 8,6,4.
This is the exact boundary between repair and enrichment.

### 2.2. A universal neutral buffer, restricted to the correct sector

**Proposition R2 (cofinal neutral stock).** For every finite nonnegative
stock request b supported on odd sources coprime to 3, there is a finite
nonnegative h>=b with V(h)=1. No such h can contain a source divisible by 3.

**Proof.** For an allowed edge u, the odd rational U(u)/u has denominator
coprime to 3. W1, derived from
[Applegate--Lagarias, Theorem 1.1](https://arxiv.org/pdf/math/0411140),
supplies a weak receipt for that reciprocal. Add e_u to obtain a neutral
packet containing u, and sum as many copies as b requests. Conversely,
`v_3(V(h))=sum h(u)v_3(u)` because every U(u) is a 3-unit. Neutrality makes
this nonnegative sum zero. □

Thus a finite signed prime-balanced rewrite schedule between nonnegative
stocks can be made stock-feasible by one common neutral buffer, provided
its intermediate stock at every 3-divisible source is already nonnegative.
Take the maximum unit-source shortfall over the finite schedule and apply
R2; the same buffer is present at both ends and may finally be removed.
This does not guarantee monotone defect decrease or grounded deficits.

The exact control uses the inherited neutral packet
`{7:2,11:1,17:1,55:1,65:1,83:1,5:1}`. Four copies dominate three separately
specified stock requests. The cross-domain resemblance to finite catalytic
context is now an actual cofinality theorem on a specified sector.

**Corollary R3 (complete neutral-extension classification).** Two weak
receipts c,d for the same n admit a common extension
`c+h_c=d+h_d`, with both h nonnegative and neutral, iff c and d have the
same multiplicity at every source divisible by3. The same holds for any
finite collection sharing this exact source profile.

**Proof.** Necessity is the last assertion of R2. For sufficiency let b be
the coordinatewise maximum of the receipts. Their common 3-source profile
means `V(b)/n` has 3-adic valuation zero. W1 supplies a nonnegative receipt
r with `V(r)=n/V(b)`. Then b+r dominates every input, has value n, and its
difference from any input is a nonnegative neutral packet. □

Thus the weak-receipt fibre of every odd 3-unit is directed under neutral
extension. For a multiple of3 the exact 3-source profile classifies its
components; all receipts anchored at their own source n lie in the one
component whose profile is one copy of n. This strengthens universal weak
coverage to **finite compatibility of all such anchored representations**.
The two-edge receipt for9 and its root path have different profiles and
cannot have a common neutral extension, although the signed repair in
section2.1 connects them. Common extension is not zero-defect repair:
its added packet can introduce new endpoint obligations.

## 3. Q2: the canonical 3-parent section and its induction prices

For odd positive v coprime to 3, define a_0(v) in {1,...,6} by
`2^a_0 v=1 mod9`, and

\[
\rho(v)=\frac{2^{a_0(v)}v-1}{3}.
\]

The six residue prices are

| v modulo 9 | 1 | 2 | 4 | 5 | 7 | 8 |
|---|---:|---:|---:|---:|---:|---:|
| a_0(v) | 6 | 5 | 4 | 1 | 2 | 3 |

**Proposition S1.** Rho(v) is the least positive odd multiple of 3 mapping
to v. Its other such predecessors are
`(2^(a_0+6t)v-1)/3=S^(3t)(rho(v))`, t>=0, where S(n)=4n+1.
Thus U bijects the normalized section
`{n positive odd:3|n,1<=v_2(3n+1)<=6}` with the positive odd 3-units.

**Proof.** The order of 2 modulo 9 is 6. Divisibility by 9 gives the factor
3 in rho(v), and its numerator divided by 3 is odd. Since v is odd,
the displayed valuation is exact. Increasing its exponent by 6 multiplies
the pre-division numerator by 64, the same effect as three S operations. □

This section is not a convergence theorem: transporting U to it by
`rho o U o U` gives a conjugate copy of U on the 3-unit states. For an
induction confined to multiples of 3, a deficit v is paid by a smaller
3-parent exactly when

\[
2^{a_0(v)}v<3n+1.                                         \tag{2}
\]

Keeping the six prices is sharper than the uniform sufficient threshold
`v<3n/64`. The small bases 3,9,15,21 have direct root paths; rho(1)=21.
For n>21 a rooted endpoint 1 is already discharged. A least unrooted
multiple of 3, if one exists, must keep every forward unit endpoint v
outside the paid region (2). The concurrent four-question note independently
recovers this section; the new step here is a guarded paid return into it.

## 4. A new paid common-future rule inside the multiplier-11 class

**Proposition J1.** For every positive n=799 mod1024, put

\[
J(n)=\frac{243n+147}{256}.
\]

Then J(n) is positive, divisible by 3, smaller than n, and

\[
U^6(n)=\frac{729n+697}{512}=U(J(n)),                         \tag{3}
\]

where the exact left valuation word is `(1,1,1,1,2,3)` and the right word
is `(1)`. Every one of the six left forward states is larger than n.

**Proof.** Ordered affine composition gives the numerator `729n+697`.
Its odd-endpoint cylinder is exactly n=799 mod1024. Every prefix slope
`3^j/2^(a_1+...+a_j)` exceeds 1, and every prefix carry is positive.
The terminal endpoint is 5 mod9, so its minimal 3-parent uses exponent 1.
Substitution gives J. Writing n=799+1024t yields J(n)=759+972t, so J is
positive, divisible by 3, 3 mod4, and `n-J(n)=40+52t>0`. Equation (3)
and the valuation 1 on its right follow immediately. □

The least example is

`799 -> 1199 -> 1799 -> 2699 -> 4049 -> 3037 -> 1139 <- 759`.

For induction on multiples of 3, the subfamily is
`n=2847+3072t -> J(n)=2703+2916t`. This is a common-future dependency,
not a claim that 759 lies on the forward orbit of 799.

Every J source is 15 mod16, outside the native starting guards of the
frozen H/G/A/B/L alphabet. This proves a new entry cell relative to that
bank, without claiming it escaped every other bank in the repository.
Its odd-source density is 1/512. It also lies in the published semigroup
cell 31 mod64, where the construction inserts multiplier 11. J handles
one sixteenth of that cell through a smaller dependency with no multiplier
in the actual common-future paths. Rooting that dependency is still required.

### 4.1. Search universe and a bank beyond the first rule

The exact search retains positive valuation words whose **every prefix
coefficient** is greater than 1. It stops a branch the first time its
terminal endpoint is uniformly 5 mod9 and its minimal 3-parent has slope
less than 1. For depth >=2 the endpoint modulo9 is independent of the
source; the test is `B/Q=5 mod9` and `2*3^r<3Q`. Each returned cell is
also checked to pay at its least positive source, proving payment on that
whole arithmetic progression. These are disjoint prefix cylinders.

Through depth 12, the first-exit counts are

| odd depth | 6 | 7 | 9 | 11 | 12 |
|---|---:|---:|---:|---:|---:|
| cells | 1 | 5 | 25 | 129 | 475 |

There are 635 cells, with total odd density `4095/524288`. The unique
depth-six cell is J; none occurs earlier in this **specified search**.
No assertion of minimality among all possible controller constructions is
made. The search deliberately targets a short inverse exponent, not all
six section prices, and only the first new exit along a branch.

### 4.2. Fractional growth credits give a terminating six-letter controller

At every native J source,

\[
\frac{J(n)+5}{n+5}\le\frac{191}{201},\qquad
\frac98\left(\frac{191}{201}\right)^3
=\frac{6967871}{7218312}<1.                                \tag{4}
\]

The first ratio decreases with n; equality is attained at 799. Use the
integer account C>=0 with

\[
\mathcal E(n,C)=(n+5)^3(9/8)^C.
\]

J earns one credit; G spends three. The inherited awards for H,A,B,L
become respectively 6,3,3,12. Cubing the proved old potential inequalities
preserves their strict bounds; (4) imports J. All guards and nonnegative
credit prefixes remain mandatory. The inherited sources are
[adaptive credit](adaptive_credit_potential_20261004.md) and
[the L payment theorem](paid_guard_budget_20261005.md).

Starting at C=0, every nonempty legal episode strictly pays its immutable
source: `x+5 <= E(x,C)^(1/3) < N+5`. It terminates: each funding action
contracts E by a fixed factor less than 1, E>=216, and the number of G
actions is bounded by the credits minted. This proves termination of
available controller actions, not termination of a search for a missing
next action or completion of the remaining root obligation.

Three J credits can fund one G. The native word JJJG has cylinder
`129055519 mod268435456`. Two J copies cannot universally fund a G by
ordinary-size payment: the limiting slope is
`(243/256)^2*(9/8)=531441/524288>1`.

### 4.3. Exact guards and translation decoding survive the new letter

Use affine letters `(3^r n+b)/2^a`:

| letter | r | a | b | native source |
|---|---:|---:|---:|---|
| H | 6 | 10 | 669 | 155 mod2048 |
| G | 2 | 3 | 5 | 11 mod16 |
| A | 4 | 7 | 85 | 187 mod256 |
| B | 4 | 7 | 73 | 7 mod256 |
| L | 2 | 4 | -3 | 219 mod256 |
| J | 5 | 8 | 147 | 799 mod1024 |

**Proposition J2 (exact native language).** A formal word in these letters
has a native positive source iff it contains neither LB nor LJ. If its
composite is `(Pn+B)/Q`, its complete source guard is

\[
Pn+B=tQ\pmod{sQ},
\quad(s,t)=\begin{cases}(16,11)&\text{last letter L},\\
(4,3)&\text{last letter J},\\(2,1)&\text{otherwise}.\end{cases}          \tag{5}
\]

**Proof.** Native H/G/A/B is equivalent to odd output; native J is
equivalent to output 3 mod4; native L to output 11 mod16. All primitive
source guards are 3 mod4. Their residues mod16 are 11 except B at 7 and
J at 15. Thus an L output can precede precisely H/G/A/L; a J output can
precede every letter. These give the two forbidden pairs. Backward
induction now works because every slope numerator is odd and every
denominator a power of 2: the final congruence imposes each exact division,
and the next native input gives the preceding required output bits.
Odd P makes (5) one residue class. Positivity follows forward from the
native positive-source formulas. □

**Proposition J3 (lossless six-letter translation).** The exact dyadic
translation beta=B/Q determines the whole formal word. Modulo9 the last
letter tags are H/J:3, G:4, A:2, B:5, L:6. For tag3, read modulo243:

\[
\beta_H-\beta_J=669/1024-147/256=81/1024.
\]

They agree modulo81 but differ modulo243. Earlier translations are
multiplied by `3^6` for H or `3^5` for J, so vanish at that resolution.
After identifying a last letter, peel it using
`beta <- (2^a beta-b)/3^r`. The exact power-of-two denominator cost drops
by a, because every nonempty composite carry is odd. The empty word has
translation zero; no nonempty word does. This proves termination and
injectivity. □

Five ternary digits distinguish H from J. This is a structural use of
deeper prime-ary digits, not an inference from a shared small number.
The old four-letter shadow and eta reader are not automatically inherited
by the enlarged alphabet. Equation (5), credit and the actual dependency
compiler remain required even though translation retains the formal word.

**Proposition J4 (repeat fuel).** For positive odd n, the maximum number
of consecutive native J operations is

\[
\max\left(0,\left\lfloor\frac{v_2(13n-147)-2}{8}\right\rfloor\right).
\]

Indeed `13J(n)-147=(243/256)(13n-147)`, while the native guard is
`13n-147=0 mod1024`. Each operation consumes eight binary digits and
requires two extra guard digits. The rational fixed point 147/13 is not
an integer, so an integer has only finite repeat fuel. Least sources for
one, two, three repeats are 799,80671,61946655.

## 5. Q3: factorization defects need a coupled prime coordinate

Recall `K_u=[u]-sum_p v_p(u)[p]+(Omega(u)-1)[1]`, taking K_1 and K_p as
zero. Actual vertex transport is `U_*[u]=[U(u)]`.

**Proposition T1 (failure of kernel invariance, exact repair).** The
prime-balanced defect module is not U_*-invariant:

\[
U_*K_9=[7]-2[5]+[1],\qquad
\Lambda(U_*K_9)=e_7-2e_5\ne0.                              \tag{6}
\]

For every positive odd u the exact transport is

\[
U_*K_u=K_{U(u)}-\sum_pv_p(u)K_{U(p)}
       +\sum_\ell\chi_\ell(u)([\ell]-[1]),                 \tag{7}
\]

where `chi(u)=v(U(u))-sum_p v_p(u)v(U(p))`.
Expand the definition of K and collect prime coefficients to prove (7).
It is an identity for every u, not a statistical approximation. In
particular, no U_*-invariant submodule containing K_9 can remain inside
the prime-balanced kernel.

The repaired information-rich coordinate is **prime charge plus composite
defect plus actual edge stock**. A finite vector remains finite after
transport; no fixed finite-dimensional invariant space is claimed.

A common-future signed packet from n to d has boundary [n]-[d] and value
n/d. Relative to [n]-[1], its error is [1]-[d], whose prime charge is
`-v(d)`. Its full decomposition is

\[
[1]-[d]=-K_d-\sum_pv_p(d)([p]-[1]).                         \tag{8}
\]

Thus a prime child still carries an obligation: for d=3 it is [1]-[3],
despite K_3=0. Composition can retain one labelled endpoint, but erasing
its prime part is invalid. Equations (7)-(8) are the precise reconciliation
with the concurrent substitution framework.

### 5.1. Change the operation: stock-supported replacement preserves defects

**Proposition T2 (defect-preserving common-future replacement).** Let finite
actual paths from n and d end at the same e, with edge stocks p and q.
If a weak receipt c_n for n contains p coordinatewise, then

\[
c_d=c_n-p+q\ge0,\qquad V(c_d)=d,\qquad D_d(c_d)=D_n(c_n).
\]

Conversely, a weak receipt c_d containing q gives
`c_n=c_d-q+p` with the same properties in reverse.

**Proof.** The removed stock is present by hypothesis. The path products
are n/e and d/e, proving the value statement. Their boundaries differ by
`Bp-Bq=[n]-[d]`; substituting into the definition of D cancels exactly. □

Any prescribed finite path can be included in a weak receipt by A1, so
these are nonempty domains at the weak level. Every specified class of
defect vectors is preserved wherever the replacement is stock-supported.
No invariant-kernel assertion about U_* is needed. The operation changes
the target and its supported path, rather than advancing every vertex label.

For J, take p to be the six-step source path and q the one edge from
J(n). T2 then reduces the integer target strictly while keeping its full
defect unchanged. In the reverse direction it transports a rooted child
certificate into a rooted source certificate, because a literal child
root path contains q. The 32 exact controls add a nonzero neutral defect
to each child certificate and verify both replacement directions: the
identity does not depend on the special case D=0.

This is a positive answer to the **stability** part of Q3 for a properly
typed operation. It does not imply that an arbitrary existing receipt
contains the next requested path, that re-anchoring preserves its defect,
or that every integer has a paid controller. These remain stock and
coverage obligations. It also explains why the prime imbalance of a signed
substitution alone disappears only after composing it with the child's
weak receipt, whose prime content is exactly the missing v(d).

## 6. Q4: finite rooted portfolios improve height but have a capacity floor

Let S be any K positive odd rooted seeds. Let N_(S,D)(X) count positive
odd n<=X that reach S in at most D actual odd steps, including depth zero.
Write L=floor(log2 X).

**Proposition C1 (finite-seed source capacity).**

\[
N_{S,D}(X)\le K\sum_{d=0}^D {L+2d\choose d}
           \le K(D+1){L+2D\choose D}.                      \tag{9}
\]

**Proof.** An actual d-step route with valuation sum A and endpoint s
satisfies

\[
2^A=\frac ns\prod_{i=0}^{d-1}(3+1/n_i)\le4^dX.
\]

Hence A<=L+2d. There are `binom(L+2d,d)` positive compositions of length d
and sum at most L+2d. For any one such word and fixed s, the affine equation
gives at most one source n. Count all words and seeds, allowing overlaps
and illegal words only to increase the upper bound. The binomial terms
increase with d, giving the last inequality. □

**Corollary C2 (precision costs depth or seed capacity).** Suppose every
residue modulo3^h has a representative in this rooted portfolio below
X_h<=exp(C h), where C is fixed, and `log K_h=o(h)`. Then D_h cannot be
o(h); in fact a uniformly successful such scheme needs D_h=Omega(h).
If D=o(h), the estimate `binom(N,D)<=(eN/D)^D` makes the logarithm of
the right side of (9) o(h), contradicting the required 3^h distinct
representatives. The same argument on a subsequence excludes D_h/h
approaching zero. This counts literal odd depth, including all steps
hidden inside a compressed controller.

This refines the inherited bounded-rank/polylogarithmic source obstruction
in [guarded pumping, sections 1-4](collatz_guarded_pumping_memory_20261004.md).
It does not preclude a growing-depth compiler or constrain arbitrary
compressed-program size without a separate depth theorem.

The exact experiment uses the first 16 odd 3-units through47 as independently
rooted seeds. For each residue modulo3^h it minimizes the one-step inverse
source `(2^a s-1)/3` over those seeds. The period `2*3^h` in a covers
all possible residue choices and larger exponents only increase sources.

| ternary precision h | worst least-source bits, seed {1} | first 16 seeds |
|---|---:|---:|
| 1 | 5 | 3 |
| 2 | 17 | 5 |
| 3 | 53 | 10 |
| 4 | 161 | 30 |
| 5 | 485 | 156 |
| 6 | 1457 | 288 |

The finite improvement is substantial. It does not turn local coverage
into exact entry at an arbitrary prescribed integer. The bound says what
has to grow in a scalable construction: rooted seed capacity, literal
certificate depth, or ordinary height.

## 7. Structural transfers recovered from apparently unrelated work

These are maps with named losses, not an assertion that the original
problems are equivalent. The first two directly produced R1 and R2.

| source and predicate | target and map | lost information / required sidecar | decisive test |
|---|---|---|---|
| [THM-3373, AMM r8 slack-kernel conformal locality](../../01-canon/theorems/THM-3373-r8-slack-kernel-conformal-locality-width-five.md): partial conformal sums retain nonnegative inequality slacks | augment the edge kernel by the endpoint-defect coordinates, then decompose sign-compatibly | need the finite edge pool, stock and defect orthant; AMM's width5 and numerical bounds do not transfer | the 9 repair succeeds; the 27 fibre stalls without its source edge |
| [THM-2281, common optimal context for finite catalytic families](../../01-canon/theorems/THM-2281-common-optimal-context-for-finite-catalytic-families.md): finitely many context requirements combine | neutral packets dominate a finite set of unit-source deficits | finiteness and the excluded 3-source sector; no knot-theoretic existence theorem is imported | neutral wild-5 packet covers three requests; v3 forbids stock at3 |
| [THM-3990, componentwise harmonic obstruction and repair quotient](../../01-canon/theorems/THM-3990-componentwise-harmonic-obstruction-and-repair-quotient.md): the actual repair operator determines its surviving quotient | test prime projection under actual U_* instead of assuming closure | prime charge must accompany composite coordinates; positivity is a separate condition | equation(6) immediately leaves the kernel |
| [THM-3357, Berggren three-branch Walsh collapse and parent circuit](../../01-canon/theorems/THM-3357-berggren-three-branch-walsh-level-collapse-and-parent-circuit.md): a level aggregate can forget branch history | compare exact controller translation with its low-resolution ternary tag | H and J need five ternary digits; exact arithmetic source guards remain essential | H/J agree modulo81 but separate modulo243 |
| [THM-4495, exact no-descent counts and Spitzer ballot](../../01-canon/theorems/THM-4495-no-descent-count-exact-order-spitzer-ballot.md): survival is an all-prefix predicate | stop the new section search at its first paid exit | terminal coefficient and all-prefix conditions are different; shortcut clock and odd-step clock have different counts | the word(2,1) contracts first and expands overall |
| [THM-3340, AMM cyclic single-donor completion](../../01-canon/theorems/THM-3340-single-donor-cyclic-rotation-proves-all-pointwise-AMM-floors.md): a donor must supply the actual local deficiency | ask whether neutral context supplies the missing source stock | a surplus at a different label does not fund an edge at n | no neutral packet can supply a missing source3^k |
| [finite-seed compiler](finite_seed_receipt_compiler_20261005.md): dependency leaves retain proofs or stay unresolved | compile J as a common-future edge and keep its smaller child | signed cancellation is not a literal forward step; root leaves and credit are separate | equation(3) certifies a join while all six source states grow |
| [guarded pumping](collatz_guarded_pumping_memory_20261004.md): bounded arithmetic rank limits ordinary source capacity | count valuation compositions with height bound | digit-surjectivity does not bound the integer represented | local coverage at h=6 still needs a 288-bit representative with16 seeds |

The tournament analogy is useful only at the level of **small reusable
certificates and quotient losses** here. There is no intrinsic pairwise
orientation in this repair problem, so no tournament is imposed on it.
Likewise, the number11 in the semigroup class and fifth ternary digit in
the decoder have explicit roles; their coexistence is not a new theorem
about the golden ratio, eta coefficients, or primes.

### 7.1. Reusable move: augment the kernel by the coordinates that must stay feasible

**Trigger/action:** when a signed kernel identity cannot be executed with
nonnegative resources, append the slack or defect coordinates before taking
conformal moves. Decompose toward a feasible target in that augmented kernel;
sign-compatible partial sums preserve stock and can decrease an integer
defect norm.
**Counterindications:** an ordinary lattice basis need not connect a
nonnegative fibre. A complete move bank in a fixed pool does not prove that
the pool contains the desired target or pay for enlarging it. Check an absent
source/resource, not only the successful small example.
**Evidence:** [THM-3373, AMM slack-kernel locality](../../01-canon/theorems/THM-3373-r8-slack-kernel-conformal-locality-width-five.md)
retains inequality slacks under partial conformal sums;
[Collatz receipt helpers R1](collatz_fusion_helpers_20261005.md)
retains edge stock and endpoint-defect signs. The receipt for9 repairs in its
seven-edge pool, while the same family for27 stalls without an edge at27.
No AMM locality constant is transferred to the Collatz problem.

The global meta-pattern index already sits at its startup-size budget; this
fully specified candidate card stays with the supporting research note.

## 8. Reproducibility, failure boundaries, and next helper questions

Run:

```bash
python3 04-computation/experiments/collatz_fusion_helpers_20261005.py
python3 -O 04-computation/experiments/collatz_fusion_helpers_20261005.py
```

Both outputs agree exactly. The 631,890 checks include: all conformal
subvectors of the displayed signed repair; power3 fibres k=1..12; section
normalization for odd unit endpoints below10002; the complete stated
first-exit word search through12; 1000 J-source samples; all9330 formal
six-letter words of lengths1..5 (7436 native,1894 forbidden), with separate
guard intersection and composed guard readers; six repeat depths; coupled
transport on501 odd labels; and direct bounded-source censuses against(9).
Finite samples supplement the all-parameter proofs; they do not establish
universal Collatz coverage. No unrestricted Graver bank was computed.

An earlier incoming Pascal note conflated endpoint coefficient payment
with all-prefix survival. The repaired
[prefix gate](collatz_pascal_prefix_gate_20261005.json) has1178 exact checks
and retains the exponential rate by a cyclic-rotation sandwich, not by
identifying the two events. The later four-question note reused that
conflation; its correction also retains this distinction. Its separate
[type audit](collatz_receipt_type_audit_20261005.json) has2,059 exact checks;
[the program](../../04-computation/experiments/collatz_receipt_type_audit_20261005.py)
checks prime-child obligations, the excluded depth-two join, the full
insertion defect and all2048 odd prefix-code residues.

The productive next questions are now narrower:

1. **Parameterized enrichment:** can one add a bounded schema of actual
   edges to repair an entire multiplier-insertion family I_m(y), while
   preserving a global rank across enrichment? R1 supplies the finite-pool
   repair mechanism; it does not pay the addition of new labels.
2. **All six section prices:** can the first-exit search use residues other
   than5 mod9, with a uniform payment account strong enough to import them?
   The 635 cells only use the cheapest inverse exponent1.
3. **Guard-aware repair compilation:** can (7)-(8) be combined with the
   lossless translation to synthesize repairs without dropping the prime
   charge or the source stock? Fixed-dimensional compression must be proved,
   not inferred from one affine obligation.
4. **Depth matched to precision:** can a recursively certified seed bank
   realize all ternary addresses at ordinary height exp(O(h)) with actual
   depth O(h)? C2 identifies this as the first scalable regime not excluded
   by capacity. Such address coverage would still need a theorem turning
   the prescribed source into the representative, rather than a nearby one.

The proposed lifted object is `(n, edge pool, stock, defect, prime charge,
guarded word, credit, rooted leaves)`, projected to n. W1 gives universal
weak fibres; R1, R2 and J add checked operations on those fibres. A global
well-founded rule that always reaches grounded leaves remains the missing
theorem. None of the local ranks above is silently extended across arbitrary
pool growth or an unproved change of representative.
