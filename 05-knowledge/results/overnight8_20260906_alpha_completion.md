# Complete alpha parity through one real-rooted carrier

Status: **PROVED ANALYTICALLY / INDEPENDENTLY AUDITED** for the
all-`A=2` completed-alpha, beta-hit theorem and the exact `A2B3`
mixed-product reduction below. The sign bank is **FINITE-EXACT**.
Universal beta-skip negativity and universal actual doubled-row
noncancellation remain **OPEN**. No scarce theorem ID is claimed.
[Independent analytic and exact audit: PASS](overnight8_20260906_alpha_independent_audit.md).

## Inheritance and the conceptual change

The closest proved input is
[THM-4440, signed duplication SOS and real-rooted Laurent return](../../01-canon/theorems/THM-4440-signed-duplication-sos-and-real-rooted-laurent-return.md),
including its strict interior-coefficient conclusion and repeated-root
scope. The complete PF path factors and virtual Hadamard carrier are
inherited from [the incoming transport theorem](nc2_hadamard_transport_overnight_hexagon_sep05.md).
The complete first-root geometry is
[THM-4436](../../01-canon/theorems/THM-4436-complete-factorial-row-simple-negative-roots-and-trinomial-phase-wall.md).
The [seventh-wave three-group identity](overnight7_20260906_laurent_midpoint_transport.md)
supplies the exact beta skip classes, with all negative indices.

The canonical hostile is still `h=1,g=5,rho=-2`:
`O star B=0`, while `O star C=15`, `O star D=10`, and `E star B=-16`.
Applying a separate zero-coefficient SOS after replacing either factor
would therefore use a false root condition. The corrected near miss is
the truncated PF factor `21+20t+5t^2`: dropping its full negative-index
ancestry destroys real-rootedness. The least-used sidecar here is the
square-root pullback of the **full** negative-root beta polynomial.

The live board has six objects:

1. Complete auxiliary support and its PF root geometry.
2. The original first-row zero coefficient, retained literally.
3. Both alpha midpoint parity classes in one binomial power.
4. Beta paths hitting versus skipping their midpoint level.
5. The two contiguous beta polynomials sharing that same path carrier.
6. The square coefficient and its signed SOS margin.

The incoming [quartic stability theorem](signed_duplication_stability_empty_core_sep06.md)
and [cubic real-core sector classification](trinomial_cubic_sector_empty_core_sep06.md)
were read before this derivation. Their exact degree and phase restrictions
are retained. Our carrier usually has much higher degree, so the sharp
quartic `7/9` margin does not apply. Our carrier is also different from
the ordinary cubic `tau+s^2+s^3`: its real-rootedness holds at every
negative phase and does not expand the cubic real-core sector itself.
Targeted searches for the quadratic pullback, alpha completion, beta-hit
response and their coefficient maps in these current notes and the
mistakes ledger found no existing statement of the theorem below.
No external priority is asserted.

## 1. All-A2 completed-alpha theorem

Take integers

```text
B>=3, h>=1, x>=0, 0<=r<B, z in {0,1},
delta=B-2, y=B*h+r, m=x+y+z, L=floor(x/delta), q=L+h.
```

Define full finite sequences, zero outside the displayed count domains,

```text
alpha_j = binom(m,z+2j),                         0<=z+2j<=m,
beta_j  = binom(x+y-2j,x+delta*j),               -L<=j<=h,
P(t)    = sum_(j=0)^h alpha_j beta_j t^j,
G(t)    = sum_(j=-L)^h beta_j t^(j+L),
alpha_double_j = binom(2m,2z+2j),               0<=2z+2j<=2m,
W(t)    = sum_j alpha_double_j (beta*beta)_j t^j.
```

Here `*` is ordinary convolution. In the definition of P, only the
nonnegative alpha indices overlap beta; alpha_double retains its genuine
negative index `j=-1` when `z=1`. The beta convolution retains every
negative auxiliary index before that overlap is taken.

**Theorem.** For every negative real `t` with `P(t)=0`,

```text
W(t)<0.                                                   (1)
```

No coprimality or first-return interpretation is needed for (1). The
complete path PF input gives G exact degree q, positive constant and
leading terms, and all roots strictly negative. Repeated roots would
cause no difficulty in this argument. First-return coprimality is a
separate consumer filter, not an analytic hypothesis.

Put `k=2h+z` and form the ordinary polynomial

```text
H_t(u) = u^(2q) (1+u)^m G(t/u^2).                          (2)
```

For each negative root lambda of G, the two roots introduced by (2)
satisfy `u^2=t/lambda>0`. Together with m copies of `-1`, these account
for every root of H_t. Thus H_t is real-rooted. More precisely,

```text
H_t(0) = beta_h t^q !=0,
deg(H_t) = m+2q,
0<k<m+2q,
(m+2q)-k = m+2L-z = x+y+2L>0.                             (3)
```

There is no hidden zero factor or missing endpoint. Literal coefficient
extraction now gives the same original root constraint:

```text
[u^k] H_t
 = sum_(j=-L)^h beta_j t^(j+L) binom(m,k-2q+2(j+L))
 = t^L sum_j beta_j t^j binom(m,z+2j)
 = t^L P(t).                                               (4)
```

The same extraction after squaring, without deleting any parity class,
is

```text
[u^(2k)] H_t^2
 = t^(2L) sum_j (beta*beta)_j t^j binom(2m,2z+2j)
 = t^(2L) W(t).                                            (5)
```

At P(t)=0, (3)--(4) meet exactly the strict interior hypothesis of
THM-4440. Therefore (5) is strictly negative. The factor `t^(2L)` is
positive, proving (1).

The inherited quantitative version is also available. Write
`H_t(u)=H_t(0) product_i(1+s_i u)`, with nonzero real s_i. Then

```text
W(t) <= - beta_h^2 t^(2h)/(4h+2z-1)
             * sum_(|I|=2h+z) (product_(i in I) s_i)^2 <0.   (6)
```

This is an explicit retained margin, not a claimed sharp constant for
these carriers. It depends on the same H_t, not independent replacements
of its factors. Positivity of all subset-square terms and the interior
degree bound make the final strictness immediate.

## 2. Exact progress on the all-h A2B3 family

Specialize `B=3,r=0,z=1,x>=1`, so `g=m=x+3h+1`, `L=x`, and `q=x+h`.
Use the full seventh-wave notation

```text
O(t)=sum_j binom(g,2j+1)t^j,
E(t)=sum_j binom(g,2j)t^j,
B(t)=beta(x,3h), C(t)=beta(x,3h-1), D(t)=beta(x,3h-2),
A_double(t)=O(t)^2+t^(-1)E(t)^2,
B_double(t)=B(t)^2+2t C(t)D(t),
P=O star B,
Q_raw=A_double star B_double,
V=O^2 star B^2,
G1=(t^(-1)E^2) star B^2,
G2=2 O^2 star(t C D),
G3=2(t^(-1)E^2) star(t C D).
```

The star is coefficientwise multiplication of the entire Laurent
sequences. The new conclusion is the uniform combined sign

```text
P(rho)=0, rho<0  ==>  W(rho)=(V+G1)(rho)<0                 (7)
```

for every h and x in scope. The old virtual carrier proved V<0 alone;
(7) absorbs the missing alpha parity in one carrier. It does **not**
separately prove G1<0. Those two statements are not equivalent.

The actual row has the exact decomposition

```text
Q_raw = W + R_skip,
R_skip = 2 A_double star(t C D) = G2+G3.                   (8)
```

Thus the sole remaining **sufficient** sign target is
`R_skip(rho)<=0` at the original P-root. Strict negativity of Q_raw is
equivalent to `R_skip(rho)<-W(rho)`, not to skip negativity. Equation
(6) supplies a quantitative payment if skip negativity proves too strong.
No uniform sign for R_skip is proved here.

The source supports are `(-6h-3,2g-6h-3,3g-6h-3)`. A primitive
first-return interpretation additionally requires `gcd(g,6h+3)=1`;
it is used in the finite first-root bank. Neither (7) nor the identity
(8) needs that filter. The true anchored second row, including its
`t^-1` lower carry, is exactly Q_raw.

## 3. One coupled mixed-product coefficient remains

Define the elementary composition polynomials

```text
F_n(t)=sum_(j=0)^floor(n/3) binom(n-2j,j)t^j.
```

Full support translation gives

```text
t^x B=F_(3q),  t^x C=F_(3q-1),  t^x D=F_(3q-2).
```

The familiar step-one/step-three recurrence
`F_n=F_(n-1)+t F_(n-3)` follows by classifying the last step of a
composition of n. The inherited midpoint split is the exact identity
`F_(2n)=F_n^2+2t F_(n-1)F_(n-2)` (here n=3q suffices).

For fixed negative t put

```text
K_B(u)=u^(2q) F_(3q)(t/u^2),
K_C(u)=u^(2q-2) F_(3q-1)(t/u^2),
K_D(u)=u^(2q-2) F_(3q-2)(t/u^2),
H_B=(1+u)^g K_B, H_C=(1+u)^g K_C, H_D=(1+u)^g K_D.
```

All three are real-rooted by the same full-path input and quadratic
pullback. They are **coupled** through q,t,g; they cannot be chosen
independently. Direct extraction yields

```text
P(t)=t^(-x)[u^(2h+1)]H_B,
W(t)=t^(-2x)[u^(4h+2)]H_B^2,
R_skip(t)=2 t^(1-2x)[u^(4h)]H_C H_D.                        (9)
```

The exact remaining sufficient statement is therefore

```text
[u^(2h+1)]H_B=0  ==>  [u^(4h)]H_C H_D>=0.                 (10)
```

Because `2t^(1-2x)<0`, the direction in (10) is correct. It is an
**OPEN coupled mixed-product inequality**, with a now-proved negative
square contribution beside it. Ordinary real-rootedness of H_C,H_D
alone does not prove (10), and their individual relevant coefficients
are not zero at the original P-root.

The contiguous relations survive as exact identities on this same
carrier. With primes denoting u derivatives,

```text
K_B' = (2/3)u[(q+2)K_C+u K_C'],
(q+1)K_D+u K_D' = 2K_C+(3/2)u K_C'.                       (11)
```

For example `(n-2Theta)F_(n-1)=(n-3Theta)F_n` with n=3q,
`Theta=t d/dt`, transforms into the first equation; replacing n by
n-1 gives the second. These transformations include the different
degrees q and q-1. Equation (11), rather than an independent-carrier
SOS, is a precise next attempt at (10).

## 4. Tests, failure boundaries, and connection contract

The standalone [companion](../../04-computation/overnight8_20260906_alpha_completion.py)
imports no repository mathematical producer. Its ordinary coefficient
convolution and Fraction quotient reduction are local implementations;
SymPy is used only to isolate real roots in the sign bank.

Its complete declared map universe is B=3,...,8; h=1,...,4;
x in {0,1,2,5}; r in {0,B-1}; z in {0,1}: 384 indexed rows,
each at t in {-1,-2,-1/3}. This intentionally includes nonprimitive
parameters and x=0, legal for the theorem. It checks both endpoints,
the interior bound and both literal coefficient extractions. These
finite controls challenge the written all-parameter identities; they
are not interpolation certificates or their proof.

Eighteen A2B3 rows, h=1,...,6 and x in {1,2,5}, separately check the
full beta split, both identities (11), and the mixed coefficient map
(9) at the same three phases.

The hostile bank starts beyond the frozen quartic endpoint: h=5,...,10,
x in {1,3,4,7,13,31}, g=x+3h+1, retaining exactly
`gcd(g,6h+3)=1`. This gives 32 indexed rows and all 240 first roots.
At every root, exact rational interval evaluation proves the three
old grouped responses, W, and R_skip strictly negative: 1,200 signs.
The three individual groups remain only **FINITE-EXACT**, even though
the combined W now has a uniform analytic proof. No h=5 symbolic
characteristic escalation was used.

The same-constraint hostile is replayed exactly. Positive-phase
quadratic pullback gives `u^2+1` already for G=1+t, so the negative
phase input cannot be silently discarded. Higher pullback degree also
does not preserve real-rootedness in general; the A=2 proof is not an
all-A theorem. The separate degree-five hostile to the incoming quartic
constant remains in force.

The connection contract is explicit:

| Item | Contract |
|---|---|
| Source | Complete A=2 beta PF path polynomial and binomial alpha power |
| Target | One ordinary real-rooted carrier H_t and its square coefficient |
| Map | G(t) to u^(2q)(1+u)^m G(t/u^2) |
| Preserved predicate | Original P(t)=0 becomes the exact same interior zero coefficient |
| Restored information | Both alpha midpoint parity classes, including the lower carry |
| Still omitted | Beta paths that skip their selected midpoint level |
| Needed sidecar | The coupled H_C,H_D mixed coefficient, or its payment against (6) |
| Cheapest decisive next test | Derive a same-carrier inequality from (11) implying (10), or find a coupled counterexample |

The conceptual change is to **complete a missing residue before taking
the square**, using a real-rooted pullback that retains the same zero
constraint. It is a structural partial closure, not evidence of a full
Laurent return theorem. The discarded step is now localized to one
mixed-product coefficient rather than hidden in an auxiliary square.

Reproduce outside the repository while this bundle is under audit:

```text
python -B 04-computation/overnight8_20260906_alpha_completion.py
python -B -O 04-computation/overnight8_20260906_alpha_completion.py
```

Both modes pass 5,638 explicit gates, with byte-identical LF output.
The semantic bank digest is
`f2b2cbf3308feb6568785245f439916d8823bb1c09dee0bbd1cbfef862be311d`.
This report and verifier are new outside-worktree files; no shared
navigation, canon, incoming proof, or Git state was edited by this lane.

**Filing:** root integrated these audited artifacts in the eighth checkpoint;
reproduction paths are relative to the repository root. Earlier outside-worktree
notes describe author provenance, not the present file location.
