# lane J1 -- the JC(2) leading-form tower: exact local law, the Euclid chain, and the search order

**status:** PROVED + VERIFIED-EXACT for sections 1-9 (including 4a', 4a) as
scoped there; section 10 collects them as necessary conditions; the search order
of section 11 is a *ranking*, not a theorem.

**headline:** the lane's "subtractive Euclidean algorithm" claim is CONFIRMED in
the divisibility form `H^(a-ib) | P_(n-i)` for `0<=i<=floor(a/b)` (PROVED), and
its *equality* form `v_L(P_(n-i)) = (a-ib)e_L` is REFUTED by an exact
automorphism witness.  Two new exact structures fall out: `P_(n-1)` is
*determined* by `Q_(m-1)` whenever `H` is not a power of a linear form, and the
tower's entire free-parameter budget is controlled by the single integer
`g_0 = deg H_0`, with a new scalar available exactly at the orders `j = 0 mod g_0`.

**script:** `04-computation/jc2_tower_depth_laneJ1.py`
**output:** `05-knowledge/results/jc2_tower_depth_laneJ1.out`
**arithmetic:** exact throughout (sympy over QQ, `fractions.Fraction` for slopes).
**script_sha256:** `fa0ee87d0c190dec58d10b455c9d75d60faa8fe21a2cad7c72e8fe2ce2eb7ad5`
**output_sha256:** `09cdb8b08f44adfcfe402945b62daec070c58d060c23668c5908dad130f27922`

**scope discipline (MISTAKE-237).**  Everything below is a *stratification of the
leading-form data of a hypothetical Jacobian pair* together with *necessary
conditions*.  Nothing here is an equivalence, a reduction, or a bridge.  In
particular: no claim that JC(2) reduces to a continued-fraction statement, none
that Zaremba's conjecture is involved beyond supplying the ordering key
`z = max partial quotient`, and none connecting this to NC2 / GMC(2) / VC(4).
The word "depth" always means the length of a continued fraction, never a
measure of difficulty of JC(2).

---

## 0. Setup and notation

`P,Q in C[x,y]`, `Jac(P,Q)=P_x Q_y - P_y Q_x = 1`, `n=deg P>=2`, `m=deg Q>=2`.
Homogeneous decomposition `P=sum_(i>=0) P_(n-i)`, `Q=sum_(k>=0) Q_(m-k)`.
Separating `Jac(P,Q)=1` by total degree gives the **graded tower**

    (G_j)   sum_(i+k=j) Jac(P_(n-i), Q_(m-k)) = 0     for 0 <= j <= n+m-3,
                                               = 1     for j = n+m-2.        (1)

`(G_0)` is the lane's `(L0)`: `Jac(P_n,Q_m)=0` forces

    P_n = c H^a,  Q_m = c' H^b,  g=gcd(n,m),  a=n/g,  b=m/g,  gcd(a,b)=1,
    deg H = g.                                                              (2)

Factor `H = c_0 H_0^d` with `H_0` **primitive** (not a proper power) and set
`g_0 = deg H_0 = g/d`; equivalently `d = gcd{e_L}` over the multiplicities `e_L`
of the distinct linear factors `L` of `H`, and `K` = number of distinct `L`.
Note `g_0 = 1` exactly when `H` is a power of a single linear form (`K=1`).
Throughout `wlog a >= b`; the opposite case is the mirror image.

`(G_(n+m-2))` collapses to `Jac(P_1,Q_1)=1`, so **`P_1` and `Q_1` are independent
linear forms** -- used repeatedly below as the tower's boundary condition.

---

## 1. Lemma J -- the local Jacobian valuation law  (PROVED + VERIFIED-EXACT)

Let `L` be a linear form, `F,G` nonzero binary forms of degrees `d_F,d_G>=1`,
`phi=v_L(F)`, `psi=v_L(G)`, `F=L^phi F'`, `G=L^psi G'`.  Then

    Jac(F,G) = L^(phi+psi-1) [ phi F' Jac(L,G') + psi G' Jac(F',L)
                               + L Jac(F',G') ],
    Jac(F,G) = L^(phi+psi-1) F'G' (phi deg G' - psi deg F')  mod L^(phi+psi).  (3)

*Proof.*  Expand the bilinear log-derivative `Jac(F,G)=FG*Jac(log F, log G)`.
For the residue use Euler: on the zero of `L`, `Jac(F',L) = -(deg F') F'` and
`Jac(L,G') = (deg G') G'`.  QED.

**Consequence (Lemma J).**

    v_L(Jac(F,G)) = v_L F + v_L G - 1,
      UNLESS  v_L(F)/deg F = v_L(G)/deg G,  in which case  >= v_L F + v_L G.   (4)

Call the exceptional case **density-matched** at `L`.

VERIFIED-EXACT: 250 random exact instances (207 generic, 43 density-matched),
**zero** violations (`PART 1`).

---

## 2. TTL -- the tropical tower law  (PROVED + VERIFIED-EXACT)

Fix `L | H` with multiplicity `e=e_L` and put

    u_i = v_L(P_(n-i)),   w_k = v_L(Q_(m-k))    in Z_(>=0) u {oo},
    u_0 = a e,            w_0 = b e.                                        (5)

**TTL.**  For every `1 <= j <= n+m-3`: the minimum of `u_i+w_k` over `i+k=j`
is attained **at least twice**, or its unique minimiser `(i,k)` satisfies the
density identity

    u_i (m-k) = w_k (n-i).                                                  (6)

*Proof.*  By (4) each term of `(G_j)` has `v_L >= u_i+w_k-1` with equality off
(6).  A strict unique minimiser would make one term's valuation strictly smaller
than every other's, contradicting `sum = 0`.  QED.

VERIFIED-EXACT on 30 genuine Jacobian pairs (random compositions of triangular
and linear automorphisms, `n+m` from 8 to 34), **every** order, every `L | H`.
The 30 graded identities `(1)` themselves were checked symbolically first.

---

## 3. T1 -- the order-one law, and its exception set  (PROVED + VERIFIED-EXACT)

At `j=1` the only pairs are `(0,1)` and `(1,0)`.

**(i) Exceptions are rare.**  `(0,1)` is density-matched iff
`w_1 = e(m-1)/g = e b - e/g`, which is an integer iff `g | e`, i.e. iff `e=g`,
i.e. iff `H` is a power of the single linear form `L`.  Same for `(1,0)`.
Hence: **if `H` has at least two distinct roots, neither order-one pair can be
exceptional.**

**(ii) The exact law.**  In the non-exceptional branch, and also in the *doubly*
exceptional branch, for every `L | H`

    u_1 - w_1 = (a-b) e_L,   equivalently   pi_H(P_(n-1)) = H^(a-b) pi_H(Q_(m-1))  (7)

where `pi_H(F)=prod_(L|H) L^(v_L(F))`.  (Doubly exceptional: `u_1=e(n-1)/g`,
`w_1=e(m-1)/g`, whose difference is `e(n-m)/g=(a-b)e` -- the law survives.)

**(iii) The divisibility survives every branch.**  In all four branches
`u_1 >= (a-b)e`, hence

    H^(a-b) | P_(n-1)      (a >= b).                                        (8)

(Singly exceptional at `(1,0)` gives `u_1 = ae - e/g >= (a-b)e` since `b>=1/g`.)

VERIFIED-EXACT: 30/30 for (8), 30/30 for (7) on the samples in its branch.

---

## 4. T2 -- the first correction of `P` is the first correction of `Q`  (PROVED)

Write `P_(n-1) = H^(a-b) A_1` (legitimate by (8)); `deg A_1 = m-1`.  Then

    Jac(P_n, Q_(m-1))  = c a H^(a-1) Jac(H, Q_(m-1)),
    Jac(P_(n-1), Q_m)  = c' b H^(b-1) Jac(H^(a-b)A_1, H)
                       = c' b H^(b-1) H^(a-b) Jac(A_1,H)
                       = c' b H^(a-1) Jac(A_1,H),

because `Jac(H^(a-b)A_1,H) = H^(a-b)Jac(A_1,H)`.  So `(G_1)` is `H^(a-1)` times

    Jac( H,  c a Q_(m-1) - c' b A_1 ) = 0,                                  (9)

so the bracket is a constant multiple of `H_0^((m-1)/g_0)` (or `0`).  But
`g_0 | g | m`, so `g_0 | m-1` forces `g_0 = 1`.  Therefore:

> **T2.**  If `H` is *not* a power of a linear form (`K>=2`, `g_0>=2`), then
>
>     P_(n-1) = (c a)/(c' b) * H^(a-b) * Q_(m-1)   exactly.                 (10)
>
> If `g_0 = 1` the only extra freedom is one scalar: `c a Q_(m-1) - c' b A_1
> = lambda H_0^(m-1)`.

Degrees check: `g(a-b) + (m-1) = ga-1 = n-1`.  Statement (10) is much stronger
than the divisibility (8): the whole component `P_(n-1)` is *determined* by
`Q_(m-1)`.

VERIFIED-EXACT: on 22 genuine Jacobian pairs (all of which have `g_0=1`, so it is
the `lambda`-branch that is tested): `H^(a-b) | P_(n-1)` in 22/22, and the residue
`c a Q_(m-1) - c' b A_1` is a scalar multiple of `H_0^(m-1)` in 22/22.  The
`g_0>=2` branch is tested through the kernel computation of section 9, which is
run at `(a,b)=(3,2),(5,2),(4,3),(5,3)` -- i.e. *inside* the counterexample regime.

### 4a'. The degenerate branch `Jac(H, Q_(m-1)) = 0`, resolved

The lane flagged this case.  It is completely settled by (9).  `Jac(H,F)=0` for a
form `F` of degree `D` holds iff `F in C * H_0^(D/g_0)` (and `F=0` if
`g_0 nmid D`).  Apply this with `D=m-1`: since `g_0 | g | m`, `g_0 | m-1` iff
`g_0=1`.  Hence

* if `K>=2` (`g_0>=2`): `Jac(H,Q_(m-1))=0` forces `Q_(m-1)=0`, and then (9) forces
  `A_1=0`, i.e. `P_(n-1)=0`.  **The degenerate branch is exactly "both first
  corrections vanish"** -- there is no intermediate degenerate configuration.
* if `K=1` (`H=c_0L^g`, `g_0=1`): `Jac(H,Q_(m-1))=0` means `Q_(m-1)=nu L^(m-1)`,
  and (9) then gives `A_1 = (ca/(c'b)) nu L^(m-1) - (lambda/(c'b)) L^(m-1)`, again
  a multiple of `L^(m-1)`.  Nothing breaks; the `lambda` scalar absorbs it.

The same argument at order `j` (with `D=j(m-1)`) shows the degenerate branch
there is governed by `g_0 | j` -- the same divisibility as everywhere else.

### 4a. The reduced tower `(R_j)` -- the next orders, explicitly

The lane asked for the degree `n+m-4`, `n+m-5`, ... conditions.  They are the
*reduced* forms of `(G_j)`.  Let `a > b`.  Divide `(G_j)` by `H^(a-(j-1)b-1)`;
the exponent is `>= 0` exactly for `j <= floor((a-1)/b)+1`.  Every term survives
as a genuine polynomial: for `i+k=j`,

    k >= 1 :  v_L >= u_i + w_k - 1 >= (a-ib)e - 1 >= (a-(j-1)b)e - 1
                  >= (a-(j-1)b-1) e,                            (uses e>=1)
    k = 0  :  Jac(P_(n-j),Q_m) = c' b H^(b-1) Jac(P_(n-j),H),  so
              v_L >= u_j + be - 1 >= (a-(j-1)b)e - 1 >= (a-(j-1)b-1) e,

using the chain (11) for `u_i` where `a-ib>=0` and `u_i>=0` otherwise (in the
last line, at `j = floor(a/b)+1`, `u_j>=0` and `be-1 >= (r-1)e` with
`r=a mod b <= b-1`).  The result:

> **(R_j).**  For `1 <= j <= floor((a-1)/b) + 1`, all terms of
>
>     sum_(i+k=j)  Jac(P_(n-i), Q_(m-k)) / H^(a-(j-1)b-1)  =  0            (R_j)
>
> are forms of degree `j(m-1) + g - 2`.

The range of validity of `(R_j)` is exactly the **first partial quotient of
`a/b`** (plus one when `b` does not divide `a`, which is always the case for a
counterexample candidate, where `b>=2` and `gcd(a,b)=1`).  The reduced tower
runs for as long as the first subtractive Euclid run does; then the divisor
exponent goes negative and the roles swap.

Written out, `(R_1)` is (9); `(R_2)`, the degree `n+m-4` condition, is

    c a H^b Jac(H, Q_(m-2))
      + Jac(P_(n-1),Q_(m-1)) / H^(a-b-1)
      + c' b Jac(P_(n-2),H) / H^(a-2b)                     = 0,           (R_2)

with `Jac(P_(n-2),H)/H^(a-2b)` read as `H^(2b-a)Jac(P_(n-2),H)` when `a<2b`.

**Substituting T2** (`g_0>=2`, so `A_1 = kappa Q_(m-1)` with
`kappa = ca/(c'b)`) and using `Q Jac(H,Q) = (1/2) Jac(H,Q^2)` and
`H^s Jac(H,F) = Jac(H,H^s F)` collapses `(R_2)` to a single Jacobian:

    Jac( H, Theta_2 ) = 0,
    Theta_2 = c a H^b Q_(m-2) + (kappa(a-b)/2) Q_(m-1)^2
              - c' b H^(2b-a) P_(n-2),      deg Theta_2 = 2(m-1).         (17)

So `Theta_2 = mu H_0^(2(m-1)/g_0)`, and `g_0 | 2(m-1)` with `g_0 | m` forces
`g_0 | 2`.  **If `g_0 >= 3` then `mu = 0` and `P_(n-2)` is completely determined
by `Q_(m-1)` and `Q_(m-2)`.**

The pattern continues: at order `j` the collapsed bracket `Theta_j` has degree
`j(m-1)`, and `g_0 | j(m-1)` with `g_0 | m` is `g_0 | j`.  So the free scalar
`lambda_j` exists **iff `g_0 | j`** -- which is exactly the `[g_0 | j]` term of
the dimension law (16), derived independently in section 9.  The two derivations
agree.

VERIFIED-EXACT: `(R_j)` holds as an identity, with the predicted degree
`j(m-1)+g-2`, on 18 genuine Jacobian pairs for every `j` in its range, including
the exactness of every division by `H^(a-(j-1)b-1)` (`PART 3c`).  Scope: those 18
pairs realise `(a,b) in {(2,1),(3,1)}` with `g in {2,3,4,6}` and `jmax_R in {2,3}`
-- so `(R_1)`,`(R_2)`,`(R_3)` are all exercised, but only on the `b=1` stratum,
because that is the only stratum where genuine pairs exist (see section 12.2).
The collapse (17) and its `g_0 | j` conclusion are proofs, not samples, and are
independently confirmed by the dimension law of section 9, which *is* computed at
`a,b >= 2`.

---

## 5. The Euclid chain -- the lane's claim, made precise  (PROVED + VERIFIED-EXACT)

> **THEOREM (Euclid chain).**  Let `a >= b`.  For every `L | H` with multiplicity
> `e` and every `0 <= i <= min( floor(a/b), n+m-3 )`,
>
>     u_i >= (a - i b) e,      equivalently     H^(a-ib) | P_(n-i).         (11)

*Proof (induction on `i`).*  True at `i=0`.  Suppose it holds up to `i` and
`u_(i+1) < (a-(i+1)b)e`.  In `(G_(i+1))` every pair `(i',k)` with `i'<=i`, `k>=1`
has `u_(i') + w_k >= (a-i'b)e + 0 >= (a-ib)e > u_(i+1)+be = u_(i+1)+w_0`, so
`(i+1,0)` is the strict unique minimiser.  TTL then forces (6):
`u_(i+1) m = w_0 (n-i-1)`, i.e. `u_(i+1) = e(ga-i-1)/g = ae - e(i+1)/g`.  Since
`1/g <= b`, this is `>= (a-(i+1)b)e`, contradicting the assumption.  QED.

**This is exactly the subtractive Euclidean algorithm on the exponent pair.**
The exponent of `H` dividing `P_(n-i)` is `>= a - ib`: one subtraction of `b` per
unit step down the tower, `(a,b) -> (a-b,b) -> (a-2b,b) -> ...`, for
`floor(a/b)` steps, i.e. exactly the first partial quotient of `a/b`.  The
lane's `(L1)` claim is **CONFIRMED in this divisibility form**.

### 5a. REFUTED: the *equality* form of the claim

The equality `u_i = (a-ib)e` is **FALSE**.  Exact witness (sample `[02]` of
`PART 2`, a genuine automorphism): `n=12`, `m=4`, `g=4`, `a=3`, `b=1`,
`H=(2x+y)^4`, `e=4`, `L=2x+y`:

    u = [12, 11, 10, 8, 7, 6, 4, 3, 3, 0, 0, 0, 0]
    w = [ 4,  3, oo, 0, 0]
    predicted (a-ib)e = [12, 8, 4, 0]                      -- u_1 = 11 != 8.

What *is* exact here is T1: `u_1 - w_1 = 11-3 = 8 = (a-b)e`.  So the correct
statement of the tower is: **the divisibility chain (11) is an inequality, and
the exact bookkeeping is the difference law (7), not a formula for `u_i` alone.**
Note also that `w` is neither monotone nor finite (`w_2 = oo`, i.e.
`Q_(m-2) = 0`) -- degenerate components occur in genuine pairs and every proof
above is written to survive them.

### 5b. The forced swap

`gcd(a,b)=1` and `b>=2` give `r = a mod b >= 1`.  Running (11) to `i=q=floor(a/b)`
leaves `u_q >= re` with `0 < r < b`: `P_(n-q)` now carries *less* `L` than `Q_m`
does, and the roles of the two towers exchange.  Section 6 shows the exchange is
forced, and section 8 shows what it costs.

---

## 6. The equal-slope law -- how the swap is controlled  (PROVED + VERIFIED-EXACT)

Put `D_i = u_0-u_i`, `E_k = w_0-w_k`, and let

    sigma_P = max_(i>=1) D_i/i,    sigma_Q = max_(k>=1) E_k/k                (12)

be the first (steepest) descent rates of the lower convex hulls of `(i,u_i)` and
`(k,w_k)`.

> **THEOREM (equal-slope law).**  Let `L | H`.  If `L` does not divide `Q_1` then
> `sigma_P <= sigma_Q`.  If `L` does not divide `P_1` then `sigma_Q <= sigma_P`.
> Since `Jac(P_1,Q_1)=1`, at most two of the `K` distinct roots of `H` divide one
> of `P_1,Q_1`; for every other `L`,
>
>     sigma_P(L) = sigma_Q(L).                                              (13)

*Proof.*  Suppose `sigma_P > sigma_Q` and let `i*` be the least index realising
the max in (12).  In scope, `u_(n-1)<=1` (from `Jac(P_1,Q_1)=1`), so
`D_(n-1)/(n-1) >= (ae-1)/(ga-1) > e/g = D_n/n`; hence `i* <= n-1 <= n+m-3`
whenever `m>=2`, and the order `j=i*` really is one of the identities (1).  At
order `j=i*`, for every `k>=1`,
`u_i+w_k >= u_0+w_0-i sigma_P-k sigma_Q > u_0+w_0-i* sigma_P = u_(i*)+w_0`,
so `(i*,0)` is the strict unique minimiser; TTL forces `u_(i*)=ae-e i*/g`, i.e.
`sigma_P=e/g`.  But `L` not dividing `Q_1` means `w_(m-1)=0`, so
`sigma_Q >= be/(m-1) = be/(gb-1) > e/g`, a contradiction.  QED.

VERIFIED-EXACT: 29 in-scope checks, all equal; and in fact all 30 samples
(including the two out-of-scope ones) satisfy (13).

**Corollary (second proof of the chain).**  `E_k <= be` always, so
`sigma_Q <= be`; by (13) `sigma_P <= be`, so `D_i <= i b e` and `u_i >= (a-ib)e`
for *all* `i` -- no induction needed, in the scope of (13).

---

## 7. Hull dichotomy  (PROVED + VERIFIED-EXACT)

Let `A`, `B` be the full lower convex hulls of `(i,u_i)` and `(k,w_k)`.

> **THEOREM.**  Let `s` be a descent rate of `A` that is not a descent rate of
> `B`, let `p, p+alpha` be the endpoints of the corresponding `A`-edge and let
> `kappa` be the `B`-vertex where `B`'s rate crosses `s`.  Then TTL forces the
> density identity (6) at `(p,kappa)` and at `(p+alpha,kappa)`.

*Proof.*  At `j=p+kappa` the neighbours `(p-1,kappa+1)` and `(p+1,kappa-1)` are
strictly larger (their rate differences are `s'-s''>0` and `s'''-s>0`), and every
other point lies above the hull; so the minimiser is unique and is a pair of hull
vertices, where `u`, `w` equal `A`, `B`.  Apply TTL.  Same at `p+alpha`.  QED.

So: **the two hulls have the same slope set, or every unmatched slope sits over
an exact density coincidence.**  VERIFIED-EXACT: 30/30 dichotomy checks pass, and
in these samples the slope sets are in fact literally equal 30/30.

---

## 8. Single-slope arithmetic: `sigma=e/t` with `t <= g-1`  (PROVED + VERIFIED-EXACT)

Suppose both `L`-hulls are single edges with the common rate `sigma` of (13), and
`L` divides neither `P_1` nor `Q_1` (so both hulls reach `0`, at `i=alpha<=n-1`
and `k=beta<=m-1`).  Then

    sigma = a e/alpha = b e/beta   =>   a beta = b alpha   =>   (gcd(a,b)=1)
    alpha = a t,  beta = b t,  sigma = e/t,  t in Z_(>=1),
    a t <= ga-1  =>  1 <= t <= g-1.                                         (14)

> **Consequence.**  The single-slope configuration is IMPOSSIBLE when `g=1`.

VERIFIED-EXACT over all `5304` triples with `g<=12`, `2<=a,b<=12` coprime,
`1<=e<=g`: the solution set of `a beta = b alpha` inside the box is *exactly*
`{(at,bt) : 1<=t<=g-1}`, and is empty for `g=1`.

---

## 9. The obstruction map and the dimension law -- the sharpest criterion
   (PROVED + VERIFIED-EXACT, and verified *inside* the counterexample regime)

Write `S_d` for the space of binary forms of degree `d` (`dim = d+1`).  `(G_j)`
is a linear equation for the two new components:

    Phi_j : S_(n-j) (+) S_(m-j) -> S_(n+m-2-j),
    Phi_j(A,B) = Jac(A,Q_m) + Jac(P_n,B)
               = c' b H^(b-1) Jac(A,H) + c a H^(a-1) Jac(H,B),             (15)
    Phi_j(P_(n-j), Q_(m-j)) = Psi_j := - sum_(i+k=j, i,k>=1) Jac(P_(n-i),Q_(m-k)).

> **THEOREM (dimension law).**  For `a >= b` and `1 <= j <= n`,
>
>     ker Phi_j = { (A,B) :  bA - a H^(a-b) B  in  C * H_0^((n-j)/g_0) },
>     dim ker  Phi_j = max(0, m-j+1) + [ g_0 | j ],
>     dim coker Phi_j = (m-2) + [ g_0 | j ].                                (16)

*Proof.*  `H^(a-b)Jac(B,H)=Jac(H^(a-b)B,H)`, so `Phi_j(A,B)=0` iff
`Jac(bA - aH^(a-b)B, H)=0` iff `bA-aH^(a-b)B` is a scalar multiple of
`H_0^((n-j)/g_0)`, which requires `g_0 | n-j`.  Since `g_0 | g | n`, that is
`g_0 | j`.  Given `B` (and `lambda`), `A` is determined, giving the kernel
dimension; `rank = (n-j+1)+(m-j+1)-dim ker`, and `dim coker = (n+m-1-j)-rank`.
QED.

VERIFIED-EXACT: `212` orders across `H in {y^2, y^3, x^2 y, x y, x(x+y),
x y(x+y), (xy)^2, x y(x+y)(x-y)}` (so `g in {2,3,4}`, `g_0 in {1,2,3,4}`) and
`(a,b) in {(3,2),(5,2),(4,3),(5,3),(2,1),(3,1)}` -- **zero** deviations, together
with an explicit check that the described kernel really annihilates `Phi_j`.

### 9a. What (16) says, in words

At each order `j` the tower has room for exactly one binary form of degree `m-j`
(that is, `Q_(m-j)`, freely chosen) **plus one extra scalar `lambda_j`, which
exists if and only if `g_0` divides `j`.**  `P_(n-j)` is then forced.  So:

* `g_0 = 1` (`H` a power of a linear form, `K=1`): a free scalar at **every**
  order -- `n+m-3` of them.  This is the maximally flexible corner, and it is
  where **every** genuine automorphism in the sample lives.
* `H` squarefree (`d=1`, `g_0=g`): free scalars only at `j = g, 2g, ...`, i.e.
  `floor((n+m-3)/g) = a+b-1` of them (for `g>=3`).  A factor `g` less freedom.

---

## 10. (b) The necessary conditions, collected

For a hypothetical `Jac=1` counterexample with leading data `(H,a,b)`,
`a>b>=2`, `gcd(a,b)=1`, `g=deg H`, `H=c_0H_0^d`, `g_0=g/d`, `K` distinct roots:

| # | condition | status |
|---|---|---|
| N1 | `min(a,b) >= 2`, `gcd(a,b)=1` | CITED.  Two separate classical facts: Jung--van der Kulk forces `a=1` or `b=1` for *automorphisms*, and `n\|m` or `m\|n` is a settled stratum of JC(2) itself.  Neither is re-derived here; `min(a,b)=1` for automorphisms is confirmed on 30/30 samples as EVIDENCE only. |
| N2 | `H^(a-ib) \| P_(n-i)` for `0<=i<=min(floor(a/b), n+m-3)` | **PROVED**, sec. 5 |
| N3 | `pi_H(P_(n-1)) = H^(a-b) pi_H(Q_(m-1))` (order-one difference law) | **PROVED**, sec. 3, unconditional when `K>=2` |
| N4 | if `g_0>=2`: `P_(n-1) = (ca/(c'b)) H^(a-b) Q_(m-1)` **exactly** | **PROVED**, sec. 4 |
| N5 | `sigma_P(L)=sigma_Q(L)` for every `L\|H` with `L` dividing neither `P_1` nor `Q_1` (at least `K-2` of them) | **PROVED**, sec. 6 |
| N6 | hull slope sets agree, or each unmatched slope carries an exact density coincidence `u_i(m-k)=w_k(n-i)` | **PROVED**, sec. 7 |
| N7 | single-slope hulls force `sigma=e_L/t_L`, `1<=t_L<=g-1`; impossible at `g=1` | **PROVED**, sec. 8 |
| N8 | free scalars in the tower occur exactly at orders `j = 0 mod g_0`; obstruction space has dimension `m-2+[g_0\|j]` at every order | **PROVED + VERIFIED-EXACT**, sec. 9 |
| N9 | the reduced tower `(R_j)` (all terms divided by `H^(a-(j-1)b-1)`) is valid for `1<=j<=floor((a-1)/b)+1`, i.e. for exactly the first partial quotient of `a/b`; if `g_0>=3` then `(R_2)` determines `P_(n-2)` from `Q_(m-1),Q_(m-2)` | **PROVED**, sec. 4a |

**The condition on the multiplicity vector of `H`** is N8 read backwards: the
multiplicity vector `(e_L)` enters the tower **only** through the single integer
`g_0 = g / gcd{e_L}`, and it controls the count of free scalars,
`floor((n+m-3)/g_0)`.  Everything else about the multiplicity vector is invisible
at this level.  The extremes are `H = c L^g` (maximal freedom, and the only
regime that any genuine automorphism realises) and `H` squarefree (freedom
divided by `g`).

**An honest non-criterion.**  Summing (16) over the whole tower gives an exact
"obstruction surplus"

    S(n,m) = sum_j (dim coker - dim ker) = n(m-2) - m(m+1)/2 + (m-1)(m-2)/2 - 1.

`S>0` does **not** exclude anything: the genuine automorphism with `(n,m)=(9,3)`
has `S=3`.  The surplus is a size statistic, not an obstruction; it is recorded
here so that nobody re-derives it and mistakes it for one.  Note also that `S` is
**independent of `g_0`** -- the `[g_0|j]` bump appears identically in `dim ker`
and `dim coker` and cancels in the difference (checked directly against the
term-by-term sum).  So the multiplicity vector changes *where* the freedom sits,
not *how much net over-determination* the tower has.

---

## 11. (c) THE SEARCH ORDER

Candidate degree pairs are `(n,m)=(ga,gb)`.  The *shape* parameter is `(a,b)`;
`g` only scales.  Order by the continued fraction of `a/b`, not by degree.

**Ordering key: `z = max partial quotient` ascending, then `a+b` ascending.**
Rationale, entirely from N2: a partial quotient `q` is a run of `q` forced
divisibilities `H^(a-ib) | P_(n-i)`, so *large* quotients are the *most*
constrained shapes and get excluded first.  The hard corner is the
all-quotients-equal-1 ray -- consecutive **Fibonacci** pairs, the
worst-approximable / golden corner named in PROBLEM-LEDGER section C.  Constant
partial quotient `q` is the `q`-th **metallic** ray, which is the same family
THM-3010 section 3 identifies as attaining maximal Newton-circuit alternation via
the Simson identity.  (That coincidence of *families* is noted; no transfer of
theorems between the two problems is claimed.)

Columns: `z`=max partial quotient (Zaremba's quantity), `dep`=CF length,
`sub`=sum of partial quotients (= total subtractive Euclid steps the tower must
run), `chain`=`floor(a/b)` (= length of the *first* forced divisibility run, N2),
`metal`=`q` if `(a,b)` is a consecutive pair of `a_k=q a_(k-1)+a_(k-2)`,
`fib`=`Y` on the golden ray.  The last three columns give the minimal degree pair
at `g=2` (the floor N7 supports on its own), at `g=9`, and at the smallest `g`
clearing Moh's `max(deg)<=100`.

    #    a    b  a/b CF          z  dep  sub chain metal fib     g=2         g=9         Moh
    -------------------------------------------------------------------------------------------
    1     3    2  [1;2]          2   2    3    1  1    Y      (  6,  4) (  27,  18) g=34 ( 102,  68)
    2     5    2  [2;2]          2   2    4    2  2    -      ( 10,  4) (  45,  18) g=21 ( 105,  42)
    3     5    3  [1;1;2]        2   3    4    1  1    Y      ( 10,  6) (  45,  27) g=21 ( 105,  63)
    4     8    3  [2;1;2]        2   3    5    2  -    -      ( 16,  6) (  72,  27) g=13 ( 104,  39)
    5     7    5  [1;2;2]        2   3    5    1  -    -      ( 14, 10) (  63,  45) g=15 ( 105,  75)
    6     8    5  [1;1;1;2]      2   4    5    1  1    Y      ( 16, 10) (  72,  45) g=13 ( 104,  65)
    7    12    5  [2;2;2]        2   3    6    2  2    -      ( 24, 10) ( 108,  45) g=9  ( 108,  45)
    8    13    5  [2;1;1;2]      2   4    6    2  -    -      ( 26, 10) ( 117,  45) g=9  ( 117,  45)
    9    11    8  [1;2;1;2]      2   4    6    1  -    -      ( 22, 16) (  99,  72) g=10 ( 110,  80)
    10   12    7  [1;1;2;2]      2   4    6    1  -    -      ( 24, 14) ( 108,  63) g=9  ( 108,  63)
    11   13    8  [1;1;1;1;2]    2   5    6    1  1    Y      ( 26, 16) ( 117,  72) g=9  ( 117,  72)
    12   19    7  [2;1;2;2]      2   4    7    2  -    -      ( 38, 14) ( 171,  63) g=9  ( 171,  63)
    13   19    8  [2;2;1;2]      2   4    7    2  -    -      ( 38, 16) ( 171,  72) g=9  ( 171,  72)
    14   17   12  [1;2;2;2]      2   4    7    1  -    -      ( 34, 24) ( 153, 108) g=9  ( 153, 108)
    15   21    8  [2;1;1;1;2]    2   5    7    2  -    -      ( 42, 16) ( 189,  72) g=9  ( 189,  72)
    16   19   11  [1;1;2;1;2]    2   5    7    1  -    -      ( 38, 22) ( 171,  99) g=9  ( 171,  99)
    17   18   13  [1;2;1;1;2]    2   5    7    1  -    -      ( 36, 26) ( 162, 117) g=9  ( 162, 117)
    18   19   12  [1;1;1;2;2]    2   5    7    1  -    -      ( 38, 24) ( 171, 108) g=9  ( 171, 108)
    19   21   13  [1;1;1;1;1;2]  2   6    7    1  1    Y      ( 42, 26) ( 189, 117) g=9  ( 189, 117)
    20   29   12  [2;2;2;2]      2   4    8    2  2    -      ( 58, 24) ( 261, 108) g=9  ( 261, 108)
    21   30   11  [2;1;2;1;2]    2   5    8    2  -    -      ( 60, 22) ( 270,  99) g=9  ( 270,  99)
    22   31   12  [2;1;1;2;2]    2   5    8    2  -    -      ( 62, 24) ( 279, 108) g=9  ( 279, 108)
    23   31   13  [2;2;1;1;2]    2   5    8    2  -    -      ( 62, 26) ( 279, 117) g=9  ( 279, 117)
    24   26   19  [1;2;1;2;2]    2   5    8    1  -    -      ( 52, 38) ( 234, 171) g=9  ( 234, 171)
    25   27   19  [1;2;2;1;2]    2   5    8    1  -    -      ( 54, 38) ( 243, 171) g=9  ( 243, 171)
    26   29   17  [1;1;2;2;2]    2   5    8    1  -    -      ( 58, 34) ( 261, 153) g=9  ( 261, 153)
    27   34   13  [2;1;1;1;1;2]  2   6    8    2  -    -      ( 68, 26) ( 306, 117) g=9  ( 306, 117)
    28   30   19  [1;1;1;2;1;2]  2   6    8    1  -    -      ( 60, 38) ( 270, 171) g=9  ( 270, 171)
    29   31   18  [1;1;2;1;1;2]  2   6    8    1  -    -      ( 62, 36) ( 279, 162) g=9  ( 279, 162)
    30   29   21  [1;2;1;1;1;2]  2   6    8    1  -    -      ( 58, 42) ( 261, 189) g=9  ( 261, 189)
    31   31   19  [1;1;1;1;2;2]  2   6    8    1  -    -      ( 62, 38) ( 279, 171) g=9  ( 279, 171)
    32   34   21  [1;1;1;1;1;1;2] 2  7    8    1  1    Y      ( 68, 42) ( 306, 189) g=9  ( 306, 189)
    33    4    3  [1;3]          3   2    4    1  -    -      (  8,  6) (  36,  27) g=26 ( 104,  78)
    34    7    2  [3;2]          3   2    5    3  -    -      ( 14,  4) (  63,  18) g=15 ( 105,  30)
    35    7    3  [2;3]          3   2    5    2  -    -      ( 14,  6) (  63,  27) g=15 ( 105,  45)
    36    7    4  [1;1;3]        3   3    5    1  -    -      ( 14,  8) (  63,  36) g=15 ( 105,  60)
    37   10    3  [3;3]          3   2    6    3  3    -      ( 20,  6) (  90,  27) g=11 ( 110,  33)
    38   11    3  [3;1;2]        3   3    6    3  -    -      ( 22,  6) (  99,  27) g=10 ( 110,  30)
    39   11    4  [2;1;3]        3   3    6    2  -    -      ( 22,  8) (  99,  36) g=10 ( 110,  40)
    40    9    7  [1;3;2]        3   3    6    1  -    -      ( 18, 14) (  81,  63) g=12 ( 108,  84)
    41   10    7  [1;2;3]        3   3    6    1  -    -      ( 20, 14) (  90,  63) g=11 ( 110,  77)
    42   11    7  [1;1;1;3]      3   4    6    1  -    -      ( 22, 14) (  99,  63) g=10 ( 110,  70)
    43   15    4  [3;1;3]        3   3    7    3  -    -      ( 30,  8) ( 135,  36) g=9  ( 135,  36)
    44   17    5  [3;2;2]        3   3    7    3  -    -      ( 34, 10) ( 153,  45) g=9  ( 153,  45)

**The golden ray** (maximal CF depth per size, every quotient `1`):

    (a,b)     CF                      dep sub   Moh-minimal (n,m)
    (5,3)     [1;1;2]                  3   4    (105,63)  at g=21
    (8,5)     [1;1;1;2]                4   5    (104,65)  at g=13
    (13,8)    [1;1;1;1;2]              5   6    (117,72)  at g=9
    (21,13)   [1;1;1;1;1;2]            6   7    (189,117) at g=9
    (34,21)   [1;1;1;1;1;1;2]          7   8    (306,189) at g=9
    (55,34)   [1;1;1;1;1;1;1;2]        8   9    (495,306) at g=9
    (89,55)   [1;1;1;1;1;1;1;1;2]      9  10    (801,495) at g=9

**Metallic rays** `a_k = q a_(k-1) + a_(k-2)` (the THM-3010 family):

    q=1 golden : (3,2) (5,3) (8,5) (13,8) (21,13) (34,21) (55,34)
    q=2 silver : (5,2) (12,5) (29,12) (70,29) (169,70)
    q=3 bronze : (10,3) (33,10) (109,33)
    q=4        : (17,4) (72,17)

### 11a. Exclusions carried by the table

* Every row already has `min(a,b)>=2` -- the `a=1`/`b=1` stratum is the settled
  `n|m` case (CITED).
* `g=1` is killed for single-slope hulls by N7 (own result, PART 4).  The `g=2`
  column is therefore the smallest degree pair this lane's own results permit;
  it is **not** claimed to be attainable.
* The `g=9` and Moh columns apply the two commonly quoted classical filters --
  `gcd(deg P, deg Q) <= 8` (Appelgate--Onishi / Nakai--Baba) and Moh's
  `max(deg P, deg Q) <= 100`.  **CITED-UNVERIFIED**: these are quoted from
  memory of the literature and were *not* re-derived here; the exact thresholds
  should be checked against the sources before any of them is used to close a
  stratum.  Nothing in sections 1-9 depends on them.
* N4 gives an immediate hard exclusion of *data*, not of `(a,b)`: any candidate
  with `K>=2` and `P_(n-1)` not proportional to `H^(a-b)Q_(m-1)` is dead
  immediately, without computing anything past order one.

---

## 12. Boundaries -- what is NOT claimed

1. **No counterexample is excluded.**  Sections 1-9 are necessary conditions on
   the leading-form data; none of them is shown to be unsatisfiable for any
   `(a,b)` with `a,b>=2`.
2. **The verification set for sections 2-8 is the `min(a,b)=1` stratum.**  No
   genuine Jacobian pair with `a,b>=2` exists to test on -- that is the
   conjecture.  Section 9 escapes this, because `Phi_j` depends only on
   `(H,a,b)`, and it was run at `(3,2),(5,2),(4,3),(5,3)`.
3. **Observations, not theorems.**  All 30 sampled automorphisms have `K=1`
   (leading forms are powers of a *single* linear form) and `min(a,b)=1`.  Both
   are consistent with Jung--van der Kulk but are recorded as EVIDENCE from a
   sample of random triangular/linear compositions with `n+m` in `[8,34]`, not
   as proofs.
4. **The ordering key is a heuristic.**  "Small `z` is the hard corner" follows
   from N2 only in the weak sense that N2's forced run is short there.  It is not
   a theorem that small-`z` shapes are harder to exclude.
5. **No bridge.**  Per MISTAKE-237: nothing here links this tower to NC2,
   GMC(2), VC(4), Jelonek data or Zaremba's conjecture.  The appearance of the
   metallic family in both this table and THM-3010 is a coincidence of *index
   families*, and is flagged as such.
6. The `g=2` column of the table records what this lane's own results permit; it
   is far below every classical filter and should not be read as a search target.

## 12a. The three sharpest open follow-ups this lane opens

1. **Is the second Euclid step forced with the same strength as the first?**
   Section 5 PROVES the first run `H^(a-ib)|P_(n-i)`, `i<=floor(a/b)`, and
   section 5b shows the swap must happen.  What is *not* proved is an analogue of
   (11) for the swapped pair `(b, r)` -- i.e. a statement `w_k >= (b-kr)e` for
   `k <= floor(b/r)`.  Section 6's equal-slope law is the natural tool; the
   missing ingredient is control of the *second* hull slope, not just the first.
   If it goes through, the whole continued fraction of `a/b` becomes a chain of
   forced divisibilities and the table of section 11 acquires real exclusion
   power instead of only a ranking.
2. **`(R_j)` collapse for `j>=3` and the `g_0 | j` scalars.**  Section 4a proves
   the collapse at `j=1,2`.  The dimension law (16) predicts the general pattern
   (`Theta_j` of degree `j(m-1)`, `lambda_j` iff `g_0|j`).  Writing `Theta_j`
   explicitly would give the recursion `P_(n-j) = F_j(Q_(m-1),...,Q_(m-j))` -- a
   formal power-series statement that `P` is a function of `Q` up to `g_0`-torsion.
3. **Can `K>=2` be excluded outright?**  Every genuine automorphism sampled has
   `K=1`.  If a counterexample were forced to have `K=1` (`H=cL^g`), then T2 and
   the reduced tower lose all their teeth simultaneously, which would say the
   tower approach cannot separate counterexamples from automorphisms.  Deciding
   `K` for a hypothetical counterexample is therefore the highest-value question
   this lane raises.

## 13. Reproduction

    python3 04-computation/jc2_tower_depth_laneJ1.py

Seven parts.  PART 1 (Lemma J), PART 2 (graded identities, TTL, T1, Euclid chain,
equal-slope law, hull dichotomy on 30 automorphisms), PART 3 (dimension law, 212
orders, run at `a,b>=2`), PART 3b (T2), PART 3c (reduced tower `(R_j)`), PART 4
(hull arithmetic, 5304 triples) each print `PASS = True`; PART 5 prints the
search-order table.  A final `=== SUMMARY ===` block prints `ALL = True`.

Headline counts from the recorded output:
`250` Lemma-J instances, `0` violations; `30` automorphism samples with all
graded identities exact; `212` obstruction orders with the dimension law exact;
`5304` hull-arithmetic triples exact; `18` reduced-tower samples exact.
