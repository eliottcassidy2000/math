---
id: THM-4042
title: "Prime-sector AP-cover exact clock and holonomic law"
status: >
  PROVED + FINITE-EXACT. For every prime sector count P, the sharp-onset
  owner formula of THM-4033 has explicit eventual winning tracks and an exact
  cyclic coefficient word. Every word period and the minimal phase Pi_P have
  closed prime-adic formulas and can be strictly smaller than
  lcm(1,...,P-1). The sharp common
  desingularizer is (P-1)!C(n,P-1); the cleared deficit is an eventual
  degree-(P-2) quasipolynomial of exact phase Pi_P, giving an explicit
  P-recursive law and a D-finite but nonalgebraic generating function. This
  concerns consecutive AP covers and proves neither AP extremality nor LRC.
source: codex sun-minimal-sturmian / prime-sector clock audit, 2026-08-24
depends_on:
  - THM-4033-prime-sector-ap-cover-eventual-owner-tail
related:
  - THM-4029-lrc14-ap-cover-twelve-owner-rational-tail
  - THM-4035-sixty-clock-separation-and-finite-kakeya-spine
  - THM-4038-ap-deficit-holonomic-sixty-phase-law
  - THM-1420-the-sixty-plateau-in-lrc
  - THM-1105-arithmetic-position-law
  - THM-563-signed-deltaw-periodicity-bound
  - THM-879-vgrid-moment-and-moebius-sinc
  - THM-1226-gcd-period-projective-charge-obstruction
  - THM-3848-rational-base-prefix-atom-tree-and-lonely-runner-separation
  - THM-4000-centered-base-split-cubic-observer-and-tripotent-crt-atlas
  - THM-4044-sixty-clock-hasse-alias-and-planar-jc-boundary-firewall
script: 04-computation/prime_sector_ap_cover_exact_clock_thm4042.py
output: 05-knowledge/results/prime_sector_ap_cover_exact_clock_thm4042.out
script_sha256: 4e4d5f98bf3746f78efa1736933f3456b143466df7c6057be923f20bdb222678
output_sha256: 1d781c506b5dea8a03b077335845c87a02ca5395bfe28c65089f92f961f74803
independent_audit_script: 04-computation/prime_sector_ap_cover_exact_clock_independent_audit_thm4042.py
independent_audit_output: 05-knowledge/results/prime_sector_ap_cover_exact_clock_independent_audit_thm4042.out
independent_audit_report: 05-knowledge/results/prime_sector_ap_cover_exact_clock_independent_audit_thm4042.md
reflection: 07-reflections/six-sixty-420-27720-master-lcm-census-2026-08-24.md
independent_audit_script_sha256: 024e2ced48a508a6903f12c1cb48663061d4d05c83b18f934585aa91782a47d0
independent_audit_output_sha256: a56beb886757a5c0b4b1fa18919326880d0a374c090322154fd910da0b287397
---

# THM-4042 -- exact clocks and holonomic laws for every prime sector count

**PROVED + FINITE-EXACT.** Fix a prime `P`. Let `D_P(m)` be the measure of
the noncover set for the consecutive AP times `0,...,m-1` on `P` sectors, as
in THM-4033, and put `n=m-1`. That theorem proves the exact geometric owner
decomposition from the sharp onset `(P^2+3)/4` for odd `P`. The result here
solves its eventual max-min selectors, determines the genuinely minimal phase,
and gives a uniform algebraic classification of the scalar tail.

## 1. Closed winning tracks

The persistent owners are exactly the reduced `a/q` with `1<=q<P`. For an
owner of denominator `q>=2`, let `a^(-1)` and `P^(-1)` denote inverses modulo
`q`. Work in the circle coordinate `theta=P x`; THM-4033's plus and minus
radii are local `theta`-lengths and are finite max-mins of terms `C/(n-c)`.
Their leading-winner data are

```text
lambda_q^+ = (P-q)/q,      s_q^+(a) = -a^(-1) mod q,
lambda_q^- = (P-q-1)/q,    s_q^-(a) = (1-P^(-1))a^(-1) mod q.           (1)
```

If

```text
E_s(n)=n-((n-s) mod q),
```

then, after a finite selector onset depending only on `P`,

\[
 \rho^+_{a/q}(n)={\lambda_q^+\over E_{s_q^+(a)}(n)},\qquad
 \rho^-_{a/q}(n)={\lambda_q^-\over E_{s_q^-(a)}(n)}.                  \tag{2}
\]

Whenever the displayed coefficient is positive, its leading track is unique.
For `q=P-1`, the minus radius is identically zero. The two `rho` values in
`(2)` remain `theta`-coordinate lengths; the deficit in the original
`x`-circle divides their sum by `P`, producing the normalization in `(3)`.
For `q=1`, the sole owner is zero,
its track is zero, and its normalized plus/minus coefficients are
`(P-1)/P` and `(P-2)/P`.

The finite selector onset need not equal the sharp *geometric* onset in
THM-4033. At `P=7`, THM-4038 separately proves that the closed coefficient
law is already valid from the sharp row `n=12`.

## 2. The cyclic word and exact phase

For `q>=2`, put

```text
U_q=(Z/qZ)^*,             t_(P,q)=P^(-1)-1 mod q,
```

and define the normalized rational word on `Z/qZ`

\[
 w_{P,q}(c)={P-q\over Pq}\mathbf1_{c\in U_q}
 +{P-q-1\over Pq}
   \#\{u\in U_q:t_{P,q}u=c\}.                         \tag{3}
\]

For `q=1`, set `w_(P,1)(0)=(2P-3)/P`. For `0<=c<=P-2`, define

\[
 A_c(r)=\sum_{q=c+1}^{P-1}w_{P,q}(c-r\bmod q).         \tag{4}
\]

Then eventually

\[
 \boxed{D_P(n+1)=\sum_{c=0}^{P-2}{A_c(n)\over n-c}.}  \tag{5}
\]

Define the least positive translation period

\[
 \delta_{P,q}=\min\{d\ge1:w_{P,q}(c+d)=w_{P,q}(c)
                    \text{ for every }c\in\mathbb Z/q\mathbb Z\},
\]

and put

\[
 \boxed{\Pi_P=\operatorname{lcm}_{1\le q<P}\delta_{P,q}.}             \tag{6}
\]

This is the exact minimal period of the coefficient vector in `(4)` and of
the rational phase functions in `(5)`. It divides
`L_P=lcm(1,...,P-1)`, but equality is not universal.

More explicitly, if `2<=q<=P-2` and `q=prod_l l^e_l`, then

\[
 \boxed{\delta_{P,q}=\prod_{\ell^{e_\ell}\parallel q}
 \ell^{\min\{e_\ell,v_\ell(P-1)+1\}}}.              \tag{7a}
\]

The boundary values are `delta_(P,1)=1` and
`delta_(P,P-1)=rad(P-1)`. Consequently, for every prime `P>=5`,

\[
 \boxed{\Pi_P=\prod_{\ell\le P-2}
 \ell^{\min\{\lfloor\log_\ell(P-2)\rfloor,
                    v_\ell(P-1)+1\}}}.              \tag{7b}
\]

The two smaller cases are `Pi_2=1` and `Pi_3=2`.

| `P` | exact `Pi_P` | `L_P` | `L_P/Pi_P` |
|---:|---:|---:|---:|
| 2 | 1 | 1 | 1 |
| 3 | 2 | 2 | 1 |
| 5 | 6 | 12 | 2 |
| 7 | 60 | 60 | 1 |
| 11 | 420 | 2520 | 6 |
| 13 | 27720 | 27720 | 1 |
| 17 | 120120 | 720720 | 6 |

Thus the naive general law `Pi_P=L_P` is **REFUTED** at `P=5`. The exact
60-clock at `P=7` is special, not a generic lcm phenomenon. For the top
denominator,

\[
 w_{P,P-1}(c)={1\over P(P-1)}\mathbf1_{\gcd(c,P-1)=1},\qquad
 \delta_{P,P-1}=\operatorname{rad}(P-1),              \tag{7}
\]

with the convention `rad(1)=1`. The radical collapse in `(7)` is the
top-boundary instance of `(7a)` and one source of the smaller clocks.

### 2a. Master-LCM form and the four-row coincidence

There is a useful equivalent form of `(7b)` which was exposed by the
repository-wide `6,60,420,27720` census.  For prime `P>=5`, put

```text
M_P=lcm(1,...,P-2),                 R_P=rad(M_P).
```

Prime by prime, `(7b)` is exactly

\[
 \boxed{\Pi_P=\gcd(M_P,(P-1)R_P)
              =R_P\gcd(M_P/R_P,P-1).}                \tag{7c}
\]

Consequently

\[
 \boxed{\frac{M_P}{\Pi_P}
 =\prod_{\ell\le P-2}\ell^{\max\{0,
       \lfloor\log_\ell(P-2)\rfloor-v_\ell(P-1)-1\}}
 ={M_P/R_P\over\gcd(M_P/R_P,P-1)}.}                  \tag{7d}
\]

In particular,

```text
Pi_P=M_P  iff  M_P/R_P divides P-1.                  (7e)
```

The scope `P>=5` matters: the separately handled row `P=3` has `Pi_3=2`,
whereas the right side of `(7c)` would be `1`.  An exact incremental audit
of `(7b)--(7e)` for every prime `P<=10000` found equality in `(7e)` only for
`P=5,7,13`; see
`04-computation/six_sixty_420_27720_master_lcm_census_20260824.py`.

At the first four odd-prime rows after `3`, another finite coincidence occurs:

```text
(Pi_5,Pi_7,Pi_11,Pi_13)=(6,60,420,27720)
                         =(lcm(1..3),lcm(1..5),
                           lcm(1..7),lcm(1..11)).       (7f)
```

Thus these four clocks are the initial-segment LCM staircase sampled at the
previous odd prime.  This is **not** a general prime-sector law.  At `P=17`,

```text
M_17=lcm(1..15)=360360,   R_17=30030,
Pi_17=R_17*gcd(12,16)=120120=M_17/3.                  (7g)
```

The previous-prime prediction would have been `lcm(1..13)=360360`.  It loses
one factor `3` because `3^2|M_17` but `3` does not divide `P-1=16`; only the
radical factor survives.  The divisibility chain also stops here:
`27720` does not divide `120120`.  Equations `(7c)--(7g)` classify the
four-row match while firewalling it from extrapolation.

## 3. Constant, sharp desingularizer, and recurrence

Summing the words gives

\[
 \kappa_P={1\over P}\sum_{q=1}^{P-1}
 {\varphi(q)(2P-2q-1)\over q}>0,\qquad
 D_P(n+1)={\kappa_P\over n}+O_P(n^{-2}).              \tag{8}
\]

For `P=7`, this is `127/35`, agreeing with THM-4029/4038. Put

\[
 Q_{P-1}(n)=\prod_{c=0}^{P-2}(n-c)
            =(P-1)!{n\choose P-1},\qquad
 Y_P(n)=Q_{P-1}(n)D_P(n+1).                           \tag{9}
\]

Then `Y_P` is eventually a quasipolynomial of degree exactly `P-2` and exact
phase period `Pi_P`. Moreover, `Q_(P-1)` is the unique monic common
desingularizer of least degree. Consequently

\[
 \boxed{\sum_{j=0}^{P-1}(-1)^{P-1-j}{P-1\choose j}
 Q_{P-1}(n+j\Pi_P)D_P(n+j\Pi_P+1)=0}                 \tag{10}
\]

eventually. This is an explicit polynomial-coefficient recurrence, so the
deficit is P-recursive.

The eventual ordinary generating function is D-finite. Equation `(8)` gives

\[
 \sum_nD_P(n+1)z^n=\kappa_P\log {1\over1-z}+O(1)
 \quad(z\uparrow1).                                  \tag{11}
\]

An algebraic function has a Newton--Puiseux, not logarithmic, local expansion.
Thus the generating function is nonalgebraic.

## 4. Proof of the winner word

Multiplication by `a` permutes the `q` owner tracks, so sort their base orbit
as `y_k=Pk/q`, `0<=k<q`. In the gap ending at `y_(k+1)`, the last missing
sector reached by positive drift has leading cost

\[
 {P\over q}-1-\{y_{k+1}\}.                            \tag{12}
\]

The fractional parts are `0,1/q,...,(q-1)/q`. The unique maximum is
`(P-q)/q`, on the track satisfying `as=-1 mod q`, proving the plus half of
`(1)`.

Under negative drift, the corresponding first missing sector has leading
cost

\[
 \{y_k\}+{P\over q}-2.                                \tag{13}
\]

Its maximum uses `{y_k}=(q-1)/q`, giving `(P-q-1)/q`; the following track is
`(1-P^(-1))a^(-1)`, proving the minus half. Whenever this coefficient is
positive, the leading winner is unique and remains the exact winner after
finitely many comparisons of rational terms. The zero case `q=P-1` was
handled separately above.

On the phase `r=n mod q`, its pole shift is `c=(r-s) mod q`. As `a` runs
through the units, the plus shifts give the first term of `(3)` and the minus
shifts give the multiplication fibre in its second term. Summing denominators
proves `(4)--(5)`, while summing their masses proves `(8)`.

## 5. Proof of minimality and sharp clearing

To prove `(7a)`, write `t=P^(-1)-1=(1-P)/P mod q`. For
`ell^e || q`, put `h=min(v_ell(P-1),e)`. The local projection of `tU_q`
consists exactly of residues of valuation `h` when `h<e`, and of zero modulo
`ell^e` when `h=e`. Its indicator therefore has least additive period

\[
 \ell^{\min\{e,h+1\}}.
\]

CRT makes the least period of `1_(tU_q)` the product in `(7a)`. The fibre
count `#{u in U_q:tu=c}` is constant on `tU_q`. If `t` is a unit, then
`tU_q=U_q` and `(3)` is a multiple of the unit indicator, whose least period
is `rad(q)`, agreeing with `(7a)`. If `t` is not a unit, `U_q` and `tU_q`
are disjoint. Their two positive values in `(3)` are distinct: writing
`s=P-q>=2` and the fibre size as `K>=2`, equality would force
`s=(s-1)K`, hence `(s,K)=(2,2)`, while `s=2` gives
`gcd(P-1,q)=gcd(P-1,P-2)=1` and makes `t` a unit. Thus every word period
preserves both level sets and has exactly the product period in `(7a)`.
For `q=P-1` the second coefficient vanishes, giving `(7)`.

Maximizing each prime exponent over `q<=P-2` proves `(7b)`; for `P>=5`,
every prime divisor of `P-1` is already at most `P-2`, so the top radical
adds no new prime. This completes the closed clock calculation.

Every denominator-`q` block is the cyclic family
`c -> w_(P,q)(c-r)`, so `(6)` is a period. Conversely, the largest pole
coordinate `c=P-2` receives only the `q=P-1` block. Any period of the full
vector must therefore be a multiple of `delta_(P,P-1)`. Subtract that
periodic block and descend through `c=P-3,P-4,...,0`; at coordinate `c=q-1`,
the denominator-`q` block is isolated. Hence every `delta_(P,q)` divides any
full period, proving `(6)`. Uniqueness of partial fractions transfers this
minimality to the phase rational functions and, after multiplication by
`Q_(P-1)`, to the phase polynomials.

For `(7)`, coprimality modulo `P-1` depends only on the residue modulo its
radical, and CRT shows no proper divisor of the radical is a period.

Finally, the plus word at `q=P-1` is positive on a translate of `U_(P-1)`.
As the phase varies, every pole `c=0,...,P-2` is therefore active on some
phase. Its residue in `(5)` is positive and cannot cancel against a distinct
pole. Any polynomial clearing every phase must vanish at all `P-1` pole
locations, hence is divisible by `Q_(P-1)`. Equation `(9)` shows that this
polynomial works. Degree, recurrence `(10)`, and the generating-function
claims follow.

## 6. Fibonacci, triangular, and scope firewalls

At `P=7`, THM-4035 proves that full Fibonacci state modulo ten and pointed
triangular state modulo thirty are reversible addresses for the exact
60-cycle. THM-4038 proves that the Fibonacci owner denominators `1,2,3,5`
produce only period 30 and that `q=4` restores period 60. Equations `(3)` and
`(6)` identify the general mechanism: cyclic owner words can have proper
periods, and their exact synchronization—not an external recurrence—sets the
clock.

The theorem requires prime `P`. At composite sector counts, owner tracks can
start on extra walls because `gcd(P,q)>1`; negative drift may delete a base
occupied sector. THM-4033's modulus-four hostile therefore blocks a naive
composite extension.

This theorem concerns consecutive AP times only. It does not show that an AP
maximizes cover measure among sparse time sets and does not prove LRC(14).

## 7. Exact companion

Run

```text
python -B 04-computation/prime_sector_ap_cover_exact_clock_thm4042.py
python -B -O 04-computation/prime_sector_ap_cover_exact_clock_thm4042.py
python -B 04-computation/prime_sector_ap_cover_exact_clock_independent_audit_thm4042.py
python -B -O 04-computation/prime_sector_ap_cover_exact_clock_independent_audit_thm4042.py
```

The exact rational audit checks `(1)` through prime `43`, reconstructs every
denominator block and its minimal period through `P=13`, checks all active
poles, and computes the closed period profile through `P=61`. Normal and
optimized runs are byte-identical. The proof above establishes the all-prime
claims; the computation is a finite hostile control. The independent script
rebuilds the selectors and cyclic words without importing the target, checks
546 winner formulas, 35 raw word blocks, and all 6,027 p-adic word periods
through prime `251`, verifies the recurrence, and retains the composite `P=4`
hostile. **QED.**
