# Independent audit: native ROOT quotient domain, boundary invoice, and ranked losses

2026-10-27. Scope: the Session 75 native ROOT packet, including Phase 1 §4a.
This is an independent proof review and independently written exact checker,
not an audit of the Session 73 fully odd map or its measures.

## Verdict

**PROOF REVIEW: PASS, with explicit convention qualifications.** The all-cutoff
quotient obstruction, first-hit mode, source-labelled initial-energy invoice,
ranked polynomial/logarithmic bounds, and marked C3 identities are correct.
The finite 627-frontier coefficient is checked independently below.

The two important qualifications are:

1. `F=0` for the signed contraction `e'=epsilon*b*e+F` is not the same condition
   as `f=0` in the inherited positive-b contraction `e'=b*e+f`. The former gives
   terminal `epsilon_path*B*A`; the latter gives `B*A`. Gauge equivalence of
   energy identities must not identify these forcing variables.
2. A global homogeneous profile with one fixed ROOT value `z` and a family of
   source-labelled transports with terminal values `epsilon_s*A_N(s)` are
   different objects. The latter invoice is correct on separate labelled paths;
   it does not construct one common ROOT boundary value for all those sources.

Neither issue invalidates the producer's displayed formulas when their stated
signed, source-labelled convention is retained. The recommended exposition is
to make both distinctions explicit before any further theorem promotion.

**FINITE-EXACT replay status: PASS.** Normal and `python3 -O` outputs are
byte-identical; their embedded checker hash matches the delivered checker.
The audit authenticates 9,000 first-hit paths (maximum depth 102, clock 176;
maximum exterior depth 67) and all 627 frontier ranks. The independent proof
review, rather than the finite scan, establishes the all-scale quantified bounds.

## 1. Direct source authentication and lineage

The audit read these blobs using `git show` from repository commit
`96e35642110104cc547516a63f37eebcf930f707`, independently of the producer's copies:

- `05-knowledge/results/paid_root_cycle_budget_20261026/README.md`, blob
  `f00fd7484371260186464a52007b856a10aa1b6a`:
  lines 44–105 define the native carry and both the positive-b and sign-gauged
  payment conventions; lines 125–162 keep the original charges and moving
  prefixes; lines 164–220 define the original dyadic debt and sidecars.
- `05-knowledge/results/root_eisenstein_first_arrival_20261018j.md`, blob
  `2da8811e64696f018835744846b7558f04f4d87e`:
  lines 44–98 give the literal original rows and Eisenstein charge; lines
  100–148 give first-arrival signs, clocks, and the marked-row energy identity.
- `05-knowledge/results/paid_root_cycle_budget_20261026/frontier.txt`, blob
  `0ab0d566e6eca73a3cbb0275182501631a79204d`, supplies the 627 labelled seeds
  and their integer ranks. Membership, closure, and exact ranks are recomputed,
  not trusted from the file's comment.

`SOURCE_PINS.json` records independent SHA-256 checks. No producer Python code
is imported. The audit's source compiler traverses all literal odd `m<4N`,
peels each to its primitive source, and accumulates its original three rows.
The native edge implementation instead uses the 2-adic valuation of `3q+1`
after initial division, then checks it against literal B-chain peeling.

The claimed historical sentence about a finite-norm all-source boundary mode
has not been located in an independently supplied source defining its Hilbert
weight and operator domain. It remains **UNLOCATED TRANSCRIPT EVIDENCE**, not a
proved premise. The verified source definitions above suffice for every scoped
result reviewed here. No S73 probability or occupation measure is substituted.

## 2. All N>=2 quotient obstruction — independently PROVED

Let `S` be the positive primitive odd seeds, namely residues 1, 3, 7 modulo 8.
For `X_s=3s+1`, the actual map has

```
a=2 for s=1 mod8, else a=1;
(3s+1)/2^a=B^j(P(s)),   B(x)=4x+1;
t=a+2j,                2^t X_(P(s))=3X_s+2^a.
```

The source charge from the pinned original source is

```
A_N(s)=sum_{B^j(s)<4N} omega^(B^j(s))
          * [(omega-omega^2) if B^j(s)<2N, else omega],
mu_N(s)=Q(A_N(s)).
```

In particular `mu_N(s)=0` for `s>=4N`. Choose

```
s_N=4N-1 when N mod4 is 0,2,3;
s_N=4N-5 when N mod4 is 1, N>=5.
```

The congruence verification is explicit:

| N mod 4 | s_N mod 8 | P(s_N) | P(s_N) mod 8 |
|---|---:|---|---:|
| 0 | 7 | 6N-1 | 7 |
| 2 | 7 | 6N-1 | 3 |
| 3 | 3 | 6N-1 | 1 |
| 1, N>=5 | 7 | 6N-7 | 7 |

Thus each source has `a=t=1`, with no B-chain stripping at its target. It lies
in `[2N,4N)`, and its next B-chain point exceeds the cutoff. Consequently
`A_N(s_N)=omega^(s_N+1)` and `mu_N(s_N)=1`. For the ordinary branch,
`P(s_N)-4N=2N-1>0`; for the exceptional branch it is `2N-7>0` because `N>=5`.
Hence the target's mass is zero.

The smallest controls are `(N,s_N,P(s_N))=(2,7,11),(3,11,17),(4,15,23),(5,15,23)`.
The two literal functions `0` and `delta_(P(s_N))` represent the same zero
class in the cutoff weighted quotient. Their pullbacks differ by 1 at `s_N`,
a coordinate of weight 1. So neither `C e=e composed with P` nor
`F e=e composed with P-c e` descends to the quotient, for either `c=b` or
`c=epsilon*b`. This is a well-definedness failure, not a divergent sum.

At N=1, `A_1(1)=-2-omega`, `A_1(3)=omega`, so positive support is exactly
`{1,3}`. Both seeds map to 1. This boundary case is forward closed, and the
threshold N>=2 is sharp.

The general null-class criterion is forward closure of positive support:
`w(s)>0 => w(P(s))>0`. If it holds, the composition norm squared is governed by
`sup_y sum_{P(s)=y}w(s)/w(y)`. Zero extension of representatives gives a valid
*killed/compressed* finite operator, but it is a chosen boundary convention,
not descent of the original ambient native pullback. Restricting to a declared
boundary extension or changing the carrier/weights is a legitimate repair;
none is supplied merely by observing that a finite sum is finite.

## 3. First-hit mode and original data — independently PROVED

For a genuine first-hit word `s=s_0,...,s_d=1`, define

```
b_i=3X_i/(3X_i+2^a_i),  c_i=epsilon_i*b_i,
B_s=product_i b_i,      epsilon_s=(-1)^sum_i t_i.
```

The affine edge identity telescopes exactly to

```
B_s=3^d X_s/(4*2^T).
```

The unique ROOT-conditioned homogeneous solution with `h(1)=z` is therefore
`h(s)=epsilon_s*z/B_s`. Defining it only on the ROOT basin, then extending by
zero off that basin, satisfies every interior equation: `P(s)` lies in the
basin iff `s` does. Such a definition cannot prove that another source arrives.

ROOT's actual self-edge is `P(1)=1`, `t=2`, `b_1=3/4`. Retaining it in the
homogeneous equation forces `z=(3/4)z`, hence `z=0`. A nonzero boundary value
requires deleting that edge or its forcing `F=z/4`. Since `d_1=7/16`, that
forcing costs `Q(z)/7` in the two-port identity. ROOT's first-hit empty path
has `B_1=1`; it must not include another factor 3/4.

For any supplied strictly positive weight w, the ROOT-supported profile has
finite squared norm iff `Q(z)*sum_{s in R} w(s)/B_s^2` is finite. No fixed
all-source native w has been established by the reviewed source text.
Counting measure already fails: the distinct primitive seeds

```
r_k=(32*64^k-5)/9,
P(r_k)=1, t(r_k)=3+6k,
b(r_k)=1-1/(16*64^k)
```

have `Q(h(r_k))=Q(z)/b(r_k)^2>Q(z)`. This proves absence of a nonzero
ROOT-boundary profile in counting ell2; it is not an assertion about arbitrary
nonroot graph components or about universal arrival.

For signed path evolution `e_(i+1)=c_i e_i+F_i`, zero F and unchanged original
input A yield terminal `epsilon_s*B_s*A`. To attain the original prescribed
terminal `epsilon_s*A`, the starting value must instead be `A/B_s`. Exact
second-port telescoping then gives

```
initial energy  = Q(A)/B_s^2,
terminal energy = Q(A),
dissipation     = Q(A)*(B_s^(-2)-1).
```

For the *literal positive-b* inherited evolution `e'=b e+f`, zero f instead
yields `B_s*A`; attaining `epsilon_s*A` requires initial `epsilon_s*A/B_s`.
The energy invoice is unchanged, but the input and forcing conventions are not.
At s=3, `b=15/16`, epsilon=-1. Preserving the physical charge -A from input A
costs `31 Q(A)` in the positive-b equation, versus `Q(A)/31` in the signed-c
or path-sign-gauged equation. These are not the same forcing variable.

Explicitly, for the same physical e,e', `f=F+(epsilon-1)b e`. A genuine path
sign gauge `g_i=sigma_i e_i`, `sigma_(i+1)=sigma_i epsilon_i`, transforms the
signed equation into `g_(i+1)=b_i g_i+sigma_(i+1)F_i`. This explains the
energy equivalence without asserting `f=F`.

Finally, one fixed global boundary z supplies terminal z for every path.
The source-dependent targets `epsilon_s*A_N(s)` require labelled copies (or an
explicit source tensor) with those terminal values. This is precisely the
setting in which the displayed sourcewise invoice is valid.

## 4. Ranked losses and frontier constants — independently PROVED/FINITE-EXACT

### 4.1 Rank formula and global depth bound

Write `q_s=1/lambda_s`, `lambda_s=(1-b_s)/(1+b_s)`. The nonroot primitive classes
have denominators

```
s=8j+1, j>=1:    q=36j+7 = 43,79,115,...;
s=4j+3, j>=0:    q=36j+31 = 31,67,103,... .
```

They are disjoint, interlace, and give
`q_k=18k+10+3(-1)^(k+1)` for k>=1. ROOT alone would contribute q=7 and is
correctly excluded. Since `(q+1)/(q-1)` decreases for q>1, selecting h distinct
nonroot sites maximizes reciprocal gain at the first h denominators. Moreover

```
log product_{k=1}^h (q_k+1)/(q_k-1)
 <= sum_{k=1}^h 2/(q_k-1)
 <= sum_{k=1}^h 2/(18k+6)
 <= integral_0^h 2/(18x+6) dx
  = (1/9) log(1+3h).
```

All inequality directions are correct: `log(1+x)<=x`, and a decreasing positive
integrand on `[k-1,k]` is at least its right-endpoint value. Squaring the
reciprocal gain proves `B_ROOT^(-2)<=(1+3d)^(2/9)` on first-hit words.

Nonrepetition is proved, not assumed, for a successful first-hit word. If a
state repeats before the hit, determinism makes the intervening segment a
cycle and prevents a later first hit. The same argument applies to successful
first entry into a forward-closed frontier. Neither argument proves that an
arbitrary source succeeds.

### 4.2 Exterior constants and arbitrary repeats

All primitive seeds below 1009 are authenticated members of V, while 1009 is
outside. The minimum exterior denominator is 4543, attained at 1009. The full
ranked denominator sequence starts its tail there at global index 252:
`4543,4567,4579,4603,...`. Its kth element is at least
`4543+18(k-1)` because this tail starts with gap 24, followed by gap 12.
Exterior denominators form a subset of this tail; including interior sites in
the comparison only weakens the upper bounds.

Consequently, for h distinct exterior sources,

```
log B_ext^(-1)
 <= sum_{k=1}^h 2/(4524+18k)
 <= (1/9) log(1+18h/4524)
  = (1/9) log(1+3h/754),

sum lambda
 <= sum_{k=1}^h 1/(4525+18k)
 <= (1/18) log(1+18h/4525).
```

This gives exactly the proposed energy exponent 2/9 and denominators 754 and
4525. Here h is **native edge count**, not the total integer braid clock T.
The h=0 empty prefix gives equality at 1 for inflation and 0 for forcing.

For an arbitrary retained prefix, group visits by their m distinct sites and
bound every visit count by M. Since all losses are nonnegative and reciprocal
gains exceed 1,

```
B_ext^(-2)<=(1+3m/754)^(2M/9),
sum lambda<=(M/18)log(1+18m/4525).
```

The convention m=M=0 handles the empty prefix. A stopped but not successful
exterior prefix need not be declared nonrepeating; retaining M is necessary.
When weighted and summed, every factor remains attached to the original
source charge mu_N(s). No arbitrary-source arrival or uniform-depth bill is
obtained by deleting M or by these depth-dependent estimates.

### 4.3 Finite internal factor

Every supplied frontier label is primitive. Independent literal traversals
verify forward closure, one-step rank decrease by exactly 1, and equality of
each file rank with its first-ROOT path length. The root rank is 0; maximum
rank is 64. The exact internal maximum is uniquely at seed 993, rank 31:

```
B_993=460166680231540515/576460752303423488,
C_V=332306998946228968225951765070086144
    /211753373595316860199691410086465225 < 157/100.
```

Thus first entry into V after h exterior edges gives
`B_s^(-2)<=C_V*(1+3h/754)^(2/9)`. The coefficient is FINITE-EXACT, while the
implication from authenticated closure and the exterior inequality applies to
any actual successful first-entry word. It is not a uniform source energy or
coherent-merging estimate.

## 5. Marked C3 carrier — independently PROVED

For `E(a+b omega)=(2a-b,-a-b,-a+2b)`, let `A2` be the integral centered plane.
The inverse coordinates are `(a,b)=((x0-x1)/3,-(x0+2x1)/3)`. Therefore

```
E(Z[omega])=L={x in A2: x0=x1=x2 mod3},
[A2:L]=3,    ||E(z)||^2=6Q(z).
```

With `R(x0,x1,x2)=(x1,x2,x0)` and `K=R-R^2`, direct calculation gives

```
E(omega z)=R E(z),
E((1+2omega)z)=K E(z),
K^2=-3I on the centered plane,
Q((1+2omega)z)=3Q(z).
```

A useful additional carrier fact is **K(L)=3A2**, hence `[L:K(L)]=3`. Indeed
`K E(a,b)=(-3b,-3a+3b,3a)` ranges over exactly 3A2. Thus `K^-1=-K/3` is a valid
real-plane inverse and integral on the transformed carrier 3A2, but K is not
an integral automorphism of the original L.

The original readout stays marked: `ell_0(z)=2a-b`, with
`ell_0(z)^2+3b^2=4Q(z)`. Rotating the charge without conjugating the observer
changes it (`ell_0(1)=2`, `ell_0(omega)=-1`). On the transformed plane charge x,
the original observer is `(x2-x1)/3`, not x0.

Bilinear polarization scales by 3 under `1+2omega`. Consequently every energy,
every signed dyadic level debt, and the positive part taken *after summing a
whole level* scale by 3. Original labels, terminals, signs, clocks, and first-
differing-bit pair owners are unchanged. Source norm scales by 3 too, so the
normalized debt is unchanged. This is a transport of the ledger, not payment.

The strengthened native hostile is independently authenticated at N=512:
`A(799)=-2-omega`, `A(1065)=omega`, seed 799 lies in V, seed 1065 does not,
and `P(1065)=799` with clock 2. Thus both actual signed charges have the same
first-entry terminal. Their Eisenstein cross term is 0, their marked cross
term is 6, and their odd-index difference is 133, so their dyadic owner is q=1.
This does not preclude cancellation by other pairs at that level.

## 6. Independent exact tests and stopping boundary

The portable checker uses only Python's standard library, integer arithmetic,
and `Fraction`; no floating-point approximation or producer-code import is
used. Validation uses explicit exceptions, not assertions removed by `-O`.

```
python3 check.py > normal.json
python3 -O check.py > optimized.json
cmp normal.json optimized.json
```

Declared finite universes:

- Literal odd inputs for every N=1..256, covering 49,408 primitive source cases;
  independently reconstruct all three source rows and compare B-chain charges.
- Every cutoff hostile N=2..20000, plus all four residue classes around bases
  `2^60`, `3^100`, and `10^150`; N=1 and N=2,3,4,5 checked separately.
- All 9,000 primitive sources s<=24000, authenticated to first ROOT, including
  every path's signed attenuation, inflated original charge, exact invoice,
  global polynomial bound, exterior bound, and combined frontier bound.
- First 32 direct-ROOT inverse-family points; first 10,000 sorted q values;
  first 512 extremal ranked products; rational log controls through h=256.
- All 627 frontier paths walked independently, exact C_V and unique maximizer.
- Synthetic repeated exterior visit multisets test the m,M algebra without
  claiming those multisets are actual native cycles.
- Integer C3 grids and all 6,561 ordered pairs in [-4,4]^2; nine authenticated
  all-the-way-to-ROOT source ledgers at N=1,2,3,4,5,17,64,129,256; and the
  separate first-entry-to-V marked hostile at N=512.

The finite ROOT ledgers intentionally stop at ROOT, not the 627-frontier, and
are used only to test transformation and pair-accounting identities. Their
numerical debt ratios are not asserted to be the inherited frontier debt
ratios or an all-scale bound.

The mathematical review has reached its coherent audit milestone. No shared
repository file was edited, no source theorem was promoted, and no commit was
made by this audit task. Publication and any required repository checkpoint
remain with the parent task. Open obligations remain universal source arrival,
a supplied all-source Hilbert carrier where needed, uniform original-source
cost, and source-specific coherent debt bounds.
