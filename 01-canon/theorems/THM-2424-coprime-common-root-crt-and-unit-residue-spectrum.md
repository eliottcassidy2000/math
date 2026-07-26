---
id: THM-2424
title: "Coprime common-root CRT realization and unit-residue spectrum"
status: >
  RESERVED / UNPROVED PROVISIONAL PROOF CANDIDATE UNDER AUDIT. For
  coprime fibre orders p,q, two profiles over the same parent have a
  canonical physical pq-root realization: evaluate them at the q- and
  p-dilates of one common inverse root. Its normalized pq-root transform
  factors exactly into the two profile transforms, and the polyphase
  identity converts every nonzero root character into global Fourier
  energy in the matching residue class modulo pq. Applied to the
  THM-2396/2397 owner-resolved clean-root cell, this gives one
  nonnegative rational two-dilation packet whose global spectrum meets
  every residue class coprime to 91, with explicit uniform and total
  energy floors; any fixed positive rational terminal word may be
  retained at all sufficiently large clocks. THM-2426 now proves that
  the primitive final-lane hypotheses supporting this specialization
  are inconsistent, so it is a counterfactual reusable carrier pending
  a transplant to a live residual. This is a derived physical packet,
  not the canonical fully masked owner endpoint, a relation current, a
  row exclusion, or a proof of LRC(14).
source: codex-2026-07-26-common-root-crt-spectrum
depends_on:
  - THM-2392-clean-toothpick-or-bounded-cross-ancestor-cage
  - THM-2396-common-core-forty-nine-orbit-word-incompatibility
  - THM-2397-clean-root-same-parent-charged-role-partition
related:
  - THM-2269-marked-expiration-root-spectrum-and-branch-state-no-go
  - THM-2305-canonical-blocker-word-handoff-hypergraph
  - THM-2380-cross-word-charged-target-correlation-and-pair-twist-gate
  - THM-2408-endpoint-prony-resultant-clock-separation-and-shared-node-boundary
  - THM-2418-alternating-base-thirteen-septimal-carry-matrix-and-rank-one-boundary
  - THM-2419-valuation-normalized-homogenization-of-affine-sideband-shells
  - THM-2426-compositional-thirteen-root-final-septimal-lane-exclusion
  - THM-2427-guard-top-thirteen-root-capacity-and-residual-types
script: 04-computation/lrc14_common_root_crt_spectrum_thm2424.py
output: 05-knowledge/results/lrc14_common_root_crt_spectrum_thm2424.out
script_sha256: a5d4aa11d0f4c8d12f500d78e0b23603ff4adc8867fcd1ffeed09215d83e3145
output_sha256: 67533d15f80cd512941856bbd5f6e8d4ca3c48bf213570ca140c09e9510a39e9
hash_basis: working-tree bytes (LF)
---

# THM-2424 -- transverse coprime fibres meet on one physical root cover

**RESERVED / UNPROVED PROVISIONAL PROOF CANDIDATE UNDER AUDIT.**

THM-2392 and THM-2397 retain two exact words over one clean parent:

```text
the unique excess on its seven siblings,

the exclusive q_* word on its thirteen predecessors.          (1)
```

They were correctly kept distinct: multiplying their separately
integrated coefficients would not create a physical current. There is,
however, a canonical operation which had not been used. A common
`91`-st inverse root simultaneously projects to one seven-sibling and
one thirteen-predecessor. Chinese remaindering then turns the pointwise
tensor, rather than an external product of integrals, into a genuine
circle function.

The resulting proof chain is

```text
same parent + coprime root orders
  -> one common physical root cover
  -> exact tensor factorization of every root character
  -> exact global Fourier residue-class energy
  -> every 91-unit residue class survives.                       (2)
```

The last arrow is spectral, not yet a canonical endpoint landing.

## 1. The coprime common-root construction

Let `p,q>=2` be coprime and put

```text
N=pq.
```

Let `Y` be a measurable subset of the circle, and let `U,V` be bounded
measurable circle functions. For a parent `y`, define their inverse-root
profiles

```text
U_y(s)=U((y+s)/p),                 s in Z/pZ,

V_y(r)=V((y+r)/q),                 r in Z/qZ.          (3)
```

On the common root coordinate `x`, put

```text
F(x)=1_Y(Nx) U(qx) V(px).                              (4)
```

For `t in Z/NZ`, the `t`-th inverse root of `y` is

```text
x_t=(y+t)/N.
```

Since

```text
q x_t=(y+(t mod p))/p                   mod 1,

p x_t=(y+(t mod q))/q                   mod 1,          (5)
```

one has the exact physical profile

```text
F(x_t)
 =1_Y(y) U_y(t mod p)V_y(t mod q).                     (6)
```

Thus the full abstract tensor of the two root profiles is realized on
one actual `N`-root fibre. Coprimality is precisely what makes every
pair of residues occur once.

## 2. Exact CRT factorization

Use normalized root transforms

```text
Uhat_y(ell)
 =(1/p) sum_s U_y(s) exp(-2 pi i ell s/p),

Vhat_y(k)
 =(1/q) sum_r V_y(r) exp(-2 pi i k r/q).              (7)
```

Let

```text
alpha=q^(-1) mod p,                 beta=p^(-1) mod q. (8)
```

For `m in Z/NZ`, define the normalized common-root transform

```text
C_m(y)
 =(1/N) sum_(t mod N)
     F((y+t)/N) exp(-2 pi i m t/N).                   (9)
```

The CRT idempotents are

```text
e_p=q alpha,                  e_q=p beta,
```

and every `t` is uniquely

```text
t=e_p s+e_q r                         mod N.           (10)
```

Substitution into (9) gives the exact factorization

```text
C_m(y)
 =1_Y(y)
   Uhat_y(m alpha mod p)
   Vhat_y(m beta mod q).                                (11)
```

No independence or averaging assumption is used. The two profiles may
be arbitrarily correlated through `y`; the product is taken before the
parent is integrated.

## 3. Polyphase converts root colour into physical frequency

For `F in L^2(T)`, use

```text
Fhat(n)=integral_T F(x)exp(-2 pi i n x)dx.             (12)
```

Expanding (9) in the Fourier basis and summing the roots gives

```text
C_m(y)
 =exp(2 pi i m y/N)
   sum_(h in Z) Fhat(m+Nh)exp(2 pi i h y).             (13)
```

Parseval on the parent circle therefore gives

```text
integral_T |C_m(y)|^2dy
 =sum_(n congruent m mod N)|Fhat(n)|^2.                (14)
```

Combining (11) and (14),

```text
sum_(n congruent m mod N)|Fhat(n)|^2
 =integral_Y
    |Uhat_y(m alpha)|^2 |Vhat_y(m beta)|^2dy.          (15)
```

This is the exact bridge. If the product on the right is nonzero on a
positive-measure parent set, some genuine physical frequency in the
class `m mod N` survives. If `m` is a unit modulo `N`, every such lift
is coprime to `N`.

There is also a bounded lift for rational step packets in every
nonzero residue class `m mod N` carrying positive energy in (15).
Suppose `F` is a finite step function and, in the progression `m+Nn`,
its oriented jumps group into `L_m` nonzero endpoint nodes. Integration
by parts gives a nonzero exponential sequence with at most `L_m`
nodes. A Vandermonde sequence cannot vanish on `L_m` consecutive
indices. Choosing the centered block therefore gives some surviving
lift with

```text
0<|n_physical|
 <=N ceil(L_m/2)-1.                                    (16)
```

The bound is finite and exact, but it supplies a physical Fourier
frequency, not a relation-current address.

## 4. The LRC clean-root packet

**POST-THM-2426 SCOPE.** The construction in Sections 4--6 is a valid
conditional consequence of the THM-2392/2396/2397 primitive final-lane
hypotheses. THM-2426 now proves that complete lane empty. The formulas
remain a reusable common-root carrier and hostile control, but do not
describe a live LRC continuation unless a separate transplant theorem
reconstructs the same parent cell in one of THM-2427's guard-top shapes
or on the `c_3<=M` side.

Use THM-2392's notation on the final LRC septimal lane. Let

```text
L_low
 =1_(E_H)
  +sum_(q_i!=q_*)1_(D_(q_i))
  +1_(D_(c_1))+1_(D_(c_2)).                           (17)
```

THM-2392 Section 3a and THM-2396 supply an owner-resolved parent cell
`Y` of mass

```text
rho=mu(Y)
 >=(1/26754)/338
 =1/9042852.                                          (18)
```

On that cell there is one fixed `d in F_7` such that

```text
J_y(s)
 :=L_low((y+s)/7)-1
 =1_({d})(s).                                         (19)
```

Let the literal exclusive top-label support be

```text
A_*(x)
 =1_(D_(q_*))(x)
   product_(i in {H} union {q_j:q_j!=q_*})
      (1-C_i(x)),                                     (20)
```

where `C_H=1_(E_H)` and the other `C_i` are the actual ordinary-role
danger indicators. THM-2397 equation (12a) proves that this is exactly
the physical single-factor deletion support. On `Y`, its thirteen-root
profile `A_y(r)=A_*((y+r)/13)` is fixed and is either

```text
singleton:  {r_0},

adjacent:   {r_0,r_0+eta},             eta!=0 mod 13. (21)
```

Apply (4) with

```text
p=7,              q=13,              N=91,

U=L_low-1,        V=A_*.
```

The resulting circle function is

```text
F_*(x)
 =1_Y(91x)[L_low(13x)-1]A_*(7x).                     (22)
```

Although `L_low-1` need not be nonnegative globally, the parent
selector in (22) restricts it to the `0/1` word (19). Hence `F_*` is a
nonnegative rational Boolean step packet, up to the usual null endpoint
set.

## 5. Every `91`-unit residue class survives

For `(p,q)=(7,13)`,

```text
alpha=13^(-1)=6 mod 7,

beta=7^(-1)=2 mod 13.                                  (23)
```

If `gcd(m,91)=1`, then

```text
ell=6m mod 7                  and
k=2m mod 13
```

are both nonzero. On `Y`,

```text
Jhat_y(ell)=(1/7)zeta_7^(-ell d).                     (24)
```

For the two possible top words,

```text
singleton:
  |Ahat_y(k)|^2=1/169;

adjacent:
  |Ahat_y(k)|^2
   =[2+2cos(2 pi k eta/13)]/169.                      (25)
```

Equations (15), (24), and (25) give, separately for every one of the
`phi(91)=72` unit residue classes,

```text
sum_(n congruent m mod 91)|Fhat_*(n)|^2
 =
  rho/8281,                                 singleton;

  rho[2+2cos(2 pi k eta/13)]/8281,          adjacent. (26)
```

Oddness of thirteen makes the adjacent quantity nonzero. More
quantitatively,

```text
min_(u!=0)|1+zeta_13^u|
 =2sin(pi/26)>2/13.                                   (27)
```

Therefore, in either status and every unit class,

```text
sum_(n congruent m mod 91)|Fhat_*(n)|^2
 >4rho/(169*8281)
 >=1/3163842975657.                                   (28)
```

Summing all unit classes gives the exact alternatives

```text
sum_(gcd(n,91)=1)|Fhat_*(n)|^2
 =
  72rho/8281,                              singleton;

 132rho/8281,                              adjacent.  (29)
```

In particular,

```text
sum_(gcd(n,91)=1)|Fhat_*(n)|^2
 >=6/6240321451.                                      (30)
```

Thus every formal packet satisfying the inherited final-lane hypotheses
would carry one owner-resolved, top-labelled, nonnegative rational
physical packet with a nonzero Fourier lift in every residue class
coprime to `91`. Post-THM-2426, this is a counterfactual carrier rather
than a live scalar row.

## 6. A delayed terminal word can be retained

Let `Q` be any fixed positive rational Boolean circle word. At a clock
`R=13^k`, its value on the thirteen-root fibre is the parent-constant
word

```text
Q(R(y+r)/13)=Q(13^(k-1)y).                            (31)
```

Define

```text
Y_R
 =Y intersection {y:Q(13^(k-1)y)=1},

rho_R=mu(Y_R).                                        (32)
```

THM-2397's two-BV estimate gives

```text
|rho_R-rho mu(Q)|
 <=Var(1_Y)Var(Q)/(12*13^(k-1)).                      (33)
```

Consequently `rho_R>0` for every sufficiently large `k`, explicitly
whenever

```text
13^(k-1)
 >Var(1_Y)Var(Q)/(12rho mu(Q)).                       (34)
```

Replacing `Y` by `Y_R` in (22) inserts the actual delayed word and
leaves every formula (26)--(29) valid with `rho_R` in place of `rho`.
Thus a fixed positive terminal word cannot erase the all-`91`-unit
common-root spectrum at large clocks.

This is the positive source--terminal correlation which the flat
histogram hostile in THM-2418 deliberately lacks: the cell `Y` fixes
the septimal excess address `d` before the root transform.

## 7. Sharp boundaries

### 7a. Coprimality is necessary

If `p=q=2`, the two residue coordinates of one common root must agree.
Take profiles

```text
U_y=(1,0),                   V_y=(0,1).
```

Their abstract tensor has a nonzero `(0,1)` cell, but no common root has
those two residues, so (4) is identically zero. The full product
factorization (11) is therefore a genuinely coprime phenomenon.

### 7b. Root energy is not a chosen endpoint phase

Equation (14) forces energy in a residue progression, not one
preselected Fourier lift. Endpoint-Prony gives (16), but the sharp
two-node hostiles of THM-2408/2416 rule out a uniform one-mode amplitude
from jump count alone.

### 7c. The packet is physical but derived

The construction (22) is stronger than the external product of
separately integrated currents in THM-2397: it is one literal circle
function, and its Fourier coefficients in (26) are ordinary physical
coefficients. It nevertheless evaluates different factors at `7x`,
`13x`, and `91x`. The parent selector remembers the owner, status,
root translate, and septimal excess address.

Canon has not identified this two-dilation packet with:

- the lawfully target-co-shifted THM-2365 current;
- the canonical fully masked THM-2305 owner endpoint;
- a same-marked-triangle Abel coefficient;
- an exact relation-neutral address in the physical lattice; or
- the same-shell residue-zero reference required by THM-2419.

A physical unit frequency is not automatically a relation current.
Moreover, THM-2426 empties the primitive final lane on which the
THM-2396/2397 specialization was built. The live transplant question is
whether one of THM-2427's seven guard-top regimes, or the `c_3<=M`
noncirculant-graft side, supplies an analogous fixed coprime parent
selector.
No row is excluded, the ledger remains `165`, and LRC(14) remains open.

## 8. Exact companion

The dependency-free companion:

- exhausts the CRT character map for all coprime `2<=p,q<=13`;
- checks the common-root factorization on thousands of integer profile
  pairs;
- verifies root orthogonality by exact integer reduction modulo the
  relevant cyclotomic polynomial and then checks the finite polyphase
  identity on rational Fourier packets;
- exhausts all `7*13*2*72=13,104` owner-address/status/unit-character
  cases in the LRC specialization by exact group-ring and cyclotomic
  reductions;
- checks the singleton/adjacent class energies, total unit energies,
  universal rational floors, and the strict chord bound;
- verifies centered endpoint-Prony controls; and
- exhibits the noncoprime `p=q=2` hostile.

Run

```bash
python3 04-computation/lrc14_common_root_crt_spectrum_thm2424.py
python3 -O 04-computation/lrc14_common_root_crt_spectrum_thm2424.py
```

Both transcripts must byte-match, after LF normalization,

```text
05-knowledge/results/lrc14_common_root_crt_spectrum_thm2424.out
```

Every truth-bearing check raises explicitly, so optimized mode executes
the same audit.
