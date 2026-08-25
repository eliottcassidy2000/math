---
id: THM-3481
title: "Rule 30 cyclic arc-norm rank and marked innovation spectrum"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  On every power-of-two
  phase cycle, integration over an arc of length m has an exact rank and
  Hasse-moment loss determined only by 2^nu_2(m)-1.  Odd arc integration is
  invertible.  Applied to the Rule 30 transition current, every odd innovative
  terminal profile retains maximal ANF degree; full Walsh support holds once
  its innovation cube has dimension at least two, while
  its physically marked value still requires the light-cone basepoint.  The
  edge-period lower bound gives an exact at-most-two-wrap trichotomy.  No Rule
  30 prize is claimed.
source: root-rule30-next-targets-20260815
audit: >
  An independent hostile audit rederived the cyclic local-ring
  factorization, ranks, image ideals, complete Hasse constraints,
  constructive odd inverse, Rule 30 marked coupling, and wrap endpoints.  It
  found and repaired the one-variable Walsh boundary before acceptance, then
  extended the operator and Boolean-table controls beyond the companion
  universe.  Ordinary and optimized replays match stored output
  byte-for-byte: ACCEPT.
depends_on:
  - THM-3458-rule30-right-edge-2-adic-odometer-and-moving-observer-boundary
  - THM-3463-rule30-mealy-section-suffix-parity-current-and-complexity-boundary
  - THM-3471-rule30-motzkin-strip-circuit-and-innovation-carry-spectrum
related:
  - THM-3468-rule30-radial-green-fold-innovation-discrepancy-and-fixed-seed-carrier-boundaries
  - THM-4050-rule30-half-arc-marked-cylinder-and-radius-nine-hostile
script: 04-computation/rule30_cyclic_arc_norm_rank_thm3481.py
output: 05-knowledge/results/rule30_cyclic_arc_norm_rank_thm3481.out
script_sha256: ccc2273232736380c0830528294517ff585cf3d12a8857b270881271b11dad95
output_sha256: e50eae451ccb594d08ce791c5fb9c92c8b4764f9f4e247e924de921d9e1ec1aa
hash_basis: raw bytes
---

# THM-3481 -- Rule 30 cyclic arc-norm rank and marked innovation spectrum

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

THM-3471 identifies the distinguished center as a moving terminal integral of
the edge transition current.  The present theorem computes exactly what that
integration operator preserves and destroys.  The answer is controlled not by
the arc length itself, but by its two-adic valuation.

## 1. Inheritance, live objects, and conventions

The closest proved mechanism is THM-3471's current identity.  For depth
`k>=2`, let

```text
T_k(h)=b_k(h+k),
T_k(h+1)+T_k(h)=Q_k(h),
xor_(h mod P_k) Q_k(h)=epsilon_k,                    (1)
```

where `P_k` is a power of two and

```text
T_k(-k)=0,       T_k(0)=c_k,
c_k=xor_(h=-k)^(-1)Q_k(h).                           (2)
```

The canonical hostile is that a phase average does not determine the marked
value at `h=-k`.  The corrected near miss is to retain both the cyclic norm
operator and the physical basepoint.  The least-used sidecar is the
`(S+1)`-adic, equivalently Hasse-moment, filtration of the phase cycle.

The concept board is:

1. a power-of-two cyclic phase module;
2. the transition current `Q_k`;
3. the arc-norm operator;
4. its `(S+1)`-adic rank loss;
5. innovation-cube ANF and Walsh spectra; and
6. the marked light-cone basepoint.

All sums in the proof are over `F_2`.  The abstract result below applies to
every function on a power-of-two cycle, independently of Rule 30.

## 2. Exact cyclic arc-norm rank

Let `p=2^d`, let

```text
V_p={q:Z/pZ -> F_2},
A_m q(h)=xor_(0<=j<m) q(h+j)                 (m>=1), (3)
```

and put `m=2^a u` with `u` odd.  Then

```text
r(m)=2^a-1.                                           (4)
```

### Theorem 2.1 (rank, kernel, and image)

For every `p=2^d` and every `m>=1`,

```text
ker_dimension(A_m)=min(r(m),p),
rank(A_m)=max(p-r(m),0).                              (5)
```

If `r(m)<p`, the image is exactly the codimension-`r(m)` cyclic ideal

```text
im(A_m)=(S+1)^r(m) V_p.                               (6)
```

If `r(m)>=p`, then `A_m=0`.

### Proof

Use the backward cyclic shift

```text
(S q)(h)=q(h-1).                                      (7)
```

Represent `q` by

```text
F_q(X)=sum_(h=0)^(p-1)q(h)X^h
```

in

```text
R_p=F_2[X]/(X^p-1)=F_2[X]/((X+1)^p).                 (8)
```

Then `S` acts as multiplication by `X`.  With

```text
N_m(X)=1+X+...+X^(m-1),                              (9)
```

the forward arc in (3) is

```text
A_m=S^(-(m-1)) N_m(S).                               (10)
```

The leading phase rotation is a unit and cannot change a kernel, rank, or
`(X+1)`-adic image ideal.

Since `u` is odd,

```text
X^u+1=(X+1)U(X),       U(1)=1.                       (11)
```

Frobenius and `(X+1)N_m=X^m+1` give the polynomial identity

```text
N_m(X)
 =(X+1)^(2^a-1) U(X)^(2^a).                          (12)
```

The last factor is a unit in the local ring (8).  Writing `Y=X+1`,
multiplication by `N_m` therefore has the same kernel and image as
multiplication by `Y^r` in `F_2[Y]/(Y^p)`.  Its image has basis
`Y^r,...,Y^(p-1)` when `r<p`, while its kernel has basis
`Y^(p-r),...,Y^(p-1)`.  When `r>=p`, it is zero.  This proves (5)--(6).

The proof includes all arcs longer than one cycle.  In particular, for any
odd integer `v`,

```text
rank(A_(vp))=1,
rank(A_(2vp))=0.                                      (13)
```

For `p<m<2p`, formula (5) remains unchanged; there is no boundary exception
at `m=p`.

## 3. Complete Hasse constraints and a constructive odd inverse

For `f in V_p`, define its cyclic Hasse moments using the standard
representatives `0,...,p-1`:

```text
M_j(f)=xor_(h=0)^(p-1) binom(h,j)f(h).                (14)
```

In the Taylor expansion at `X=1`,

```text
F_f(1+Y)=sum_(j=0)^(p-1) M_j(f)Y^j.                  (15)
```

Consequently, when `r=r(m)<p`,

```text
A in im(A_m)
 iff M_0(A)=...=M_(r-1)(A)=0.                        (16)
```

These are the complete universal linear constraints: there are exactly `r`
independent conditions, and (5) leaves no further one.  When `r>=p`, all `p`
Taylor coefficients vanish and `A=0`.  Thus:

```text
m odd:          no Hasse loss;
m=2 times odd: M_0=0;
m=4 times odd: M_0=M_1=M_2=0;
m=8 times odd: M_0,...,M_6=0;                        (17)
```

and so on.

Odd arc integration is not merely full rank; it has an explicit inverse.
If `m` is odd, choose `n` with

```text
mn=1 mod 2p.                                         (18)
```

The norm-composition identity gives

```text
N_m(S)N_n(S^m)=N_(mn)(S)=1.                          (19)
```

Indeed, `mn=1+2tp`; in the last norm, every residue class occurs `2t`
times except residue zero, which occurs once more.  Restoring the harmless
rotation in (10) reconstructs `q` from `A_mq`.  Therefore:

```text
boxed: every odd-length cyclic arc profile is a lossless recoding
       of the complete phase current.                                (20)
```

This is an operator theorem, not a claim that the recoding is sparse in the
innovation coordinates or cheap in a particular machine model.

## 4. Rule 30 terminal profiles

Return to (1), put `p=P_k`, and define for any `m>=1`

```text
A_(k,m)(h)=xor_(0<=j<m) Q_k(h+j).                    (21)
```

THM-3471 gives the physical marked evaluation

```text
boxed: c_k=A_(k,k)(-k).                               (22)
```

Theorem 2.1 immediately yields, with `a=nu_2(m)`,

```text
rank(q |-> A_(k,m))=max(P_k-2^a+1,0),                (23)
```

and the first `min(2^a-1,P_k)` Hasse conditions in (16).  In particular,
for odd `k` the actual length-`k` terminal integration loses no phase-current
information at all.

Summing (21) over every phase counts each current bit `m` times, so

```text
xor_(h mod p) A_(k,m)(h)
 =(m mod 2) epsilon_k.                                (24)
```

Equation (24) is the bridge from cyclic rank to innovation spectrum.

### 4.1 Odd innovative terminal integrals remain spectrally full

Suppose `epsilon_k=1` and `m` is odd.  Let

```text
Gamma_k:Z/P_k Z -> F_2^d,
d=log_2 P_k,                                          (25)
```

be THM-3463/THM-3471's innovation-coordinate bijection below depth `k`, and
pull `A_(k,m)` through `Gamma_k`.  By (24), its truth table has odd support.
The coefficient of the full monomial in the Boolean ANF is the parity of the
truth table.  Hence

```text
deg_ANF(A_(k,m) o Gamma_k^(-1))=d,                   (26)
```

and every innovation variable is essential.

For `d>=2`, every unnormalized real Walsh coefficient is congruent to two
modulo four.  At the zero character it is `2^d-2w` with `w` odd; at every
nonzero character the pure character sum vanishes and the result is minus
twice a sum of an odd number of signs.  Therefore

```text
boxed: every Walsh coefficient is 2 mod 4 and nonzero. (27)
```

The one-variable cube is the sharp boundary: at `k=3`, one Walsh coefficient
vanishes.  Taking `m=k` proves (26) for every **odd innovative physical
terminal profile**, and proves (27) when `d>=2`.  The depth `k=3` profile on
the one-variable cube is the explicit exception.  This is strictly about the
whole spatial phase profile; it does not determine the marked sample (22).

THM-4050 later replaces the physical length `k` by the light-cone length
`floor(k/2)+1`, applies the odd-length criterion there, and identifies the
nearest-zero adaptive endpoint with a marked spatial cylinder. It does not
turn this unpointed rank theorem into a temporal balance result.

## 5. The at-most-two-wrap physical trichotomy

THM-3458 proves

```text
ceil(k/2)<=P_k.                                       (28)
```

Thus the physical length-`k` arc never exceeds two current cycles.  Write
`p=P_k`.  Since one cycle has XOR `epsilon_k`, there are exactly four cases:

```text
k<p:          A_(k,k) is the raw length-k arc;
k=p:          A_(k,k) is the constant epsilon_k;
p<k<2p:       A_(k,k)=epsilon_k+A_(k,k-p);
k=2p:         A_(k,k)=0.                              (29)
```

In the third case, `k-p<k/2`, so the sole wrap removes more than half of the
arc after recording one holonomy bit.  At the two-cycle endpoint,

```text
c_k=0.                                                (30)
```

Moreover that endpoint is necessarily innovative.  If `k=2p`, then at time
`p=P_k` the packed row `R_p` has its unique extreme-left top bit at index
`2p=k`.  Hence

```text
epsilon_k=bit_k(R_(P_k))=1,
P_(k+1)=2p.                                           (31)
```

At the one-cycle endpoint, (22) and (29) give simply

```text
c_k=epsilon_k       when k=P_k.                       (32)
```

These are all-`k` consequences.  They do not assume that either endpoint case
occurs infinitely often.

## 6. Decisive hostile and prize scope

At an odd innovative depth with `P_k>1`, the terminal profile has odd support
inside an even-sized phase cycle.  It is therefore neither zero nor one
everywhere.  Rotating the mark realizes both possible sampled bits while
preserving the unpointed cyclic current, its orbit statistics, and its arc
rank.  Even maximal ANF degree, lossless odd arc integration, and (when
`log_2 P_k>=2`) full Walsh support therefore fail to choose `c_k`; the
boundary-calibrated phase `h=-k` is load-bearing.

Equivalently, in spin notation

```text
sigma_k(h)=(-1)^T_k(h),
(-1)^A_(k,m)(h)=sigma_k(h+m)sigma_k(h).               (33)
```

The center sign is one marked open-line product, not an additive current
average.  A small additive discrepancy for the current would not by itself
control its parity.  The live Prize-2 obligation is a **pointed diagonal
product correlation across depths**, not another Haar phase density.

The preservation ledger is:

| map | preserves | destroys | required sidecar |
|---|---|---|---|
| `q -> A_mq`, `m` odd | the full phase current, constructively | no linear phase information | marked evaluation remains separate |
| `q -> A_mq`, `nu_2(m)=a` | the quotient of rank `p-2^a+1` | exactly a `2^a-1`-dimensional repeated-root sector | Hasse moments / a preimage representative |
| terminal profile -> support parity | innovation holonomy for odd `m` | every ordered phase value | physical basepoint `-k` |
| phase cycle -> innovation cube | all phase values and Haar measure | the native cyclic law | nonlinear odometer action |
| unpointed profile -> marked center | orbit statistics | the selected coefficient | light-cone calibration |

This theorem proves no center nonperiodicity, no temporal density `1/2`, and
no computational lower bound.  It narrows the marked-origin target and gives
an exact obstruction to treating arc integration as a smoothing or
complexity-reducing operation at odd depths.

## 7. Exact verification

Reproduce the finite audit with

```bash
python3 04-computation/rule30_cyclic_arc_norm_rank_thm3481.py
python3 -O 04-computation/rule30_cyclic_arc_norm_rank_thm3481.py
```

The operator audit constructs every direct circulant arc matrix for

```text
p=2^d,       0<=d<=7,
1<=m<=2p,                                            (34)
```

and independently checks all `510` ranks, the `1,546` forced Hasse moments,
all `255` constructive odd inverses, the rank-one one-cycle boundary, and the
zero two-cycle boundary.

The physical audit computes the Rule 30 edge periods through depth `30`,
checks every phase modulo every `P_k`, and compares the cyclic profile with an
independent centered light-cone evolution.  It checks (22), (24), every
available wrap case, and the complete ANF/Walsh tables at all innovative
depths through `29`.  Universal statements are the proofs above; the finite
universes audit the implementation and supply explicitly labelled exact
controls only.
