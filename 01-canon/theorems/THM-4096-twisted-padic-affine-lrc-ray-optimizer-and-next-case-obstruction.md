---
id: THM-4096
title: "Twisted p-adic affine margin optimizer and next-ray obstruction"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY VERIFIED-EXACT. On the canonical Cover14 ray
  V_m={1,...,12,182m}, the normalized THM-2057 witness margin has an exact
  rational affine expansion against the correctly typed twisted values
  L_p(-1,omega_p^2)=(p-1)/12. The first point m=1 is the unique positive-side
  point and has a complete nonnegative optimizer; every m>=2 lies strictly
  below zeta(-1) and admits no such nonnegative expansion. Two exact Dedekind
  observers bracket the carrier for m>=2. This is a rational affine shadow,
  not a common p-adic-valued sum, an f_14 identification, or LRC(14).
source: codex-padic-zeta-tournament-20260825
depends_on:
  - THM-2057-scaled-zeta-core-one-tail-closure
  - THM-3706-lrc-dedekind-checksum-unbounded-bockstein-depth
related:
  - HYP-3779
  - MISTAKE-497
  - MISTAKE-503
  - MISTAKE-504
script: 04-computation/lrc_twisted_padic_affine_margin_optimizer_thm4096.py
output: 05-knowledge/results/lrc_twisted_padic_affine_margin_optimizer_thm4096.out
independent_audit_script: 04-computation/lrc_twisted_padic_affine_margin_optimizer_thm4096_independent_audit.py
independent_audit_output: 05-knowledge/results/lrc_twisted_padic_affine_margin_optimizer_thm4096_independent_audit.out
script_sha256: 277a9dca809654f0af575ca1812a6983e6163a29dea79717ce60601014afffbd
output_sha256: 41d9d95a654b13c04e623e15e1ea278d69cdb5ae2252ead67b13a59073157318
independent_audit_script_sha256: 392ad6f67d3be5edd4c43d9366c775cd2ac9a2111c2da494706f11455e8dddef
independent_audit_output_sha256: da5cb1670357369dbbdb6a366448e634fcbce09452c76eb70484d6b9bbe884bf
hash_basis: raw LF bytes
---

# THM-4096 -- the twisted p-adic affine optimizer stops after the first ray point

**PROVED + VERIFIED-EXACT + INDEPENDENTLY VERIFIED-EXACT.** This theorem
repairs the character error in HYP-3779, extracts the exact optimization
problem that survives it, and then tests the next point rather than
extrapolating from the exceptional first one.

## 1. The Cover14 witness and its normalized carrier

For an integer `m>=1`, put

```text
V_m={1,2,...,12,182m},       D_m=182m+1.                    (1)
```

THM-2057 proves that `V_m` is `1/14`-lonely at

```text
t_m=14m/D_m.                                                (2)
```

Measure the phase past the endpoint `1/14`, and normalize it on the
Bernoulli scale:

```text
Delta_m=t_m-1/14=(14m-1)/(14D_m),
C_m=-(14^2/12)Delta_m=-7(14m-1)/(6D_m),
z_infinity=zeta(-1)=-1/12.                                 (3)
```

Then

```text
boxed: C_m-z_infinity=(15-14m)/(12D_m).                     (4)
```

Consequently `m=1` is the unique positive integer on the positive side of
`z_infinity`:

```text
C_1=-91/1098=z_infinity+1/2196,
C_m<z_infinity                    for every m>=2.           (5)
```

This sign change is the first hostile control. The special identity at
`m=1` does not continue down the ray.

## 2. Correctly typed p-adic vertices

For an odd prime `p`, let `omega_p` be the Teichmuller character. The
Kubota--Leopoldt interpolation formula at `n=2`, on the character branch
`chi=omega_p^2`, gives

```text
Z_p^(2):=L_p(-1,omega_p^2)
        =(1-p)zeta(-1)
        =(p-1)/12
        =z_infinity+p/12.                                  (6)
```

Thus `Z_7^(2)=1/2`. This is the value that HYP-3779 incorrectly called the
trivial-character value `zeta_7(-1)`. Equation `(6)` follows directly from
equation (1.1) of [Luochen Zhao, *Sum Expressions for Kubota--Leopoldt
p-adic L-functions*](https://arxiv.org/abs/2201.08870): with
`chi=omega_p^2`, its classical character `chi omega_p^(-2)` is trivial.

For the level-14 two-prime optimization below, also define the rational
Euler-factor vertex

```text
Z_2^(2):=(2-1)/12=1/12.                                    (7)
```

No 2-adic character convention is needed for the theorem: `(7)` is used
only as a rational vertex. More generally, all vertices in `(6)--(7)` are
rational numbers. Combining vertices attached to different primes below is
therefore an optimization in their common rational affine line. It is not an
addition operation between distinct fields `Q_p`.

The standard trivial-character branch is genuinely different at `p=7`.
For every odd prime `p>=5` and every integer `j>=1`, set

```text
k_j=2+(p-3)p^j.                                             (8)
```

Then `k_j -> 2` p-adically, `p-1` divides `k_j`, and `k_j=2 mod p`.
Interpolation on the trivial branch and von Staudt--Clausen give

```text
zeta_p(1-k_j)=(1-p^(k_j-1))(-B_(k_j)/k_j),
v_p(zeta_p(1-k_j))=-1,
p zeta_p(1-k_j)=1/2 mod p.                                 (9)
```

Continuity at `-1` preserves the last two statements, so

```text
boxed:
v_p(zeta_p(-1))=-1,
p zeta_p(-1)=1/2 mod p.                                   (10)
```

In particular `zeta_7(-1)` lies in `4/7+Z_7` and cannot equal the unit
`1/2`. The primary referee checks six finite interpolation controls; the
all-`j` conclusion is the exact argument `(8)--(10)`, not a finite
extrapolation.

## 3. The complete nonnegative affine optimizer

Let `P` be a finite nonempty set of primes. For coefficients
`alpha_p>=0`, put

```text
alpha_infinity=1-sum_(p in P) alpha_p,
H_P(alpha)=alpha_infinity z_infinity
           +sum_(p in P) alpha_p Z_p^(2).                  (11)
```

Here `(7)` supplies the rational vertex when `2 in P`. By `(4)` and `(6)`,

```text
boxed:
C_m=H_P(alpha)
iff
sum_(p in P) p alpha_p=(15-14m)/D_m.                       (12)
```

The right side is positive only at `m=1`. Therefore:

```text
* m=1 admits nonnegative affine representations;
* every m>=2 admits none.                                  (13)
```

At `m=1`, any solution automatically has
`sum alpha_p<=1/(2*183)<1`, so `alpha_infinity>0`; the
representations are honest convex combinations. For the level primes
`P={2,7}`, write `a=alpha_2` and `b=alpha_7`. The entire feasible segment is

```text
2a+7b=1/183,
0<=b<=1/1281,
a=1/366-(7/2)b.                                            (14)
```

For the equal-cost objective `a+b`, the value decreases strictly with `b`.
Its unique minimizer is

```text
boxed: (a,b)=(0,1/1281),       min(a+b)=1/1281.             (15)
```

If both level-prime coefficients are required to be positive, `1/1281` is
only an infimum and is not attained. More generally, for any finite `P`, the
minimum of `sum alpha_p` is `1/(183 max(P))`, uniquely attained by putting
all mass at `max(P)`. Over arbitrary finite prime support the infimum is zero
and is not attained. Hence the optimizer itself does not select the level
prime `7` unless the allowed local lanes are fixed in advance.

## 4. Two Dedekind observers expose the next-case geometry

Let `s(h,k)` be the classical sawtooth Dedekind sum. THM-3706 proves

```text
F_m=s(14,D_m)=91m(13m-14)/(6D_m),
A_m=s(14m,D_m)=91m(13-14m)/(6D_m).                          (16)
```

Direct subtraction from `(3)` gives the factorizations

```text
F_m-C_m=(m-1)(1183m+7)/(6D_m),
C_m-A_m=(m-1)(1274m-7)/(6D_m).                             (17)
```

The remaining comparison in the full order has the independent positive
factorization

```text
F_m-1/2=((m-2)(1183m+546)+1089)/(6D_m)>0       (m>=2).     (17a)
```

At `m=1` there is a triple collision

```text
F_1=C_1=A_1=-91/1098.                                      (18)
```

For every `m>=2`, the factors in `(4)`, `(17)`, and `(17a)` are strict,
producing the complete order

```text
boxed: A_m<C_m<-1/12<1/2<F_m.                              (19)
```

The first new case is already separated by visible rational gaps:

```text
m=2:  -91/73 < -63/730 < -1/12 < 1/2 < 364/365,
C_2-z_infinity=-13/4380.                                  (20)
```

Thus neither natural Dedekind observer equals the normalized carrier after
the first point. The collision at `m=1` is a boundary event, not evidence of
a ray-wide identity.

Finally, in the orientation used by Atiyah's formula (4.25), the signature
eta invariant is

```text
eta_sig(L(k,h))=-4s(h,k).                                  (21)
```

Hence `eta_sig(L(183,14))=182/549`, not `-91/1098`; see
[Michael Atiyah, *The Logarithm of the Dedekind
Eta-Function*](https://webhomes.maths.ed.ac.uk/~v1ranick/papers/atiyahlg.pdf).

## 5. Scope and reproduction

The exact connection ledger is:

```text
source:       THM-2057 witness phase t_m on the 182m Cover14 ray
map:          Delta_m |-> C_m=-(14^2/12)Delta_m
target:       rational affine line through zeta(-1) and Z_p^(2)
preserved:    exact rational value and side of zeta(-1)
destroyed:    p-adic topology, character field, and local-place identity
sidecar:      the character omega_p^2 and the signed moment sum p alpha_p
hostile:      m=2, where the required moment is negative
consequence:  a complete first-point optimizer and a ray-wide no-go.       (22)
```

This theorem does not identify a common adelic invariant, the `f_14` cusp
residual, or an LRC owner/intertwiner. It does not strengthen the already
proved loneliness of `(1)`. It converts a false analogy into one exact
positive boundary case and an infinite obstruction beyond it.

Reproduce the exact arithmetic audit with

```bash
python3 -B 04-computation/lrc_twisted_padic_affine_margin_optimizer_thm4096.py
python3 -B -O 04-computation/lrc_twisted_padic_affine_margin_optimizer_thm4096.py
```

Independent reproduction:

```bash
python3 -B 04-computation/lrc_twisted_padic_affine_margin_optimizer_thm4096_independent_audit.py
python3 -B -O 04-computation/lrc_twisted_padic_affine_margin_optimizer_thm4096_independent_audit.py
```

The referee checks `2,000` ray identities, all `1,999` strict interval
certificates, twenty direct Dedekind sums, the complete `{2,7}` optimizer,
seven finite-support optimizer gates, six rational Euler-factor vertices,
six Bernoulli interpolation controls through index `198`, the eta
normalization, and the exact `m=2` hostile. It uses only exact rational
arithmetic and explicit `require` gates, so optimized execution retains all
checks.

The independent referee reconstructs the phase from folded residues rather
than importing `(2)`, evaluates `1,024` Dedekind sums by the Euclidean
reciprocity algorithm rather than sawtooth enumeration, computes Bernoulli
numbers by the Akiyama--Tanigawa transform, and proves optimizer uniqueness
with exchange identities. It checks `4,096` ray points, 101 points of the
complete `{2,7}` segment, twenty-five nested prime banks through `97`, and
three hostile versions of the retracted claims. Both transcripts match their
frozen outputs under normal and optimized execution.

**QED.**
