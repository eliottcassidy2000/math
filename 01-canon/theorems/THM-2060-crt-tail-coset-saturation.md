---
id: THM-2060
title: "Sharp CRT tail-coset saturation and the dyadic two-tail seam"
status: >
  PROVED COROLLARY / SHARP REPACKAGING. A translated q-grid meets the open
  1/14-danger arc in at most ceil(q/7) points. Consequently every THM-2059
  tail-histogram bin has at least q-ceil(q/7) points, where
  q=h/gcd(N,h)=a/gcd(a,w) is independent of the clock. This gives an exact
  packet lower bound, a finite-clock exceptional-residue atlas, and reduces a
  primitive scaled two-tail LRC(14) row to the sole dyadic/odd seam. The
  qualitative one-tail sheet dodge and hereditary-primitivity consequence
  were already proved in THM-760/761/765; the new content is the sharp
  per-coset count, its CRT formulation, and the resulting atlas composition.
source: codex-2026-07-21-LRC-CRT-coset-saturation
script: 04-computation/lrc_crt_tail_coset_saturation_codex_20260721.py
result: 05-knowledge/results/lrc_crt_tail_coset_saturation_codex_20260721.out
script_sha256: 8e425115997dc1e0fadbd0ba56a880e0e27b384cfe4709e2dfe800b2a8ec7366
result_sha256: a022b5e7bb4b6b9528365836c5546c7a977e3b7b6d0e7dc6019614e6fcc8df58
hash_basis: normalized repository blobs (LF)
depends_on:
  - THM-2059
  - THM-761
related:
  - THM-760
  - THM-765
  - THM-769
  - THM-775
  - THM-2057
  - THM-2058
  - HYP-8846
---

# THM-2060 -- Sharp CRT tail-coset saturation

Put

```text
|x|_m=min(x mod m,m-(x mod m)).
```

The weak LRC(14) condition is `14|x|_m>=m`; its complement is the
**open** danger arc. Keeping that endpoint convention is essential below.

## A. Sharp translated-grid count

Let `0<lambda<1/2`. Any translate of the `q`-grid in `R/Z` has at most

```text
ceil(2 lambda q)                                           (1)
```

points in the open arc `{z:||z||<lambda}`. The bound is sharp over all
translates. Hence every translate has at least

```text
kappa_lambda(q)=q-ceil(2 lambda q)                         (2)
```

points in the closed safe complement. In particular, every translate is
guaranteed a safe point exactly when

```text
lambda <= (q-1)/(2q).                                     (3)
```

At `lambda=1/14`,

```text
kappa(q)=q-ceil(q/7),                                     (4)
```

which is positive exactly for `q>=2`.

### Proof

If `m` grid points lie in one lift of the open arc, the first and last are
separated by `(m-1)/q<2lambda`. Thus `m-1<2lambda q`, which gives (1),
including the integral endpoint because the arc is open. Conversely,
`ceil(2lambda q)` consecutive grid points have span strictly below
`2lambda` (with the same endpoint observation), so a translate puts them
inside the arc. This proves sharpness and (2). The inequality
`ceil(2lambda q)<=q-1` is equivalent to (3). Equation (4) follows. QED.

## B. Exact sheet multiplicity for one tail

Let `C` be a finite set of positive integers, `a,w>0`, and suppose `theta`
is `lambda`-safe for `C`. Put

```text
c=gcd(a,w),             q=a/c,
t_k=(theta+k)/a,        k in Z/aZ.                        (5)
```

Every one of the `a` lifts has exactly the same core clearance:

```text
||(au)t_k||=||u theta||             for u in C.            (6)
```

The tail phases form a translated `q`-grid, each point repeated `c` times.
Consequently at least

```text
c kappa_lambda(q)=c(q-ceil(2lambda q))                    (7)
```

of the lifts are simultaneously `lambda`-safe for `aC union {w}`.

At `lambda=1/14`, `a` not dividing `w` is equivalent to `q>=2`, so every
chosen weak-safe core phase has a weak-safe full-row lift. When `a|w`, (7)
only becomes zero: it does **not** say that the full row is unsafe.

### Proof

Equation (6) is immediate. As `k` runs modulo `a`, the residues `wk mod a`
form the subgroup `c Z/aZ`; its `q` values each occur `c` times. Translate
by `w theta/a` and apply Part A. QED.

This is the same sheet geometry used qualitatively in THM-760/761 and in the
deck obstruction THM-765. Part B records the sharp safe-lift multiplicity;
it is not a new hereditary-primitivity proof.

## C. THM-2059 histograms have full tail support

Use THM-2059's notation

```text
Q=Na,       g=gcd(w,Q),       h=Q/g,       u=w/g,
d=gcd(N,h),
A_N(C)={r mod N:14|cr|_N>=N for every c in C},
beta_j=#{s mod h:14|us|_h>=h and s=j mod d}.              (8)
```

If `c=gcd(a,w)` and `q=a/c`, then exactly

```text
h/d=q,       lcm(N,h)=Nq,       Q/lcm(N,h)=c.             (9)
```

Every tail bin obeys the sharp uniform lower bound

```text
beta_j >= q-ceil(q/7)=kappa(q)          for every j mod d. (10)
```

Thus, writing `P_N` for THM-2059's histogram dot product,

```text
P_N(C;a,w) >= kappa(q)|A_N(C)|,                           (11)
# safe k mod Na >= c kappa(q)|A_N(C)|.                    (12)
```

The tail histogram has full support on `Z/dZ` if and only if `q>=2`, or
equivalently `a` does not divide `w`. For `q=1`, the bin `j=0` consists only
of the tail-dangerous residue `s=0` and is empty.

### Proof

For every prime `p`, write `n=v_p(N)`, `A=v_p(a)`, and `W=v_p(w)`. Then

```text
v_p(h)=max(n+A-W,0),
v_p(d)=min(n,v_p(h)),
v_p(h/d)=max(A-W,0)=v_p(a/gcd(a,w)).                       (13)
```

This proves the first identity in (9); the other two follow from
`lcm(N,h)=Nh/d` and `Q=Na`.

Each congruence class modulo `d` contains `q=h/d` residues modulo `h`.
Multiplication by the unit `u` maps it bijectively to another such class.
After division by `h`, its image is a translated `q`-grid. Part A shows that
at most `ceil(q/7)` of those residues are tail-dangerous, proving (10).
Now THM-2059 gives

```text
P_N=sum_j alpha_j beta_j >= kappa(q)sum_j alpha_j,
```

which is (11), and its lift factor is `Q/lcm(N,h)=c`, proving (12). The
full-support statement follows from (4); its `q=1` converse is the displayed
zero-bin observation. QED.

## D. Several tails and the primitive two-tail seam

For `S=aC union {w_1,...,w_r}`, put

```text
c_i=gcd(a,w_i),       q_i=a/c_i.
```

At any `1/14`-safe core phase, tail `i` makes at most
`c_i ceil(q_i/7)` of the `a` lifts dangerous. Therefore

```text
sum_i ceil(q_i/7)/q_i < 1                                 (14)
```

certifies a common weak-safe lift. This is the exact open-arc sharpening of
THM-761's multi-exception union bound; no independence between tails is used.

For `q>=2`,

```text
ceil(q/7)/q <= 1/2,
```

with equality only at `q=2`. Consequently a primitive thirteen-speed row

```text
S=aC union {w_1,w_2},       |C|=11,       a>=2             (15)
```

is weakly lonely at threshold `1/14`, except possibly when

```text
a=2,       w_1 and w_2 are odd.                            (16)
```

Indeed, if neither tail is divisible by `a`, lower-dimensional LRC gives a
core phase with margin at least `1/12`, and (14) fails to certify only when
`q_1=q_2=2`. Then `a/2` divides the core and both tails, so primitivity forces
`a=2`. If one tail is divisible by `a`, absorb it into the scaled core;
settled LRC for the resulting twelve-speed core gives margin at least `1/13`,
and Part B lifts past the remaining nondivisible tail. Both tails cannot be
divisible by `a` in a primitive row.

This is the LRC(14) entrance to the dyadic two-tail seam. THM-769 proves the
parallel capacity statement at threshold `1/13`, and THM-775 studies its
later dyadic deletion descent. Nothing here proves that (16) is empty.

## E. Finite-clock exceptional-residue atlas

For a core clock `N`, define

```text
E_N(C)={b mod N:
          14|br|_N<N for every r in A_N(C)}.               (17)
```

Assume `M(C)>=1/14`. If `aC union {w}` has no weak-safe phase, Part B forces
`a|w`. Writing `b=w/a`, the dilation invariance

```text
M(aC union {ab})=M(C union {b})                            (18)
```

then forces `b mod N in E_N(C)` for every clock `N`.

For a finite clock family `F`, put `L=lcm(F)` and

```text
E_F(C)={b mod L: b mod N in E_N(C) for every N in F}.      (19)
```

This compresses the infinite tail line to finitely many periodic residue
rays. If `E_F(C)` is empty, the whole one-tail family is closed. THM-2057's
small full-unit clocks are the special case where the clock forces the
zero residue, and intersections of clocks become its lcm divisibility tax.

If `M(C)=1/14` and `F` contains the reduced denominator of every point in the
finite weak-safe set of `C`, then (19) is exact:

```text
M(C union {b})<1/14       iff       b mod L in E_F(C).      (20)
```

THM-2058's boundary packet decomposition is an exact way to list those
denominators. In the strict bulk `M(C)>1/14`, a finite clock family need not
sample every safe component, so (19) is only a necessary sieve unless a
separate completeness argument is supplied.

## Guardrails and prior-art boundary

1. Do not reverse Part B. The condition `a|w` means only that universal tail
   saturation is unavailable. For `a=2`, `C={1,...,12}`, `w=26`, and `N=14`,
   the row is `2{1,...,13}` and `t=1/28` is safe, although `q=1`.
2. Do not promote a finite-clock sieve to an iff in the bulk. For `C={1}` and
   `F={2}`, the residue `b=2` kills the chosen clock, but `{1,2}` is safe at
   `t=1/3`.
3. The one-tail existence consequence and hereditary-primitivity reduction
   predate this theorem: see THM-760, THM-761, and THM-765. The precise new
   bridge is from their sheet orbit to THM-2059's individual CRT histogram
   bins, plus the sharp count and the finite-clock composition (17)--(20).

## Assumption challenge and tournament analysis

The natural vertices are the `a` lift sheets, not runners, arcs, or pair-sum
rulers. Each tail contributes a danger hyperedge on this common vertex set.
The quotient to the orders `q_i` preserves capacity but destroys translations
and overlaps; the THM-2059 histograms restore exactly the reduction-class
sidecar needed for one tail. A tournament orientation is artificial because
danger incidence is symmetric and several tails may hit the same sheet. Only
at the saturated two-sheet seam does a genuine binary ownership gauge emerge:
the two odd tails must own opposite sheets. The useful fingerprint is the
covering deficit `a-|union_i D_i|`, not a runner score sequence.

## Computational audit

The companion script uses integer arithmetic only. It checks the valuation
identities, every tail-histogram bin over a bounded exhaustive box, the sharp
danger count, the packet lower bound against direct enumeration, the two-tail
capacity edge, the exact boundary atlas for `{1,...,13}`, and both guardrail
examples. It uses explicit runtime checks under both ordinary Python and
`python -O`; the frozen output ends in `PASS`.
