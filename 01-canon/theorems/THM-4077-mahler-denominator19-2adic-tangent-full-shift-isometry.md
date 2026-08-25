---
id: THM-4077
title: "Mahler denominator-19 2-adic tangent full-shift isometry"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY VERIFIED-EXACT. For every fixed
  reset depth s, THM-4074's normalized denominator-19 programmer extends
  uniquely to an onto isometry of Z_2. Composed with the THM-2228 carry
  coordinate, it is an onto isometry from odd 2-adics to all infinite carry
  words beginning in one. Positive odd parameter points are exactly the words
  realized after the runway by one fixed launch in this family; general
  2-adic points require moving starts and moving observation times converging
  to a tangent chart at -9/19. Parameter termination and output-state
  termination are distinct dense loci. Their separate intersections with the
  Haar-null strict safe set remain open, so no Z-number is produced or ruled
  out.
source: codex-frontier-synthesis-creative-20260825d / Mahler inverse-limit lane
audit: >
  PASS. The primary affine-recurrence path checks 84 complete residue
  permutations for 0<=s<=6 and heights 1<=h<=12, including 57,330 parity
  gates, 28,658 one-bit isometry gates, and all 7,161 odd-start carry
  cylinders through h=10. The independent binary-geometric/Hensel path lifts
  three coherent hostile targets through 64 bits for every s, performs 2,646
  candidate gates and 448 safe-logarithm identity gates, then directly replays
  36 ordinary denominator-19 launches through 240,804 orbit steps. Normal and
  optimized outputs byte-match; both scripts have zero assert nodes and zero
  floating literals.
depends_on:
  - THM-2228-mahler-three-halves-carry-tail-and-integral-stabilization
  - THM-4072-mahler-safe-terminal-fibre-product-and-finite-state-obstruction
  - THM-4074-mahler-denominator19-postterminal-arbitrary-delay
related:
  - THM-3848-rational-base-prefix-atom-tree-and-lonely-runner-separation
script: 04-computation/mahler_denominator19_2adic_tangent_isometry_thm4077.py
output: 05-knowledge/results/mahler_denominator19_2adic_tangent_isometry_thm4077.out
script_sha256: 594ba68d0b47e5ebaa650cfceac679eabafa27f516f3ba7fb934c2932c5502ae
output_sha256: 39a093d7ad26d73829b9ebf8af549be8838d2266f4ae43dc26387e9f806848d0
independent_audit_script: 04-computation/mahler_denominator19_2adic_tangent_isometry_thm4077_independent_audit.py
independent_audit_output: 05-knowledge/results/mahler_denominator19_2adic_tangent_isometry_thm4077_independent_audit.out
independent_audit_script_sha256: 07aeb83e5032913f3c20aa1a3e2764fa2bfcce272e3b799b949e259eb758f023
independent_audit_output_sha256: 4d794c4579f4f6da92cfb3a12d91a34321b5b8a9afab43a5f9835cae6129c8b4
hash_basis: raw LF bytes
---

# THM-4077 -- the finite programmer closes to a 2-adic tangent chart

**PROVED + VERIFIED-EXACT + INDEPENDENTLY VERIFIED-EXACT.** THM-4074
programs every finite odd-start carry cylinder, but changes the launch with
the requested horizon. The compatible finite permutations do have an exact
inverse limit. This theorem identifies that limit and proves why it is a
moving-start/moving-time chart rather than one infinite orbit.

## 1. The normalized map is an onto isometry of `Z_2`

Fix `s>=0` and put

```text
L=18*2^s,                 k=s+3,
g=3^L,
U=9*3^k*(g-1)/(19*2^k).                                  (1)
```

The order of `2 mod 19` and LTE, as used in THM-4074, show that `U` is a
positive odd integer and `g=1 mod 4`. For `n>=0`, define

```text
F_s(n)=U*[n]_g=U*(g^n-1)/(g-1),       F_s(0)=0.          (2)
```

If `n>n'`, then

```text
F_s(n)-F_s(n')
 =U*g^n'*(g^(n-n')-1)/(g-1).                             (3)
```

LTE for the odd base `g=1 mod 4` gives

```text
boxed: v_2(F_s(n)-F_s(n'))=v_2(n-n').                   (4)
```

Thus `F_s` is an isometry on the dense subset of nonnegative integers in
`Z_2` and has a unique continuous extension

```text
Fhat_s: Z_2 -> Z_2.                                      (5)
```

For every `h`, `(4)` makes the reduction of `F_s` an injection on the
`2^h` residue classes, hence a permutation. The compatible finite
permutations show that `(5)` is onto. Passing `(4)` to the inverse limit shows
that it remains an isometry, and in particular preserves parity.

The analytic notation

```text
Fhat_s(t)=U*(g^t-1)/(g-1)                                (6)
```

is legitimate for `t in Z_2`: expand `g^t` by the convergent binomial series
around `g=1`. It also gives the continuous affine equation

```text
Fhat_s(t+1)=g*Fhat_s(t)+U.                               (7)
```

## 2. Every infinite odd-start carry word has one parameter

Let `Phi:{0,1}^N -> Z_2` be the THM-2228 carry homeomorphism. Its finite
cylinder theorem says more precisely that two words share their first `h`
digits exactly when their `Phi` values agree modulo `2^h`; also
`c_0=Phi(c) mod 2`. Hence `Phi` is an isometry for the binary prefix metric.

Write

```text
[1]={c in {0,1}^N:c_0=1},
Z_2^odd={t in Z_2:t=1 mod 2}.                             (8)
```

Then

```text
boxed: Psi_s=Phi^(-1) o Fhat_s : Z_2^odd -> [1]          (9)
```

is an onto isometry. This is a rooted-cylinder coding. It is not a conjugacy
between addition or the binary shift and the Mahler map.

For `c in [1]`, let `t_h` be the unique odd representative

```text
1<=t_h<2^h,               F_s(t_h)=Phi(c) mod 2^h.      (10)
```

The finite permutations give

```text
t_(h+1) congruent t_h (mod 2^h),
t_h -> t=Psi_s^(-1)(c).                                    (11)
```

These are exactly the compatible THM-4074 programming classes.

## 3. Fixed launch versus tangent blow-up

For each representative in `(10)`, set

```text
m_h=L*t_h,
A_(s,t_h)=9*(2^m_h-1)/19.                                (12)
```

This is a positive integer. THM-4074's exact launch calculation gives

```text
T^(m_h+k)(A_(s,t_h))=F_s(t_h),                           (13)
```

where `T(a)=ceil(3a/2)`. Therefore the next `h` carries at the moving
observation state `(13)` are the first `h` digits of `c`.

For a fixed positive odd `t`, the inherited identities give the exact full
word factorization

```text
c(A_(s,t))=(100)^(L*t/3) 0^k Psi_s(t).                     (13a)
```

Each displayed `100` returns the strict-safe follower graph to state zero,
as does each following zero. Thus the infinite tail begins at the same rooted
safe state, not merely at a carry word with the right residue.

Because canonical representatives can only stay fixed or increase by `2^h`,
the following are equivalent:

```text
t_h is eventually constant;
t is a positive odd ordinary integer;
the LSB-first parameter bits of t are eventually zero;
one fixed launch A_(s,t) realizes the whole post-runway word c.             (14)
```

If `(14)` fails, `t_h` is nondecreasing and unbounded, hence `m_h -> infinity`.
Yet

```text
v_2(A_(s,t_h)+9/19)=m_h,                                  (15)
F_s(t_h) -> Fhat_s(t).                                    (16)
```

Thus all moving launches converge to the same formal state `-9/19`, their
observation times diverge, and their shifted states converge to the chosen
chart value. The full shift in `(9)` is precisely this renormalized tangent
chart. It is not the 2-adic continuation of the original launch map and does
not describe one fixed orbit when `(14)` fails.

Let

```text
N_odd={1,3,5,...},
D_s=Psi_s(N_odd).                                          (17)
```

Then `D_s` is exactly the fixed denominator-19 post-runway word locus in
`[1]`. It is countable, dense, and co-dense: positive odd integers are dense
in odd `Z_2`, while nonpositive and nonintegral odd points meet every
cylinder. It is not the locus of all generic fixed integer launches; Section
4 identifies that larger locus.

## 4. Parameter termination is not output termination

There is a second eventual-zero coordinate:

```text
G_s=Fhat_s^(-1)(N_odd).                                    (18)
```

A parameter lies in `G_s` exactly when the **output/native** binary digits of
`Fhat_s(t)` eventually vanish, so that the coded carry word is the itinerary
of a genuine positive odd integer. Since `F_s(n)` is a positive odd integer
for every positive odd `n`,

```text
N_odd is a proper subset of G_s.                            (19)
```

For strictness, let

```text
tau_s=Fhat_s^(-1)(1).                                      (20)
```

This point is odd but is not a rational integer. Indeed `F_s(0)=0`, while
`F_s(n)>=U>1` for `n>=1`. For `n=-r<0`, the continuation of `(6)` is the
negative rational `U(g^(-r)-1)/(g-1)`, so it cannot equal `1`. Nevertheless

```text
Psi_s(tau_s)=Phi^(-1)(1),                                  (21)
```

the genuine integer itinerary beginning `1011`. THM-2228 shows that this
word is unsafe. Thus output termination does not repair the strict tail gate.
It also gives the strict word-locus inclusion

```text
D_s proper subset Psi_s(G_s)=Phi^(-1)(N_odd),               (21a)
```

with `tau_s` witnessing the difference.

Conversely, the strict safe periodic word

```text
c_*=(100)^infinity,             Phi(c_*)=-9/19,             (22)
```

has the explicit odd parameter

```text
tau_*s=Fhat_s^(-1)(-9/19)
      =log(1-(2/3)^k)/log(g) in Z_2.                         (23)
```

These are 2-adic logarithms. Both numerator and denominator have valuation
`k`, so the quotient is an odd 2-adic integer. It is not a rational integer.
Nonnegative inputs have nonnegative `F_s` values. If an input were `-r<0`,
`(23)` would force

```text
g^(-r)=1-(2/3)^k,
```

but the left side is at most `3^(-18)` and the right side is at least
`19/27`. The safe word therefore exists in the tangent chart but its native
state is still the nonordinary formal point `-9/19`.

## 5. The two exact open intersections

Let `S` be the THM-4072 strict safe carry set and define

```text
Acal_s=Psi_s^(-1)(S intersect [1]).                         (24)
```

This set is nonempty by `(22)`--`(23)`. The two termination questions are
now exact and visibly different:

```text
the full carry itinerary of some A_(s,t) in the fixed-s denominator-19
family lies in S, and hence supports a Mahler Z-orbit
  iff Acal_s intersect N_odd is nonempty;                   (25)

some positive odd integer starts a Mahler Z-orbit
  iff Acal_s intersect G_s is nonempty.                     (26)
```

Equation `(25)` uses the periodic safe launch prefix and the fixed-launch
boundary `(14)`. Equation `(26)` is the THM-2228 strict-tail plus ordinary-
state characterization transported through `(9)` and `(18)`. Both
intersections remain **OPEN**, with the first contained in the second by
`(19)`. The hostile points cross the coordinates:

```text
tau_s  in G_s but not Acal_s,
tau_*s in Acal_s but not G_s.                               (27)
```

Density of the termination loci cannot settle either intersection:
THM-3848 proves that `S` is Haar-null.

## 6. Exact audits and scope

The primary audit constructs `(2)` by the affine recurrence `(7)`, checks
every finite permutation and its one-bit lift, then composes it with a direct
implementation of the THM-2228 carry coordinate. The independent audit uses
binary-doubling geometric sums instead, lifts the safe periodic, state-one,
and all-one targets one bit at a time through height 64, checks the logarithm
identity in `(23)`, and directly runs the large ordinary launches in `(12)`.

Replay from the repository root:

```bash
python3 -B 04-computation/mahler_denominator19_2adic_tangent_isometry_thm4077.py
python3 -B -O 04-computation/mahler_denominator19_2adic_tangent_isometry_thm4077.py
python3 -B 04-computation/mahler_denominator19_2adic_tangent_isometry_thm4077_independent_audit.py
python3 -B -O 04-computation/mahler_denominator19_2adic_tangent_isometry_thm4077_independent_audit.py
```

Both outputs byte-match their frozen companions. The theorem closes the
inverse limit of the finite programmer; it proves neither intersection
`(25)` nor `(26)`, supplies no positive Z-number, and gives no universal
rejection theorem.
