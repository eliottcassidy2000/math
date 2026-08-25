---
id: THM-4062
title: "LRC(14) affine-intercept obstruction to static divisor-star complement compression"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. Every divided
  exception mask factors through its exact denominator, primitive step, and
  a label-gauged affine intercept. THM-4056/4059 retain the first two but
  erase the intercept. Two primitive typed d=4 rows have identical complete
  static labelled divisor/depth packets and incidence (11,9,8), yet one
  covers all four selected lifts while the other leaks two. The same depth
  sign transfers exactly to THM-4042's AP owner word; at P=11 it expands the
  formal coefficient clock from 420 to 2520. This is an information-loss
  theorem and signed contrast, not an LRC(14) proof or cover improvement.
source: codex-frontier-synthesis-breakthrough-20260825 / LRC complement-compression lane
audit: >
  PASS after repairing the phase-representative gauge. The primary exact
  companion checks the typed rows, strict masks, label covariance, packet
  aliases, signed AP words, and all 2,520 phase vectors. An independent
  Euclidean-depth implementation rebuilds both rows and both signed/unsigned
  P=11 clocks directly from THM-4042 winner tracks, including 54 formal
  Fourier identities. Both normal/optimized pairs byte-match; both scripts
  have zero assert nodes and zero float literals.
depends_on:
  - THM-4024-lrc14-complete-divisor-incidence-envelope
  - THM-4030-lrc14-d4-affine-defect-lattice-boundary
  - THM-4042-prime-sector-ap-cover-exact-clock-and-holonomic-law
  - THM-4056-divisor-phase-compiler-and-duffin-schaeffer-lcm-clock
  - THM-4059-stern-brocot-depth-packet-character-and-divisor-star-convolution
related:
  - THM-4004-lrc14-three-detuned-divisor-comb-profile
  - THM-4052-lrc14-affine-component-width-escape-cones
  - THM-4057-stern-brocot-depth-pullback-and-rational-edge-tournament-gauge
script: 04-computation/lrc14_affine_intercept_divisor_star_obstruction_thm4062.py
output: 05-knowledge/results/lrc14_affine_intercept_divisor_star_obstruction_thm4062.out
independent_audit_script: 04-computation/lrc14_affine_intercept_divisor_star_obstruction_thm4062_independent_audit.py
independent_audit_output: 05-knowledge/results/lrc14_affine_intercept_divisor_star_obstruction_thm4062_independent_audit.out
script_sha256: aab25543b47d72307891b12d87c0cb74ecdb88f13c444e95ed7a3610807acfc1
output_sha256: d352a8ba923b30cbb0227718120207bf63d6def886b36096704a3fb2ba0cbb48
independent_audit_script_sha256: 4d29ccb89cbfef22b8aede289aa9dfcd33a47fdd96ccf9f3e2cbff361e86bde6
independent_audit_output_sha256: 29167744830c288ad103e84e64a6a3d941d91627936795a465c7d7362835b7c8
hash_basis: raw LF bytes
---

# THM-4062 -- the missing complement coordinate is an affine intercept

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.** THM-4024
compresses a hypothetical `11+2` row by divisor incidence, THM-4056 retains
each primitive exact-denominator packet, and THM-4059 adds its full Stern
depth sign. This theorem determines the precise coordinate that all three
still lose when one asks whether exceptions cover a selected family of
divided LRC lifts.

## 1. Label-gauged affine-intercept factorization

Fix `d>=1`, a positive integer exception `delta`, and a **chosen real
representative** `y` of a divided-pack phase. The labelled lifts are

```text
x_j=(y+j)/d,                 j in C_d.                 (1)
```

The representative and the labels form one gauge: replacing `y` by `y+1`
rotates every label by minus one. Put

```text
g=gcd(d,delta),
q=d/g,
a=delta/g mod q,
tau_delta(y)=delta*y/d mod 1.                         (2)
```

Use THM-4056's convention `U_1={0}` when `q=1`; otherwise `a in U_q`.
For every label,

```text
delta*x_j=tau_delta(y)+a*(j mod q)/q mod 1.           (3)
```

Indeed, write `delta/g=a+qL`. Then
`delta*j/d=(a+qL)j/q`, and `Lj` is integral. Thus the strict danger mask is
exactly

```text
D_delta(y)=pi_q^(-1){r in C_q:
  ||tau_delta(y)+a*r/q||<1/14},                       (4)
```

where `pi_q:C_d->C_q` is reduction.

The phase gauge is load-bearing. The scalar `delta/d` does not define a map
on `R/Z` unless it is integral. Instead `(4)` is covariant:

```text
D_delta(y+1)=D_delta(y)-1 in C_d.                     (5)
```

Equivalently, `tau_delta(y+1)=tau_delta(y)+a/q mod 1`, which is exactly the
same relabelling in `(4)`.

For an exception family `E`, the exact survivor count is

```text
N_(E,d)(y)=sum_(j in C_d) product_(delta in E)
 (1-1[||tau_delta(y)+a_delta*(j mod q_delta)/q_delta||<1/14]).        (6)
```

The selected phase is fully spoiled iff `N_(E,d)(y)=0`. Formula `(6)` is the
exact nonlinear replacement for a scalar packet-star test.

THM-4056's labelled packet records `(q,a)`, and THM-4059 determines
`epsilon(a,q)` and every twisted packet statistic from that static data.
Neither records `tau_delta(y)`. For one labelled phase, the intercepts are
sufficient; across a pack-safe component one must retain the covariant affine
sections `y->delta*y/d`, equivalently the exception heights relative to the
selected phase gauge.

## 2. A same-packet, opposite-coverage `d=4` hostile

Take `d=4`, divided pack `H={1,...,10}`, and representative `y=1/11`. Its
four lifts are

```text
x_j=(1+11j)/44,                 j=0,1,2,3.            (7)
```

Compare the exception triples

```text
E =(2,9,11),
E'=(2,1,3).                                            (8)
```

Their owner-labelled residues modulo four agree term by term:

```text
(2,9,11) mod 4=(2,1,3)=(2,1,3) mod 4.                (9)
```

Hence both have the identical inverse THM-4056 packet list

```text
((2,1),(4,1),(4,3))                                   (10)
```

and identical THM-4059 depth signs

```text
(-1,-1,-1).                                           (11)
```

More strongly, the complete static residue labels agree, so every
denominator weight, Ramanujan mode, depth weight, or twisted-depth mode that
factors through this labelled `C_4` compiler agrees.

Attach either triple to the same eight body owners

```text
U_0=(8,12,20,24,28,32,36,40)                         (12)
```

and pair speeds `(4,16)`, corresponding to
`s=1,t=4,(p_0,q_0)=(1,4)`. Both thirteen-speed rows are primitive and
pairwise distinct. Their eleven-owner bodies have THM-4024 incidence

```text
(c_1,c_2,c_4)=(11,9,8),                              (13)
```

and the eight divided body owners together with the divided pair are exactly
`H`.

Direct strict-mask arithmetic gives

```text
D_2(1/11) ={0,2},
D_9(1/11) ={3},
D_11(1/11)={1},

D_2(1/11)={0,2},
D_1(1/11)={0},
D_3(1/11)={0}.                                       (14)
```

The first union is all of `C_4`; the second is `{0,2}` and leaks labels
`1,3`. The close distances are

```text
delta=2:  labels 0,2 have distance 1/22;
delta=9:  label 3 has distance 1/22;
delta=11: label 1 has distance 0;
delta=1:  label 0 has distance 1/44;
delta=3:  label 0 has distance 3/44<1/14.             (15)
```

All omitted distances exceed `1/14`. Every pack speed `4h` has distance
`||h/11||>=1/11` on every lift, so labels `1,3` are genuinely safe for the
whole second row. The first row is THM-4030's selector hostile, not an LRC
counterexample: at `x=21/22` its full clearance is `1/11`.

The lost coordinates appear explicitly as

```text
tau_9=9/44 versus tau_1=1/44,
tau_11=1/4 versus tau_3=3/44.                         (16)
```

Therefore no selected-phase complement-compression map factoring only
through the complete static labelled divisor/depth packet can preserve
coverage. This does not rule out compression retaining the affine intercept
field, the pack-safe components, or the finite THM-3818 producer.

## 3. The depth-signed AP owner word

There is also a positive transfer from THM-4059 into THM-4042's eventual AP
coefficient controller. Fix a prime `P` and `2<=q<P`, and put

```text
t=P^(-1)-1 mod q.                                     (17)
```

Define the signed word

```text
w^epsilon_(P,q)(c)
 =(P-q)/(Pq) epsilon(c,q)1[c in U_q]
 +(P-q-1)/(Pq) sum_(u in U_q,tu=c)epsilon(u,q).       (18)
```

For `q=1`, separately set

```text
w^epsilon_(P,1)(0)=(2P-3)/P.                         (19)
```

This is THM-4059's declared `+1` zero-phase bookkeeping weight, not a Stern
depth assigned to `0/1`.

THM-4042's plus owner map is `a->-a^(-1)`. Modular inversion and the subtree
mirror preserve full depth, so the first pushforward in `(18)` has weight
`epsilon(c,q)`. The minus winner gives the multiplication fibre `tu=c`,
which proves the second term.

Use the positive-exponent convention

```text
hat(f)(k)=sum_(c in C_q)f(c)exp(2*pi*i*k*c/q),         (20)
E_q(k)=sum_(a in U_q)epsilon(a,q)exp(2*pi*i*k*a/q).   (21)
```

Then `(18)` has the exact Fourier transform

```text
hat(w^epsilon_(P,q))(k)
 =(P-q)/(Pq)E_q(k)+(P-q-1)/(Pq)E_q(kt).               (22)
```

This is a formal signed contrast of the eventual THM-4042 coefficient word.
It is not a nonnegative cover measure.

## 4. The signed `P=11` clock restores the full LCM

For `0<=c<=9`, define the eventual formal coefficient controller

```text
A^epsilon_c(r)=sum_(q=c+1)^10
  w^epsilon_(11,q)(c-r mod q).                        (23)
```

Every block is `q`-periodic, so `lcm(1,...,10)=2520` is a period. Four exact
blocks force the reverse divisibility:

```text
q=5: (0,6,-6,-6,6)/55,                    period 5;
q=7: (0,1,1,1,1,1,1)/11,                 period 7;
q=8: (0,-3,0,3,0,3,0,-3)/88,             period 8;
q=9: (0,1,-3,0,-1,-1,0,-3,1)/99,         period 9.   (24)
```

If `T` is a period of the full vector, coordinate `c=9` isolates the `q=10`
block. Descending in `c` and subtracting already isolated higher-`q` blocks
shows that `T` is a period of every word. Therefore

```text
lcm(5,7,8,9)=2520 divides T.                          (25)
```

Thus the signed controller has exact period `2520`. THM-4042 proves that the
unsigned `P=11` controller has period `420`; the depth sign restores the
missing 2- and 3-power depth through `q=8,9` and expands the clock by six.
This does not improve AP coverage or give an LRC certificate.

## 5. Audit and boundary

The primary companion verifies all 2,520 signed phase vectors and the strict
`d=4` masks. The independent companion computes depth from the Euclidean
algorithm rather than THM-4059's inverse formula, rebuilds THM-4042's owner
tracks, proves `(22)` in the rational group ring at every frequency, and
independently obtains signed period `2520` and unsigned period `420`.

Reproduce both normal/optimized pairs from the repository root with

```text
python3 -B 04-computation/lrc14_affine_intercept_divisor_star_obstruction_thm4062.py
python3 -B -O 04-computation/lrc14_affine_intercept_divisor_star_obstruction_thm4062.py
python3 -B 04-computation/lrc14_affine_intercept_divisor_star_obstruction_thm4062_independent_audit.py
python3 -B -O 04-computation/lrc14_affine_intercept_divisor_star_obstruction_thm4062_independent_audit.py
```

The result is a sharp no-go for **static** complement compression. It does
not close the physical finite producer, classify every `d=4` row, compare
whole pack-safe components, prove AP extremality, or prove LRC(14). **QED.**
