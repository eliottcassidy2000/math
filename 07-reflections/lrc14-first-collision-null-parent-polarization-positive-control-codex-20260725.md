# LRC14 first-collision null-parent polarization positive control

**Status:** FINITE-EXACT + INDEPENDENTLY HOSTILE-AUDITED.

This note records a sharp positive control at the first repeated septimal
layer. It is not a scalar-row exclusion and is not proved canon. Its purpose
is to distinguish two possible failures:

```text
immediate deletion/word target supports fail to meet

versus

they meet exactly, but the common-right prefix does not transport to a
canonical terminal word.
```

The exact control rules out the first failure on the reduced THM-2377
same-layer tensor. The second remains open.

## Inheritance

The closest proved mechanisms are:

- THM-2377, `THM-2377-septimal-valuation-collision-and-bockstein-carry-gate`,
  whose same-layer control uses speeds `(13,1,7,1274)`;
- THM-2370,
  `THM-2370-deletion-martingale-drift-conservation-and-sharp-clone-hostile`,
  whose Section 6 gives linear deletion decomposition at a common bare right
  endpoint; and
- THM-2380,
  `THM-2380-cross-word-charged-target-correlation-and-pair-twist-gate`,
  whose target correlation detects exact same-target support overlap.

THM-2380's sharp disjoint-support hostile shows that separate drift floors do
not force a shared target. The present calculation is the opposite control:
the pre-collision parent is target-null, so its two children must have equal
and opposite charged currents.

## The reduced common-right split

Put

```text
d(y)=1_(||y||<1/14),             g=1-d.
```

Use `sigma` for the first target-shift coordinate, to avoid confusing it
with the speed `s_0=1`. Define

```text
B(r,sigma,t)
 =integral_T
   d(13x)
   d(1274x-r/13) g(1274x-t/13)
   g(7x+t/13) dx,

W(r,sigma,t)
 =integral_T
   d(13x) g(x)
   d(1274x-r/13) g(1274x-t/13)
   g(7x+t/13) dx,

A(r,sigma,t)
 =B(r,sigma,t)-W(r,sigma,t)

 =integral_T
   d(13x) d(x)
   d(1274x-r/13) g(1274x-t/13)
   g(7x+t/13) dx.                                  (1)
```

Thus `B=A+W` is a positive split at one common right endpoint. No factor
receives `sigma`; the reduced tensor is blind to the first target
coordinate.

Use THM-2365's lawful action

```text
T_(u,v):(r,sigma,t)->(r+v,sigma+u,t+v)
```

from base point `(1,0,0)`. Write

```text
B_(u,v)=B(1+v,u,v),

W_(u,v)=W(1+v,u,v),

A_(u,v)=A(1+v,u,v).                                (2)
```

Before the speed-one factor is inserted, the active septimal valuations
are strictly ordered:

```text
nu_7(13)=0 < nu_7(7)=1 < nu_7(1274)=2.
```

THM-2377's distinct-layer theorem therefore makes `B` target-null.
The direct interval count on the common endpoint grid `231868` gives

```text
B_(u,v)=78/8281,

W_(u,v)=N_v/8281,

A_(u,v)=M_v/8281,                                  (3)
```

where

```text
N=(78,74,71,71,71,71,71,71,71,71,71,71,74),

M=78-N=(0,4,7,7,7,7,7,7,7,7,7,7,4).               (4)
```

The integer cell counts are respectively

```text
B: (2184)^13,

W: (2184,2072,1988^10,2072),

A: (0,112,196^10,112).                              (5)
```

## Exact opposite-phase target overlap

Use THM-2380's normalized target transform. All coefficients with nonzero
first target coordinate vanish because (1) is `sigma`-blind. For
`q!=0` in the second coordinate,

```text
Nhat(q)
 =7+3(zeta^q+zeta^(-q))
 =7+6 cos(2 pi q/13)
 >0,

Mhat(q)=-Nhat(q).                                   (6)
```

The lower bound is strict because `7-6=1`. Hence

```text
what(0,q)=Nhat(q)/107653,

ahat(0,q)=-what(0,q),

ahat(0,q)conj(what(0,q))
 =-(7+6 cos(2 pi q/13))^2/107653^2<0.              (7)
```

All twelve nonzero target colours therefore overlap, in one strict
negative real half-plane.

The operator explanation is simpler than the trigonometry. After removing
target means,

```text
A^circ=-W^circ
```

because the parent is constant. The centered THM-441
correlation-adjoint identity becomes

```text
R^circ=-(1/169) W^circ star W^circ,

Rhat(q)=-|what(q)|^2.                               (8)
```

At this common-endpoint prefix, the target-null parent canonically fixes
the two children's relative sign, so no extra phase quadrature is needed
there. This is not an absolute terminal phase reference and does not
physically realize THM-2380's endpoint-matched pair twist after later
factors are restored.

## Exact correlation and energy

The normalized translated correlation is independent of its first
coordinate and equals

```text
R(delta_u,delta_v)
 =K_(delta_v)/(13*8281^2),                          (9)
```

with

```text
K=(5562,5587,5620,5629,5629,5629,5629,
   5629,5629,5629,5629,5620,5587).                 (10)
```

Thus

```text
R(0)=5562/891474493,

mean_delta R(delta)=5616/891474493,

R(0)-mean R=-54/891474493.                         (11)
```

Parseval gives

```text
sum_(q!=0)|Nhat(q)|^2=702,

sum_(q!=0)|Nhat(q)|^4=77766.                       (12)
```

Consequently the absolute nonzero-target Gram mass and THM-2380 complete
Cayley energy are

```text
54/891474493,

5982/10331448031704891637,                         (13)
```

respectively.

## Deep-colour refinement

For a fixed deep colour `a`, define

```text
J_a^W(sigma,t)
 =13^(-1) sum_r W(r,sigma,t) zeta^(ar),

calW_a(t)=zeta^(-at)J_a^W(0,t).                    (14)
```

If `j(k)` is THM-2377's deep difference profile and

```text
C_a=sum_k j(k)zeta^(ak),
```

then

```text
calW_a(t)=C_a N_t/8281,

calA_a(t)=C_a M_t/8281.                            (15)
```

Here

```text
C_0=144/91,

C_a=-(13+12 cos(2 pi a/13))/91<0       for a!=0.   (16)
```

The same strict anti-alignment therefore holds for all

```text
12*12=144
```

nonzero deep-colour/target pairs. This is still a finite colour response,
not one exact `(m,X)` marked triangle passing through both packets. These
are actual finite Fourier coefficients of the reduced physical tensor,
but the bank is derived from THM-2377's rank-one factorization; it is not
an independently realized owner-typed or terminal-word bank.

## The first failed implication and the surviving target

The retained child `W` is only the immediate post-collision THM-2377
tensor slice. It has no full `A_0` safe core, canonical THM-2305 owner,
terminal word, delayed word, or complete blocker-status packet. Equal
septimal valuations also do not canonically orient which tied factor is
inserted second.

THM-2370's common-right decomposition is exactly the scope in which
`B=A+W` is linear. Restoring later factors or replacing the bare right
endpoint by the canonical fully masked endpoint introduces cross-layer
terms. THM-2378 separately warns that a reduced lawful carrier need not
extend to a scalar cover.

Therefore the strongest honest conclusion is:

```text
on the immediate first-collision prefix:
  target support overlap and phase alignment are exact;

to reach the canonical terminal word:
  a canonical, owner-typed, tie-oriented, common-endpoint restoration
  filtration is still required, with target-shift covariance.         (17)
```

The cheapest next theorem should show that the common-right prefix is a
canonical terminal-word filtration stage, or that every later restoration
adds only target-null current. Another static Gram identity does not address
that transport. This is the first missing transport sidecar on this route,
not a claim that it is the unique or universally minimal sidecar.

No scalar row is excluded. The LRC(14) ledger remains `165`.

## Exact reproduction

Run

```bash
python3 04-computation/lrc14_first_collision_null_parent_positive_control.py
python3 -O 04-computation/lrc14_first_collision_null_parent_positive_control.py
```

and compare both transcripts, after LF normalization, with

```text
05-knowledge/results/lrc14_first_collision_null_parent_positive_control.out
```

The exact companion independently recounts the strict interval cells,
checks `B=A+W`, reconstructs (4)--(13), derives the twelve exact
cyclotomic deep transforms, and checks all `144` nonzero
deep-colour/target products. LF SHA-256:

```text
script  998c4a06314678c51dad66f03bf8cabd52297c1db20b4da207e75199ffd44b02
output  faa0a21d2ff62b2eb84e694c5a52cb7c9c88f894028453e5cce85482e87a1eac
```

Two independent hostile audits accepted the exact interval recount,
normal/optimized/stored replay, every displayed fraction, the repaired
THM-441 normalization in (8), all twelve cyclotomic deep transforms, all
`144` deep/target products, and the common-endpoint transport boundary.
