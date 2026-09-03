# LRC(14): the anchor strip is a Euclidean-remainder operator

**Session status:** one analytic theorem, one exact restricted-covariance
normal form, one sharp current-only cubic envelope, and one decisive
second-moment hostile. LRC(14) remains open.

## Inheritance pass and concept board

The closest proved mechanism was THM-4228's periodic-observable endpoint
cocycle. Its BV estimate says a fast pullback is approximately independent
of a finite-component region. The present anchor strip is special enough to
compute the error exactly. THM-4233 already contains a period-seven
quasipolynomial, but for a different moving primitive two-comb observable;
THM-859 carries exact common-dilation conjugacy and ramification, but not the
intersection with an undilated anchor strip.

The hostile was already visible in the current session: common odd dilation
leaves the full half-turn current energy unchanged, while the anchor-relative
core energy moves. The corrected near miss was therefore

```text
core energy = (6/7) full energy.
```

That is exact only on one Euclidean residue class. The least-used sidecar was
the address of the thin strip after pushing it through the common dilation.

The live board was

```text
Brownian residue current | anchor strip | Euclidean quotient/remainder
q-cubic nonlinear rebate | absolute-current histogram | addressed walls.
```

## 1. Exact polyphase/grid identity

For `A(x)=1_{||x||<1/14}`, an odd multiplier `m`, an even anchor `2h`, and
any half-periodic `f`, put

```text
c=gcd(m,2h)=gcd(m,h),   D=m/c,   H=2h/c,
D=7q+r,                 0<=r<=6.
```

Pushing Haar measure through `x=mt` turns the anchor indicator into the
average of `A` on a shifted `D`-grid. That average is not merely close to
`1/7`; it is exactly

```text
[q+1_{||Hx+q/2||<r/14}]/D.                            (1)
```

Thus, on ordinary `dt` over `[0,1/2)`,

```text
integral f(mt)A(2ht)
 =(q/D) integral f +(1/D) integral_R f,               (2)

R={||Hx+q/2||<r/14},        |R|=r/14.                 (3)
```

The parity of `q` is load-bearing: it decides whether the residual comb is
centered at zero or at the half phase. Formula `(1)` is the exact aliasing
law of a duty-cycle `1/7` square wave. In signal-processing language, the
anchor restriction is a polyphase decimator with a single quantization
residual; this is a literal identity, not just an analogy.

More generally, a wall with radius numerator `s/14` and half-phase bit
`epsilon` remains one wall after a degree-`D` average. If `sD=7Q+R`, its
state updates by

```text
(s,epsilon) -> (R,D*epsilon+Q mod 2).                 (3a)
```

Sequential updates compose exactly, so nested restrictions form a finite
fourteen-state Euclidean transducer rather than a growing wall tree.

Three controls expose the geometry:

```text
D=7q:             exact independence, no residual;
D=7q+1, q even:  the original narrow anchor wall returns;
D=7q+6, q odd:   the residual is the complement of that narrow wall.
```

## 2. Every restricted pair covariance reduces to one primitive wall

For actual odd speeds `a,b`, write `m=gcd(a,b)`, `a=mu`, `b=mv`. The
observable `sigma_u sigma_v` is half-periodic. Formula `(2)` gives

```text
integral_strip sigma_a sigma_b
 =q K(u,v)/(2D)+(1/D) integral_R sigma_u sigma_v,      (4)
```

where `K` is the already known full-circle Brownian residue kernel. Summing
`(4)` over ordered pairs gives the exact anchor-strip current energy for an
arbitrary twelve-tail family.

This is a better statement of the remaining debt than “a triple Fourier
correlation”: each pair is reduced to its primitive ratio plus one labelled
nested wall. The labels cannot generally be dropped, because different
pairs have different `D,H,q,r`.

The stochastic-process interpretation is now precise. The rank-three
Brownian kernel is the unconditional covariance. Conditioning on anchor
danger adds the conditional-bias term `integral_R sigma_u sigma_v/D`. Exact
independence occurs when the reduced fibre degree is divisible by seven;
otherwise the failure of independence is a single residual event.

## 3. Sharp bounds and what they do not prove

For a nonnegative `f` with half-base mass `F`, `(2)` gives

```text
qF/D <= strip(f) <= (q+1)F/D,                          (5)

|strip(f)-F/7|
 <=max(r,7-r)F/(7D)                    (r!=0).          (6)
```

If `0<=f<=M`, the residual-set mass `r/14` yields the sharper range-aware
upper term `min(F,Mr/14)/D`. These bounds are sharp for arbitrary measurable
observables with the given range and mass; no physical current equality is
claimed.

For `f=1_{min(a,b)=0}`, the same formula and the LRC-through-thirteen
Lipschitz floor give the valid cofinal estimate `D>91 max(B)`. This is not a
new LRC family. After primitive normalization the row is `DB union {H}` with
`gcd(D,H)=1`, and THM-616 already proves the much stronger statement for
every `D>1`. The q-zero calculation is useful only as a normalization and
operator control.

## 4. The second moment is the wrong stopping point

Let

```text
p(n)=1-n+binom(n,2)-binom(n,3),
C=a-b.
```

The nonlinear part of `p(min(a,b))` is `|p(a)-p(b)|/2`. Minimizing that
quantity over the forgotten sum `a+b`, at fixed `d=|C|`, gives the exact
current-only envelope

```text
g(0)=g(1)=g(2)=0,   g(3)=1/2,
g(d)=d(d^2-6d+11)/12                    (d>=4).         (7)
```

This is sharp pointwise. It strictly refines the variance rebate
`max(C^2-4,0)/12`, while still forgetting total depth.

On the primitive denominator-complete control

```text
h=420,
B=(1,3,5,7,9,11,13,15,17,19,21,45),                  (8)
```

the exact variance certificate is negative:

```text
B3+(V_core-12/7)/12 = -71225183/162954792.
```

Even integrating the positive part exactly remains negative:

```text
B3+(1/12) integral max(C^2-4,0)
 =-3392354543/9777287520.
```

But the sharp current histogram succeeds:

```text
B3+integral g(|C|)=2572403/33948915>0,
q-cubic exact        =526429/5819814>0.                (9)
```

So a perfect anchor-strip energy theorem would still not complete the
variance bridge on this canonical row. The next statistic is not necessarily
the full joint `(a+b,a-b)` law: the one-dimensional absolute-current
histogram already contains enough signal here.

## 5. Exact dilation hostiles and controls

For the base `(8)`, the full half-current energy is

```text
V=35666/15015,
```

and the `m=1` strip energy is `S=15361/45220`. The exact remainder operator
gives

```text
strip(13B) =(2V-S)/13,
strip(49B) =V/7        after gcd reduction D=7,
strip(127B)=(18V+S)/127.                               (10)
```

Thus the full pair energy is identical on the whole common-dilation ray, but
the anchor strip is not. The `D=7` equality is a positive exact-independence
control; `D=13` and `D=127` are complementary and nested-wall controls.

## 6. Next sharp problems

1. **Residual-bank compression.** Group the 144 ordered tail pairs by their
   exact `(D,H,r,q mod 2)` wall. Test whether the Brownian three-feature
   factorization survives within each group, or whether owner signs are the
   minimal additional sidecar.

2. **Signed residual cancellation.** Absolute bounds on `(4)` throw away the
   same signal that earlier approaches lost. Compute the residual pair Gram
   matrices for the canonical hostiles and ask whether their negative modes
   cancel by residue class or addressed tooth.

3. **Exploit the proved transducer.** The radius/phase update `(3a)` is now
   exact. Compute the orbit types actually generated by the pair-gcd bank and
   ask whether the third-tooth handoff words factor through those fourteen
   states. Calling it a continued-fraction mechanism still requires a map
   from addressed tooth transitions to the transducer.

4. **Current-histogram forcing.** Find a low-order representation of
   `integral g(|C|)` using tail probabilities `P(|C|>=j)` or discrete convex
   increments. The target is an addressed renewal inequality for those tail
   probabilities, not another second-moment estimate.

5. **Hostile discipline.** Any proposed compression must retain the three
   exact controls `(10)` and the negative variance margins `(9)`. A scalar
   full-energy or independence surrogate is already refuted.

The theorem artifact is
`01-canon/theorems/THM-4345-lrc14-halfperiodic-anchor-strip-euclidean-remainder-and-current-envelope.md`;
the exact wall audit is
`04-computation/lrc14_anchor_current_strip_probe_anchor_current_20260902.py`.
