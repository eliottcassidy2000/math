# Independent referee report: proposed THM-4089

## Verdict

**PASS**, provided the theorem is stated as a global optimization and
four-cell obstruction for the **exact pinned modular-template margin formula**.
It is not a theorem about the truth or falsity of any p-adic zeta
irrationality statement, and it is not a no-go for other continuation
templates or irrationality methods.

The audited external source is the detached checkout at commit
`b46a1770901551961710e155d775aae7c5ea39e7`.  The public branch still points
to that commit at audit time.  The proposed result uses only the formulas
displayed there; it does not validate the draft's BGG, Hasse, torsor, Bost,
continuation, or Jensen interfaces.

## Symbolic calculus audit

Write

```text
m=s+1,  c=(p+1)/12,  K=floor(s/xi),
I(xi)=(s+1/2) H_K-K xi+(s+1-(K+1)xi)_+^2/(2(K+1)),
J(xi)=integral_0^xi (1-cx)_+ dx,
F(xi)=s^2 xi-(s-1)J(xi)+s I(xi).
```

Since `tau=2F/(s+1)^2`, minimizing `tau` is equivalent to minimizing
`F`.  On the `K=k` cell,

```text
I'(xi)=(k+1)xi-(s+k+1)   while the positive part is active,
I'(xi)=-k                 after it switches off.
```

At `xi=(s+1)/(k+1)` the values and derivatives join, with derivative
`-k`.  At a floor seam `xi=s/k`, where `K` changes from `k` to `k-1`,
the two values are

```text
(s+1/2)H_k-s,
(s+1/2)H_(k-1)-s+s/k+1/(2k),
```

and are equal because `(s+1/2)/k=s/k+1/(2k)`.  Both derivatives are
`-k`.  This includes the `K=1` to `K=0` seam at `xi=s`; on the last cell
`I=(s+1-xi)^2/2`.  Thus `I` is `C^1` and convex on `1<xi<s+1`.

The function `J` is `C^1` and concave: its derivative decreases linearly
from `1` to `0` and then stays zero.  Hence `F` is convex.  In the first
floor cell (`K=s-1`, positive part active, Hasse integral unsaturated),

```text
F'(xi)=[s^2+(s-1)c]xi-(s^2+s-1).
```

Its zero is

```text
xi*=(s^2+s-1)/(s^2+(s-1)c)
   =12(s(s+1)-1)/(12s^2+(s-1)(p+1)).
```

For `p in {2,3,5,7}`, `c<1`, and the three required cell inequalities reduce
to transparent positive quantities:

```text
xi*>1                    iff (s-1)(1-c)>0,
xi*<(s+1)/s              since (s+1)den-s*num=s+(s^2-1)c>0,
c xi*<1                  iff c s^2<s^2.
```

The derivative of a convex function vanishes there.  It cannot vanish on a
second interval: on an active piece its slope is positive, while on a flat
piece of `I` the residual derivative is already positive.  Therefore `xi*`
is the unique global minimizer of `tau` on the full source domain.

For the `Y` variable put `f(x)=acos(x)-x sqrt(1-x^2)` and
`G(Y)=s Lambda_p(Y)-C_p(Y)`.  Since `f'(x)=-2sqrt(1-x^2)`,

```text
G'(Y)=2 pi[-s+4 sum_(p|q, qY<1) phi(q)/q sqrt(1-q^2Y^2)].
```

At every collision boundary `Y=1/q`, the departing layer has both value and
first derivative zero.  Thus `G'` is continuous.  The always-active `q=p`
term is strictly decreasing on `(0,1/p)`, while every other active term is
nonincreasing, so `G'` is strictly decreasing.  Its right limit is `-2pi s`.
Its left limit is `+infinity`: the subfamily `q=p*2^j` already contributes a
fixed positive amount per active index with `qY<=1/2`.  Hence `G` has a unique
global maximizer.

If `G'(L)>0>G'(R)`, monotonicity gives `G(Y)<=G(L)` to the left and
`G(Y)<=G(R)` to the right.  Concavity gives on `[L,R]`, and also at `R`,

```text
G(Y) <= G(L)+G'(L)(R-L).
```

Therefore the verifier's tangent expression is a valid global upper bound,
not a sampled maximum.

## Numerical and implementation audit

The candidate verifier was read line by line.  Its rational floor tests,
active-layer tests, fixed-point square roots, Machin enclosure for `pi`,
positive atanh enclosure for logarithms, alternating atan enclosure, interval
sign handling, and upward decimal rounding are correct in the exercised
domains.  Normal, optimized, and frozen line transcripts agree.

Candidate artifacts:

```text
6a4063e4a041fe30f419d34a594ea8e56b32fcf7315bd4137ca4f0eae71c4a70  next_case_margin_audit.py
cbe529a95384e41943adca5c45e546dc8484431e3e78349b18b151ece35c3f54  next_case_margin_audit.out
```

I also wrote a genuinely independent verifier.  It imports neither the
candidate verifier nor the pinned certificate.  It uses exact `Fraction`
interval arithmetic, a different atan range reduction, exact alternating and
positive-series remainder bounds, and independent integer-square-root
enclosures.  It checks every symbolic seam for each `s`, both derivative-sign
brackets, and the four tangent bounds.  Normal and `-O` transcripts match its
frozen output.

```text
8a65db82de9911565b18ca1b522cd87566112a33f4798875afbbf34f544f89c5  independent_margin_referee.py
6199c1ef483ab385d79f63746226c9e45c6b0265b6657f9e21dfb375c1c5c2e7  independent_margin_referee.out
```

Its exact interval endpoints imply:

```text
(2,31): max M < -1.943953720741976712...
(3,13): max M < -1.655957196706988949...
(5,7):  max M < -3.841753001089363046...
(7,5):  max M < -3.609764970278347647...
```

These agree with the candidate's outward-rounded 15-digit upper bounds.

Reproduction:

```text
python -B .scratch/padic_repo_audit_20260825/next_case_margin_audit.py
python -B -O .scratch/padic_repo_audit_20260825/next_case_margin_audit.py
python -B .scratch/padic_margin_referee_20260825/independent_margin_referee.py
python -B -O .scratch/padic_margin_referee_20260825/independent_margin_referee.py
```

## Required canon wording

1. Quantify the optimization over exactly
   `p in {2,3,5,7}`, odd `s>=3`, `1<xi<s+1`, and `0<Y<1/p`.
2. Define the full saturated integral `J_p`, even though `xi*` lies before its
   saturation wall.  Do not present the unsaturated quadratic as the global
   definition.
3. State that `xi*` is the unique global minimizer of `tau`, and that each of
   the four printed constants is an outward-rounded upper bound for the
   global maximum over both `xi` and `Y`.
4. Call the conclusion an obstruction to obtaining a **positive margin from
   this pinned modular-template formula by tuning `xi,Y`**.  It does not prove
   rationality, does not refute irrationality, and does not rule out another
   admissible continuation template, a changed arithmetic cost, or another
   method.
5. Keep the external draft's 22 singleton irrationalities
   **AUTHOR-CLAIMED / UNREFEREED**.  This elementary theorem does not promote
   them and should not list them as proved dependencies.

With those scope clauses, promotion is mathematically sound.
