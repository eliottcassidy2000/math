# Collatz blueprint: energy obstructions and repaired analytic scope

**Date:** 2026-09-21. **Status:** PROVED elementary obstructions and repairs;
FINITE-EXACT controls; floating displays are diagnostics only. No Collatz
convergence theorem or classification of all signed cycles is claimed.

**Source being audited:** the attached `Pasted text.txt`, headed
“Grand Unified Geometric Collatz Proof,” sections 2 and 3. The positive
odd accelerated map throughout this note is
`T(n)=(3n+1)/2^{v_2(3n+1)}`. Its integer domain is not silently replaced by
a real differential equation or a rational dynamical system.

## Inheritance and live concept board

The closest proved mechanism is the arbitrary finite halving-word cylinder
and its squarefree refinement in
[the earlier synthesis](arithmetic_braids2_20260917_synthesis.md) and
[squarefree symmetry, §§2–3](arithmetic_braids2_20260917_squarefree_symmetry.md).
The canonical hostile is an arbitrarily long all-`k=1` growing prefix.
The corrected near miss is that a natural density of sources is not a
time average on each orbit. The useful sidecar is the exact exponent word
and its source congruence, retained alongside real height.

| Concept | Predicate and representation | Operation / decisive hostile |
|---|---|---|
| Proposed energy | Real function sampled on positive odds | Differentiate inside one dyadic band; inspect its boundary |
| Growth cylinders | `n=-1 mod 2^(L+1)` | Translate close to a phase minimum without changing the first L exponents |
| Finite-time descent | One decrease within a globally fixed number of steps | Long growth followed from the last energy minimum |
| Squarefree density | Source frequency versus orbit frequency | Fixed point 1; growing paths squarefree at every node |
| Mechanical words | Binary block balance | The supplied word contains both 00 and 11 |
| Cycle closure | Rational affine equation plus integrality | The exact rational 5/7,11/7 cycle |

The anchor is the energy claim; the niche is the two general obstruction
theorems; the wildcard is the finite-word/mechanical-word mismatch. Every
connection below keeps its domain, measure or lost coordinate explicit.

## 1. The actual derivative and dyadic discontinuities

The integer part in the outer cosine contributes an integer multiple of
`2*pi` and therefore disappears. Put

```text
u={log_2 x},
a(u)=pi*cos^2(pi*u/2),
A(u)=3/2+(1/2)*cos(a(u)),
V(x)=ln(x)*A(u),                      0<=u<1.
```

**PROVED.** `A(0)=1`, `lim_(u->1-) A(u)=2`, and `1<=A(u)<2`. On `0<u<1`,

```text
A'(u)=(pi^2/4)*sin(a(u))*sin(pi*u)>0,
V'(x)=[A(u)+(ln(x)/ln(2))*A'(u)]/x.                         (E1)
```

Thus V is strictly increasing on every open dyadic band in `x>1`.
The derivative printed in the blueprint is not the derivative of its
displayed V. Already at `x=sqrt(2)`, the printed derivative is `1/sqrt(2)`,
whereas (E1) is `(3/2+pi^2/8)/sqrt(2)`.

At `x=2^m`, for every integer `m>=1`,

```text
lim_(x->2^m-) V(x)=2m*ln(2),
V(2^m)=lim_(x->2^m+) V(x)=m*ln(2).                          (E2)
```

The jump is `-m*ln(2)`. The proposed V is consequently not continuous and
has no derivative there. Its envelope is log-periodic, but V itself is
not scale invariant: `V(2x)=V(x)+ln(2)*A({log_2 x})`.

A spatial derivative is in any event not a Collatz time increment. An
exact hostile inside one band is `19 -> 29`: its halving exponent is one,
and both inputs lie in `(16,32)`, so `V(29)>V(19)`. More generally,

```text
20*2^m-1 -> 30*2^m-1,       m>=0,                          (E3)
```

has exponent one and both endpoints in `(16*2^m,32*2^m)`.
This refutes one-step descent at arbitrarily large integers, independently
of any floating evaluation.

## 2. Arbitrarily long failure of descent for this specific V

**PROVED.** For every integer `L>=1` and **every** integer `m>=2L+4`, set

```text
c=2^(L+1)-1,                 n=2^m+c.
```

Then the first L halving exponents of n are all one, and

```text
V(T^j(n))>V(n)               for every 1<=j<=L.             (E4)
```

**Exact path.** Since `n=-1 mod 2^(L+1)`, induction gives

```text
n_j=T^j(n)=3^j*(2^(m-j)+2^(L+1-j))-1,       0<=j<=L.
```

For `j<L`, `n_j=3 mod4`, so the next exponent is exactly one.
Also `n_j+1=(3/2)^j*(n+1)`, whence `n_j/n>(3/2)^j` for `j>=1`.

**The starting phase is close to the envelope minimum.** The elementary
inequalities `1-cos(t)<=t^2/2` and `|sin(t)|<=|t|` give

```text
A(u)-1 = [1-cos(pi*sin^2(pi*u/2))]/2 <= pi^6*u^4/64.
```

Here `u=log_2(1+c/2^m)<=c/(2^m*ln(2))`, while `ln(n)<(m+1)*ln(2)`.
Using only `pi<4` and `ln(2)>1/2`, the initial excess satisfies

```text
E=V(n)-ln(n)
 <= pi^6/(64*ln(2)^3)*(m+1)*(c/2^m)^4
 < 512*(m+1)*(c/2^m)^4
 < 512*(2L+5)/16^(L+3)
 <= 7/128.                                                   (E5)
```

For the penultimate bound, `(m+1)/16^m` decreases with m; substitute
`m=2L+4` and `c<2^(L+1)`. The last expression decreases for `L>=1`,
and its value at `L=1` is `7/128`. Finally,
`7/128<1/3<ln(3/2)`, the second strict inequality following from
`ln(1+t)>t/(1+t)` for `t>0`.

Since `V(n_j)>=ln(n_j)`, for every `1<=j<=L` we obtain

```text
V(n_j)-V(n) >= ln(n_j/n)-E > ln(3/2)-7/128 > 0.
```

This proves (E4). It rules out a globally fixed lookahead bound for descent
of this V, including after removing any finite set of small starting inputs.
It does not rule out descent after an input-dependent, unbounded number
of steps, and does not supply a divergent Collatz orbit.

## 3. General obstruction for bounded corrections to logarithmic height

**PROVED.** Let `a>0` and let b be any bounded real function on the positive
odd integers. Define `U(n)=a*ln(n)+b(n)`. For every `L>=1` and `N>=1`, there
is an odd `n>N` for which the first L exponents are one and

```text
U(T^j(n))>U(n)               for every 1<=j<=L.             (E6)
```

No regularity or periodicity of b is required. Put
`D=sup b-inf b`, and choose a positive integer r such that
`r*a*ln(3/2)>D`. Start an all-one prefix of length `R=r+L` at
`n_0=2^(R+1)q-1>N`. For all `r<=j<=R`,

```text
U(n_j)-U(n_0)>j*a*ln(3/2)-D>0.
```

Choose the **last** index i attaining the minimum of the finite list
`U(n_0),...,U(n_R)`. Necessarily `i<r`; by the last-minimum choice, every
subsequent list entry is strictly greater than `U(n_i)`. In particular,
the next L entries are, and `n_i>=n_0>N`. This proves (E6).

The source-to-target map here is “a long growth cylinder, followed by its
last energy minimum.” It preserves an exact suffix of the exponent word
and strict height growth; it does not determine where that suffix later
terminates. A bounded log-periodic additive correction is a special case.
The proposed blueprint instead multiplies ln(n) by A, so the next theorem
is needed to cover that distinct class.

## 4. General obstruction for a positive log-periodic envelope

**PROVED.** Let `A:[0,1)->R` have a global minimum `a>0`, attained at a
phase `u_*`. Suppose for some `C,epsilon,beta>0`,

```text
0<=A(u_*+delta)-a<=C*delta^beta      for 0<=delta<epsilon,
```

where epsilon can be reduced so that `u_*+epsilon<1`. No continuity across
the phase seam is required. For `U_A(n)=ln(n)*A({log_2 n})`, for every
`L,N>=1` there exists `n>N` with the first L exponents one and

```text
U_A(T^j(n))>U_A(n)           for all 1<=j<=L.               (E7)
```

To prove it, put `M=2^(L+1)` and choose

```text
n_m=M*ceil((2^(m+u_*)+1)/M)-1.
```

Then `n_m=-1 mod M` and `0<=n_m-2^(m+u_*)<M`. For all sufficiently large
m its fractional log phase is `u_*+delta_m`, with
`delta_m=O(2^(-m))` on the indicated right side. Hence

```text
ln(n_m)*(A({log_2 n_m})-a)=O(m*2^(-beta*m)) -> 0.
```

For `1<=j<=L`, the exact growing prefix gives

```text
U_A(T^j(n_m))-U_A(n_m)
 >= a*ln(T^j(n_m)/n_m)-ln(n_m)*(A({log_2 n_m})-a)
 > a*ln(3/2)-o(1)>0.
```

This includes any positive continuous Lipschitz log-periodic profile that
attains its minimum, and also the discontinuous blueprint profile, whose
right-side bound at zero has `beta=4`. The positive minimum and local
regularity are hypotheses of this theorem, not consequences of boundedness.

## 5. Squarefree density does not supply an orbit average

The inherited source-density theorem gives relative squarefree densities
`9/pi^2,6/pi^2,9/pi^2` in the odd rows `1,3,5 mod6`, respectively. Across
all odd integers the density is `8/pi^2`. The unrestricted integer density
`6/pi^2` therefore cannot be inserted unchanged even as an odd-source
sampling model.

More decisively, the universal orbit-time version is false: `T(1)=1`, so
the squarefree-indicator average on this orbit is exactly one. If every
positive odd orbit converges to 1, every such limiting time average is one.
For the ordinary unaccelerated `1,4,2` cycle the frequency is `2/3`.
If the intended assertion applies only to hypothetical divergent orbits,
these examples do not disprove that restricted assertion, but neither a
source density nor the displayed formula establishes it.

The stronger inherited hostile is **PROVED for each fixed L**: a positive
density of q make all `L+1` nodes

```text
n_j=3^j*2^(L+1-j)*q-1,       0<=j<=L,
```

squarefree while every one of the first L steps grows with exponent one.
The companion script independently verifies the earlier `L=10,q=1`
example. This refutes a fixed-length growth obstruction based only on
squarefreeness of those nodes. The next subsection supplies an actual
intersection argument with the energy obstruction. It still does not
assert that the particular single-point sequence in §2 is squarefree.

### 5a. A short-interval sieve gives the simultaneous obstruction

**PROVED.** Fix `L>=1`, and let `delta_L>0` be the inherited density of q
for which all forms `3^j*2^(L+1-j)*q-1`, `0<=j<=L`, are squarefree.
For real `X->infinity`, and any positive interval length H satisfying
`H=O(X)` and `H/sqrt(X)->infinity`,

```text
#{integer q in [X,X+H]: all L+1 forms squarefree}
                    = (delta_L+o(1))*H.                    (E9)
```

All constants and the convergence in this statement are for **fixed L**.
For a fixed prime cutoff z, CRT counts the q surviving those primes as
`delta_(L,z)*H+O_(L,z)(1)`. An omitted square divisor has
`p<=sqrt(A_L*(X+H))`, where `A_L=2*3^L`. At each p the at most L+1 forms
exclude at most `(L+1)*(H/p^2+1)` interval points. Thus the omitted tail
is at most

```text
(L+1)*H/z + O_L(sqrt(X+H)).
```

Divide by H, first let X tend to infinity with z fixed, and then let z
tend to infinity. The finite products tend to the positive `delta_L`;
the tail vanishes in this order. This proves (E9), without assuming
independence of squarefreeness at successive nodes.

Take `M=2^(L+1)` and the interval of q defined by

```text
2^m <= M*q-1 <= 2^m*(1+1/m).
```

Its left endpoint is `X_m=(2^m+1)/M` and its length is
`H_m=2^m/(M*m)`, so (E9) applies. There are consequently
`(delta_L+o(1))*2^m/(M*m)` such q with every node squarefree. For all
these sources `n=M*q-1`, the initial phase is at most `1/(m*ln(2))`,
and the same bound as in §2 gives

```text
0<=V(n)-ln(n)<=pi^6/(64*ln(2)^3)*(m+1)/m^4 -> 0.
```

Thus, for every fixed L and all sufficiently large m, these are
**simultaneously squarefree at every node, strictly height-growing for
L steps, and satisfy `V(T^j(n))>V(n)` for every `1<=j<=L`.**
This is a short-interval abundance theorem, not a claim that the specific
source `2^m+2^(L+1)-1` is squarefree.

Both general obstructions also survive the same restriction. In §3,
start the last-minimum argument with an all-squarefree prefix of its
required fixed length R; the selected suffix remains all-squarefree.
In §4, let `t_m=2^(m+u_*)` and choose all-squarefree sources in
`[t_m,t_m*(1+m^(-s))]`, where `s*beta>1`. Their q-interval length still
dominates `sqrt(X_m)`, while the initial excess is
`O(m^(1-s*beta))->0`. This supplies the missing intersection, so (E6)
and (E7) hold with all L+1 displayed nodes required to be squarefree.

### 5b. Source averages still do not imply deterministic drift

For odd inputs sampled by initial congruence classes, the first exponent
has frequency `P(k)=2^(-k)`, and a specified finite word has relative
frequency `2^(-sum k_i)`. This gives the formal asymptotic log multiplier

```text
sum_(k>=1) 2^(-k)*ln(3/2^k)=ln(3/4)<0.
```

The actual height increment is instead
`ln(T(n)/n)=ln(3+1/n)-k*ln(2)`. At n=1 it is zero with k=2.
Neither expression is automatically an increment of V. A probability
space, its sampling law, and a theorem connecting it to each fixed
integer's orbit are missing from the blueprint's deterministic conclusion.
Even an almost-everywhere statement on the odd 2-adic integers would not
by itself include all positive integers: that subset is countable and
has Haar measure zero.

## 6. Exact finite-word and rational-cycle hostiles

The supplied `2,4,2,4,4,2,4,2,2` becomes `0,1,0,1,1,0,1,0,0` after mapping
2 to 0 and 4 to 1. It contains both `11` and `00`. A binary mechanical
word `b_j=floor((j+1)*alpha+rho)-floor(j*alpha+rho)`, `0<alpha<1`, has
the number of ones in each block of length r equal to either
`floor(r*alpha)` or `ceil(r*alpha)`: this follows by telescoping. For
r=2 it cannot contain both block sums zero and two. The displayed finite
word is therefore not a factor of a Sturmian word.

The lower Wythoff sequence `W(k)=floor(k*phi)` is an increasing integer
sequence, not itself a binary word. Its first nine successive gaps,
multiplied by two and beginning with k=0, are
`2,4,2,4,4,2,4,2,4`; the ninth entry differs from the attachment.
An exact integer computation uses `W(k)=(k+isqrt(5*k*k))//2`.
No observable for “spaces between structural trajectory shifts” is
defined in the attachment, so repairing this digit alone supplies no map
to a Collatz word, phase, or XOR operation.

For a proposed cycle word `(k_1,...,k_L)`, define `K_0=0`,
`K_j=sum_(i=1)^j k_i`, `K=K_L` and

```text
B=sum_(j=0)^(L-1) 3^(L-1-j)*2^K_j.
```

The affine closure equation is

```text
n_0=B/(2^K-3^L).                                          (E8)
```

The bare expression `2^S-3^m` is not an equation. In (E8), rational
solvability does not imply integer solvability. The literal “only rational
solution” claim has the exact positive hostile

```text
5/7 --k=1--> 11/7 --k=3--> 5/7.
```

Both numerators and denominators are odd; the displayed exponents are the
exact 2-adic valuations of `3n+1`. This is a nontrivial rational cycle of
the natural odd-denominator extension, not an integer Collatz
counterexample. Invoking a theorem on linear forms in logarithms without
its hypotheses, an effective bound, and an integrality argument cannot
remove this missing step. A finite enumeration of the three familiar
negative cycles likewise cannot establish that they exhaust all negative
integer cycles; the earlier census was explicitly finite.

## 7. Strongest useful repair and the remaining target

The corrected derivative, jump formula and all-one cylinders are exact.
The first useful boundary is now sharp: for the proposed V, **no globally
fixed number of steps guarantees even one decrease below its initial
energy**. The same obstruction holds for the two broader classes in §§3–4.

A logically sufficient replacement for the positive Collatz conjecture
is: for every odd `n>1`, there exists a finite, input-dependent `j>=1`
such that `V(T^j(n))<V(n)`. Because `V(n)>=ln(n)` and `V(1)=0`, its
sublevel sets on positive integers are finite; repeated strict decreases
cannot continue inside a fixed sublevel set without reaching 1. Conversely,
reaching 1 proves such a decrease. Thus this replacement is equivalent
to positive Collatz convergence, not a new proof of it. The work needed
is a per-input stopping argument; a mean drift, dyadic picture, or fixed
finite window does not supply one.

## Reproduction and audit boundary

```text
python 04-computation/experiments/collatz_blueprint_20260921_energy.py
python -O 04-computation/experiments/collatz_blueprint_20260921_energy.py --output energy-optimized.json
```

[Source](../../04-computation/experiments/collatz_blueprint_20260921_energy.py)
and [JSON](../../04-computation/experiments/collatz_blueprint_20260921_energy.json)
declare every finite universe. Exact controls cover `1<=L<=64`, three
specified choices of m per L, all associated steps and rational error
bounds; the same-band family for `0<=m<=64`; the word-balance witness;
the rational cycle; and the squarefree length-ten path by factorization
and independent trial-square testing. A separate exact search records
squarefree witnesses in the short dyadic intervals of §5a, with a rational
upper bound below `1/3` for the initial energy excess. These witnesses are
finite controls, not the proof of the asymptotic count (E9).
Floating derivative/energy samples
are separated under `floating_diagnostics_NOT_PROOFS`. Universal proofs
are the arguments above, not extrapolations from those computations.

**Independent read-only audit:** the quadratic/geometry lane checked the
derivative, the exact constant `7/128`, both general energy obstructions,
and the finite-sublevel equivalence in §7; it then separately checked the
short-interval limit order and both squarefree intersections in §5a. All
passed. Ordinary and optimized Python replays produce identical JSON.
These are proof review and executable controls, not a Lean formalization.
