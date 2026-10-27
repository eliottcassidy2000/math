# Session 74: all-period golden fibres and the source-preserving lift law

**Status: PROVED + INDEPENDENTLY AUDITED (repository all-period gap closure,
not a historical-novelty claim).** Root explicitly authorized this scoped
promotion after the final independent audit. All claims are elementary
and self-contained below. This closes a bounded repository gap; it is not a
claim of historical novelty (see PRIMARY_LITERATURE.md). The integer-cycle
census and the finite-valuation obstruction are inherited, not newly proved.
No Collatz convergence or integer-cycle classification is claimed.

Source pin: `4959b4933cba7e6561b9f19b1d68278ecf6cb2f5`.

## 1. The recovered eleven is three different objects

* THM-4507, `01-canon/theorems/THM-4507-finite-valuation-collatz-potential-obstruction.md`, already excludes arbitrary nonlinear functions of finitely many polynomial valuations at any finite prime set, including `{2,3,11}`, as a globally nonincreasing correction to `c log n`, even with any remainder bounded on each valuation fibre. Its positive witnesses change with block length. Its rational-shadow construction is not a divergent positive orbit.
* THM-4591, `01-canon/theorems/THM-4591-integer-cycles-of-the-terras-map-to-period-301993-are-the-five-pythagorean-cycles.md`, gives the periods `{1,2,3,11}` only as a FINITE-EXACT census through Terras period 301,993. Here `T(x)=x/2` for even `x` and `(3x+1)/2` for odd `x`. The cycle through `-17` has 11 steps and 7 odd steps; its standard-map period is **18**, not 11.
* THM-4592, `01-canon/theorems/THM-4592-the-eleven-squares-fibonacci-torus-is-the-period-five-collatz-cycle-set.md`, uses the **standard** map `C(x)=x/2` or `3x+1`. Its eleven is `L_5=11` and its clock is five. Its existing all-odd-length bijection language is explicitly restricted to FINITE-EXACT `n<=22` by MISTAKE-584. The new proof below targets precisely that bounded gap, not the integer-cycle census.

Concept board: rooted cyclic words / quadratic-integer carries / ideal quotients /
period extension / transported source mass / positive shadow precision.

## 2. All-period golden fibre theorem

Let `phi^2=phi+1`, `O=Z[phi]`, `R_n=O/(phi^n-1)`, and let `W_n` be the set
of **rooted** binary words of length `n>=1` with no cyclically adjacent `11`,
including across the last/first seam. Put

```
F(w) = sum_{j=0}^{n-1} w_j phi^(n-1-j),
kappa_n(w) = F(w) mod (phi^n-1).
```

**Theorem.** For every odd `n>=1`, `kappa_n:W_n -> R_n` is bijective. For
every even `n>=2`, it is surjective and its only non-singleton fibre is

```
{0^n, (10)^(n/2), (01)^(n/2)} -> 0.
```

This is a theorem on rooted words, not necklaces; rotation is retained as a
dynamical action. For all `n`, left shift corresponds to multiplication by
`phi` in `R_n`.

### 2A. Every collision gives a bounded integral periodic carry

Suppose `kappa_n(u)=kappa_n(v)`. Then

```
z_0 = (F(u)-F(v))/(phi^n-1) in O,
z_{i+1} = phi z_i - (u_i-v_i),      z_n=z_0.
```

All carries are integral because `phi` and the digits are integral. Under the
real embedding, each `F(w)/(phi^n-1)` is the infinite periodic expansion
`sum_{j>=1} w_{(j-1) mod n} phi^(-j)`, lying in `[0,1]`: in each disjoint
pair of places at most one digit is one, so the sum is at most
`sum_{k>=0} phi^(-(2k+1))=1`. The same holds at every rotated word. Hence
`|z_i|<=1`.

Under conjugation `psi=1-phi=-1/phi`, periodicity and the recurrence give
`M<=|psi|M+1` for `M=max_i |z_i'|`; therefore `|z_i'|<=phi^2`.
For `z=a+b phi`,

```
|b| sqrt(5)=|z-z'| <= 1+phi^2 < 2 sqrt(5).
```

Thus `b in {-1,0,1}`. Solving the two elementary interval inequalities gives
exactly the seven possible carries

```
0, +/-1, +/-(phi-1), +/-(phi-2).
```

The weaker generic bound `|z|<=phi` instead gives exactly nine carries, adding
`+/-phi`; these are retained in the independent finite graph certificate.
The stronger word-specific bound removes them without enumeration.

### 2B. Seven-state hand proof of all possible periodic collisions

Write `B=phi-1` and `A=phi-2`. For `d=u_i-v_i`, the complete table of
transitions remaining within the seven-carry set is

| current carry | `d=-1` | `d=0` | `d=1` |
|---|---|---|---|
| `-1` | `-B` | — | — |
| `-B` | `0` | `-1` | — |
| `A` | `-A` | `-B` | — |
| `0` | `1` | `0` | `-1` |
| `-A` | — | `B` | `A` |
| `B` | — | `1` | `0` |
| `1` | — | — | `B` |

* If a periodic carry orbit contains `1`, the next digit pair is `(1,0)`
  and the next carry is `B`. Returning immediately from `B` to zero would
  require another `(1,0)`, violating the no-`11` rule for `u`. Therefore the
  next pair is `(0,0)` and the carry returns to `1`. The two-cycle is forced
  forever: one word is alternating and the other zero. The negative case
  interchanges the two words.
* Each of `A` and `-A` has only the other as a predecessor. Any periodic
  orbit containing either is therefore their two-cycle, with digit pairs
  `(0,1),(1,0)` alternately. These are the two alternating words.
* If a periodic orbit contains neither `+/-1` nor `+/-A`, it cannot contain
  `+/-B` either, since their next carries are zero or `+/-1`; periodicity
  would require re-entry from `+/-1` or `+/-A`. The carry is identically
  zero, so the two words coincide.

Thus every nontrivial collision is one of the three pairings among the
three stated even words. All three do have zero image: for even `n`, the
alternating sums are `(phi^n-1)/(phi^2-1)` and `phi` times that; since
`phi^2-1=phi` is a unit, they lie in `(phi^n-1)`.

### 2C. Surjectivity is an index argument, not a sampled count

The number of rooted legal cyclic words is
`trace([[1,1],[1,0]]^n)=L_n`, including `n=1`. Multiplication by
`g=a+b phi` on the integral basis `(1,phi)` has matrix

```
[[a,b],[b,a+b]],  determinant a^2+ab-b^2.
```

For `g=phi^n-1`, the absolute determinant is

```
|R_n| = |(phi^n-1)(psi^n-1)| = L_n-1-(-1)^n.
```

It is nonzero for every `n>=1`. The proven fibre list therefore has exactly
`|R_n|` distinct images: `L_n` for odd `n`, `L_n-2` for even `n`. This proves
surjectivity in every length. Left-shift equivariance follows from
`F(shift w)=phi F(w)-w_0(phi^n-1)`.

## 3. Rational standard-Collatz sources, proved separately

For `w in W_n`, compose the affine branches dictated by its bits. If `k` is
the number of ones and `h=n-k`, the composition has the form

```
x -> (3^k x + B_w)/2^h,
x_w = B_w/(2^h-3^k).
```

The denominator is nonzero. If `k=0`, the fixed point is zero. If `k>0`,
`h>=1`, and `2^h-3^k` is odd. The affine numerator `B_w` is a sum, one term
per one-bit, with a factor `2^(number of preceding zero-bits)` and an odd
power of three. If `w_0=1`, precisely the first one-term is odd, since
cyclic no-`11` forces a zero before every later one; if `w_0=0`, every term
is even. Thus `x_w` has parity `w_0` in `Q intersect Z_2`.

The image under the first branch is the fixed point for the rotated word
(by cyclic composition), so the same argument proves **every** prescribed
parity, and `C^n(x_w)=x_w`. Conversely a rational 2-adic integer fixed by
`C^n` supplies its unique parity word in `W_n` and satisfies this same affine
fixed-point equation. Hence `W_n` and `Fix(C^n;Q intersect Z_2)` are in
bijection, independently of the golden quotient. The three exceptional
even words are precisely the sources `0,-1,-2`.

Consequently the all-period fibre theorem also applies to these rational
periodic points. It does not classify which points are ordinary integers or
positive integers; those arithmetic predicates still require the actual
word and `B_w/(2^h-3^k)`.

## 4. The source-preserving changing-modulus law

Let `m,r>=1`, `n=rm`, and

```
Q_{m,r}=1+phi^m+...+phi^((r-1)m).
J_{m,r}: R_m -> R_rm,     a -> Q_{m,r} a.
rho_{m,r}: R_rm -> R_m,   a -> a mod (phi^m-1).
```

**Theorem.** `J` is an injective additive homomorphism, commutes with
multiplication by `phi`, and represents period extension of the **same**
rational source, because `F(w repeated r times)=Q_{m,r} F(w)`. Furthermore

```
rho_{m,r} o J_{m,r} = multiplication by r on R_m.
```

Proof: `(phi^(rm)-1)=Q_{m,r}(phi^m-1)`. If `Qa` lies in the product ideal,
cancellation in the integral domain `O` gives `a in (phi^m-1)`; this is
injectivity, not merely an equal-cardinality observation. Modulo
`phi^m-1`, every summand of `Q` is one. Generally `J` is **not a unital
ring map**. The ordinary quotient map `rho` is a ring map and is not an
inverse to `J`.

**Integer-clock sidecar.** Commutation with the one-step action is not a
claim that a repeated traversal has no cost. The underlying rational source
and its minimal period are unchanged, while the chosen word presentation
has length `m -> rm`. Literally traversing that presentation repeats the
original length-`m` presentation `r` times: its odd count `k` becomes `rk`, its even count `h` becomes
`rh`, and its standard-map step count `m` becomes `rm`. Additive traversal
costs and accumulated drift must therefore be multiplied by `r` as well.

### Multiplication twist and the square-zero eleven-cell ideal

The independent auditor supplied a further exact guard:

```
Q_{m,r}^2 = r Q_{m,r} in R_rm,
J(a)J(b)=r J(ab).
```

Indeed, writing `z=phi^m`, multiplication of the two sums
`1+z+...+z^(r-1)` has exactly `r` terms in every exponent class modulo `r`,
and `z^r=1` in `R_rm`. The image of `J` is the ideal generated by `Q`.
For `(m,r)=(5,11)`, all products in its eleven-element image vanish,
because `11 R_5=0`. Thus the lifted source cells form a **square-zero
ideal**, not an embedded copy of the field `F_11` as rings. The additive
embedding and one-step dynamics remain faithful. The whole ring `R_55`
is not being identified with the dual numbers or with this eleven-element
ideal. This addition is credited to the independent session auditor.

### Exact return-loss index for every modulus

Write `phi^m-1=a_m+b_m phi`, so `a_m=F_(m-1)-1`,
`b_m=F_m`, and put `d_m=gcd(a_m,b_m)`, `N_m=|R_m|`.
The multiplication matrix above has entry-gcd `d_m`, and determinant of
absolute value `N_m`; its Smith factors are therefore

```
R_m (additively) = Z/d_m Z + Z/(N_m/d_m) Z.
```

Consequently the exact kernel size of `rho o J = [r]` is

```
K(m,r)=gcd(r,d_m) * gcd(r,N_m/d_m).
```

Every nonempty return fibre has precisely `K(m,r)` labels, and the image
has `N_m/K(m,r)` labels. Return is lossless iff `gcd(r,N_m)=1`.
This is a label-fibre cardinality, not Shannon entropy without a measure.
Source weights on each merged fibre must be summed, not equalized.
For `m=5`, the Smith factors are `(1,11)`; for `m=10`, they are `(11,11)`.
Thus repetition eleven destroys all base label information in both cases,
while the actual ideal lift remains injective.

### Elevenfold hostile and the old factor-two twist

At `m=5`, `R_5=F_11`, `phi=4`. The eleven rational sources are zero, the
five points of the `1/13` cycle, and the five points of the `-5` cycle.
For `r=2`, projecting their lifts back multiplies labels by two, exactly
THM-4592's corrected colour swap. More generally nonzero labels preserve or
swap the two Paley colours according as `r mod 11` is a square or nonsquare.

At `r=11`, **all eleven projections are zero**, although their lifts in
`R_55` remain eleven distinct cells. The ambient ring has
`L_55=312119004989` cells. Thus even this exceptionally small exact source
set becomes completely invisible to the naive return projection. This is
an explicit all-period mechanism behind a changing-modulus loss, not a
numerological identification of 11s.

For arbitrary source weights `mu(w)`, use `J_* mu` on the lifted eleven
cells. Its total mass and source distinctions are preserved. Applying
`rho_*` at `r=11` puts the entire mass at zero. Replacing the transported
law by uniform mass on all `312119004989` cells is a different measure and
is not authorized by this correspondence. The exact certificate assigns
weights `1,...,11`, retains mass 66, and detects its collapse to zero.
When `m` is even, the original three-way zero fibre must itself be retained
as a source sidecar; the ring coordinate never recovered those three sources.

**Connection contract:** rooted periodic source -> repeated-word ideal lift;
map `J` above; preserves source, one-step equivariance, and source mass (not free traversal cost); loses nothing
at the additive-label level, but even zero already has a three-source fibre;
projection `rho` destroys `r`-torsion information; retain `Q`, the lifted
cell, and source weights. Cheapest hostile: the elevenfold lift from 5 to 55.

## 5. A denominator-free eleven-step shadow (control, not new obstruction)

The exact Terras cycle is

```
-17,-25,-37,-55,-82,-41,-61,-91,-136,-68,-34 -> -17.
```

For every positive integer `M` and `k>=1` with `M*2^(11k)>17`,

```
n=M*2^(11k)-17,
T^(11k)(n)=M*3^(7k)-17=y.
```

All intermediate parity bits are those of the repeated cycle: initial
agreement has `11k` binary digits and one bit pays each shortcut step.
The endpoint expansion is the inherited multiplier `(2187/2048)^k`.
For `M=2^5*3^3*11^2`, exact endpoint valuations of
`X,X+1,X-1,3X+5,2X+1` coincide at `2,3,11`. The tested lengths are
`k=1,2,3,11,22` (up to 242 actual steps), not an infinite positive orbit.

The hostile `P(X)=X+17` breaks this particular fixed shadow: its 2-adic
valuation drops by `11k` and its 3-adic valuation rises by `7k`. This is
why a single integer cycle does not replace THM-4507's polynomial-root
avoidance across infinitely many rational shadows. No new universal
potential obstruction is claimed here.

## 6. Reproduction and audit surface

`verify_golden.py` uses only exact integers/Fractions and explicit exceptions,
not assertions. Run from any working directory:

```
python3 -B PATH/verify_golden.py > receipt.json
python3 -B -O PATH/verify_golden.py > receipt.optimized.json
cmp receipt.json receipt.optimized.json
```

The two retained receipts agree byte-for-byte, SHA256
`8be63c669fae6bf6c5dd8f5f13b52b80b109cb28bfb1b669cf16bc3470260958`.
They record 95,755 checks, the complete 36-vertex/59-edge coarse graph,
its four recurrent SCCs, all legal words through length 22, 5,775 direct
rational-source checks through length 16, 22 source-preserving lift tests
with all eleven sources, and the integer shadows. Invalid internal-`11`
and invalid cyclic-seam words are deliberate collision-producing controls.
The finite automaton is a certificate of an all-period implication because
the analytic carry bounds prove that **every** collision enters its finite
state universe; the word census itself remains explicitly bounded.

The separate exact additive-quotient control `verify_return_loss.py` tests
all `1<=m<=18`, `1<=r<=22` (396 Smith loss cases), checks the multiplication
twist on an integral basis, and checks all 121 products in the eleven-cell
square-zero ideal. Run normally and under `-O`; both retained
`return_loss*.json` receipts agree. Its computation is a control on the
explicit Smith proof, not a substitute for that all-integer proof.

## Independent audit acceptance

The separate session audit packet `session74-golden-audit/AUDIT.md` accepts
the all-period collision proof, rational source legality, quotient index,
rooted/primitive-period conventions, ideal lift, multiplication twist,
integer-clock sidecar, and Smith return-loss formula. Its independent code
imports no producer code; its final normal/optimized receipts agree on
105,461 checks, including 264 direct-quotient Smith tests. Final independent
receipt SHA256:
`0b191472543d422b8ad7a65e360b9c7365af468d7a81af798e6df346b0cd9a71`.
The auditor's sharp graph has 28 vertices/41 edges; the producer's coarse
graph has 36 vertices/59 edges. These are intentionally different complete
carry regions, not conflicting counts. Root authorized the scoped
PROVED + INDEPENDENTLY AUDITED promotion after reading the final report.
