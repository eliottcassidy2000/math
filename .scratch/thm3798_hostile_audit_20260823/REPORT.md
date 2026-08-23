# Independent hostile audit -- THM-3798 common-AP step-three peel

## Verdict and promotion status

**PASS.  No mathematical repair is needed.**  The complete step-three sign
table, endpoint UFD power laws, four scalar-bucket divisor arguments, both
adjacent ODE systems, all constant and support-shrink seams, output swap, and
the step-four hostile table have independent proofs below.

The current canonical file
`01-canon/theorems/THM-3798-cubic-pseudoplane-common-ap-step-three-support-peel.md`
correctly remains **RESERVED / UNPROVED EMPTY STUB** until promotion.  This
audit recommends promotion to

```text
PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
```

The scope must remain exact common-AP `2 x 4` or `4 x 2` support of step
three.  Common step at least four, two-disjoint-pair and three-chain-plus-
isolated grammars, arbitrary `2 x 4`, arbitrary Darboux pairs, and `JC(2)`
remain **OPEN**.

Two harmless Markdown delimiters should be repaired when the proof is moved
to canon:

```text
factor `p)       -> factor `p`
`Delta^2).       -> `Delta^2`.
```

## Verification repair lineage

The first candidate executable was not acceptable as a divisibility
certificate.  It formed `cancel(S/Delta^m)` in the fraction field and then
checked `S=Delta^m cancel(S/Delta^m)`, which is tautological even when the
quotient has a denominator.  It also declared both integration constants
`nu,kap` nonzero.  That verification artifact is **SUPERSEDED**:

```text
old script_sha256 = 2dcfde2f039a9a4740209f973157517792134c7fcb652b5e141c68ee0a781038
old output_sha256 = 2cc0a1c62f401f48685b766dedc73299b08e9435a19c3f22287af7672cdf54fd
old active checks = 17
```

The repaired current candidate supplies explicit denominator-free
quotients, tests `nu=0` and the `kap=0` boundary, contains no executable
Python assertion, and executes all `24` gates under normal and optimized
Python:

```text
current script_sha256 = 6b6a1867ec9311bf108dea2147c39504e66fec78b92c62f0f302043c2db8dad4
current output_sha256 = 9c65319ccf272948767aa4083435299ad072ed660c1b7b62d592dcf464f0c522
current semantic_sha256 = 597dc0e9819318b6d6b1ce6beec3a9821f7d0be6cd42840e43bb7acb47615cbc
current active checks = 24
```

Normal and `python -O` stdout are identical and both LF-normalize exactly to
the repaired frozen transcript.  The repair changes verification, not the
theorem or its proof.

## Inheritance pass and audit board

- **Closest proved mechanism:** THM-3796,
  `01-canon/theorems/THM-3796-cubic-pseudoplane-first-two-by-four-support-peel.md`,
  supplies the collision grammar, exact-support sign/zero law, and the
  step-three rows as an open control.
- **Underlying proved sidecar:** THM-3785,
  `01-canon/theorems/THM-3785-linear-higher-pole-russell-pseudoplane-maximal-observable.md`,
  supplies the exact homogeneous modules and bracket formula.
- **Canonical hostile:** the seven step-four sign survivors.  They show that
  the endpoint-power mechanism is step-specific and cannot be promoted to an
  arbitrary-step claim.
- **Corrected near miss:** the superseded fraction-field cancellation gates
  above.  Exact polynomial division and an actual `Delta`-remainder path
  replace them.
- **Least-used sidecar:** the upper adjacent ODE in the `(-1,-7)` row.  It is
  not needed by the candidate proof, but it forces `g_2=mu p` and strengthens
  that row from a guaranteed `Delta^2` factor to `Delta^3`.

The audit concept board was: contribution addresses; sign/zero endpoints;
UFD exponent lattice; the squarefree arm divisor; adjacent first integrals;
profile-module/support seams; and the step-four boundary.  The anchor was the
step-three proof, the niche was the unused upper ODE, and the wildcard was the
step-four hostile table.

## 1. Exhaustive step-three sign rows

Put the `F` offsets at `0,3` and the `G` offsets at `0,3,6,9`.  Their eight
contributions have relative addresses and multiplicities

```text
0,3,6,9,12                  with multiplicities 1,2,2,2,1.
```

A scalar at a singleton address would be a homogeneous Darboux pair, which
THM-3785 excludes.  Thus its relative address is
`q in {3,6,9}` and the scalar weight equation is

```text
a+b=-q-2.                                                   (1)
```

The bottom singleton commutes.  Here is the complete sign/zero law used at a
singleton:

- if nonzero weights have opposite signs, their commutation equation makes a
  positive product of their nonzero polynomial profiles constant; the
  negative profile contains `Delta` and this is impossible;
- if exactly one weight is zero, commutation differentiates the zero-weight
  profile and makes it a scalar, which is translated away and is not active;
- `(0,0)` commutes identically;
- therefore an active singleton can commute only at equal strict signs or at
  `(0,0)`.

Since the right side of (1) is negative, the bottom cannot be `(0,0)` and
must have `a,b<0`.  Hence

```text
-q-1 <= a <= -1,                 b=-q-2-a,             (2)
```

which is an exhaustive finite interval, not a search cutoff.  The top
singleton has weights `(a+3,b+9)` and sum `10-q`.

- At `q=3`, this sum is `7`, so both top weights must be positive.  Equation
  (2) leaves `a=-2,-1`, hence `(-2,-3),(-1,-4)`.
- At `q=6`, the top sum is `4`; the same strict inequalities leave
  `a=-2,-1`, hence `(-2,-6),(-1,-7)`.
- At `q=9`, the top sum is `1`.  Two positive integers have sum at least
  two, two negative integers have sum at most minus two, and `(0,0)` has sum
  zero.  No row exists.

Thus the exact table is

```text
q=3: (-2,-3), (-1,-4)
q=6: (-2,-6), (-1,-7)
q=9: none.                                               (3)
```

An independent wide census over `-80<=a,b<=80` returns exactly (3), while
the argument above proves exhaustion for all integers.

## 2. Endpoint UFD power law, including units

For nonzero same-sign integers `u,v`, set `U=|u|`, `V=|v|` and
`d=gcd(U,V)`.  The equation

```text
[u,f;v,g]=u f g'-v f'g=0
```

has the same sign-adjusted content

```text
U f g'=V f'g.
```

Therefore in `k(w)`,

```text
(g^U/f^V)'=0,                    g^U=lambda f^V,       (4)
```

where `lambda in k*`; characteristic zero makes the constant field of
`k(w)` exactly `k`.  Write `U=dU_0`, `V=dV_0` with
`gcd(U_0,V_0)=1`.  At every irreducible polynomial `pi`, unique
factorization gives

```text
U ord_pi(g)=V ord_pi(f),
ord_pi(f)=U_0 n_pi,              ord_pi(g)=V_0 n_pi.  (5)
```

Consequently `f=alpha h^U_0` and `g=beta h^V_0` for nonzero field units
`alpha,beta`.  Algebraic closedness is used exactly here: choose an
`U_0`-th root of `alpha`, absorb it into `h`, and obtain

```text
f=t^(U/d),                       g=lambda_0 t^(V/d),  (6)
```

with `lambda_0 in k*`.  No unit or root choice was discarded.

For `Delta=w^3-c^3`,

```text
gcd(Delta,Delta')=1,             disc_w(Delta)=-27c^6 != 0.           (7)
```

This uses `char(k)=0` and `c!=0`.  Equivalently, all three roots are simple
over the algebraically closed field.  Thus `Delta` is radical: if
`Delta|t^m` for a positive integer `m`, every irreducible factor of `Delta`
has positive valuation in `t`, so `Delta|t`.  This validates every passage
from a negative endpoint power to an arm factor in its root `t`.

The independent executable checks the logarithmic-derivative identity for
all endpoint exponent pairs used here and checks `1,940` exact valuation
lattice cases as a hostile control.  The unbounded proof is (5), not that
finite check.

## 3. Scalar address three

Write the two `F` profiles from bottom to top as `p,q` and the four `G`
profiles as `g_0,...,g_3`.

### Row `(-2,-3)`

Bottom commutation and (6) give

```text
p=t^2,                          g_0=lambda t^3,        lambda in k*.  (8)
```

Because `p` has negative weight, `Delta|p`; squarefreeness gives
`Delta|t`.  The scalar bucket factors before any division:

```text
S=[-2,t^2;0,g_1]+[1,q;-3,lambda t^3]
 =t^2[-2g_1'+3lambda q t'+3lambda q't].               (9)
```

Hence `Delta^2|S`.  This remains true if the weight-zero profile `g_1` is
constant and therefore removable; no hidden active-zero assumption is used.

### Row `(-1,-4)`

Bottom commutation gives `g_0=lambda p^4`.  Both `p` and `g_1` have weight
`-1`, so write `p=Delta P`, `g_1=Delta N`.  The derivative-of-`Delta` terms
cancel exactly:

```text
[-1,p;-1,g_1]=Delta^2(P'N-PN').                      (10)
```

The other scalar summand is

```text
[2,q;-4,lambda p^4]=8lambda q p^3p'+4lambda q'p^4,   (11)
```

which is divisible by at least `Delta^3`.  Thus `Delta^2|S` for every
specialization of the profiles and constants.

## 4. Scalar address six and both adjacent systems

### Row `(-2,-6)`

The two endpoints give

```text
g_0=lambda p^3,                   g_3=mu q^3,
lambda,mu in k*.                                         (12)
```

The lower adjacent equation has the exact residual form

```text
g_1=3lambda p^2q+H,              [-2,p;-3,H]=0.       (13)
```

For an arbitrary upper profile `z`, the upper adjacent collision is exactly

```text
[-2,p;3,mu q^3]+[1,q;0,z]
 =q(z-3mu p q^2)'.                                  (14)
```

Since `q` is nonzero, the full polynomial solution is

```text
g_2=3mu p q^2+nu,                 nu in k.            (15)
```

In particular `nu=0` is included.  Direct substitution gives

```text
S=6(lambda-mu)p(pq^2)'+[1,q;-3,H].                   (16)
```

If `H=0`, (16) has the factor `p`, hence `Delta`; if
`lambda=mu` it is zero, not a nonzero scalar.  If `H!=0`, (13) and the UFD
law give

```text
p=t^2,                           H=kappa t^3,          kappa in k*.  (17)
```

Now `Delta|t`.  The first term of (16) has a factor `t^2`, and

```text
[1,q;-3,kappa t^3]=3kappa t^2(qt'+q't),              (18)
```

so `Delta^2|S`.  The formal identity remains polynomial at `kappa=0`, where
it meets the `H=0` branch; the independent remainder path tests both
`kappa=0` and `kappa!=0`, and both `nu=0` and `nu!=0`.

### Row `(-1,-7)`

Bottom commutation gives `g_0=lambda p^7`.  Put

```text
K=g_1-7lambda p^6q.
```

The lower adjacent equation is exactly

```text
-pK'+4p'K=-p^5(K/p^4)'=0.                           (19)
```

The constant field is `k`, so its full polynomial solution is

```text
g_1=p^4(7lambda p^2q+nu),          nu in k.           (20)
```

For the scalar bucket, write the independent negative profiles as
`p=Delta P` and `g_2=Delta N`.  The first summand is the Wronskian (10).
The second differentiates `p^4` at most once and therefore retains at least
`p^3`, hence `Delta^3`.  Thus `Delta^2|S`, including `nu=0`.

There is also a useful redundant upper-seam audit.  Top commutation gives
`g_3=mu q`, since both top weights are `2`.  With `L=g_2-mu p`, the unused
upper collision becomes

```text
2qL'+q'L=0,                       (qL^2)'=0.           (21)
```

The exact weight-`2` module has `q=w^2Q(w^3)`, so active `q` is nonconstant.
If `qL^2` is a nonzero field constant, `q` would be a unit in `k[w]`; if it
is zero, the domain property gives `L=0`.  Therefore `L=0` and
`g_2=mu p`.  The first scalar Wronskian then vanishes identically, sharpening
this row to a `Delta^3` factor.  The candidate did not need this stronger
constraint, but it closes the omitted adjacent-profile seam independently.

In every row, `S` is the profile of a weight-zero bracket bucket and hence
lies in `B_0=k[w^3]`.  Since `Delta=w^3-c^3` is a nonunit, divisibility by
`Delta` makes `S` either zero or nonconstant.  It cannot be a nonzero scalar.

## 5. Support shrink, zero weights, and profile modules

THM-3785 gives the exact module

```text
B_u=x^u w^rho(u) Delta^m(u) k[w^3],
rho(u)=u mod 3 in {0,1,2},
m(u)=max(0,ceil(-u/3)).                              (22)
```

The audit recomputed (22) on every weight appearing in the four rows.
Consequences used above are exact:

- every negative profile contains `Delta` and is nonconstant;
- the upper `F` profile has weight `1` or `2`, hence respectively contains
  `w` or `w^2` and is nonconstant;
- the only zero-weight profiles are `g_1` in `(-2,-3)` and `g_2` in
  `(-2,-6)`, both inside collision buckets rather than singleton endpoints;
- a scalar part of either zero-weight profile can be translated away without
  changing the bracket or another support component.

All endpoint constants `lambda,mu` are nonzero because the exact endpoint
profiles are nonzero.  The integration constant `nu` is arbitrary.  In the
`H!=0` branch, `kappa` is nonzero; its zero specialization is the already
audited `H=0` seam.

No constant specialization opens a support loophole.  In `(-2,-6)`,
`g_2=3mu p q^2+nu` cannot vanish or become scalar because `p q^2` is a
nonzero nonconstant multiple of `Delta`.  In `(-1,-7)`, the factor
`7lambda p^2q+nu` cannot vanish for the same reason, and the upper seam gives
the nonzero profile `g_2=mu p`.  A cancellation that makes another displayed
profile zero only leaves the exact-support hypothesis; the divisor proof was
performed in the larger class and remains valid.  Independently, the smaller
`2 x 3` boundary is already closed by THM-3787,
`01-canon/theorems/THM-3787-cubic-pseudoplane-complete-low-support-darboux-nonentry.md`.

## 6. Exact output swap

For arbitrary weights and profiles,

```text
[v,g;u,f]=-[u,f;v,g].                                (23)
```

Swapping the outputs preserves the exact two supports and their common step,
and changes a nonzero scalar bracket only by sign.  The weight-zero
translation convention is symmetric.  Thus every `4 x 2` cell is exactly
the swapped `2 x 4` cell, with no new module, zero, or support-shrink seam.

## 7. Step-four hostile boundary

At step four the relative addresses are `0,4,8,12,16`.  Repeating the exact
negative-sum interval and singleton sign/zero test gives

```text
q=4:  (-3,-3), (-2,-4), (-1,-5)
q=8:  (-3,-7), (-2,-8), (-1,-9)
q=12: (-3,-11).                                      (24)
```

These seven rows survive only the necessary sign/zero gate.  The audit does
not construct their coefficient systems and does not exclude them.  They are
positive **OPEN** controls and the sharp advertised failure boundary of the
present proof.  The THM-3796 disjoint-pair control
`r=4,V={0,1,4,5}` also remains open.

## Independent exact artifacts

```text
.scratch/thm3798_hostile_audit_20260823/
  audit_thm3798_common_ap_step_three.py
  audit_thm3798_common_ap_step_three.out

audit_script_sha256 = 9c31a03c6ebbe430c3a2b24c6f89a6660e6a7901de0603cf8ad910857f0bc74d
audit_output_sha256 = 4ef837e89e778f4ed452cde26c6c2c3c75956e362c5f0014733c28a81a6469d8
audit active checks = 6105
audit AST bare asserts = 0
```

The independent program uses two redundant paths:

1. formal first-jet polynomial algebra with exact polynomial division and
   reconstruction of every claimed `Delta` factor;
2. actual polynomial representatives in the exact weight modules, with
   Euclidean remainders modulo `Delta` or `Delta^2`.

Reproduction:

```text
python .scratch/thm3798_hostile_audit_20260823/audit_thm3798_common_ap_step_three.py
python -O .scratch/thm3798_hostile_audit_20260823/audit_thm3798_common_ap_step_three.py
```

Both executions LF-normalize exactly to the frozen independent output.  The
research cards used were “Search the statement before the method,” “Verify
the consequence, not the model's own assumptions,” “Use redundant paths as
detectors,” and “Correct the object before sharpening the technique.”  No new
meta-pattern promotion is warranted.
