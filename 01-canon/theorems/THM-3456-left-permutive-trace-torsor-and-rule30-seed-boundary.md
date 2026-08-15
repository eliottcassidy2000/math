---
id: THM-3456
title: "Left-permutive trace torsors and the Rule 30 seed boundary"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  A binary
  left-permutive radius-one cellular automaton has an exact trace/right-half
  coordinate system and uniform finite trace fibres.  For Rule 30 this gives
  an explicitly language-dependent SOP2 witness, a ternary mask-intersection
  compiler, and exact factorial and polynomial-lift controls, but none retains
  the distinguished single-seed boundary or proves an open Rule 30, LRC,
  factorial, or Jacobian statement.  The corollary SOP2 => SOP3 is CITED from
  Chernikov's preprint, not reproved here.
source: root-rule30-260813291-20260815
audit: >
  independent all-rule trace, global-mask, cyclic-gluing, factorial compiler,
  finite-width Jacobian, scope, source-typing, hostile-mutation, dependency,
  and ordinary/optimized/stored byte-replay audit; documentation gates clean
depends_on:
  - THM-2810-factorial-hankel-faithfulness-and-bounded-radial-carrier-no-go
  - THM-3395-small-sheet-typed-cover-star-cochain
related:
  - THM-2050-period14-top-germs-do-not-determine-global-loneliness
  - THM-3359-modular-c-finite-supports-harmonic-density-and-periodic-scar
  - THM-3367-berggren-spinor-pencil-hessian-gauge-and-affine-line-keller-closure
  - THM-3457-rule30-right-edge-2-adic-odometer-and-moving-observer-boundary
script: 04-computation/rule30_left_permutive_trace_torsor_thm3456.py
output: 05-knowledge/results/rule30_left_permutive_trace_torsor_thm3456.out
script_sha256: 9b965c0e02c28843fadf14a3522620ea37866a12f79ce15a0d106502b81f7bd9
output_sha256: 4a080d603ffda54673b0059a9c30621545ebc986de42a1314ac92ea7d1e18f19
hash_basis: raw bytes
external: "Artem Chernikov, SOP2=SOP3, arXiv:2608.13291v1 (2026-08-13): CITED PREPRINT theorem SOP2 => SOP3"
---

# THM-3456 -- left-permutive trace torsors and the Rule 30 seed boundary

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

The elementary cellular-automaton statements below are proved for all
horizons.  The exact companion audits bounded universes independently.  The
model-theoretic corollary uses Chernikov's **CITED PREPRINT** theorem
`SOP2 => SOP3`; it is deliberately scoped to an explicit two-sorted expansion.

## 1. The trace/right-half coordinate theorem

Work over `F_2`.  Let `F:{0,1}^Z->{0,1}^Z` be a binary radius-one cellular
automaton with local rule

```text
f(l,c,r)=l xor h(c,r),                                  (1)
```

where `h:{0,1}^2->{0,1}` is arbitrary.  Thus `F` is left-permutive.  For a
configuration `x`, write

```text
tau_F(x)=((F^t x)_0)_(t>=0),
rho_+(x)=(x_1,x_2,...).                                 (2)
```

Then

```text
Theta_F:x |-> (tau_F(x),rho_+(x))                       (3)
```

is a bijection

```text
{0,1}^Z  ->  {0,1}^{N_0} x {0,1}^{N_+}.               (4)
```

More precisely, for every `t>=0` there is a Boolean function `H_t` such that

```text
(F^t x)_0=x_(-t) xor H_t(x_(-t+1),...,x_t).             (5)
```

Consequently a desired trace and an arbitrary positive half-line reconstruct
`x_0,x_(-1),x_(-2),...` successively and uniquely.

### Proof

At `t=0`, (5) is immediate.  In the radius-one light cone for `(F^t x)_0`,
the only path from the extreme-left input `x_(-t)` uses the left input at
every gate.  Each such gate has Boolean derivative one in that input by (1),
so `x_(-t)` occurs with coefficient one and no other occurrence can cancel
it.  Equivalently, induction on `t` gives (5) directly.  Once the trace and
positive half-line are fixed, its equation at time zero fixes `x_0`; the
equation at time `t` then fixes `x_(-t)` after all arguments of `H_t` have
already been fixed.  This constructs the inverse of (3) coordinatewise.  QED.

## 2. Uniform finite fibres and the ensemble boundary

Fix a horizon `T`.  Restricting (3) gives a finite bijection

```text
(x_(-T),...,x_T)
  |-> ((F^t x)_0)_(0<=t<=T), (x_1,...,x_T).             (6)
```

Therefore every center-trace word of length `T+1` has exactly

```text
2^T                                                        (7)
```

initial cones of length `2T+1` above it.  Under independent fair initial
bits, the vertical trace is itself an independent fair-bit process and is
independent of the positive initial half-line.

There is a spatial version.  A one-step open block map from `m+2` input bits
to `m` output bits has exactly four preimages of every output: choose the two
right boundary bits and solve the remaining inputs from right to left using
left permutivity.  Iterating through depth `T` gives exactly

```text
4^T                                                        (8)
```

preimages from an input block of length `m+2T` to every output block of
length `m`.

Equations (7)--(8) concern a free input ensemble.  They say nothing about a
fixed initial configuration.  For Rule 30, the single-seed problem fixes the
positive half-line and the desired trace simultaneously; it removes exactly
the `T` free bits counted in (7).  This is the first load-bearing loss.

The loss is genuine for the whole left-permutive class.  Rule 90 has
single-seed center trace `1,0,0,...`, whereas Rules 60 and 150 have constant
single-seed center trace `1,1,1,...`, even though all have the same bijection
(3).  Hence neither trace surjectivity nor uniform ensemble balance implies
single-seed aperiodicity or density one half.

For completeness, these infinite controls follow from Laurent-polynomial
evolution over `F_2`.  Rule 60 has kernel `1+X`, whose constant coefficient in
every power is one.  Rule 90 has kernel `X^(-1)+X`; its nonzero-time constant
coefficient is either absent or the even central binomial coefficient.
For Rule 150, write `t` in binary and factor
`(X^(-1)+1+X)^t` into its Frobenius powers.  The largest selected signed power
of two exceeds the sum of all smaller ones, so only the all-zero selection
contributes to the constant term, with coefficient one.

## 3. The typed SOP corollary and its exact scope

Add two sorts:

```text
C = configurations,        W = finite binary words,       (9)
```

and a primitive relation

```text
E(x,w)  iff  w is an initial segment of tau_F(x).          (10)
```

The formula `E(x;w)` has `SOP2`.  Indeed, index parameters by the binary tree
`2^{<omega}` using the corresponding words.  Along every infinite branch,
the prefix conditions are consistent by the surjectivity of (3); two
incomparable words define disjoint trace cylinders and hence give an
inconsistent pair.

Chernikov's theorem `SOP2 => SOP3` therefore implies that the complete theory
of this explicit incidence expansion has `SOP3`.  The cited theorem constructs
an `SOP3` formula after retaining typed witness--parameter pairs; it does not
identify causal cells or CA states with model-theoretic vertices.

This corollary is intentionally language-dependent.  Rule 90 has the same
trace-prefix incidence tree, although its distinguished-seed trace is
eventually zero.

A separate hostile applies to any future **causal-path-indexed** witness; it
is not a claim about the binary prefix parameters in (10).  The depth-two
causal addresses `LR`, `CC`, and `RL` all terminate at the same physical cell
`(2,0)`.  Thus a family required to distinguish those paths cannot factor
through the endpoint quotient without a full path/address sidecar.

## 4. Rule 30 as a ternary intersection compiler

Rule 30 is

```text
f(l,c,r)=l xor c xor r xor c r
        =l xor (c or r).                                (11)
```

For three subsets `A,B,C` of a finite ground set, apply (11) pointwise and put

```text
W_A=A xor (B union C),       P=A xor B xor C.           (12)
```

Then

```text
W_A xor P=B intersect C.                                (13)
```

Cyclically rotating the three inputs recovers all three pair intersections.
Thus three Rule 30 gates plus XOR certify that `A,B,C` are pairwise disjoint.

For `r>=1` owner masks, give the owners distinct ternary codewords of length

```text
d=ceil(log_3 r).                                        (14)
```

At digit `j`, union masks whose owner code has digit zero, one, or two, and
apply the three rotated gates.  Some owner pair overlaps iff some digit
reports a group-pair overlap: choose a differing digit for an overlapping
pair in the forward direction, and expand any reported group intersection in
the reverse direction.  A report need not identify which owner pair caused
it.  Hence `3d` Rule 30 gates give an exact global pairwise-disjointness
certificate.  Adding equality of the total union with the ground set
certifies an exact partition.  The code depth (14) is sharp within this
nonadaptive three-group coding scheme, since `3^d` distinct owner labels are
necessary.

Equivalently, the union of all reported defects is exactly

```text
union_(i<j) (mask_i intersect mask_j).                  (14a)
```

This is a symmetric defect-colour compiler, not a tournament: intersections
have no intrinsic orientation.  It is also only an LRC prefilter.  In the
`q=6` typed-sheet controls of THM-3395, both the positive row with speeds
`(2,8,10)` and the hostile row `(2,8,14)` have the identical partition masks

```text
{0,3}, {1,4}, {2,5},                                    (15)
```

but only the first has the required compatible star cochain.  The mask
compiler destroys affine gaps and common source time, exactly the data needed
for the cochain equations.  It therefore proves no LRC row safe.

Periodic gluing supplies a second sharp warning.  On an `N`-cycle, the
all-one output word has a Rule 30 predecessor iff `3|N`, and in that case has
exactly the three rotations of

```text
(100)^(N/3).                                             (16)
```

To prove this, choose a one in a predecessor.  Requiring the next two output
bits to equal one forces two zeros followed by another one, and the pattern
continues around the cycle.  Thus the open fibres (8) do not survive an
untyped cyclic quotient.  In particular, one-step all-safe reachability on a
13-cycle is empty.  This refutes only that naive one-step cyclic encoding; it
says nothing about reachability or safety in the actual LRC phase problem.

## 5. Exact factorial near miss and sequence compiler

Let the factorial functional be

```text
L(l^a c^b r^d)=a! b! d!.                               (17)
```

The unique real multilinear polynomial agreeing with Rule 30 on the Boolean
cube is

```text
g=l+c+r-c r-2l c-2l r+2l c r.                          (18)
```

Directly,

```text
L(g)=0,                 L(g^2)=6.                      (19)
```

So `g` is an exact FC(3) one-moment near miss, detected at the next power; it
is not a factorial-conjecture counterexample.

There is an expansion-free moment compiler.  Write `D_j` for the derangement
number `!j` and put

```text
M_k=sum_(j=0)^k (-1)^j binom(k,j) D_j^2.               (20)
```

Then for every `m>=0`,

```text
L(g^m)=sum_(q=0)^m binom(m,q) q!
          sum_(s=0)^q binom(q,s)(-2)^s M_(m-q+s).      (21)
```

Proof: for independent unit exponentials `C,R`, put
`o=C+R-CR=1-(1-C)(1-R)`.  Since
`E[(1-C)^j]=(-1)^j D_j`, equation (20) is `E[o^k]`.  Now
`g=o+l(1-2o)`; expanding first in `l` and then in `o` gives (21).

For comparison, the integer ANF lift

```text
p=l+c+r+c r                                               (22)
```

has only positive factorial moments.  If

```text
u_n=L(p^n)/n!,                                           (23)
```

then

```text
sum_(n>=0) u_n z^n
  =(1-z)^(-3) Phi(z/(1-z)^2),
Phi(w)=sum_(d>=0) d! w^d,                               (24)
```

as a formal power series.  Equivalently,

```text
u_n=sum_(d=0)^n d! binom(n+d+2,n-d),                    (25)
```

and, for `n>=4`,

```text
u_n=(n+2)u_(n-1)-(2n-1)u_(n-2)
       +(n-3)u_(n-3)+u_(n-4),                          (26)
u_0,u_1,u_2,u_3=1,4,13,45.                             (27)
```

To derive (24)--(25), expand `p^n` and let `d` be the number of selected
bilinear factors `cr`.  After the factorial readout, summing the remaining
numbers of `l,c,r` factors contributes

```text
d! z^d (1-z)^(-1)(1-z)^(-(d+1))(1-z)^(-(d+1)).
```

Summing over `d` gives (24), and coefficient extraction gives (25).
The formal series `Phi` satisfies

```text
w^2 Phi'(w)+(w-1)Phi(w)+1=0.
```

After the substitution in (24), this becomes

```text
z^2(1-z)^2 U'(z)+(-1+3z-3z^2+z^4)U(z)+(1+z)=0.
```

Its coefficient of `z^n`, for `n>=4`, is exactly (26).

This is an efficient exact sequence compiler, not a zero-moment mechanism.

The Boolean quotient and factorial carrier cannot be silently identified.
Boolean truth imposes `x_i^2=x_i`, while the characteristic-two factorial
carrier kills every exponent at least two and so imposes `x_i^2=0` in its
faithful quotient.  A common quotient collapses `x_i`; characteristic two is
also outside the good-prime range used on the live FC frontier.

## 6. Polynomial-lift and Jacobian boundaries

For each finite width `w>=1`, the open right-edge ANF lift with
`0<=k<w`,

```text
K_k=x_k+x_(k-1)+x_(k-2)+x_(k-1)x_(k-2),                (28)
```

with exterior lower indices zero, is unit triangular over characteristic
zero.  Its Jacobian determinant is

```text
product_(k=0)^(w-1) partial K_k/partial x_k = 1,        (28a)
```

and it has a recursive polynomial inverse, so every finite truncation is a
tame Keller automorphism.  But it is not a faithful real Boolean lift:
at `(x_k,x_(k-1),x_(k-2))=(0,1,1)` its value is `3`, rather than the Rule 30
bit `1`.

Conversely, the unique multilinear real truth lift uses

```text
o_k=x_(k-1)+x_(k-2)-x_(k-1)x_(k-2),
B_k=x_k+o_k-2x_k o_k.                                   (29)
```

For the same finite width, it is faithful on `{0,1}` but its triangular
Jacobian determinant is

```text
product_(k=0)^(w-1) (1-2o_k),                           (30)
```

which is nonconstant from width two onward.  Thus the two canonical lifts
preserve complementary predicates: ANF preserves triangular Keller form,
whereas the multilinear lift preserves the real truth table.  This supplies
only a tame control and a no-go for the naive lift; it gives no reduction of
the planar Jacobian conjecture.

Over `F_2` the representative gauge is even more severe.  If `Phi` is the
periodic Rule 30 ANF map, then

```text
Psi_i=Phi_i^2+x_i^2+x_i                                (31)
```

has the same Boolean truth table and formal Jacobian `I`, yet remains
noninjective because the all-zero and all-one states have the same image.  At
period two this is

```text
Psi(x,y)=(x+x^2 y^2, y+x^2 y^2).                       (32)
```

Its characteristic-zero Jacobian determinant is
`1+2xy^2+2x^2y`, not a constant.  This is a
characteristic/representative boundary, not a counterexample to `JC(2)` over
characteristic zero.

## 7. Preservation and loss ledger

```text
source:     free-input left-permutive CA / Rule 30 Boolean masks
target:     trace x right-half coordinates; typed prefix cylinders;
            pair-overlap masks; factorial and polynomial lifts
map:        (3), (12)--(13), (18), (22), (28)--(29)
preserves:  exact trace prefixes and free-input fibres; Boolean truth only
            for the explicitly faithful lifts; pairwise mask intersections
destroys:   distinguished-seed boundary, path after endpoint quotient,
            affine LRC gaps/common time, factorial carrier under Boolean
            quotient, and either truth or Keller form under the two lifts
sidecar:    fixed seed/right half; full causal address; labelled LRC cochain;
            characteristic and polynomial representative
consequence: exact compilers and sharp hostiles; no Rule 30 prize, LRC(14),
             FC(3), or JC(2) conclusion
```

## 8. Exact audit

The companion verifies the Rule 30 truth table, trace fibres through horizon
eight, codec round trips through six, open blocks and iterates, the three
single-seed controls, `299,592` small mask families, the `q=6` masks, both
factorial compilers, cyclic all-one predecessors through length twelve, and
all displayed prefixes and recurrences.

```bash
python3 04-computation/rule30_left_permutive_trace_torsor_thm3456.py
python3 -O 04-computation/rule30_left_permutive_trace_torsor_thm3456.py
```

Both runs reproduce the stored output byte for byte.  The finite ranges audit
the implementation; the universal quantifiers are discharged by the proofs
above.
