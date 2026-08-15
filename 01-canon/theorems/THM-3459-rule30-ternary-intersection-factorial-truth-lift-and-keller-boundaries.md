---
id: THM-3459
title: "Rule 30 ternary intersection, factorial truth lift, and Keller boundaries"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  Three rotated Rule 30
  gates recover all pair-overlap defects and give a sharp ternary-code global
  disjointness compiler.  The faithful real truth lift has factorial moments
  L(g)=0 and L(g^2)=6 with an exact derangement compiler; the positive ANF
  lift has a closed formal OGF and recurrence.  The tame ANF, faithful
  multilinear, periodic, and characteristic-two polynomial representatives
  preserve incompatible predicates, yielding exact LRC/FC/JC stopping
  boundaries but no open-problem conclusion.
source: root-rule30-260813291-20260815
audit: >
  independent all-rule trace controls, global-mask identity, cyclic gluing,
  factorial compilers, finite-width Jacobians, hostile optimized mutation,
  dependency/scope, hash, and ordinary/optimized/stored byte-replay audit
depends_on:
  - THM-2810-factorial-hankel-faithfulness-and-bounded-radial-carrier-no-go
  - THM-3395-small-sheet-typed-cover-star-cochain
  - THM-3456-left-permutive-trace-bijection-and-rule30-seed-boundary
related:
  - THM-2050-period14-top-germs-do-not-determine-global-loneliness
  - THM-3359-modular-c-finite-supports-harmonic-density-and-periodic-scar
  - THM-3367-berggren-spinor-pencil-hessian-gauge-and-affine-line-keller-closure
  - THM-3458-rule30-right-edge-2-adic-odometer-and-moving-observer-boundary
script: 04-computation/rule30_frontier_ports_thm3459.py
output: 05-knowledge/results/rule30_frontier_ports_thm3459.out
script_sha256: 06df8bb5834853b3623ffa04dfa034b5a7189363f216940794eaf67a786756de
output_sha256: 98f3f8ec6776c8e1d015ec9ca7a09ca162d0f65c0cecb9e85faf16eed20f8fba
hash_basis: raw bytes
---

# THM-3459 -- Rule 30 ternary intersection, factorial truth lift, and Keller boundaries

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

THM-3456 already proves the finite-alphabet left-permutive trace bijection,
its explicit single-seed boundary sidecar, and the enriched trace-incidence
`SOP`/`SOP2`/`SOP3` statements.  This theorem does not duplicate that result.
It starts from the Rule 30 local identity and asks which LRC, factorial,
sequence, and polynomial predicates survive four exact representations.

## 1. Rule 30 as a ternary intersection compiler

Over `F_2`, Rule 30 is

```text
f(l,c,r)=l xor c xor r xor c r
        =l xor (c or r).                                (1)
```

For three subsets `A,B,C` of a finite ground set, apply (1) pointwise and put

```text
W_A=A xor (B union C),       P=A xor B xor C.           (2)
```

Then

```text
W_A xor P=B intersect C.                                (3)
```

Cyclically rotating the inputs recovers all three pair intersections.

For `r>=1` labelled masks `M_0,...,M_(r-1)`, give the labels distinct ternary
codewords of length

```text
d=ceil(log_3 r).                                        (4)
```

At digit `j`, union the masks whose label digit is zero, one, or two and run
the three rotations.  The union of every reported defect is exactly

```text
union_(i<k) (M_i intersect M_k).                        (5)
```

Indeed, an overlapping pair differs in some code digit and is then separated
into two groups, while every reported group intersection expands into one or
more owner-pair intersections.  A report certifies existence but need not
identify the pair responsible.  Hence `3d` Rule 30 gates give an exact global
pairwise-disjointness test.  Adding

```text
union_i M_i = ground set                                (6)
```

certifies a partition.  The depth (4) is sharp within this nonadaptive
three-group coding scheme because `3^d` distinct labels are necessary.

This is a symmetric defect-colour compiler, not a tournament: intersections
have no intrinsic orientation.

### LRC boundary

The compiler is a lawful early sieve and not an LRC exit.  In THM-3395's
`q=6` controls, both speeds `(2,8,10)` and the hostile `(2,8,14)` have the
same partition masks

```text
{0,3}, {1,4}, {2,5},                                    (7)
```

but only the first has a compatible affine star cochain.  Passing to masks
destroys the gaps and common source time.  The lawful pipeline is therefore

```text
mask-disjointness sieve -> labelled star solver -> physical clock.     (8)
```

No arrow may be omitted.

## 2. Periodic gluing is a separate gate

For comparison, a one-step open spatial block map from `m+2` input bits to
`m` Rule 30 output bits has exactly four preimages of every output: choose the
two right boundary bits and solve the remaining inputs right-to-left by left
permutivity.  Iterating through depth `T` gives exactly `4^T` preimages from
length `m+2T` to length `m`.

On an `N`-cycle, the all-one output word has a Rule 30 predecessor iff `3|N`.
When it does, the predecessors are exactly the three rotations of

```text
(100)^(N/3).                                             (9)
```

To prove this, choose a one in a predecessor.  The adjacent all-one output
conditions force two zeros and then another one; propagation closes only in
period three.  Thus the open spatial fibres do not survive an untyped cyclic
quotient.  In particular the all-one word on a 13-cycle is
unreachable in this naive one-step encoding.  This says nothing about the
actual LRC phase reachability or safety predicate.

## 3. Faithful truth lift: an exact factorial near miss

Let

```text
L(l^a c^b r^e)=a! b! e!.                               (10)
```

The unique real multilinear polynomial agreeing with Rule 30 on the Boolean
cube is

```text
g=l+c+r-c r-2l c-2l r+2l c r.                          (11)
```

Directly,

```text
L(g)=0,                 L(g^2)=6.                      (12)
```

Thus `g` is an FC(3) one-moment near miss detected at the next power, not a
factorial-conjecture counterexample.

All its moments have an expansion-free compiler.  Let `D_j=!j` be the
derangement number and define

```text
M_k=sum_(j=0)^k (-1)^j binom(k,j)D_j^2.                (13)
```

Then, for every `m>=0`,

```text
L(g^m)=sum_(q=0)^m binom(m,q)q!
          sum_(s=0)^q binom(q,s)(-2)^s M_(m-q+s).      (14)
```

Proof: for independent unit exponentials `C,R`, put
`o=C+R-CR=1-(1-C)(1-R)`.  Since
`E[(1-C)^j]=(-1)^jD_j`, equation (13) is `E[o^k]`.  Now
`g=o+l(1-2o)`; expansion first in `l` and then in `o` gives (14).

The Boolean and factorial carriers cannot be silently identified.  Boolean
truth imposes `x_i^2=x_i`; the characteristic-two factorial quotient has
`x_i^2=0`.  Their common quotient collapses `x_i`, and `p=2` is outside the
good-prime range of the live FC frontier.

## 4. The ANF lift has a closed moment sequence

The integer ANF representative

```text
p=l+c+r+c r                                               (15)
```

has positive factorial moments.  Put

```text
u_n=L(p^n)/n!.                                           (16)
```

Then, as formal power series,

```text
U(z)=sum_(n>=0)u_nz^n
    =(1-z)^(-3)Phi(z/(1-z)^2),
Phi(w)=sum_(d>=0)d!w^d.                                 (17)
```

Equivalently,

```text
u_n=sum_(d=0)^n d! binom(n+d+2,n-d),                    (18)
```

and, for `n>=4`,

```text
u_n=(n+2)u_(n-1)-(2n-1)u_(n-2)
       +(n-3)u_(n-3)+u_(n-4),                          (19)
u_0,u_1,u_2,u_3=1,4,13,45.                             (20)
```

To prove (17)--(18), let `d` count selected bilinear factors `cr` in `p^n`.
After the factorial readout, the remaining `l,c,r` counts have generating
term

```text
d! z^d(1-z)^(-1)(1-z)^(-(d+1))(1-z)^(-(d+1)).          (21)
```

Summing over `d` gives (17).  Moreover

```text
w^2Phi'(w)+(w-1)Phi(w)+1=0,                            (22)
```

which becomes

```text
z^2(1-z)^2U'(z)+(-1+3z-3z^2+z^4)U(z)+(1+z)=0.         (23)
```

The coefficient of `z^n` is (19).  This is an efficient exact sequence
compiler, not a zero-moment mechanism.

## 5. Four polynomial representatives and the Keller boundary

Fix a finite width `w>=1`, take `0<=k<w`, and set exterior negative indices
to zero.  The open-edge ANF lift is

```text
K_k=x_k+x_(k-1)+x_(k-2)+x_(k-1)x_(k-2).                (24)
```

It is unit triangular over characteristic zero, has

```text
product_(k=0)^(w-1) partial K_k/partial x_k=1,          (25)
```

and a recursive polynomial inverse.  Every finite truncation is a tame
Keller automorphism.  Width one is the faithful identity.  For every `w>=2`
the lift is not truth-faithful: at width two, `x_1=x_0=1` gives integer output
`K_1=2` instead of Boolean zero; for `w>=3`, the local input `(0,1,1)` gives
integer `3` instead of the Rule 30 bit `1`.

Conversely, let

```text
o_k=x_(k-1)+x_(k-2)-x_(k-1)x_(k-2),
B_k=x_k+o_k-2x_k o_k.                                   (26)
```

This is the faithful real multilinear lift, but

```text
det Jac(B)=product_(k=0)^(w-1)(1-2o_k)                 (27)
```

is nonconstant from width two onward.

Closing the ANF lift on an `N`-cycle gives

```text
F_i=x_(i-1)+x_i+x_(i+1)+x_i x_(i+1).                  (28)
```

At zero, `Jac(F)=I+S+S^(-1)`; its determinant vanishes for `3|N` and has
absolute value three otherwise.  At all minus one the Jacobian is `S^(-1)`
and has determinant of absolute value one.  Thus the periodic lift is not
Keller for any `N>=3`.

Finally over `F_2`,

```text
Psi_i=F_i^2+x_i^2+x_i                                  (29)
```

has the same Boolean truth table and formal Jacobian `I`, yet remains
noninjective because the all-zero and all-one states collide.  At period two,

```text
Psi(x,y)=(x+x^2y^2,y+x^2y^2),                          (30)
```

whose characteristic-zero determinant is
`1+2xy^2+2x^2y`, not a constant.

The four representatives preserve complementary predicates:

| boundary / lift | truth faithful | Jacobian status |
|---|---:|---:|
| open ANF over characteristic zero | only `w=1` | tame Keller |
| open real multilinear | yes | identity at `w=1`; nonconstant for `w>=2` |
| periodic ANF over characteristic zero | no | nonconstant/vanishing |
| periodic Frobenius gauge over `F_2` | yes | formal `I`, noninjective |

This is a characteristic/representative/boundary ledger.  It supplies tame
controls and exact no-gos for naive lifts, not a reduction of `JC(2)`.

## 6. Preservation and loss ledger

```text
source:      Rule 30 Boolean masks and truth table
targets:     overlap defects; factorial sequences; polynomial maps
maps:        (2)--(5), (11), (15), (24), (26), (28)--(29)
preserves:   global pairwise-disjointness; Boolean truth only for faithful
             lifts; every declared factorial moment; finite Keller form only
             for the open ANF representative
destroys:    owner identity in an aggregate defect, affine LRC common time,
             factorial carrier under Boolean quotient, and either truth or
             Keller form under the canonical characteristic-zero lifts
sidecars:    labelled star cochain/physical clock; all powers; boundary;
             characteristic and polynomial representative
consequence: exact prefilters, moment/sequence compilers and stopping
             boundaries; no Rule 30 prize, LRC(14), FC(3), or JC(2) result
```

## 7. Exact audit

The optimization-stable companion rechecks THM-3456's trace controls, all
Rule 30 Boolean inputs, open blocks/iterates, `299,592` small mask families,
the `q=6` masks, both factorial compilers, and cyclic all-one predecessors
through length twelve.  Independent audits extended all sixteen binary
left-permutive rules, the moment identities, the recurrence, and the
finite-width Jacobians.  Deliberately corrupted optimized runs fail closed.

```bash
python3 04-computation/rule30_frontier_ports_thm3459.py
python3 -O 04-computation/rule30_frontier_ports_thm3459.py
```

Normal, optimized, and stored outputs are byte-identical.  Finite ranges audit
the implementation; equations (3)--(5), (9), (14), (17)--(23), and
(25)--(30) prove the universal statements.
