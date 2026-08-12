---
id: THM-3353
title: "Split-prime parabolic branch transplant and unary transducer compiler"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  For every rational
  split prime p=a^2+b^2, an explicit odd-cusp Euclidean compiler produces two
  p^2-spaced U-spine subsequences on which the p-coordinate edge of the
  fixed-hypotenuse parent torsor has an exact Berggren address
  D^((p-1)/2) U^(s-1) S, modulo unit/global-conjugation gauge.  The formulas,
  valuation-one invoice, two complementary local roots, unary transducers,
  and arbitrary-rank dispersion are proved below, exact-tested, and
  independently hostile-audited.
source: codex-2026-08-12-split-prime-parabolic-transplants
depends_on:
  - THM-3334-berggren-parabolic-spine-gaussian-collision-torsor
  - THM-3336-primitive-gaussian-multiplication-content-curved-farey-triangulation
related:
  - THM-3333-gaussian-square-farey-pythagorean-triangular-light-cone
  - THM-3334-berggren-parabolic-spine-gaussian-collision-torsor
  - THM-3336-primitive-gaussian-multiplication-content-curved-farey-triangulation
  - THM-3345-prime-xor-ancestry-path-groupoid-and-source-dependent-berggren-cost
script: 04-computation/split_prime_parabolic_transplants_thm3353.py
output: 05-knowledge/results/split_prime_parabolic_transplants_thm3353.out
script_sha256: 49c1da9fa3fdda4686dda9000b1fb2c68d33e50f07263f159695de362df25dd8
output_sha256: 0f314f40163d1703b4e03fdf8217eca4843dd847f05bc4df158a9fbfb8e00667
hash_basis: LF-normalized bytes
---

# THM-3353 -- split-prime parabolic branch transplants

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

No literature-priority claim is made.  The theorem assembles elementary
Gaussian conjugation, the Euclidean parameter form of the Berggren tree, and
finite-state unary transduction into one explicit compiler.  The quotient
typing below is essential: the literal action on signed Gaussian lifts can
differ from a one-bit toggle by unit and global-conjugation gauges.
[MISTAKE-373](../MISTAKES.md)
records the minimal raw-lift witness and repaired scope.

## 1. Inheritance and the live question

[THM-3334](THM-3334-berggren-parabolic-spine-gaussian-collision-torsor.md)
identifies the consecutive-parameter `U`-spine and its fixed-hypotenuse
Gaussian factor-choice fibres.  [THM-3336](THM-3336-primitive-gaussian-multiplication-content-curved-farey-triangulation.md)
proves that changing one split-prime allocation is a lawful Gaussian operation
with a content invoice.  [THM-3345](THM-3345-prime-xor-ancestry-path-groupoid-and-source-dependent-berggren-cost.md)
then finds one infinite branch formula,

```text
U^(25s)  |-->  DD U^(s-1) ADD,                           (1)
```

for the prime `5`, while leaving analogous split-prime branch transducers
open.  The closest hostile is also THM-3345: no prime colour, folded weight,
fixed descendant word, or rooted automorphism determines the ancestry move on
a whole Boolean fibre.  The least-used sidecar is the odd--odd cusp fixed by
the two otherwise distinct parameter gates `A` and `D`.

This candidate proves that (1) is the first row of an effective two-branch
compiler for **every** rational prime `p=1 mod 4`.  It does not produce one
finite automaton uniform in `p` or on all sources of all Gaussian fibres.

## 2. Parameter matrices and the odd-cusp compiler

Use the Euclid spinor map

```text
Phi(m,n)=(m^2-n^2,2mn,m^2+n^2)
```

and the parameter matrices

```text
B_U=[2 -1; 1 0],       B_A=[2 1; 1 0],
B_D=[1  2; 0 1].                                         (2)
```

Direct expansion gives `B Phi(v)=Phi(Bv)` for the matching three-dimensional
Berggren matrix.  For a root-to-child word `w=w_1...w_l`, write

```text
B(w)=B_(w_l)...B_(w_1).                                  (3)
```

Thus `B(w)(2,1)` is the Euclid parameter pair of the node with address `w`.

Let `p=1 mod 4` be prime and choose its unique ordered representation

```text
p=a^2+b^2,            a>b>0.                             (4)
```

Put

```text
d=a^2-b^2,        e=2ab,       kappa=sign(d-e),
H=[d e; -kappa e kappa d],
h=H(1,1)=(d+e,|d-e|).                                   (5)
```

Then

```text
H^T H=p^2 I,             det(H)=kappa p^2.               (6)
```

The coordinates of `h` are positive, odd, coprime, and unequal.  There is a
unique word `R_p` satisfying

```text
B(R_p)(3,1)=h.                                           (7)
```

Here is a constructive proof.  For a coprime odd pair `(x,y)` with `x>y>0`,
the unique inverse step is

```text
x<2y:       (x,y) -> (y,2y-x)       [U^-1],
2y<x<3y:    (x,y) -> (y,x-2y)       [A^-1],
x>3y:       (x,y) -> (x-2y,y)       [D^-1].              (8)
```

The equalities `x=2y` and `x=3y` are impossible except for the terminal
primitive pair `(3,1)`.  Every nonterminal step preserves positive coprime odd
coordinates and strictly lowers their sum.  It therefore terminates at
`(3,1)`, and the disjoint ratio intervals make the recovered word unique.

The collision that drives the compiler is

```text
B_A(1,1)=B_D(1,1)=(3,1).                                (9)
```

Consequently the two suffixes

```text
S_A=A R_p,                    S_D=D R_p                  (10)
```

have matrices `M_X=B(S_X)` satisfying

```text
M_X(1,1)=h,                   det(M_X)=eta_X in {+1,-1},
eta_D=-eta_A.                                             (11)
```

Equations (7)--(11) are a finite Euclidean algorithm that compiles `p` into
two explicit suffix words.  No search in the Berggren tree is required.

## 3. Determinant matching forces a consecutive intercept

Let `P=[0 1;1 0]`.  For each `X in {A,D}`, choose the unique matrix

```text
H_X in {H,HP}               with det(H_X)=eta_X p^2.     (12)
```

Right multiplication by `P` preserves `H_X(1,1)=h` and reverses the
determinant.  Let `c_X` be the second column of `M_X`, and define

```text
u_X=p(1,1)-H_X^T c_X.                                   (13)
```

The determinant gauge in (12) forces

```text
u_X=(r_X+1,r_X)                                          (14)
```

for an integer `r_X`.  To see this, write the columns of `H_X` as `k_1,k_2`
and let `J(x,y)=(-y,x)`.  Since `H_X^T H_X=p^2I` and its orientation is
`eta_X`,

```text
k_1-k_2=-eta_X J h.
```

Also `det(h,c_X)=det(M_X)=eta_X`.  Therefore

```text
(u_X)_1-(u_X)_2
 =-(k_1-k_2).c_X
 =eta_X det(h,c_X)=1,                                   (15)
```

which proves (14).  This is the hidden role of the `A/D` cusp collision: it
turns a rational Gaussian rotation back into a consecutive integral spinor.

## 4. The universal branch identity

Put `q=(p-1)/2`.  From (2),

```text
B_D^q(2,1)=(p+1,1),
B_U^(s-1)B_D^q(2,1)=(ps+1,p(s-1)+1).                   (16)
```

Since `M_X(1,1)=h`, equations (6), (11), and (13) give

```text
H_X u_X
 =p h-p^2 c_X
 =p M_X(1,1-p).                                         (17)
```

Combining (14), (16), and (17) yields the exact identity, for every integer
`s>=1`,

```text
B(D^q U^(s-1) S_X)(2,1)
 = (1/p) H_X (p^2s+r_X+1, p^2s+r_X).                   (18)
```

Define

```text
t_X(s)=p^2s+r_X.                                        (19)
```

Whenever `t_X(s)>=1`, the right input in (18) is exactly the `U`-spine
spinor `(t+1,t)`.  A uniform explicit lawful threshold is

```text
s_X,0=max(1,ceil((1-r_X)/p^2)).                          (20)
```

For path formulas below we use the harmless stronger threshold `t_X(s)>=2`.

The matrix `H/p` is multiplication by `(a-ib)/(a+ib)` in real coordinates,
up to the output conjugation chosen in (5).  The optional right factor `P` in
(12) is an input unit/conjugation gauge.  On the parent quotient

```text
X_(C_t) ~= F_2^omega(C_t)/<1>,
```

`H_X/p` represents the `p`-coordinate toggle.  On raw signed Gaussian lifts
it may additionally apply unit and global-conjugation gauges, which complement
every allocation bit.  Thus it is not an arbitrary equal-norm rotation, but
the one-coordinate statement is a statement on the parent quotient rather
than on a chosen raw lift.  Equation (6) gives

```text
N(B(D^q U^(s-1)S_X)(2,1))=C_(t_X(s)),
C_t=2t^2+2t+1.                                          (21)
```

The left side of (18) is a Berggren word from the root, so its Euclid pair is
primitive, opposite-parity, and in the positive chamber.

## 5. The toggled prime occurs exactly once

The residue `r_X mod p` is one of the two roots of `C_t=0 mod p`, and in fact

```text
v_p(C_(t_X(s)))=1                  for every s>=1.       (22)
```

The exact valuation is not an experimental pattern.  From (13),

```text
N(u_X)/p = 2p-2 h.c_X+p N(c_X),
N(u_X)/p = -2 h.c_X mod p.                              (23)
```

The vector `h mod p` is nonzero and isotropic.  Its orthogonal line over
`F_p` is therefore its own span.  If `h.c_X=0 mod p`, then `c_X` would be a
multiple of `h`, contradicting

```text
det(h,c_X)=eta_X=+-1.                                   (24)
```

Thus the quotient in (23) is nonzero modulo `p`.  Replacing `u_X` by
`p^2s(1,1)+u_X` does not change its norm divided by `p` modulo `p`, proving
(22).  The branch therefore realizes the `p`-coordinate edge on the parent
torsor quotient, modulo the unit/global-conjugation gauge just stated, with
intrinsic folded weight

```text
{p,C_(t_X(s))/p}.                                       (25)
```

This valuation-one invoice is why the modulus in the compiler is `p^2`, not
merely `p`.

## 6. The two gates are the two local roots

Write `T=B(R_p)` with columns `x,y`.  Then

```text
h=3x+y,             c_A=x,             c_D=2x+y=h-x.   (26)
```

The determinant choices imply, after possibly naming one matrix first,

```text
H_D=H_A P.
```

Substitution in (13) gives the stronger integral relation

```text
u_D=(2p-p^2)(1,1)-P u_A,
r_A+r_D+1=2p-p^2.                                      (27)
```

In particular,

```text
r_D=-1-r_A mod p.                                       (28)
```

The two classes are distinct and exhaust the two roots of `C_t=0 mod p`.
They select one non-Hensel lift modulo `p^2` above each local root; they do not
claim to list all `p-1` valuation-one lifts above that root.

For the first controls, the compiler returns

```text
p=5:   R_5=DD,
       r_A=1,   S_A=ADD;       r_D=-17, S_D=DDD,

p=13:  R_13=AA,
       r_A=-37, S_A=AAA;       r_D=-107,S_D=DAA.        (29)
```

The `A` row at `p=5` is exactly THM-3345's family (1).  The companion family
is

```text
U^(25s-18) |--> DD U^(s-1) DDD.                         (30)
```

## 7. Exact ancestry costs and fixed-prime transducers

Take `s` large enough that `t=t_X(s)>=2`.  The source and target words are

```text
source: U^(t-1),
target: D^q U^(s-1) S_X.                                (31)
```

Their first letters differ, so their longest common prefix is empty.  If
`ell_X=|S_X|`, their unique Berggren-tree path cost and absolute depth jump
are therefore

```text
L_X(s)=(p^2+1)s+r_X+q+ell_X-2,
J_X(s)=|(p^2-1)s+r_X-q-ell_X|.                          (32)
```

Both are unbounded and the compression ratio in the moving `U` block is
exactly `p^2:1`.

For each **fixed** pair `(p,X)`, (31) is a deterministic subsequential
transduction on the regular unary domain

```text
{U^(p^2s+r_X-1):s>=s_X,0}.                              (33)
```

A finite controller checks the exponent modulo `p^2`, emits the fixed prefix
`D^q`, emits one `U` per subsequent block of `p^2` source letters after the
fixed offset, and appends `S_X`.  Thus THM-3345's one unary transducer is not
sporadic: every split prime supplies two effectively compiled transducers.

The number of states, prefix, offset, and suffix depend on `p`.  No claim is
made that their union is one rational relation or that a single finite-state
machine reads arbitrary primes and arbitrary Boolean-fibre sources.

## 8. Every split-prime branch has arbitrary Boolean rank and dispersion

Fix `p`, one gate `X`, and integers `h,B>=1`.  Choose distinct primes

```text
q_1,...,q_h =1 mod 4,                 q_i!=p.
```

Each congruence `C_t=0 mod q_i` has two roots.  Chinese remaindering one root
at every `q_i` together with

```text
t=r_X mod p^2                                           (34)
```

gives an infinite arithmetic progression contained in the compiled branch
(19).  For every such `t`, (22) retains exactly one `p` factor and all the
`q_i` divide `C_t`.  Hence

```text
omega(C_t)>=h+1.                                        (35)
```

THM-3334's fixed-hypotenuse factor-choice theorem makes the corresponding
parent fibre have affine Boolean dimension at least `h`.  Taking a sufficiently
large CRT representative also makes both quantities in (32) exceed `B`.

Therefore **each fixed split prime**, not only `5`, supports simultaneous
arbitrarily high prime-toggle rank and arbitrarily large ancestry path/depth
dispersion along an explicit regular U-spine sublanguage.

## 9. Boundaries, hostiles, and exact evidence

For `p=3 mod 4`, `-1` is not a square modulo `p`; equivalently `C_t` has no
root.  This is the sharp inert-prime hostile.  The two compiled arrows are
still source-dependent operations, not fixed descendant words on a whole
fibre.  Their symmetric prime label supplies no intrinsic tournament
orientation.  The ambient Berggren tree still has zero graph `H^1`, so the
compiler supplies no LRC current, JC flux, or nontrivial ancestry cohomology.

The exact companion checks every split prime below `5,000` (`329` primes),
both branches and `24` rows per branch (`15,792` rows).  It compares the
two-dimensional parameter action with an independent three-dimensional
Berggren evaluation, checks orthogonality, determinants, cusp descent,
complementary roots, valuation one, primitivity, norm preservation, path
costs, and CRT rank controls.  All `50` inert primes below `500` are rootless.
The finite sweep additionally observes `s_X,0=1` throughout its universe, but
the universal theorem uses the explicit safe threshold (20).

The independent hostile audit also checked `1,282` further split primes in
`5,000<p<30,000` by Gaussian division (`2,564` branches and `5,128` rows),
and an independent cusp breadth-first search found every one of the `8,156`
coprime odd pairs with coordinate sum at most `401` exactly once.  Its
minimal raw-lift typing witness is `p=5,X=A,s=1`:

```text
(2+i)(16+5i)=27+26i  |-->  (2+i)(16-5i)=37+6i,
```

whereas literal conjugation of only the `5`-factor gives `37-6i`, the global
conjugate of the displayed target.  Both represent the same `5`-coordinate
edge in the parent quotient, with folded weight `{5,281}`.

Reproduce with

```bash
python3 04-computation/split_prime_parabolic_transplants_thm3353.py
python3 -O 04-computation/split_prime_parabolic_transplants_thm3353.py
```

The two transcripts must match the stored output after LF normalization.
The universal claim rests on equations (4)--(35), not on the finite range.
