---
id: THM-3131
title: "C3 local resolvent splitting and matching-Newton gate"
status: >
  PROVED + VERIFIED-EXACT.  Let a separable depressed quartic with global
  S4 splitting field have a valuation whose decomposition and inertia groups
  are the same C3.  The completion of its ramified moved-sheet quartic field
  is exactly the completion of the full S4/V4 resolvent field and of the
  splitting field.  Hence all three squared-pair-sum Kummer roots are already
  squares locally: a pure C3 divisor carries no local V4 character.  If
  n=v(q), their common splitting-field value is 2n.  When 3 does not divide
  n, both lower resolvent traces cancel strictly and v(Disc)=4n; when 3
  divides n, the two trace values are exact and v(Disc)>4n.  The matching
  clutch is always diagonal with zero ternary checksum.  Although the local
  V4 plane splits, the leading unit of q retains exactly the tame C3
  cubeclass raised to n.  For an actual fixed-plus-escaping-C3 graph quartic,
  this cubeclass test is automatically saturated: if the primitive graph
  coordinate has pole order m, then n=-m and its leading depressed coefficient
  is c_m A_X^3 tau^m, where c_m is -1 off 3|m and 1/8 on 3|m.  In THM-3081's
  coordinate-line Keller scope this becomes q_0=c_m H_X(theta)^3 K(theta)^m;
  the decoder fixes K/L but not the individual terminal-prefactor cubeclass
  [L].  Thus q_0 alone cannot exclude the branch.  This is a local C3 test and
  information-loss theorem, not a C3, S4, JC(2), or globalization theorem.
source: jc-c3-resolvent-2026-08-02
depends_on:
  - THM-2598-quartic-v4-resolvent-torsor-and-universal-cusp-boundary
  - THM-2685-equivariant-kummer-boundary-parity-completion-and-divisor-residue-gate
  - THM-3046-quartic-resolvent-root-valuation-binary-ternary-clutch
  - THM-3081-terminal-toric-residue-parameter-mobius-rigidity-and-autonomous-decoder
related:
  - THM-3057-tame-quartic-inertia-clutch-index-resonance
  - THM-3059-quartic-twojet-even-jelonek-c3-escape-counterexample
script: 04-computation/jc_c3_local_resolvent_matching_gate_thm3131.py
output: 05-knowledge/results/jc_c3_local_resolvent_matching_gate_thm3131.out
script_sha256: 8ba868cad34507c4b32f8112ff541330236a0419336b232ae7c5e813f07f6554
output_sha256: 540191f75468c4ff76124d73474772a5d139fd1da16fc750c354767337160f16
hash_basis: LF-normalized bytes
---

# THM-3131 -- a pure `C3` place splits the local `V4` packet

**PROVED + VERIFIED-EXACT.**

## 1. Exact setup and statement

Let `(K,v)` be a discretely valued field of residue characteristic
different from two and three, and assume that `K` contains the cube roots of
unity.  Let

```text
f(T)=T^4+pT^2+qT+r,                  q!=0,               (1)
```

be separable with splitting field `M/K` and global Galois group `S4`.  Fix an
extension `w` of `v` to `M`.  Assume its decomposition and inertia groups are

```text
D_w=I_w=<sigma> isomorphic to C3.                         (2)
```

Thus `sigma` acts as a three-cycle on the four quartic roots.  Normalize

```text
w(M*)=Z,                    w restricted to K=3v.        (3)
```

Let `V=V4 normal S4`, put `E=M^V`, and let `L_a=M^(H_a)`, where `H_a` is the
stabilizer of a root moved by `sigma`.  At the primes induced by `w`,

```text
M_w = E_w = (L_a)_w.                                     (4)
```

For the three complementary pairings of the depressed roots, let

```text
beta_1=z_0+z_1=-(z_2+z_3),
beta_2=z_0+z_2=-(z_1+z_3),
beta_3=z_0+z_3=-(z_1+z_2),

U_i=beta_i^2.                                            (5)
```

The `U_i` lie in `E`, are cycled by `sigma`, and are the roots of

```text
S(U)=U^3+2pU^2+(p^2-4r)U-q^2.                           (6)
```

The equality `(4)` gives the stronger-than-unramified conclusion

```text
[U_1]=[U_2]=[U_3]=0 in E_w^*/E_w^(*)2.                  (7)
```

Put `n=v(q)`.  Then

```text
w(U_1)=w(U_2)=w(U_3)=2n.                                (8)
```

There is one nonbinary phase which survives the local split.  Choose
uniformizers `t` of `K_v` and `pi` of `M_w` so that

```text
sigma(pi)=zeta pi,                 pi^3=epsilon t,
q=t^n q_0,                         epsilon,q_0 units.    (8a)
```

Then in the common residue field `k` of `(4)`,

```text
[bar(q_0)]=[bar(epsilon)]^n in k^*/k^(*)3.              (8b)
```

Thus `bar(q_0)` is a cube when `3|n`; when `3` does not divide `n`, it
recovers the tame cubic extension class, up to the invertible exponent `n`.

The two possible matching-Newton patterns are exactly

```text
3 does not divide n:
  3v(p)>2n,
  3v(p^2-4r)>4n,
  v(Disc(f))=4n;                                        (9)

3 divides n:
  3v(p)=2n,
  3v(p^2-4r)=4n,
  v(Disc(f))>4n.                                       (10)
```

Here `v(0)=+infinity`.  The inequalities and equality remain valid for
negative valuations.  Thus
`(9)--(10)` are a four-scalar necessary test for a claimed pure `C3` place.

Finally, let `s_i=w(U_j-U_k)` for the three unordered pairs.  The `U_i` are
a common translate of THM-3046's matching-product resolvent roots, so these
are exactly the same three root-difference values.  They are all equal.
Consequently THM-3046's integral matching clutch is

```text
kappa=(s,s,s) mod 2,                  tau=3s mod 3=0.    (11)
```

This is not merely equality of quartic and resolvent discriminants: it fixes
the local Kummer specialization, the two lower Newton coefficients, and the
complete binary/ternary matching quotient.

Condition `(2)` is the clean residue-degree-one formulation.  For a general
tame `C3` inertia place, the same conclusions hold after strict
henselization, which removes the unramified part of the decomposition group.
They need not hold over the unstrictified field if nontrivial residue gluing
is deliberately retained.

## 2. The source and resolvent completions really coincide

For a finite Galois extension with decomposition group `D` and subgroup
`A`, the local degree above the corresponding prime of `M^A` is

```text
[M_w:(M^A)_w]=|D intersect A|.                          (12)
```

Here an order-three subgroup of `S4` intersects `V4` trivially.  It also
intersects the stabilizer of any moved sheet trivially.  Therefore

```text
|D intersect V|=|D intersect H_a|=1,                   (13)
```

and `(12)` proves `(4)`.  By contrast, if `H_b` stabilizes the unique fixed
sheet, then `D subset H_b`; that source branch has local degree one.  The
choice of the moved sheet is therefore load bearing.

This supplies an actual map behind the quartic/resolvent comparison:

```text
ramified quartic source completion -> full resolvent completion
```

is the identity subfield of `M_w`.  Passing to `S4/V4` loses the global
choice among the four `V4` translates, but it loses nothing in this one
local field because those translates become four split local factors.

## 3. The `V4` Kummer plane specializes to zero

Each `U_i` is fixed by `V4`: a double transposition can exchange the two
entries of a pair or exchange the complementary pairs, changing `beta_i` by
at most a sign.  Thus `U_i in E`.  Equation `(4)` places all quartic roots,
and hence every `beta_i`, in `E_w`.  Formula `(5)` immediately proves `(7)`.

There is also an intrinsic squareclass audit which does not hide the role of
valuation normalization.  Since `sigma` cycles the `U_i`, their values are
equal.  Their product is `q^2`, so `(3)` gives

```text
3w(U_i)=w(q^2)=6v(q),                                  (14)
```

which proves `(8)` and in particular makes the common value even.  Choose a
uniformizer `pi` with `sigma(pi)=zeta pi`.  The ratios `U_i/U_j` have value
zero and residues which differ by powers of `zeta^(2n)`.  These constants
are squares, and Hensel's lemma makes each ratio a square.  The three local
squareclasses are therefore equal.  But

```text
[U_1]+[U_2]+[U_3]=[q^2]=0,                              (15)
```

and multiplication by three is the identity on an `F2`-space.  Hence their
common class is zero, again proving `(7)`.

Thus THM-2685's boundary word is not only `000`: the entire standard Kummer
plane restricts trivially in the local squareclass group.  A pure `C3`
divisor cannot itself carry either of the two nonzero `V4` characters needed
by THM-2655.  Any global carrier must be detected by gluing between places,
global units, or `Cl[2]`.

The cubic character has not disappeared.  Because `M_w/K_v` is totally
tamely ramified of degree three, choose `pi` as in `(8a)`.  Write

```text
beta_1=pi^n(b+O(pi)),                    b in k^*.       (15a)
```

The other two leading coefficients are `zeta^n b` and `zeta^(2n)b`.
The definitions `(5)` give the exact sign identity

```text
beta_1 beta_2 beta_3=-q.                                (15b)
```

Taking leading residues in `(15b)` and using
`zeta^(n(0+1+2))=1` gives

```text
bar(q_0)=-bar(epsilon)^n b^3.                           (15c)
```

Since `-1` is a cube, `(15c)` proves `(8b)`.  The local information ledger is
therefore sharp: the binary `V4` packet splits, while the ternary tame phase
is retained by one residue cubeclass.  This cubeclass, not the common
discriminant, is the least sidecar needed to remember the local `C3` cover.

## 4. The exact `3|n` Newton boundary

Let `a` be the leading residue of `pi^(-2n)U_1`.  Inertia acts trivially on
the residue field.  After cyclically ordering the roots,

```text
in(U_1,U_2,U_3)
 =a pi^(2n)(1,zeta^(2n),zeta^(4n)).                     (16)
```

If `3` does not divide `n`, the three displayed residues are distinct and
their sum is zero.  Their pairwise-product sum is also zero.  Comparing with
the first two elementary symmetric functions in `(6)` gives

```text
w(2p)>2n,                     w(p^2-4r)>4n.             (17)
```

Because `w` restricts as `3v`, these are the first two lines of `(9)`.
Moreover each difference `U_i-U_j` has value exactly `2n`.  The root formula
for the cubic discriminant gives

```text
w(Disc(S))=2 sum_(i<j)w(U_i-U_j)=12n.                  (18)
```

For the normalized squared-pair-sum resolvent `(6)`, direct root
factorization gives the exact identity

```text
Disc(S)=Disc(f).                                        (19)
```

Dividing `(18)` by the ramification index three proves `v(Disc(f))=4n`.

If `3` divides `n`, all three leading residues in `(16)` equal `a`.  Their
sum and pairwise-product sum have leading coefficients `3a` and `3a^2`,
which are nonzero in residue characteristic different from three.  Hence

```text
w(2p)=2n,                       w(p^2-4r)=4n.           (20)
```

Every difference now has value strictly greater than `2n`.  Inertia cycles
the three differences, so their values are one common integer `s>2n`.
Equations `(18)--(19)` become

```text
v(Disc(f))=2s>4n,                                      (21)
```

proving `(10)`.  The same cyclic equality of difference values proves
`(11)` in both cases.  This also shows why the boundary is exactly
`3|v(q)`: it is the point where the leading `C3` character becomes trivial.

## 5. Contact with the terminal Mobius decoder

Assume now the additional hypotheses of THM-3081: `(1)` is a primitive graph
quartic for a planar polynomial Keller map, the chosen moved-sheet source
branch lies over a target coordinate line, has tame index three and residue
degree one, and its toric key tower has reached the primitive terminal stage.
THM-3081 gives a terminal value-zero parameter `theta` with

```text
C(u)=C(theta),                theta is Mobius in u.      (22)
```

By `(4)`, this is also the residue field of the full resolvent completion.
Thus the coefficient `a` in `(16)` is a rational function of `theta`, and
the local quartic matching packet has a concrete terminal expression

```text
in(U_i)=a(theta) pi^(2n) zeta^(2ni).                    (23)
```

Moreover `(15a)` has `b in C(theta)^*`, and `(15c)` identifies the surviving
ternary phase as the cubeclass of `q_0`.  The terminal decoder can therefore
record the local `C3` Kummer parameter, but not either nonzero character of
the split `V4` kernel.

Because `C(theta)` is a rational function field over `C`, its cubeclasses are
detected exactly by divisor multiplicities modulo three.  Equation `(8b)`
therefore has the immediately executable residue-divisor form

```text
div_theta(q_0)=n div_theta(epsilon) mod 3.               (23a)
```

In particular, if `3|n`, every zero and pole of `q_0(theta)` has multiplicity
divisible by three.  A single simple zero or pole excludes that local branch.
If `3` does not divide `n`, `(23a)` reconstructs the divisor class of the tame
cubic unit `epsilon` because `n` is invertible modulo three.  This is cheaper
than running the Laurent-key tower: depress the graph quartic, read `q`, and
factor only its leading residue along the proposed component.

This is a useful hostile test for an abstract proposed local packet.  For an
actual fixed-plus-escaping graph quartic it is never an independent
obstruction: Section 6 proves that the graph roots saturate `(23a)`
identically.

But `(7)` says that every displayed matching root represents the zero local
Kummer class.  Therefore THM-3081's square

```text
K(T)/L(T)=constant*T^(A-1)*(linear(T))^2                (24)
```

cannot, at this one place, recover a **nonzero** `V4` character.  The two
squares live on the same local field, but the `V4` origin has split away.
The missing sidecar is necessarily global or multi-place: for example a
marked sheet/half-face, transition data between conjugate divisors, or the
unit/class-group carrier of THM-2655/2685.

This is the precise stopping reason for a direct decoder-to-resolvent
identification.  It replaces a tempting analogy by an actual local map and
an exact account of the information that map destroys.

## 6. Automatic graph-quartic saturation and the missing cubeclass

Retain the coordinate-line Keller scope of Section 5 and assume the actual
fixed-plus-escaping-`C3` graph anatomy.  Let `X` be a primitive affine source
coordinate, let `X_f` be its finite fixed-sheet value, and choose the Kummer
uniformizer from `(8a)`.  Put

```text
tau=epsilon^(-1),                  t=tau pi^3,

X_0=A_X(u)pi^(-m)+O(pi^(-m+1)),   m>=1,
X_i=sigma^i(X_0)
   =A_X(u)zeta^(-im)pi^(-m)+O(pi^(-m+1)),               (25)
```

where `A_X in C(u)^*` and `X_f` has nonnegative value.  Let `z_f,z_i` be the
four roots after formal depression and use the three moved pair sums

```text
beta_i=z_f+z_i=X_f+X_i-(X_f+X_0+X_1+X_2)/2.            (26)
```

The exact identity `(15b)` remains `prod_i beta_i=-q`.

If `3` does not divide `m`, the leading moved trace vanishes because

```text
1+zeta^(-m)+zeta^(-2m)=0.
```

The trace correction in `(26)` has strictly larger value, so

```text
beta_i=A_X zeta^(-im)pi^(-m)+O(pi^(-m+1)),
q=-A_X^3 pi^(-3m)+O(pi^(-3m+1)).                       (27)
```

If `3|m`, the three moved initials are all `A_X pi^(-m)`.  Depression gives

```text
z_f=-(3/4)A_X pi^(-m)+O(pi^(-m+1)),
z_i= (1/4)A_X pi^(-m)+O(pi^(-m+1)),
beta_i=-(1/2)A_X pi^(-m)+O(pi^(-m+1)),

q=(1/8)A_X^3 pi^(-3m)+O(pi^(-3m+1)).                   (28)
```

Since `pi^3=t/tau`, both lanes have the uniform form

```text
n=v_t(q)=-m,
q=t^(-m)(q_0+O(t)),

q_0=c_m A_X^3 tau^m,
c_m=-1 if 3 does not divide m,
c_m=1/8 if 3 divides m.                                (29)
```

Both constants are cubes in `C`.  The expression is independent of the
normal-parameter gauge: under `pi'=h(u)pi`, one has

```text
A_X'=A_X h^m,                  tau'=tau h^(-3),
(A_X')^3(tau')^m=A_X^3 tau^m.                           (30)
```

Thus the horizontal divisor law is not merely a congruence:

```text
div(q_0)=3 div(A_X)+m div(tau),
div(q_0)=n div(epsilon) mod 3.                          (31)
```

When `3|m`, equation `(29)` explicitly reads

```text
q_0=(A_X tau^(m/3)/2)^3.                               (32)
```

When `3` does not divide `m`, a simple zero or pole of `q_0` is accompanied
by the exactly required nontrivial order of `tau` modulo three.  Hence no
zero/pole pattern of `q_0` by itself can violate the THM-3131 cubeclass gate
for this graph branch.

The calculation also occurs directly in the actual, nondepressed graph
coefficients.  If

```text
f(T)=T^4+a_3T^3+a_2T^2+a_1T+a_0,
R(T)=cT^4+r_3T^3+r_2T^2+r_1T+r_0=c f(T),               (33)
```

then formal depression gives exactly

```text
q=a_1-a_2a_3/2+a_3^3/8
 =(8c^2r_1-4cr_2r_3+r_3^3)/(8c^3).                    (34)
```

Off `3|m`, the leading traces in `(34)` cancel and the escaping triple in
`a_1` supplies `-A_X^3`.  On `3|m`, the three tied terms have respective
leading constants `-1`, `9/2`, and `-27/8`, whose sum is `1/8`.  This is the
graph-coefficient version of `(27)--(28)` and explains rather than merely
detects the Newton boundary.

There is a sharper contact with THM-3081.  At its terminal stage choose
integers `A_*,B_*` with `A_*g+B_*e=1`.  Write its Mobius parameter as

```text
theta=(a u+b)/(c u+d),                 Delta=ad-bc,
lead(X)=rho^(-m)H_X(theta),
tau=rho^3 K(theta),
lead(U)=rho^(3-e)L(theta).                              (35)
```

Substitution in `(29)` cancels the entire value-one scale:

```text
q_0=c_m H_X(theta)^3 K(theta)^m,
[q_0]=[K]^m in C(theta)^*/C(theta)^(*)3.                (36)
```

THM-3081's autonomous decoder determines only

```text
K/L= kappa/(3 Delta) theta^(A_*-1)(a-c theta)^2.        (37)
```

Since every nonzero complex constant is a cube, `(36)--(37)` give the exact
transport law

```text
[q_0]
 =[L]^m [theta]^(m(A_*-1)) [a-c theta]^(2m).            (38)
```

If `3` does not divide `m`, the exponent `m` is invertible on cubeclasses,
so `q_0` transports the class of `K` and, after the displayed Mobius
correction, the class of `L`.  If `3|m`, it transports no cubeclass at all.
The decoder controls `K/L`; it does **not** independently control `[L]`.
That terminal two-form-prefactor class is the precise missing coordinate in
a `q_0`-only argument.  Equivalently, before taking the norm, one may retain
one marked leading amplitude of `beta_i`; the product `prod beta_i=-q`
forgets that `C3` phase.

This separates the complementary primes cleanly.  Squaring the `beta_i`
trivializes the local `2`-primary `V4` packet by `(7)`, while multiplying
them retains only the automatically compatible `3`-primary tame class by
`(29)`.  A further exclusion must therefore couple a marked pair sum, the
companion `b` and Keller congruence of THM-2621, or genuinely global/multi-
place algebraization of `L` or `tau`.  No such coupling is proved here.

## 7. Exact controls and boundaries

Run

```bash
python3 04-computation/jc_c3_local_resolvent_matching_gate_thm3131.py
python3 -O 04-computation/jc_c3_local_resolvent_matching_gate_thm3131.py
```

Both executions match the stored transcript.  The companion verifies:

1. `Disc(f)=Disc(S)` symbolically for independent `p,q,r`;
2. the three intersection sizes in `(13)` and the fixed-sheet boundary;
3. the nondivisible family `T^4-tT`, with `n=1` and discriminant value four;
4. the divisible family obtained from
   `beta_i=s^3(1+zeta^i s)`, whose values are
   `(n,v(p),v(p^2-4r),v(Disc))=(3,2,4,14)`;
5. THM-3059's generic `S4` hostile, which has
   `(-2,-1,-2,-8)` and lands exactly on `(9)`.
6. the exact depression formula `(34)`, the fixed-plus-cubic hostile
   `q=-t^(-1)-u^3/8`, both graph leading constants in `(29)`, normal-gauge
   invariance, and the cancellation of `rho` in `(36)`.

The theorem does **not** assert that every `C3` place has trivial global
`V4` torsor, that a larger decomposition group is irrelevant before strict
henselization, that `(24)` or `(37)` is a physical matching root, that `[L]`
is independently controlled, or that an arbitrary
Jelonek component can be straightened to a coordinate line.  It excludes a
claimed local `C3` model only when one of `(7)--(11)`, or the
completion intersection ledger fails.  It does not exclude all `C3`, `S4`,
degree-four Keller maps, or `JC(2)`.
