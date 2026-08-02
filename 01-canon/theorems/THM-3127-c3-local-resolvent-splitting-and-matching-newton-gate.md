---
id: THM-3127
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
  cubeclass raised to n.  In THM-3081's
  additional coordinate-line Keller scope, its terminal Mobius coordinate
  therefore parametrizes the same local resolvent residue field but cannot
  recover a nonzero V4 character without a second-place or marked gluing
  sidecar.  This is a local C3 exclusion test, not a C3, S4, JC(2), or
  globalization theorem.
source: jc-c3-resolvent-2026-08-02
depends_on:
  - THM-2598-quartic-v4-resolvent-torsor-and-universal-cusp-boundary
  - THM-2685-equivariant-kummer-boundary-parity-completion-and-divisor-residue-gate
  - THM-3046-quartic-resolvent-root-valuation-binary-ternary-clutch
related:
  - THM-3057-tame-quartic-inertia-clutch-index-resonance
  - THM-3059-quartic-twojet-even-jelonek-c3-escape-counterexample
  - THM-3081-terminal-toric-residue-parameter-mobius-rigidity-and-autonomous-decoder
script: 04-computation/jc_c3_local_resolvent_matching_gate_thm3127.py
output: 05-knowledge/results/jc_c3_local_resolvent_matching_gate_thm3127.out
script_sha256: 0f4abf813ae0012e014a474f5d5f06651d4601c31bc5341c1a5faf70cf2be24a
output_sha256: 1993bab698ee5984fe976b8e84d3cb53b1dc4b4ac740c097023ca1313db0240c
hash_basis: LF-normalized bytes
---

# THM-3127 -- a pure `C3` place splits the local `V4` packet

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

## 6. Exact controls and boundaries

Run

```bash
python3 04-computation/jc_c3_local_resolvent_matching_gate_thm3127.py
python3 -O 04-computation/jc_c3_local_resolvent_matching_gate_thm3127.py
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

The theorem does **not** assert that every `C3` place has trivial global
`V4` torsor, that a larger decomposition group is irrelevant before strict
henselization, that `(24)` is a physical matching root, or that an arbitrary
Jelonek component can be straightened to a coordinate line.  It excludes a
claimed local `C3` model only when one of `(7)--(11)`, or the
completion intersection ledger fails.  It does not exclude all `C3`, `S4`,
degree-four Keller maps, or `JC(2)`.
