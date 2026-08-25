---
id: THM-4059
title: "Stern--Brocot depth packet character and divisor-star convolution"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. Modular inversion inside every unit packet
  preserves full Stern--Brocot depth and admits an explicit parity formula
  from the mod-two Farey-flank hexagon.  The depth-character packet sums
  satisfy an exact binary-Farey functional equation and recurrence, a complete
  mod-four classification, a divisor convolution with the maximal-vertex
  imbalance of the THM-4057 tournament, and arbitrary weighted/twisted-Fourier
  THM-4056 clock identities.  Their primitive opposite-parity columns are the
  fixed-height sums of the Berggren A-Walsh character.  Exact Legendre/Paley
  pockets occur at 5 and 13.  These results prove no asymptotic, named-constant
  irrationality, Duffin--Schaeffer correlation estimate, LRC(14), or global
  classification of character packets.
source: codex-khinchin-ds-rational-tournament-20260824
audit: >
  Independent hostile referees returned promotion verdicts after checking
  every proof boundary and replaying three normal/optimized companions. The primary
  checks canonical Farey flanks, full-depth modular inversion,
  and the inverse parity law through q=1200; all 7,600,457 unit pairs through
  q=5000; an independent binary-Farey recurrence and the mod-four law; direct
  tournament stars and Mobius inversion; arbitrary rational and signed
  Duffin--Schaeffer weights on six clocks; the complete 202,861-node primitive
  Berggren spinor slice of height at most 1000; ternary rows; Paley pockets;
  and hostile nonmultiplicativity controls in normal and optimized modes.
depends_on:
  - THM-4056-divisor-phase-compiler-and-duffin-schaeffer-lcm-clock
  - THM-4057-stern-brocot-depth-pullback-and-rational-edge-tournament-gauge
related:
  - THM-873-ramanujan-fourier-expansion-of-interval-core-good-sets
  - THM-2632-farey-v4-theta-channel-and-hurwitz-crt-parity-sidecar
  - THM-3357-berggren-three-branch-walsh-level-collapse-and-parent-circuit
  - THM-3509-reduced-fraction-harmonic-k4-face-and-fibonacci-unit-cassini-ray
  - THM-3756-odd-square-ordinal-berggren-affine-descent
script: 04-computation/stern_brocot_depth_packet_character_thm4059.py
output: 05-knowledge/results/stern_brocot_depth_packet_character_thm4059.out
script_sha256: 4279d31e869764ea80febfdfaad90689cb5e679847929c04cb7ed5768efd56f4
output_sha256: 313f79dc4f21b69a62847796211a55230264d5d4ccf7abb8796c6ac702fdd3d8
secondary_script: 04-computation/stern_depth_packet_divisor_star_thm4059.py
secondary_output: 05-knowledge/results/stern_depth_packet_divisor_star_thm4059.out
secondary_script_sha256: 875ab2ea14eca894b34d754d106b74ad963057926cb20833515442e675cfa275
secondary_output_sha256: a0c39269d766aad2aaa522b86c4ae366d8f12370bcbb5fc3b6a05886333d86be
independent_audit_script: 04-computation/stern_depth_packet_divisor_star_thm4059_independent_audit.py
independent_audit_output: 05-knowledge/results/stern_depth_packet_divisor_star_thm4059_independent_audit.out
independent_audit_script_sha256: c792be30b3bfa5efc0d0d1b821ad99cb09adddef080aa758bab53209fafd72bd
independent_audit_output_sha256: 65c39db3cda9994cc492246984d9ab86ea3682b19e45339dac684995dc775219
hash_basis: raw bytes
---

# THM-4059 -- the Stern--Brocot checkerboard across denominator packets

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

**Terminology warning.** Here “depth character” means the `{+1,-1}`-valued
Stern checkerboard weight.  It is generally not a group character of `U_q`;
the exact character pockets and the first obstruction are classified only in
the finite scope stated in Section 6.

THM-4056 grades reduced rationals by exact denominator and THM-4057 grades
them by Stern--Brocot depth.  This theorem computes their first exact
interaction.  It also separates three involutions which are easy to conflate:

```text
a/q -> q/a                 rational reciprocal,
a/q -> (q-a)/q             subtree mirror,
a in U_q -> a^(-1) in U_q  modular inversion.           (1)
```

The first changes denominator, the second and third remain in one denominator
packet, and only the first is endpoint reversal in the user's rational-edge
dictionary.

## 1. Full depth and its modular-inverse parity formula

For coprime positive `p,q`, let

```text
D(p/q)=sum of the canonical continued-fraction digits - 1,
epsilon(p,q)=(-1)^D(p/q).                               (2)
```

THM-4057 proves that `D` is Stern--Brocot edge depth and that both rational
reciprocal and, on `(0,1)`, the mirror `p/q -> (q-p)/q` preserve it.

Fix `q>1` and `a in U_q`, represented by `1<=a<q`.  Let

```text
1<=u<q,             au=1+kq.                            (3)
```

Thus `u` is the internal modular inverse of `a`, and `0<=k<a`.

### Theorem 1.1 (modular inversion preserves full depth)

For every such pair,

```text
D(u/q)=D(a/q).                                         (4)
```

If

```text
a/q=[0;c_1,...,c_r]                                    (5)
```

is canonical and `v` is the denominator of its penultimate convergent, the
continuant determinant gives

```text
av=+1 mod q  or  av=-1 mod q.                          (6)
```

Reversing the continuant word gives `v/q` and preserves the sum of the
digits.  If the reversed word ends in `1`, canonicalizing it replaces
`[...,c,1]` by `[...,c+1]`, again preserving the digit sum.  Hence
`D(v/q)=D(a/q)`.  Equation `(6)` says `u=v` or `u=q-v`; the latter is the
depth-preserving subtree mirror.  This proves `(4)`.

### Theorem 1.2 (closed inverse-parity formula)

The depth sign is

```text
epsilon(a,q)=
  (-1)^(a+u),       q odd,
  (-1)^(k+1),       q even.                            (7)
```

To prove `(7)` without assuming the continued-fraction word, form the
canonical Farey-parent flank

```text
M(a,q)=[ k    a-k ]
       [ u    q-u ],             det M=-1.             (8)
```

Its columns are nonnegative, sum to `(a,q)^T`, and have determinant `-1`, so
they are the unique oriented parents of `a/q`.  The root flank is

```text
J=[0 1;1 0].                                           (9)
```

Each left or right Stern--Brocot descent right-multiplies the flank by a
column-addition shear.  Modulo two, the two shears are the adjacent
transpositions of the three nonzero vectors of `F_2^2`.  Therefore the sign
of the resulting permutation relative to `J` is exactly `(-1)^D`.  This is
the bipartite character of THM-2632's Cayley-graph `C_6` on `S_3`, not one
of its three theta channels.

The six residue cases finish the proof.  When `q` is odd,
`k=1+au mod 2`, and the table is

| `(a,u) mod 2` | `M mod 2` | relative sign |
|---|---|---:|
| `(0,0)` | `[1 1; 0 1]` | `+1` |
| `(0,1)` | `[1 1; 1 0]` | `-1` |
| `(1,0)` | `[1 0; 0 1]` | `-1` |
| `(1,1)` | `[0 1; 1 0]` | `+1` |

This is `(-1)^(a+u)`.  When `q` is even, both `a,u` are odd.  For even and
odd `k`, respectively, the two flanks are

```text
[0 1;1 1] -> -1,             [1 0;1 1] -> +1,          (10)
```

which is `(-1)^(k+1)`.  This includes the boundary `q=2`.

The same flank proof gives a global positive-rational form. For arbitrary
coprime `p,q>0` with `q>1`, let `r` be the least positive inverse of `p`
modulo `q` and write `pr-qh=1`. Then

```text
D(p/q) ≡ pq+pr+rh (mod 2).                              (10a)
```

Indeed the parent flank is `M=[h,p-h; r,q-r]`. Set
`G=JM=[r,q-r; h,p-h]`. The sign of `[a,b;c,d] in GL_2(F_2)` is
`(-1)^(ac+bc+bd)`; substitution gives `(10a)`. For `q=1`, separately,
`D(p/1)=p-1`. Formula `(7)` is exactly the specialization `1<=p<q`.

When `q` is even, the carry `k` is also the inverse-lift sheet bit:

```text
a^(-1) mod 2q = u       if k is even,
                  u+q   if k is odd.                    (10b)
```

Thus the sign is negative on the stationary inverse lift and positive on the
switched lift. This connects the checkerboard to the dyadic clock extension
without making it multiplicative.

## 2. Packet sums and a binary-Farey recurrence

Define the exact-denominator depth weight, using THM-4056's total cyclic
compiler convention,

```text
epsilon_1(0)=S(1)=1            (declared bookkeeping weight),
S(q)=sum_(a in U_q) epsilon(a,q),       q>1.             (11)
```

This does not assign a Stern depth to `0/1`. Every exact-denominator `S`-
packet sum below excludes the bookkeeping `q=1` term; in particular `(15)`
starts at `q=2`. Also,

```text
S(2)=-1.                                                    (12)
```

Let the locally finite formal series over all positive primitive pairs be

```text
F(X,Y)=sum_(p,q>=1, gcd(p,q)=1) epsilon(p,q)X^pY^q.     (13)
```

### Theorem 2.1 (binary-Farey functional equation)

One has

```text
F(X,Y)=XY-F(XY,Y)-F(X,XY).                              (14)
```

Every primitive pair other than `(1,1)` has a unique parent obtained by
subtracting its smaller coordinate from its larger one.  The two children of
`(p,q)` are `(p,p+q)` and `(p+q,q)`, and each increases `D` by one.  Summing
the two sign-reversing children gives `(14)`.

The map `(a,q-a)->(a,q)` is one of those children, so

```text
sum_(q>=2) S(q)t^q=-F(t,t).                             (15)
```

Using the symmetry `F(X,Y)=F(Y,X)` in `(14)` yields the exact recurrence

```text
S(q)=2 sum_(2a+b=q, gcd(a,b)=1) epsilon(a,b),   q>2.   (16)
```

### Corollary 2.2 (complete mod-four law)

For `q>2`,

```text
S(q)=2 mod 4
```

if and only if

```text
q=4,
or q=p^e,
or q=2p^e,
```

where `p` is an odd prime with `p=3 mod 4`; in every other case

```text
S(q)=0 mod 4.                                           (17)
```

Modulo two the signs in `(16)` are all one.  Its index set is

```text
1<=a<q/2,          gcd(a,q)=1,
```

which has `phi(q)/2` elements.  Hence `S(q)=phi(q) mod 4`.  Factoring the
totient shows `phi(q)=2 mod 4` exactly in the three cases displayed above.
This also proves that every `S(q)` with `q>2` is even; `q=2` is the unique
mirror-fixed exception.

## 3. The tournament star is the divisor transform of the packet character

Let `sigma_D` be THM-4057's depth tournament gauge.  At the maximal vertex
`q`, define its lower-star imbalance

```text
B(q)=sum_(1<=a<q) sigma_D(a,q).                          (18)
```

Thus `B(q)` is incoming-from-smaller minus outgoing-to-smaller.  Reducing
`a/q` by `g=gcd(a,q)` gives exact denominator `d=q/g`.  THM-4056's compiler
groups these lower vertices into one copy of every `U_d`, `d|q`, so

```text
A(q)=sum_(d|q) S(d),
B(q)=sum_(d|q,d>1) S(d)=A(q)-1,
S(q)=sum_(d|q) mu(q/d)A(d),                             (19)
S(q)=sum_(d|q) mu(q/d)B(d),       q>1.                 (19a)
```

Here `A` is the total cyclic compiler including its declared zero phase;
`B(1)=0` and `B` is the positive lower star. Equation `(19a)` is the Möbius
inverse of the zero-phase-excluded divisor sum, not a second value for `S(1)`.

Consequently the indegree and outdegree of the maximal vertex in the
tournament on `{1,...,q}` are

```text
indeg(q)=(q-1+B(q))/2,
outdeg(q)=(q-1-B(q))/2.                                 (20)
```

The first values are

```text
S(2..20)=
-1,2,-2,0,-2,6,0,-2,-4,2,-4,0,2,8,-4,-4,2,10,0,

B(2..20)=
-1,2,-3,0,-1,6,-3,0,-5,2,-7,0,7,10,-7,-4,-1,10,-7.   (21)
```

The sequence `S` is not multiplicative:

```text
S(2)S(5)=0 != -4=S(10),
S(3)S(5)=0 != 8=S(15).                                 (22)
```

Since `|S(q)|<=phi(q)<=q`, absolute convergence for `Re(s)>2` gives

```text
sum_(N>=1) B(N)/N^s = zeta(s) sum_(q>=2) S(q)/q^s.      (22a)
```

This is a divisor-convolution factorization, not an Euler product: `S` is
nonmultiplicative. Also, `(18)--(20)` describe the lower star, equivalently
the full star only when `q` is the apex of the initial tournament. They do
not give the degree of vertex `q` inside a larger finite tournament.

## 4. Weighted LCM clocks and twisted Ramanujan packets

For `x!=0 in C_N`, reduce

```text
x/N=a/d_N(x)
```

and put

```text
chi_N(x)=epsilon(a,d_N(x));             chi_N(0)=0.     (23)
```

Here `chi_N` is the positive-packet character, extended by zero only for
notation on all of `C_N`.  Separately, the total compiler character
`tilde chi_N=chi_N+1_{0}` assigns the declared bookkeeping weight `+1` at
zero, without assigning a Stern depth to `0/1`.  Thus the unweighted sum of
`chi_N` is `B(N)`, while the unweighted sum of `tilde chi_N` is `A(N)`.

For every commutative-ring-valued denominator weight `W`, THM-4056's labelled
compiler gives

```text
sum_(x in C_N, x!=0) W(d_N(x))chi_N(x)
  =sum_(d|N, d>1) W(d)S(d).                            (24)
```

The unweighted left side is `B(N)`.  On the four LCM clocks it is

```text
N:       6,   60,   420,   27720,
B(N):   -1,  -11,   -39,    -151.                       (25)
```

For a complex-valued `W`, define the twisted primitive packet

```text
E_d(k)=sum_(a in U_d) epsilon(a,d)exp(2*pi*i*k*a/d).    (26)
```

Then every clock Fourier mode has the exact refinement

```text
sum_(x!=0) W(d_N(x))chi_N(x)exp(2*pi*i*k*x/N)
  =sum_(d|N,d>1) W(d)E_d(k).                            (27)
```

At `k=0`, `(27)` is `(24)`.  In general `E_d` is not a Ramanujan sum because
the depth sign is numerator-specific already at `d=5`.  Equation `(27)` is a
twisted packet transform, not the denominator-only transform of THM-4056.

Taking

```text
W(d)=2psi(d)/d                                          (28)
```

gives the exact formal **signed two-radius/multiplicity mass** of the
divisor-complete Duffin--Schaeffer block on `d|N`, `d>1`; each summand is one
exact-denominator layer. It is Fourier/multiplicity data, not a signed
Lebesgue union measure once arcs can saturate, and never an infinite
correlation estimate.
For the prefix `d<=Q` in `C_(L_Q)`, the exact-denominator filter remains
mandatory.

## 5. Berggren rows and LCM-height columns

Restrict `U_q` to primitive Euclid spinors `(a,q)` of opposite parity and put

```text
G_B(q)=sum epsilon(a,q) over that subpacket.             (29)
```

### Theorem 5.1 (fixed-height Berggren column)

For every `q>1`,

```text
G_B(q)=S(q),       q even,
G_B(q)=S(q)/2,     q odd.                               (30)
```

When `q` is even, every unit `a` is odd.  When `q` is odd, the mirror
`a->q-a` pairs an even and an odd unit and preserves full depth.  Each pair
therefore contributes twice the retained even-numerator sign, proving `(30)`.

Under the standard Berggren code, THM-4057 proves

```text
epsilon(w(1/2))=-(-1)^(#A(w)).                          (31)
```

Thus every complete ternary-depth row has signed mass

```text
-sum_(|w|=h)(-1)^(#A(w))=-(2-1)^h=-1,                  (32)
```

whereas `(30)` gives the orthogonal fixed-denominator-height columns.  This is
a row/column relation between ternary ancestry and LCM height, not an
identification of the two gradings.

For an even divisor clock, the signed primitive-Pythagorean phase mass is

```text
P(N)=sum_(q|N,q>1) G_B(q).                              (33)
```

Together with the unsigned THM-4056 count, the master-clock rows are

| `N` | primitive phases | signed depth mass `P(N)` |
|---:|---:|---:|
| `6` | `4` | `-2` |
| `60` | `52` | `-16` |
| `420` | `367` | `-45` |
| `27720` | `25987` | `-201` |

## 6. Exact Paley pockets and a Khinchin hostile

The numerator function `a->epsilon(a,q)` is usually not a group character.
Even denominators already fail at the identity. The first odd failure is

```text
epsilon(2,9)=epsilon(4,9)=-1 != epsilon(2,9)^2.         (34)
```

The denominator-eleven row supplies another hostile:
`epsilon(2,11)=+1` but `epsilon(4,11)=-1`.

There are nevertheless exact small pockets:

- it is the trivial character on `U_3`, `U_7`, and `U_15`;
- on `U_5` and `U_13`, it is exactly the Legendre symbol.

The negative sets in the latter two rows are

```text
q=5:   {2,3},
q=13:  {2,5,6,7,8,11}.                                 (35)
```

They are precisely the quadratic nonresidues.  Therefore `(26)` is the
ordinary Ramanujan sum at `q=3,7,15`, while the classical quadratic Gauss
evaluation gives, for `p=5,13`,

```text
E_p(0)=0,
E_p(k)=(k/p)sqrt(p)          when p does not divide k.  (36)
```

Since `5,13=1 mod 4`, these are Paley **graph** packets, not Paley
tournaments.  No global classification of character-valued depth packets is
claimed.

The first Paley pocket is also a sharp finite-Khinchin hostile:

```text
2/5=[0;2,2],        4/5=[0;1,4].                        (37)
```

Both have two positive digits, digit product `4`, and geometric mean `2`, but

```text
epsilon(2,5)=-1,             epsilon(4,5)=+1.           (38)
```

So even the exact finite digit length and product do not recover the depth
checkerboard.  This strengthens MISTAKE-231's warning without saying anything
about irrationality of Khinchin's constant.

## 7. Transfer and loss ledger

| source | target / map | preserved | destroyed / required sidecar |
|---|---|---|---|
| `a in U_q` | modular inverse `u` | denominator, full depth/sign, and reconstructible numerator | none: this is a fixed-`q` involution |
| `a in U_q` | scalar `epsilon(a,q)` | checkerboard sign | numerator and full depth magnitude |
| Farey flank | `GL_2(F_2)=S_3` | depth parity / bipartite sign | full digit word and depth magnitude |
| packet signs | scalar `S(q)` | signed imbalance | numerator labels and distribution |
| labelled lower star at `q` | reduced packets indexed by `d|q` | vertex, exact denominator, and sign | none |
| packet totals `{S(d):d|q}` | scalar `B(q)` | total lower-star imbalance `(19)` | individual divisor contributions; inversion needs the family `{B(d):d|q}` |
| positive packet family | `C_N\{0}` | every positive labelled sign and twisted Fourier mode | zero is excluded |
| total compiler extension | `C_N` via `tilde chi_N` | the positive packet family plus its `W(1)`/`A(N)` bookkeeping term | zero has declared weight `+1`, but no Stern depth |
| signed DS layers | `(24),(27)` | formal two-radius multiplicity/Fourier data | union geometry, saturation, and infinite correlations |
| Berggren tree | ternary rows / height columns | A-Walsh sign and exact column mass | branch order if only row/column totals retained |
| digit word | length and product | finite Khinchin-type scalar | checkerboard sign by `(37)--(38)` |

Internal modular inversion in `(3)` is not rational reciprocal, and the
Legendre row at thirteen is not the Paley tournament used elsewhere in the
repo.  These type boundaries are part of the theorem.

## 8. Reproduction and scope

Run

```bash
python3 -B 04-computation/stern_brocot_depth_packet_character_thm4059.py
python3 -B -O 04-computation/stern_brocot_depth_packet_character_thm4059.py
python3 -B 04-computation/stern_depth_packet_divisor_star_thm4059.py
python3 -B -O 04-computation/stern_depth_packet_divisor_star_thm4059.py
python3 -B 04-computation/stern_depth_packet_divisor_star_thm4059_independent_audit.py
python3 -B -O 04-computation/stern_depth_packet_divisor_star_thm4059_independent_audit.py
```

Each normal/optimized transcript pair is byte-identical.  The primary exact
companion includes:

- `437,785` primitive flank cases through denominator `1200`;
- all `7,600,457` unit pairs through denominator `5000`, with an independently
  evaluated binary-Farey recurrence and prefix hash;
- direct tournament stars through `1000`, weighted rational clock controls,
  and the master clocks;
- the complete `202,861` primitive Berggren spinors with height at most
  `1000`, compared as a set with the generated tree;
- nonmultiplicativity, noncharacter, q=1/q=2, prefix, and digit-product
  hostile controls.

The secondary audit independently checks 298,035 arbitrary coprime positive
pairs and every packet/lower star through 3000. The third path checks
1,216,587 modular inverses through denominator 2000, compiler/Möbius recovery
through 1000, the mod-four and half-shell laws, degree directions, and the
minimal q=9 and `(2,5,10)` hostiles.

The functional equation, inverse-depth theorem, mod-four law, divisor
convolution, clock transform, and all-height Berggren row/column formulas are
proved symbolically above; finite computation is only their audit.

Nothing here proves an asymptotic for `S(q)` or `B(q)`, a global character
classification, the Duffin--Schaeffer theorem, irrationality or transcendence
of Khinchin's constant or `e+pi`, LRC(14), or any transfer to planar `JC(2)`.
