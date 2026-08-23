---
id: THM-3819
title: "Chebyshev Pell evaluation, Berggren return cocycle, and conductor staircase"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  After evaluation
  of THM-3757's canonical indexed tower at z=2, a lossless integer coordinate
  change identifies the evaluated pair with the positive square-triangular
  Pell state.  Consecutive selected states on THM-3756's
  Berggren L-spine are separated by the state-dependent negative-Pell return
  times 7,41,239,... .  The same return appends exactly that many cyclic
  summands to THM-3745's standardized conductor-defect staircase and acts by
  two iterations of one Mobius map on THM-3744's scalar loneliness maxima.
  Evaluation loses the JC polynomial geometry, the forest inner root is not
  the Pell square-root sidecar, actual conductor rings do not append, and the
  scalar observer loses LRC carries and owners.  The exact companion and an
  independent hostile audit both pass.
source: root + square_operation_wildcard / incoming-signal extension session, 2026-08-23
audit: >
  INTERNAL EXACT AUDIT + INDEPENDENT HOSTILE AUDIT PASS.  The
  assertion-free companion checks 6,982 polynomial, integer, forest,
  conductor-shape, LRC, mod-thirteen, support-two, and hostile gates.  Normal
  and optimized outputs byte-match the frozen transcript.  The independent
  audit rederived every all-depth identity, all fourteen modular states, the
  finite support-two intersection, and the hostiles.  It repaired the wording
  so only the post-evaluation coordinate change, not evaluation of polynomial
  geometry, is called lossless.
depends_on:
  - THM-3742-square-triangular-pell-mod13-central-sign-projective-cycle
  - THM-3744-pell-prefix-loneliness-constant-carry-exact-formula
  - THM-3745-monomial-plane-branch-conductor-triangular-pell-selector
  - THM-3756-odd-square-ordinal-berggren-affine-descent
  - THM-3757-pell-chebyshev-three-charge-hyperelliptic-obstruction-tower
  - THM-3818-scaled-inert-cubeclass-support-two-pair-packet
related:
  - THM-3335-square-triangular-pell-markov-pythagorean-selector
  - THM-3793-inert-prime-sum-all-scale-two-cube-singleton
  - THM-3811-ramification-class-unit-criterion-and-nonlinear-cubic-packet
script: 04-computation/chebyshev_pell_berggren_conductor_staircase_thm3819.py
output: 05-knowledge/results/chebyshev_pell_berggren_conductor_staircase_thm3819.out
script_sha256: 3cc0a69a126711c325a03cf107b3bda2ea740018a88a1e4130ac821c6977f997
output_sha256: 7cb9237be31b773bb82620db608faeb4b7027e52be3684d2ab429d0fa7ba23b4
semantic_sha256: 7012073dd5917c62003333e0f600a84ec78a9255d00d851b3865db6c7155069f
hash_basis: raw LF bytes
---

# THM-3819 -- one Pell state, three lawful operations, and three losses

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  The proof below
is all-depth integer and polynomial algebra.  The exact companion checks its
conventions, finite controls, and consequence objects; an independent audit
rederives the maps and sharp loss boundaries.

No literature-priority or global-novelty claim is made.  This is a typed
composition of the proved repository mechanisms cited below.  It proves no instance of
the planar Jacobian conjecture or LRC(14).

## 1. Inheritance pass and notation

[THM-3757](THM-3757-pell-chebyshev-three-charge-hyperelliptic-obstruction-tower.md)
uses a polynomial Pell identity to build smooth planar noncoordinates.
[THM-3742](THM-3742-square-triangular-pell-mod13-central-sign-projective-cycle.md)
types the corresponding integral Pell orbit and signed clock modulo thirteen.
[THM-3756](THM-3756-odd-square-ordinal-berggren-affine-descent.md) turns the
Berggren tree into an affine odd-root forest, while
[THM-3745](THM-3745-monomial-plane-branch-conductor-triangular-pell-selector.md)
gives the triangular normalization defect and
[THM-3744](THM-3744-pell-prefix-loneliness-constant-carry-exact-formula.md)
gives the exact Pell-prefix loneliness packet.

The canonical near miss is the shared seed

```text
(x,q)=(3,1) in the Pell conic,
(x,d)=(3,1) in the odd-root forest.                         (1)
```

The two second coordinates have different meanings and already separate at
the next state.  The least-used sidecar is the **return time**, not another
square or triangular scalar.

Let

```text
P_0=0,       P_1=1,       P_(j+2)=2P_(j+1)+P_j             (2)
```

be the Pell numbers.  For `n>=1`, retain THM-3757's normalization

```text
tau=2z-1,
G_n(z)=C_n(tau)+(tau-1)U_(n-1)(tau),
H_n(z)=C_n(tau)+(tau+1)U_(n-1)(tau),                       (3)
```

where `C_n,U_n` are Chebyshev polynomials and

```text
z G_n(z)^2-(z-1)H_n(z)^2=1.                               (4)
```

Define

```text
g_n=G_n(2),                         h_n=H_n(2),
x_n=2g_n-h_n,                       q_n=(h_n-g_n)/2,
m_n=(x_n+1)/2,
a_n=3g_n-2h_n,                      b_n=g_n.               (5)
```

Here `h_n` is the Berggren return time below, not THM-3756's generic
centered-difference parameter.

## 2. After evaluation, the two integer states are losslessly equivalent

### Theorem 2.1 (exact evaluation and indexing)

For every `n>=1`, all quantities in `(5)` are positive integers.  Initially,

```text
(g_1,h_1,x_1,q_1,m_1,a_1,b_1)=(5,7,3,1,2,1,5).          (6)
```

Their exact indexing is

```text
a_n=P_(2n-1),
q_n=P_(2n)/2,
x_n=P_(2n)+P_(2n-1),
b_n=g_n=P_(2n+1),
h_n=P_(2n+1)+P_(2n).                                     (7)
```

They satisfy

```text
h_n^2-2g_n^2=-1,                   x_n^2-8q_n^2=1,       (8)
```

and

```text
[g_(n+1)]   [3 2][g_n]       [x_(n+1)]   [3 8][x_n]
[h_(n+1)] = [4 3][h_n],      [q_(n+1)] = [1 3][q_n],   (9)

m_(n+1)=m_n+h_n,
(a_(n+1),b_(n+1))=(b_n,6b_n-a_n).                       (10)
```

*Proof.*  At `z=2`, equation `(3)` first gives

```text
g_n=C_n(3)+2U_(n-1)(3),             h_n=C_n(3)+4U_(n-1)(3).
```

Thus the inverse definitions in `(5)` give
`x_n=C_n(3)` and `q_n=U_(n-1)(3)`.  The standard Chebyshev unit formula is

```text
(3+sqrt(8))^n=C_n(3)+sqrt(8)U_(n-1)(3).                 (11)
```

It identifies `(x_n,q_n)` with THM-3742's positive norm-one orbit and gives

```text
g_n=x_n+2q_n,                       h_n=x_n+4q_n.         (12)
```

The displayed relations prove integrality, positivity, and losslessness on
the evaluated state.  Evaluation of `(4)` gives `2g_n^2-h_n^2=1`.
Multiplication by `3+sqrt(8)` proves the second recurrence in `(9)`;
conjugating by `(12)` gives the first.  Standard Pell identities give `(7)`.
Finally,

```text
x_(n+1)=3x_n+8q_n=x_n+2(x_n+4q_n)=x_n+2h_n,            (13)
```

which proves the first formula in `(10)`; direct substitution proves the
Markov mutation.  QED

The first rows are

```text
n    m_n   q_n    x_n    g_n    h_n
1      2     1      3      5      7
2      9     6     17     29     41
3     50    35     99    169    239
4    289   204    577    985   1393.                    (14)
```

## 3. The negative-Pell coordinate is the Berggren return cocycle

Let `v_n=(m_n,1)` in THM-3756's ordinal chart.  Its odd-root pair is
`(2m_n-1,1)=(x_n,1)`.

### Theorem 3.1 (state-dependent return law)

For every `n>=1`,

```text
v_(n+1)=L^(h_n)v_n,
L^(h_n)(x_n,1)=(x_(n+1),1).                             (15)
```

*Proof.*  THM-3756 gives `L(r,1)=(r+1,1)` and, in odd-root
coordinates, `L(x,1)=(x+2,1)`.  Equations `(10)` and `(13)` prove `(15)`.
QED

The return times are

```text
h_n=7,41,239,1393,8119,... .                            (16)
```

This is an operation-compatible **cocycle**, not a fixed tree substitution.
The square-triangular vertices form a return set on the boundary ray, not a
Berggren subtree.  Their primitive triples are

```text
Xi(m_n,1)=(x_n,4q_n^2,4q_n^2+1).                       (17)
```

Indeed THM-3756 gives
`Xi(m,1)=(2m-1,2m(m-1),2m(m-1)+1)`, while `(8)` gives
`m_n(m_n-1)/2=q_n^2`.

## 4. The return time counts conductor-defect summands

For an abstract polynomial variable `X`, put

```text
D_m=direct_sum_(i=1)^(m-1) k[X]/(X^i).                 (18)
```

This is THM-3745's standardized normalization-quotient shape.

### Theorem 4.1 (operation-compatible staircase)

For every field `k`, every `n>=1`, and every monic polynomial of degree
`m_n`, THM-3745 gives

```text
delta_n=T_(m_n-1)=q_n^2,
conductor exponent=m_n-1=depth(v_n)+1.                 (19)
```

Moreover,

```text
D_(m_(n+1))
 =D_(m_n) direct_sum
  direct_sum_(i=m_n)^(m_n+h_n-1) k[X]/(X^i).           (20)
```

Thus one forest `L`-edge appends one cyclic summand and one Pell return
appends exactly `h_n` summands, of total length

```text
q_(n+1)^2-q_n^2
 =sum_(i=m_n)^(m_n+h_n-1)i
 =h_n(2m_n+h_n-1)/2.                                  (21)
```

*Proof.*  Equation `(19)` is THM-3745 with `(8)`.  The boundary word of
`(m_n,1)` is `L^(m_n-2)`, so its depth is `m_n-2`.  Equation `(20)` directly
compares the two standardized sums using `m_(n+1)-m_n=h_n`; summing their
lengths gives `(21)`.  QED

The word **standardized** is load-bearing.  Equations `(20)--(21)` compare
module shapes after renaming their polynomial variables.  They assert no
inclusion of branch rings or normalization quotients of different degrees.

## 5. The scalar loneliness observer is one Mobius map

For THM-3744's prefixes `S_N={P_1,...,P_N}`, put

```text
u_n^-=M(S_(4n-1)),                 u_n^+=M(S_(4n+1)).  (22)
```

### Theorem 5.1 (alternating scalar semiconjugacy)

For every `n>=1`,

```text
u_n^- =a_n/x_n=(3g_n-2h_n)/(2g_n-h_n),
u_n^+ =x_n/(2b_n)=(2g_n-h_n)/(2g_n).                   (23)
```

If `f(u)=1/(4-2u)`, then

```text
f(u_n^-)=u_n^+,                    f(u_n^+)=u_(n+1)^-. (24)
```

Hence one Pell/forest return is observed on the optimum scalar by `f^2`.

*Proof.*  Equation `(23)` is THM-3744's factorization together with `(5)`.
The identity `4x_n-2a_n=2b_n` proves the first half-step.  Equations
`a_(n+1)=b_n` and `x_(n+1)=2g_n+h_n` prove the second.  QED

The first rows are

```text
n       u_n^-       u_n^+
1        1/3         3/10
2        5/17       17/58
3       29/99       99/338.                            (25)
```

This is operation-compatible for the scalar maximum.  It does not transport
THM-3744's runner-labelled residue and carry packet.

## 6. Modulo thirteen the forest realization remains a cocycle

The coordinate change

```text
(x,q) |-> (g,h)=(x+2q,x+4q)                           (26)
```

has determinant two and is invertible modulo thirteen.  It conjugates
THM-3742's matrix `[3,8;1,3]` to `[3,2;4,3]`; hence `(g,h)` also has one
signed fourteen-cycle on the relevant nonzero norm fibre.  But the forest
section replaces `q` by the constant inner root one, and its update is

```text
x_(j+1)=x_j+2h_j.                                     (27)
```

Since `M^7=-I` and `M^14=I`, telescoping gives, from every norm-one state,

```text
sum_(i=0)^6 h_(j+i) =-x_j       (mod 13),
sum_(i=0)^13 h_(j+i)=0          (mod 13).              (28)
```

The companion verifies `(26)--(28)` on all fourteen states.  The half-period
displacement depends on the initial `x_j`, and the boundary `L`-translation
has order thirteen.  Thus the signed Pell `C14` is not an autonomous
Berggren `C14`.

## 7. Exact intersection with the incoming support-two address

The Markov factors obey `a_n+b_n=2x_n`.  Inside THM-3743's cap
`a+b<=356`, the values `x_1,x_2,x_3,x_4=3,17,99,577` prove that the complete
Pell intersection is

```text
(a,b)=(1,5),(5,29),(29,169).                          (29)
```

Only

```text
(5,29),              5+29=34=2*17                    (30)
```

satisfies THM-3793's inert cube-free primitive-sum hypothesis.  Its reversible
address is

```text
5^3+29^3=24514,                                       (31)
```

with exactly one positive distinct two-cube representation.  This gives
`78` labelled placements but preserves no other speed, owner, phase, arrival,
or loneliness predicate.

At this one merged row the complete arithmetic packet is

```text
(a,b)=(5,29),       (x,q)=(17,6),       m=9,
(g,h)=(29,41),      Xi(9,1)=(17,144,145),
delta=36,           (u^-,u^+)=(5/17,17/58),
cube address=24514.                                      (31a)
```

THM-3818 supplies one further finite residue reduction.  For a selected labelled copy of
the pair `(5,29)`, the primitive pair-sum modulus is `D=34`, and its safe set
on the pair-sum lattice is exactly

```text
G_(5,29)=Z/34Z minus {0,7,14,20,27}.                    (31b)
```

Indeed safety is `14 d_34(5j)>=34`, so the unsafe residues have
`5j in {0,+/-1,+/-2}`.  Since `5^(-1)=7 mod 34`, these are precisely the five
listed residues.  Thus the other eleven runners in any hypothetical
counterexample having this selected minimal pair must cover the remaining
`29` residue classes by THM-3818's danger sets.  Under a common global dilation
of the entire labelled row by `c`, the full residue schedule repeats `c` times.
If only the selected pair is scaled, only its safe-set pattern repeats; the
eleven danger sets still require the ambient residues modulo `34c`.  This is a
concrete finite covering obligation, not a row exclusion and not an off-lattice
statement.

## 8. Operation and loss ledger

| source | target and map | preserved predicate / operation | destroyed information | needed sidecar |
|---|---|---|---|---|
| THM-3757 `(G_n,H_n)` | evaluate at `z=2`, invert `(12)` | Pell state, norm, depth update | polynomial identity, roots, critical resultant, generic-fibre differential | full profiles and `psi_n` |
| Pell `(x_n,q_n)` | forest `(m_n,1)` | selected triple and return order | `q_n` is not the inner root; native edge time | `q_n` or cocycle `h_n` |
| boundary `(m,1)` | standardized `D_m` | degree, triangular length, conductor exponent, one-summand `L` update | `F`, multiplication, tangent slopes, actual cross-degree map | `F` and coefficient-ring identification |
| marked state plus channel | scalar `u_n^-/u_n^+` | exact optimum and Mobius update | speeds, residues, carries, owner, phase | prefix/carry packet and channel bit |
| inert marked row `(5,29)` | pair-sum lattice modulo `34` | exact 29-class safe set and cover obligation | off-lattice phases and eleven-runner arrival | eleven labelled ambient residues |
| Markov `(a_n,b_n)` | THM-3793 cube address | coprime ratio and mutation | placement and remaining LRC semantics | complete relation row |

No row supplies the next row's sidecars automatically.

## 9. Sharp hostile boundaries

### 9.1 Selector nonclosure and no fixed word

At `m=2`, `T_(m-1)=1`; one `L`-child has `m=3` and `T_2=3`, not a
square.  The first return is `L^7`.  More generally, starting on `s=1`, an
`M` or `R` child sets the second coordinate greater than one; later children
never restore one.  Every positive word preserving the boundary is therefore
`L^c`.  Return lengths `7` and `41` exclude any fixed word.

### 9.2 The second coordinates separate immediately

```text
M(3,1)=(17,6)               in Pell coordinates,
L^7(3,1)=(17,1)             in forest odd-root coordinates.       (32)
```

Six is a triangular square root; one is `sqrt(C-B)`.  They cannot be
identified.

### 9.3 Evaluation does not retain JC geometry

At depth one, `G_1=4z-3,H_1=4z-1`.  Replacing `G_1` by
`G_tilde=G_1+(z-2)=5z-5` retains the evaluation `(5,7)`, but

```text
z G_tilde^2-(z-1)H_1^2=(z-1)(9z^2-17z-1),             (33)
```

not one.  The evaluated state cannot certify THM-3757's absorbed roots,
smoothness, or generic-fibre cohomology.

### 9.4 Actual conductor rings do not append

For `F(T)=T^m`, let `A_m=k[t^m,t^(m+1)]`.  The consecutive marked rings

```text
A_9=k[t^9,t^10],                 A_50=k[t^50,t^51]     (34)
```

are incomparable.  Although `t^50` lies in `A_9`, `t^51` does not: an
equality `51=9a+10b` forces `b=6 mod 9`, whose smallest possibility gives
`60`.  Conversely `t^9` does not lie in `A_50`.  The standardized operation
does not lift to ring inclusion.

### 9.5 The scalar maximum is not the carry packet

```text
M({1,2})=M({1,2,5})=1/3.                               (35)
```

At phase `1/3`, the labelled distance packets have lengths two and three.
The scalar observer cannot reconstruct runner ownership or carry
multiplicity.

### 9.6 Constant Pell data cannot pay THM-3811's class debt

An element `p+q_1 omega+q_2 theta` whose coefficients are constants extracted
from `(g_n,h_n)` is inside THM-3811's affine-profile search.  Its exponent-one
norm equation is excluded by that theorem's exact Groebner gate; at exponent
at least two, degree zero violates its `D>=2e-1` floor.  This kills only the
constant-profile insertion, not THM-3811's open class problem.

## 10. Reproduction and scope

Run

```bash
python3 -B 04-computation/chebyshev_pell_berggren_conductor_staircase_thm3819.py
python3 -B -O 04-computation/chebyshev_pell_berggren_conductor_staircase_thm3819.py
```

Both modes byte-match the stored transcript.  The assertion-free companion
performs `6,982` gates:

- polynomial Chebyshev/Pell identities through depth `18`;
- integer Pell, parity, Pythagorean, forest-return, conductor-depth, and
  defect-summand invoices through depth `240`;
- both LRC factorizations and Mobius half-steps through depth `60`;
- the conjugacy, sharp orbit, central sign, and cocycle sums on all `14`
  norm-one states modulo thirteen;
- the complete support-two intersection, inert row, and two-cube fibre; and
- every computational hostile in Sections 9.1--9.5.

The semantic digest is
`7012073dd5917c62003333e0f600a84ec78a9255d00d851b3865db6c7155069f`.
The checks audit constants and consequence objects; the all-depth quantifiers
follow from `(11)--(13)`, the affine `L` law, the direct-sum decomposition,
and the displayed fraction identities.

The proved gain is a lossless coordinate chart on the evaluated Pell state, a
state-dependent Berggren return cocycle, an operation-compatible standardized
conductor staircase, and a scalar LRC semiconjugacy.  The planar Jacobian
conjecture, LRC(14), actual conductor-ring transport, and LRC carry/owner
transport remain open.  The independent audit confirms that evaluation itself
still loses the polynomial geometry exhibited in Section 9.3.
