# One dihedral presentation, three different dynamics

**Status: VERIFIED-EXACT for the Berggren/Farey and Fibonacci matrix
factorizations and for the `C4/C6` completion census; VERIFIED-NUMERICAL only
for the Keller old-`L` reflection sections through the audited depth-five
window.  No identification of the three actions, current, or LRC conclusion.**

## 1. Inheritance pass

THM-3334 makes the consecutive-parameter Berggren `U`-spine a parabolic
Farey ray and proves its quadratic `Q`-sequence.  THM-3509 identifies all
reduced fractions with primitive four-term Fibonacci windows and locates the
Cassini sign on an oriented two-edge sidecar of `K4`.  THM-3536 gives the
signed `C4` angle-language class and its two missing antipodals.  THM-3537
and its numerical sequel expose the old-`L` binary fold inside the fixed
ternary inverse tower.

The common operation had not been isolated: in all three lanes, the moving
element is a product of two reflections.  This is an exact common **source
presentation**, not a target-preserving map between the resulting objects.

## 2. The universal reflection--rotation grammar

Let

```text
D_infinity=<r,s | r^2=s^2=1>,       c=sr.              (1)
```

Then

```text
r c r=c^(-1).                                             (2)
```

Every pair of involutions therefore supplies a representation of the
infinite dihedral group, with `c` as its rotation.  The image of `c` can be
finite, parabolic, or hyperbolic; equation `(2)` alone does not choose among
those possibilities.

## 3. Berggren's `U` is an exact product of integral Lorentz reflections

Use the Pythagorean Lorentz form `J=diag(1,1,-1)` and THM-3334's matrix

```text
U=[1 -2 2; 2 -1 2; 2 -2 3].                            (3)
```

Put

```text
R_B=diag(-1,1,1),
S_B=U R_B=[-1 -2 2; -2 -1 2; -2 -2 3].                 (4)
```

Direct multiplication gives

```text
R_B^2=S_B^2=I,       S_B R_B=U,
R_B U R_B=U^(-1),
R_B^T J R_B=S_B^T J S_B=U^T J U=J.                     (5)
```

This realizes `(1)` inside `O(2,1;Z)`.  It also explains the word
“parabolic” without metaphor:

```text
(U-I)^3=0,       (U-I)^2!=0.                            (6)
```

On the invariant parabola

```text
P(a)=(a,(a^2-1)/2,(a^2+1)/2),
```

`R_B` reverses `a` and `U` translates `a` by two.  Starting at `(3,4,5)`,
the exact scalar heights are

```text
Q_n=2c_n+1=(2n+3)^2+2
   =11,27,51,83,123,171,227,... .                       (7)
```

The reduced Euclid parameters form the Farey ray

```text
1/2,2/3,3/4,4/5,... .                                   (8)
```

Indeed the parameter action

```text
T=[0 1; -1 2],       T(k,k+1)^T=(k+1,k+2)^T            (9)
```

has its own integral factorization

```text
R_T=[0 1;1 0],       S_T=[1 0;2 -1],
R_T^2=S_T^2=I,       S_T R_T=T.                        (10)
```

Consecutive fractions in `(8)` have cross-determinant one.  Thus this is a
literal parabolic ray in the Farey graph, not only a list of rational labels.

## 4. The Fibonacci ray is the hyperbolic representation of the same grammar

Let

```text
Q=[1 1;1 0],       M=Q^2=[2 1;1 1].                    (11)
```

Define

```text
R_F=[0 -1;1 0],       S_F=[-1 2;-1 1].                 (12)
```

In `GL_2(Z)` one has

```text
R_F^2=S_F^2=-I,       S_F R_F=M.                       (13)
```

Hence their images are involutions in `PGL_2(Z)` and again represent
`D_infinity`.  This time the rotation is hyperbolic:

```text
det M=1,       tr M=3,       chi_M(t)=t^2-3t+1.         (14)
```

Its primitive orbit is

```text
M^n(1,0)^T=(F_(2n+1),F_(2n))^T,                        (15)
```

so the fractions in `(0,1)` are

```text
1/2,3/5,8/13,21/34,55/89,... .                         (16)
```

Cassini's identity says consecutive terms of `(16)` also have cross-
determinant one.  Equations `(8)` and `(16)` are therefore two rays in the
same Farey graph beginning at `1/2`: the Berggren spine has linear
denominators and parabolic growth, while the Fibonacci ray has exponential
denominators and hyperbolic growth.  Sharing a root fraction does not merge
their ancestry.

## 5. The Keller old-`L` rotation is a finite quotient, presently numerical

The verified numerical old-`L` atlas has

```text
g_(r+1)=(A_r,B_r,id)(01),
A_r^2=B_r^2=id,       C_r=B_r A_r.                     (17)
```

Thus every observed depth also receives a dihedral representation, now in a
finite symmetric group on `3^r` sheets.  The observed rotation orders are

```text
ord(C_1),...,ord(C_4)=2,4,12,36.                        (18)
```

The global old-`L` permutation doubles the `C_r` cycle lengths and appends
`3^r` fixed leaves.  This is the exact mechanism within the accepted
numerical window.  Unlike `(5)` and `(13)`, an all-level algebraic
factorization of the local inverse-chart monodromy has not been proved.

The useful common diagram is therefore

| representation of `(1)` | rotation image | dynamical type | retained data |
|---|---|---|---|
| `O(2,1;Z)` | Berggren `U` | parabolic, quadratic | Lorentz form, Euclid parameters |
| `PGL_2(Z)` | Fibonacci `Q^2` | hyperbolic, exponential | reduced fraction, Cassini orientation |
| `Sym(3^r)` | Keller `C_r` | finite periodic | ancestry blocks, inertia cycles |

Only the presentation and reversal law are shared.  Metric growth, arithmetic
content, branch ancestry, and coefficient fields are destroyed by passing to
the abstract dihedral source.

## 6. Why sizes four and six still do not give tournaments

A directed rotation orbit of size `m` supplies only its `m` successor pairs.
The exact completion census is

```text
m=4: 4 observed pairs, 2 missing pairs,   4 tournament completions;
m=6: 6 observed pairs, 9 missing pairs, 512 tournament completions. (19)
```

None of these completions is invariant under the cyclic rotation.  More
generally, a rotation-invariant tournament on an even number of vertices
would have constant outdegree `(m-1)/2`, which is impossible.  The
reflection colours in `(1)` also disappear after retaining only the directed
cycle.

Thus a `C4` or `C6` orbit is an honest cyclic carrier with missing edges, not
“essentially” a tournament.  THM-3536's nonzero signed-`C4` class requires an
edge cochain; THM-3509's Cassini sign requires an oriented antipodal fibre.
Neither sidecar is selected by the bare dihedral presentation.

## 7. Subsets of the harmonic series

Both rays define subsets of reduced fractions and, through their denominator
sets, subsets of the natural numbers:

```text
D_B={2,3,4,...},             D_F={F_(2n+1):n>=1}.       (20)
```

Consequently

```text
sum_(d in D_B) 1/d diverges,
sum_(d in D_F) 1/d converges.                           (21)
```

The first is a cofinite harmonic tail; the second converges because Fibonacci
denominators grow exponentially.  This is a concrete warning about “a subset
of the harmonic series”: the Boolean support is primary, while its reciprocal
mass can have completely different analytic behaviour even for two canonical
Farey rays generated by the same reflection--rotation grammar.

## Connection contract

| field | exact answer |
|---|---|
| common source | infinite dihedral presentation `(1)` |
| Berggren map | `r -> R_B`, `s -> S_B`, `c -> U` |
| Fibonacci map | projective `r -> R_F`, `s -> S_F`, `c -> Q^2` |
| Keller map | numerical `r -> A_r`, `s -> B_r`, `c -> C_r` through `r<=4` |
| preserved | product-of-reflections law and reversal `(2)` |
| destroyed | parabolic/hyperbolic/finite type, metric growth, ancestry, amplitudes, coordinate fields |
| tournament loss | two missing pairs at size four; nine at size six |
| harmonic loss | support phase and ancestry are not recovered from convergence or density |

Reproduce the exact matrix, Farey-neighbour, Fibonacci, and tournament rows
with

```text
python -B 04-computation/berggren_fibonacci_dihedral_reflection_cospan_20260816.py
python -B -O 04-computation/berggren_fibonacci_dihedral_reflection_cospan_20260816.py
```

The Keller leg is intentionally not imported into that exact companion.  It
retains its separate VERIFIED-NUMERICAL status and forward-residual gates.
