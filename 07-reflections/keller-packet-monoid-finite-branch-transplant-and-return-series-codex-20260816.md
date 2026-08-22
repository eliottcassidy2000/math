# The finite survivor defines a transplant monoid, but complete packets never return

**Status: PROVED STRUCTURAL COROLLARY OF INDEPENDENTLY AUDITED THM-3528;
VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  The fixed Keller packet family now has the monoid structure suggested
at the start of the session.  Complete packets multiply by adding their two
grades, and the cleared cubic norm is multiplicative.  The resulting operator
shifts any returned old-boundary factor through the canonical ancestry

```text
L -> H -> J -> G -> R_5 -> R_6 -> R_7 -> R_8 -> ... . (1)
```

This is the exact branch transplant.  THM-3529 subsequently proves that its
positive-defect antecedent cannot occur for any complete packet: the finite
survivor divisor is `V(F^*L)`, and the x-free minimum-beta face excludes it.
Thus the transplant is a sharp conditional monoid law and hostile mechanism,
while the fixed positive-level raw orbit has no old-`L` return at all.

## 1. Why the local `1+2` anatomy is load-bearing

At the generic old boundary `L=0`, the cubic inverse fibre decomposes into

```text
one regular finite survivor + two divergent conjugate branches.           (2)
```

For a packet `A(e,m)`, the two divergent values have exact orders `-e/2` and
`-e/2`.  Their sum fixes the pole order `-e`.  The finite survivor contributes
the only variable term

```text
s(P)=v_L(P(q_fin))>=0.                                (3)
```

Thus

```text
v_L(N(P))=-e+s(P),
ord_L(L^eN(P))=s(P).                                  (4)
```

The folded pair controls clearing; the distinguished singleton controls
return.  This is an exact use of the `1+2` cubic anatomy, unlike the earlier
speculative comparison with a ternary combinatorial tree.

The old finite-sheet computations showed `s=0` through input `R_7` and output
`R_8`.  THM-3528 proves that a hypothetical later `s>0` would not stop the
packet tower but would instead mark a divisor return.  THM-3529 now proves
the stronger statement `s(P)=0` for every complete packet, so every
positive-level raw output is old-`L`-coprime.

## 2. The complete-packet monoid

Let `[P]` denote a nonzero complete-packet polynomial modulo multiplication
by `Q^*`.  Products obey

```text
A(e,m)A(f,n)=A(e+f,m+n).                              (5)
```

Every one of the five complete initial forms multiplies to the stated output
form, so this is not merely an addition law for two exposed weights.  Define

```text
T(P)=L^eN(P),       grade(P)=(e,m).                   (6)
```

Then

```text
T(PQ)=T(P)T(Q),
grade(T(P))=[[7,-2],[3,-2]] grade(P).                 (7)
```

Constants satisfy `T(cP)=c^3T(P)`, so (7) is literal after quotienting by
nonzero rational scalars.  The raw orbit is `P_n=T^n(L)`.

This gives three distinct gradings:

```text
inverse fibre degree:       3^n,
packet/Pell norm grade:     (-8)^n,
factor-ancestry shift:      t^n.                      (8)
```

They share an index but retain different information.  Neither the cubic
degree nor the Pell norm records divisor returns.

## 3. Exact transplant law

Suppose the first old-boundary factor at some raw rung has multiplicity `s`:

```text
P_n=L^sR,       s>0.                                  (9)
```

Dividing the five complete faces by those of `L^s` gives

```text
R has A(e_n-s,m_n),       s<=e_n-m_n.                 (10)
```

The bound comes from the minimum-beta factor `z^(e-m)` and is sharp at the
packet level.  Multiplicativity now gives, for every `k>=0`,

```text
P_(n+k)=P_k^s T^k(R).                                 (11)
```

So one finite-branch return draws a full diagonal through the factor table:

```text
level n      contains L^s,
level n+1    contains H^s,
level n+2    contains J^s,
level n+3    contains G^s,
...                                                               (12)
```

Equation (11) is stronger than a recurrence coincidence: it is an exact
polynomial divisibility statement.  It is also one-way.  If
`s_j=s(P_j)>0`, then the factor occurs in the **next** rung
`P_(j+1)=L^sR`; for every genuine later output, `m_(j+1)>0`, so `R` is
nonconstant and reducibility propagates.  The seed hostile `P_0=L=L*1`
shows why positive multiplicity alone does not imply reducibility when the
cofactor is a unit.  Defect zero by itself does not prove irreducibility,
squarefreeness, or image status.

## 4. A subset of the naturals is a return word, not a harmonic scalar

The complete finite-sheet ledger is

```text
D(t)=sum_(n>=0) s(P_n)t^n in N[[t]].                  (13)
```

For the fixed raw orbit, THM-3529 gives `D(t)=0`.  In a general return-word
formalism its support is a subset of `N`, multiplicities are retained, and
transplant is the unilateral shift `D(t)->tD(t)`.  This is the faithful
realization of the “subset of the harmonic series” intuition: retain the
indexed word.

Scalar evaluation at reciprocals is not faithful.  The exact `1..13` census
has

```text
8192 subsets -> 3712 sums,
2944 collision values, maximum multiplicity 3,                        (14)
```

and the three-way fibre

```text
1/2=1/3+1/6=1/4+1/6+1/12.                           (15)
```

Thus a harmonic sum can be a statistic of the return word, but cannot be its
address.  A formal power series or indicator sequence is the lawful carrier.

## 5. Ternary tree, Fibonacci-style recurrence, and Berggren ancestry

The inverse fibre still branches cubically away from the Jelonek spine, and
the packet coordinates still obey

```text
u_(n+2)=5u_(n+1)+8u_n.                                (16)
```

These statements live on different objects.  The ternary tree counts inverse
points; (16) evolves Newton-face grades; (11) evolves divisor factors.  The
finite survivor in (2), rather than one of the two escaping branches, is what
couples the first and third objects.

The attachment's distinguished Berggren `U`-spine is parabolic: its matrix is
unipotent, its coordinates are quadratic in depth, and its scalar tail has
constant second difference eight.  The Keller packet ray is hyperbolic in
`Q(sqrt(57))`; its induced primitive-triple map is a Lorentz similitude of
multiplier `64` and determinant `-512`, not a Berggren isometry.  Both inhabit
the primitive-Pythagorean universe, but their ancestry words are different:

```text
Berggren U-spine:       one unimodular tree edge repeated;
Keller packet ray:      one nonunimodular Pell jump repeated;
Keller factor return:   one monoid factor shifted by T.                (17)
```

Calling all three “the same Fibonacci tree” would discard exactly the
operation that distinguishes them.

## 6. The four-object tournament window

The first four packet objects

```text
L,H,J,G                                                    (18)
```

have a total chronological order under iteration.  Orienting every earlier
object toward every later one gives a transitive tournament `T4` with six
pair comparisons.  This is a legitimate tournament because reachability is
intrinsic after the operator `T` is fixed.

It is not the LRC pointed-six carrier.  That carrier is a bidirected `P4` with
pair census `(3,0,3)`, while (18) has `(0,6,0)`.  The XOR sidecar shows how a
chosen transitive `T4` gauge can turn the P4 edge set into a reversal mask,
but the conversion loses amplitudes, roots, return multiplicities, and the
`F13` seam.  The common count six is therefore a useful coordinate atlas,
not a current or flux map.

## Connection contract

| field | exact answer |
|---|---|
| source | complete five-face packet polynomials for the fixed inverse chart |
| operation | `T(P)=L^eN(P)` |
| monoid law | products add `(e,m)` and `T(PQ)=T(P)T(Q)` |
| return observable | `s(P)=v_L(P(q_fin))` |
| transplant | `L^s|P_n` implies `P_k^s|P_(n+k)` |
| preserved | packet completeness, factor multiplicity, iteration offset |
| separate grades | fibre `3^n`, Pell `(-8)^n`, return shift `t^n` |
| faithful subset carrier | formal defect series in `N[[t]]` |
| lossy carrier | scalar harmonic subseries |
| raw finite-sheet ledger | identically zero by THM-3529 |
| still open at this reflection state | factorization, images, separability, components, general JC |

## Reproduction

```text
python -B 04-computation/keller_packet_monoid_branch_transplant_audit_20260816.py
python -B -O 04-computation/keller_packet_monoid_branch_transplant_audit_20260816.py
```

The synthetic defects in the companion are hostile arithmetic controls, not
geometrically realizable defects of complete packets after THM-3529.  Normal
and optimized outputs match.  Script/output/
semantic LF SHA-256 are
`8256c7179c415e8588a0612608f7c253baf026d8a237cd8cbf01720758e8b5dc`,
`b934628ab80fefe5fdc6662d10f12bffa4f6f1fb25a3006b9254f3f4b7d204d6`,
and `f004cd7643933e81a2fbce73a3df6d72c9e4943702cca7ee17b32d6690be3874`.
