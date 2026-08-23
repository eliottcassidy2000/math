---
id: THM-3716
title: "Monomial Broughton Hamiltonian obstruction family"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  Over every
  characteristic-zero field, for all
  m>=2 and n>=1 the nonsingular noncoordinate polynomial
  Q=x+x^m y^n has unimodular gradient, but the canonical complementary
  volume-form class is not in the image of J(-,Q).  Equivalently Q has no
  polynomial Jacobian mate.  The obstruction is an isolated two-edge
  monomial chain with infinitely many forced nonzero coefficients.  This is
  a counterexample-search no-go family, not a proof of JC(2).
source: root / 2026-08-22
audit: >
  PASS.  An independent derivation checked the Bezout complement, curl sign,
  unique low/high predecessors, transport coefficient and recurrence, the
  characteristic-zero boundary, the principal syzygy argument, and the
  reducible-fibre noncoordinate conclusion.  Normal, optimized, and frozen
  exact transcripts agree with all three recorded hashes.
depends_on: []
related:
  - THM-3709-cohn-alternating-two-by-two-elementary-decoration-nonentry
script: 04-computation/jc2_monomial_broughton_hamiltonian_obstruction_thm3716.py
output: 05-knowledge/results/jc2_monomial_broughton_hamiltonian_obstruction_thm3716.out
script_sha256: 6d3819f25dbb73a194ad579814610ed9fd0997460f98a936148db086bf9c5950
output_sha256: 89470553e8bf0e3972e20b4ddfc72aefcd7081a51e94bca0191c764359435d27
semantic_sha256: 15230359a117c8d77dbcbbc6b4a251a600e85b742f7331b6c86ec18baa3da688
hash_basis: raw LF bytes
---

# THM-3716 -- every monomial Broughton component has nonzero Hamiltonian debt

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  This turns the
Hamiltonian-cokernel lens into an explicit infinite family.  It includes the
normalized `X+X^2T` resonance that closes the last constant cell in THM-3709.

Let `k` be a characteristic-zero field, let `m>=2`, `n>=1`, and put

```text
Q=x+x^m y^n,                 z=x^(m-1)y^n.              (1)
```

Use `J(F,G)=F_xG_y-F_yG_x` and
`curl(P dx+R dy)=P_y-R_x`.

## 1. A nonsingular noncoordinate component

The derivatives are

```text
Q_x=1+mz,                    Q_y=n x^m y^(n-1).          (2)
```

They satisfy the polynomial Bezout identity

```text
(1-mz)Q_x +(m^2/n)x^(m-2)y^(n+1)Q_y=1.                (3)
```

Thus `Q` has no critical point over any field extension.  It is nevertheless
not a coordinate: its zero fibre is the reducible curve

```text
Q=x(1+x^(m-1)y^n)=0.                                  (4)
```

Indeed a coordinate has prime principal ideal, whereas `(4)` has two
nonunit factors.

Define the complementary polynomial one-form

```text
alpha=(m^2/n)x^(m-2)y^(n+1) dx +(mz-1)dy.             (5)
```

Equation `(3)` says exactly

```text
alpha wedge dQ=dx wedge dy.                            (6)
```

Its curl is the nonzero monomial

```text
tau=curl(alpha)
   =[m(m+n)/n]x^(m-2)y^n.                              (7)
```

For a prospective Jacobian mate, a Hamiltonian shear `alpha+f dQ` is closed
precisely when

```text
J(f,Q)=tau.                                            (8)
```

We prove that `(8)` has no polynomial solution.

## 2. The Hamiltonian derivation is a two-edge lattice walk

For every monomial `x^i y^j`, direct differentiation gives

```text
J(x^i y^j,Q)
 =-j x^i y^(j-1)
  +(ni-mj)x^(i+m-1)y^(j+n-1).                         (9)
```

The two output exponents differ by the fixed vector

```text
v=(m-1,n).                                             (10)
```

Consequently the target monomial in `(7)` belongs to an isolated bidiagonal
chain.  Its only nonnegative low-edge preimage is

```text
a_0=(m-2,n+1).                                         (11)
```

There is no high-edge preimage, because that would have `x`-exponent `-1`.
After the first forced coefficient, cancellation of the high edge forces
successive input exponents

```text
a_k=(i_k,j_k)
    =(m-2+k(m-1),(k+1)n+1),             k>=0.          (12)
```

No monomial outside this chain can contribute to one of its output sites:
the low and high preimages are unique.

## 3. The chain never terminates

Write `c_k` for the coefficient of `x^i_k y^j_k` in a hypothetical `f`.
At the first output site, `(7)` and the low edge in `(9)` force

```text
c_0=-m(m+n)/[n(n+1)] !=0.                              (13)
```

The high-edge multiplier at `a_k` is

```text
n i_k-m j_k=-(m+(k+2)n),                               (14)
```

which is nonzero in characteristic zero.  At every later output site, the
previous high edge and the current low edge must cancel, so

```text
c_k=-[m+(k+1)n]/[(k+1)n+1] c_(k-1),       k>=1.        (15)
```

Every numerator and denominator in `(15)` is a positive integer.  Hence all
`c_k` are nonzero.  Equation `(8)` would require infinitely many monomials,
contradicting `f in k[x,y]`.  Therefore `(8)` has no polynomial solution.

If a polynomial `P` satisfied `J(P,Q)=1`, then `dP` would be a closed
complement obeying `(6)`.  The difference `dP-alpha` is necessarily a
multiple `f dQ`: from

```text
(dP-alpha) wedge dQ=0                                  (16)
```

and the Bezout identity `(3)`, the syzygy module of `(Q_y,-Q_x)` is generated
by `(Q_x,Q_y)`.  Closure would then give `(8)`, a contradiction.  Thus `Q`
has no polynomial Jacobian mate.

## 4. Meaning, reproduction, and boundary

The theorem separates two properties that a counterexample search must not
conflate:

```text
gradient ideal (Q_x,Q_y)=k[x,y]       -- no affine critical point,
[tau] in coker J(-,Q) is nonzero      -- no volume-preserving mate. (17)
```

For `m=2,n=1`, `(7)--(8)` become exactly the `6y` Broughton obstruction in
THM-3709.  Formula `(15)` is the all-parameter version of Cohn's factorial
repair tail: the local Bezout identity opens the path, but no coefficient is
zero at a finite seam.

Reproduce the exact identities with

```bash
python3 -B 04-computation/jc2_monomial_broughton_hamiltonian_obstruction_thm3716.py
python3 -B -O 04-computation/jc2_monomial_broughton_hamiltonian_obstruction_thm3716.py
```

The assertion-free companion checks `(1)--(15)` for `m=2..9`, `n=1..8`,
all monomials in a `7 by 7` hostile bank, and twelve forced tail layers.
The finite grid guards the formulas; the lattice-chain proof is universal.
This theorem does not classify all nonsingular polynomials, all Hamiltonian
cokernels, Keller pairs, or JC(2).  **QED.**
