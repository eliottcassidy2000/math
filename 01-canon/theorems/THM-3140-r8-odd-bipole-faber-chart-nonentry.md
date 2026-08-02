---
id: THM-3140
title: "R=8 explicit odd-bipole Faber-chart nonentry"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  The explicit
  d=11 odd-bipole response realizing THM-3133's first sharp abstract boundary
  cannot enter a normalized R=8 nonsplit polynomial exact-square-prefix Faber
  source chart.  Divisor support forces A_src=c(x^2-1)^6 and
  (x^2-1)^7|B_src; the exact infinity fan then forces s_F constant and d_F of
  degree at most two, making ord_infinity(K_Q)>=2, whereas the response
  identity requires ord_infinity(K_Q)=-1.  THM-3151 subsequently supersedes
  this scope boundary by excluding every response in the (11,11) cell and
  all later resonant equality cells; neither theorem proves the planar
  Jacobian conjecture.
audit: >
  An independent hostile derivation rebuilt the Faber rows, divisor argument,
  R=8 pole orders, all three infinity chambers, zero-polynomial boundaries,
  wall resultant 65536, and final K-order contradiction.  Its separate script
  passed 1,554 gates and byte-matched under normal and optimized Python.  The
  canonical companion independently agrees by direct multinomial extraction
  and a logarithmic-derivative recurrence; normal, optimized, and stored
  transcripts match exactly.
source: root/frontier-synthesis-2026-08-02
depends_on:
  - THM-2827-uniform-pole-order-faber-nonresonance-atlas-and-double-pole-exclusion
  - THM-3133-common-simple-zero-faber-exclusion-and-odd-bipole-boundary
related:
  - MISTAKE-317
  - THM-3151-resonant-odd-bipole-equality-cell-nonentry-and-degree-floor
script: 04-computation/jc_r8_odd_bipole_faber_chart_nonentry_thm3140.py
output: 05-knowledge/results/jc_r8_odd_bipole_faber_chart_nonentry_thm3140.out
script_sha256: 077b1f00db153af66e2666fb22e1ff448cf423756588a712fb1f0dfc6e5827da
output_sha256: 67728a70fa4631202ef59af50965a2222725e8950dfed106fe3ebc507b5c07c3
hash_basis: LF-normalized bytes
---

# THM-3140 -- R=8 explicit odd-bipole Faber-chart nonentry

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-3133 proves that the first possible balanced response cell in the
normalized Faber bank is

```text
R=8,                 D=11,                 N=22,
pole passport=(11,11).                                      (1)
```

It also constructs an explicit odd-bipole response attaining this abstract
boundary.  The present theorem shows that this particular response does not
lift through the source/Faber equations.  It does not classify all responses
with passport `(11,11)`.

## 1. The target response and chart equations

Put

```text
T_p=x^2-1,
E=E_11=x^11-(11/2)x^9+(99/8)x^7-(231/16)x^5
       +(1155/128)x^3-(693/256)x,                            (2)
V=T_p^13,             G=E/T_p^12,
F=E^2/T_p^11,         M=E T_p.                              (3)
```

The THM-3133 first integral is

```text
2T_p E'-22xE=693/128.                                       (4)
```

It implies that `E` is squarefree and coprime to `T_p`; consequently `M` is
squarefree.  Suppose, for contradiction, that `(3)` enters a normalized
nonsplit polynomial exact-square-prefix Faber chart

```text
Q=E_30+sum_(j=1)^6 a_j E_(4j-2).                            (5)
```

Use the source coordinates

```text
T=A_src^2/V,
d_F=C_0-B_src^2/(4V),
s_F=A_src B_src/(2V)-E_0.                                  (6)
```

The exact flux equations and the THM-2827 polynomial-ring sidecar are

```text
Phi_Q=0,              Psi_Q in C,
K_Q=T H_Q,            A_src K_Q=lambda M,       lambda!=0,
H_Q in Q[T,s_F,Td_F].                                      (7)
```

## 2. Divisor support fixes the source orders

Let `alpha` be a zero of `A_src` away from `V`, of order `a>=1`.  All three
arguments of `H_Q` in `(7)` are regular there and

```text
ord_alpha(T)=2a,       ord_alpha(K_Q)>=2a,
ord_alpha(A_src K_Q)>=3a>1.                                 (8)
```

This contradicts the squarefreeness of `M` and the last identity in `(7)`.
Thus every zero of `A_src` lies over `x=1` or `x=-1`.

At either of these two points

```text
ord(V)=13,                         ord(M)=1.                  (9)
```

The exact `R=3k+2` resonance formula of THM-2827, with `R=8`, `k=2`, and
`delta=1`, forces

```text
ord(A_src)=1+(2k+1)delta=6.                                (10)
```

The accompanying pure-`q` inequality gives `ord(B_src)>=7`; order six is
the planted hostile and remains polar.  There are no other zeros of
`A_src`, so for a nonzero constant `c`,

```text
A_src=cT_p^6,                  T=c^2/T_p,
T_p^7 | B_src.                                              (11)
```

It follows directly from `(6)` that `d_F` and `s_F` are polynomials.

## 3. The exact infinity fan

Normalize the valuation at infinity by `ord_infinity(x)=-1`.  Equation
`(11)` gives

```text
ord_infinity(q)=1,             ord_infinity(T)=2,            (12)
```

where `q^2=T`.  For nonzero `d_F,s_F`, write
`D=deg(d_F)` and `S=deg(s_F)`, with leading coefficients `d_0,s_0`; zero
polynomials are treated separately by deleting their monomials.  Exact
coefficient extraction from

```text
(1+2d_F t^2+q t^3+(d_F^2-s_F)t^4)^(j-1/2)                  (13)
```

gives three exhaustive chambers for the top row `j=8`.

### Chamber I: `D<2S+2`

The unique top faces are

```text
in(Phi_8)=-(6435/1024)s_0^7,
in(Psi_8)= (6435/8192)s_0^8.                                (14)
```

Their orders are `-7S` and `-8S`.  Every retained lower row `j<=6` has
orders `-(j-1)S` and `-jS`, respectively.  When `S>=1`, the top row is
strictly lower and cannot be cancelled.

### Chamber II: `D=2S+2`

After removal of common nonzero factors and the substitution

```text
z=d_0q_0^2/s_0^2,                                           (15)
```

simultaneous cancellation of the two top faces would require

```text
f(z)=(z-1)(z^2-6z+1)=0,
g(z)=z^4-28z^3+70z^2-28z+1=0.                              (16)
```

But

```text
Res_z(f,g)=65536!=0.                                        (17)
```

At least one top flux therefore remains uncancelled, again strictly below
all retained lower rows.

### Chamber III: `D>2S+2`

The unique `Psi_8` face is

```text
(6435/8192)q_0^8d_0^4                                      (18)
```

of order `8-4D`, strictly below every `Psi_j`, `j<=6`.

The three chambers show that `s_F` cannot be nonconstant.  Thus `s_F` is a
constant, possibly zero.  With `S=0`, chamber III also excludes every
nonzero `d_F` of degree at least three.  Hence

```text
s_F is constant (possibly zero),
d_F=0 or deg(d_F)<=2.                                       (19)
```

The exact `s_F=0`, `d_F=0`, and `deg(d_F)=0,1,2,3` boundaries are included
in the independent hostile audit; no convention `deg(0)=0` is being used.

## 4. Final order contradiction

By `(12)` and `(19)`, each of `T`, `s_F`, and `Td_F` is regular at infinity.
The ring statement in `(7)` therefore gives

```text
ord_infinity(H_Q)>=0,        ord_infinity(K_Q)>=2.           (20)
```

The alternative `H_Q=0` is already incompatible with the nonzero right side
of `(7)`.  But `(3)`, `(7)`, and `(11)` require

```text
K_Q=lambda M/A_src=(lambda/c)E/T_p^5,
ord_infinity(K_Q)=10-11=-1.                                 (21)
```

Equations `(20)` and `(21)` contradict one another.  The explicit odd-bipole
response `(2)--(3)` cannot enter the chart `(5)`.  QED.

## 5. Exact controls and independent audit

The canonical companion builds rows `j=1,...,8` twice: directly by
multinomial extraction and independently by a logarithmic-derivative
recurrence.  It checks `H_j in Q[T,s_F,Td_F]`, every chamber and retained
lower row through `D,S<=80`, the wall factorization and resultant, the
odd-bipole ODE and clean factors, both divisor gates, and the final orders.
Run

```text
python3 04-computation/jc_r8_odd_bipole_faber_chart_nonentry_thm3140.py
python3 -O 04-computation/jc_r8_odd_bipole_faber_chart_nonentry_thm3140.py
```

and compare byte-for-byte with the declared output.

The independent hostile audit does not import the canonical companion.  It
passes `1,554` gates, including the zero-polynomial cases, under both normal
and optimized Python.  Its script and output hashes are, respectively,

```text
e1c973978e6e270c2916b62bd0869700dd069822a90b6a2c3bafafa760e850db
1af56f4c6ee0fcd5caf756acdab6960890a7d81cb5bce418b8883238295282fe
```

## 6. Connection contract and scope

The source is THM-3133's explicit sharp response, the target is the
normalized `R=8` source chart, and the map is source lifting followed by the
divisor and infinity fans.  It preserves the response divisor identity and
flux cancellation.  Valuation alone forgets coefficients, the missing
`R-1` row, and zero-polynomial boundaries; the exact wall resultant,
polynomial-ring sidecar, and squarefree response are load-bearing.

This theorem itself excludes one canonical realization of the first abstract
passport `(11,11)`.  THM-3151 later excludes every response with that
passport by a response-independent infinity fan.  Neither result excludes
another Faber branch or arbitrary chart entry, constructs a Keller map, or
proves `JC(2)`.

**End of proof.**
