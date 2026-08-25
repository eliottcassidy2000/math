---
id: THM-4091
title: "Integral coordinate-change boundary for LCM and coefficientwise depth"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED. Integral formal coordinate
  changes preserve every finite LCM-cleared denominator staircase. They also
  preserve the literal coefficientwise depth-one lattice, but for every depth
  e>=2 universal preservation already fails at output degree three. An exact
  scaled composition-matrix criterion characterizes the coordinate changes
  that preserve a fixed higher-depth lattice. This is an arithmetic transport
  theorem only; it supplies no decay, nonvanishing, or irrationality result.
source: codex-padic-zeta-irrationality-20260825 / denominator-transport niche
audit: >
  PASS. The symbolic proof was independently checked line by line. The exact
  companion exhausts all 625 integral coordinate changes with four
  coefficients in [-2,2] through output degree six: 13,125 derivative/depth-one
  gates, 52,500 LCM-matrix gates, and 7,500 degree-one/two minimality gates.
  It finds the first possible depth-two failure at (n,k)=(3,2) and replays the
  sharp hostile for e=2,...,8. Normal and optimized outputs byte-match; the
  script contains no assert statements or floating literals.
depends_on: []
related:
  - THM-4056-divisor-phase-compiler-and-duffin-schaeffer-lcm-clock
  - THM-4057-stern-brocot-depth-pullback-and-rational-edge-tournament-gauge
script: 04-computation/integral_coordinate_change_depth_thm4091.py
output: 05-knowledge/results/integral_coordinate_change_depth_thm4091.out
script_sha256: c0125199819821e051f07c2e1418e6f1f09628acd22b1146b7e90806d4726a4c
output_sha256: a86fae7a2586cf41f34bba7a3a22ffa9d8193ab0f4bb43243ffc217d03fcf006
hash_basis: raw working-tree bytes (LF)
---

# THM-4091 -- LCM depth survives every integral coordinate change; literal depth does not

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.** There are two natural
ways to encode a denominator staircase. The cumulative `lcm(1,...,N)` version
is functorial under every integral formal change of variable. The sharper
coefficientwise condition `n^e a_n in Z` has an exact boundary: it is still
functorial at depth one, but not at any depth at least two.

Throughout, let

```text
R(q)=sum_(k>=0) a_k q^k in Q[[q]],
h(f) in f Z[[f]],
L_N=lcm(1,...,N),                    N>=1, e>=1.       (1)
```

Define the depth-`e` lattice

```text
D_e={R: a_0 in Z and k^e a_k in Z for every k>=1}.    (2)
```

## 1. Cumulative LCM clearing is always preserved

If `R in D_e`, then for every `N>=1`,

```text
boxed: L_N^e R(h(f)) mod f^(N+1) lies in Z[f]/(f^(N+1)).  (3)
```

Only terms with `k<=n<=N` contribute to `[f^n]R(h)`. For every such `k`,

```text
L_N^e a_k=(L_N/k)^e (k^e a_k) in Z,                  (4)
```

and `[f^n]h^k` is integral. Summing proves (3); the constant coefficient is
integral by definition of `D_e`.

This is the exact sense in which an LCM denominator bound survives an
integral coordinate change. It gives no coefficientwise `n^e` statement.

## 2. Depth one is exceptionally functorial

For every `R in D_1` and every `h in f Z[[f]]`,

```text
boxed: R(h) in D_1;
       n [f^n]R(h) in Z for every n>=1.               (5)
```

The mechanism is the formal derivative identity

```text
(n/k)[f^n]h^k=[f^(n-1)]h^(k-1)h' in Z,                (6)
```

obtained by comparing the coefficient of `f^(n-1)` in
`(h^k)'=k h^(k-1)h'`. Therefore

```text
n[f^n]R(h)
 =sum_(k=1)^n (k a_k) ((n/k)[f^n]h^k) in Z.           (7)
```

Equation (6), not a divisibility accident in small degrees, is the positive
depth-one mechanism.

## 3. Every higher depth fails at the first possible degree

For each `e>=2`, take

```text
R_e(q)=q^2/2^e in D_e,             h(f)=f+f^2.        (8)
```

Since `h^2=f^2+2f^3+f^4`,

```text
3^e [f^3]R_e(h)=3^e/2^(e-1),                         (9)
```

which is not an integer. Hence `R_e(h) notin D_e` for every `e>=2`.

Output degree three is universally minimal. At `n=1`, only `k=1` can occur.
At `n=2`, the `k=1` contribution is integral because `a_1` and `h` are
integral, while the `k=2` contribution is cleared exactly by the factor
`2^e`. Thus no member of `D_e` and no integral `h` can fail before `n=3`.
The hostile also uses the smallest possible nonlinear input index, `k=2`.

## 4. Exact preservation criterion for a fixed coordinate change

Fix `e>=1` and `h in f Z[[f]]`. Composition by `h` preserves `D_e` if and
only if

```text
boxed: (n/k)^e [f^n]h^k in Z
       for every n,k>=1.                              (10)
```

The entry vanishes when `k>n`, so only `1<=k<=n` matters. Necessity follows
by testing the basis element `R=q^k/k^e`. For sufficiency, write

```text
n^e[f^n]R(h)
 =sum_(k=1)^n (k^e a_k)
                  ((n/k)^e[f^n]h^k),                  (11)
```

a sum of integers. Thus (10) is the exact scaled composition-matrix test;
it is neither merely necessary nor merely sufficient.

## 5. Relation to the irrationality frontier

THM-4056's LCM clock gives denominator addresses but explicitly does not give
irrationality. The same firewall applies here. Equations (3)--(11) supply a
transport law for arithmetic lattices; an irrationality proof still needs a
nonzero integer linear form tending to zero after all denominator clearing,
as isolated in THM-4057.

The external August 2026 p-adic-zeta draft motivating this audit explicitly
avoids rowwise transport of literal `n^e` denominators after changing modular
coordinates. THM-4091 proves the elementary boundary behind that caution:
cumulative LCM clearing is safe, literal depth one is safe for the derivative
reason (6), and literal higher depth is unsafe without the matrix sidecar
(10). It does **not** verify any geometric, Frobenius, slope, or irrationality
claim in that draft.

Reproduce the finite controls with

```bash
python -B 04-computation/integral_coordinate_change_depth_thm4091.py
python -B -O 04-computation/integral_coordinate_change_depth_thm4091.py
```
