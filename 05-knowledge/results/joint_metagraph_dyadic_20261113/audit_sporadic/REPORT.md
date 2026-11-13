# Independent audit of signed pair lifts and eta depth

**Verdict: PASS for Propositions1–4 and the stated primary eta inputs.**
No mathematical repair is required. One optional title precision is noted
below. This audit is independent of the author checker: it imports none of
its code. Normal and optimized audit receipts agree byte for byte.

The audited files and their exact hashes are recorded in [receipt.json](receipt.json).
The scope is [BRIDGES.md](../sporadic/BRIDGES.md), especially signed orbital
transport, the pair lift, integral dyadic roots, coefficient valuations, and
nonrationality at a modular cusp. The inherited Golay, Ogg, and Monster
claims are not newly audited computationally here.

## Proof audit

1. **Signed stabilizers.** For a transported arc, an endpoint swap and the
   global reversal sign must cancel on every setwise pair stabilizer.
   This gives both necessity and sufficiency; choosing one orientation per
   unordered-pair orbit gives exactly two independent choices. The sign is
   a homomorphism, as assumed. The theorem is about a specified action,
   including its kernel, not merely the abstract group's order.

2. **Dihedral pair lift.** For odd composite m as well as odd prime m, a
   nonidentity rotation fixes no unordered pair of original vertices: a
   swap would force a nontrivial element of even order. The induced cyclic
   action is therefore free on the pair-vertices and their unordered pairs.
   Its orbit count is `binom(m(m-1)/2,2)/m`, exactly the asserted exponent.
   A negative reflection fixes `(m-1)/2` pair-vertices; two fixed vertices
   give an immediate impossible reversed arc. Thus m>=5 fails, while m=3
   is equivalent to the original signed action. The A5 argument is sound:
   its sign homomorphism is trivial, then simplicity would force any
   nontrivial image to embed A5 in an odd-order automorphism group.

3. **Logarithmic derivative and roots.** The forbidden-divisor schedule is
   `N does not divide d`, not `gcd(N,d)=1`; the recurrence correctly retains
   this difference. The coefficient of q is -r for every N>1. It both
   obstructs any integral 2^s-th root when `2^s` does not divide r and proves
   equality of the valuation bound at n=1. Conversely dividing Euler
   exponents constructs every allowed root. Frobenius modulo2 gives first
   odd exponent exactly `2^v2(r)`. The logarithmic derivative gives
   `2^v2(r) | n a(n)`, including zero coefficients, so the stated bound is
   correct. No implication about all coefficients being nonzero is made.

4. **Cusp nonrationality.** The inverse Möbius transform is exactly
   `gamma^-1(1/N+iy)=-1/N+i/(N^2 y)`. Trivial character permits equality of
   function values, and the leading q^-1 term gives exponential growth
   `exp(2*pi/(N^2 y))`. Multiplication by q changes this by a factor tending
   to one. Rational functions have only finite-order poles, so they cannot
   have that radial growth. A fixed finite matrix output, including any
   finite alteration of initial coefficients, has rational generating
   function. The proof does not exclude nonlinear updates, variable
   matrices, unbounded integers, growing state, or all notions of memory.

**Optional wording precision:** the title's phrase “eta recurrences require
 growing memory” is broader than Proposition4 by itself. The body already
states the proper boundary. A title such as “eta coefficients escape fixed
finite-matrix recurrences” would express exactly the proved obstruction.
This does not block the theorem or its stated consequences.

## Primary-source verification

Reopened 2026-11-13:

* [Aricheta, arXiv:1809.07421v2, AppendixB](https://arxiv.org/html/1809.07421):
  explicitly lists the eight required eta ratios for levels2,3,4,5,7,9,13,25
  as principal moduli of X_0(N). These are the curves before taking the
  Fricke quotient. The listed composite powers are8,3,1 at4,9,25.
* [Savitt, Theorem1.1 and Remark1.2](https://math.jhu.edu/~savitt/papers/pdfs/newman.pdf):
  the eta exponents +r at1 and -r atN satisfy the two congruences with sums
  -24 and24, with weight0. The character is trivial because N^r is a
  square for these eight levels. The cusp sign formula leaves only the
  infinity pole; at the intermediate cusp denominators2,3,5 in levels4,9,25
  the order is zero. The audit computes the width-normalized orders
  `r*(N-c^2)/(24*gcd(c^2,N))` for every c dividing N, giving order-1 at
  infinity, order+1 at denominator1, and zero at the intermediate cusps.

These inputs justify modular invariance and the simple pole used in
Proposition4. The elementary root and derivative proofs require neither
principal-modulus classification nor moonshine.

## Independent exact computation

[check.py](check.py) uses Python's standard library and explicit exceptions.
Its universe and controls are:

* For levels2,3,4,5,7,9,13,25, all coefficients through degree80. One path
  computes the full Euler series E(q), then forms `E(q)^r/E(q^N)^r` by
  formal polynomial division. Another path independently computes the
  filtered-divisor logarithmic recurrence. The first13 coefficients agree
  with the author receipt. The audit tests every permitted integral
  dyadic root by polynomial powering, the first forbidden root at q,
  exact parity support depth, and all coefficient valuation inequalities.
* For odd m=3,5,7,9, the signed stabilizer criterion is evaluated literally
  on every unordered pair and every group element, including composite9.
  It checks original vertices, pair-vertices, cyclic survivors, and the
  negative reflection's fixed points. At m=3 and5 every natural labelled
  tournament is tested directly for signed invariance, producing2 and4.
  At m=3 all pair-vertex tournaments give the positive count2. At m=5 all
  512 cyclic-invariant ten-vertex orientations are constructed, checked
  literally under every rotation, and shown to fail the signed reflection.
  They are labelled orientations, not512 isomorphism classes.
* A separate negative involution on three points with one fixed point is
  a positive control; a negative involution on four points with two fixed
  points is the hostile control. Exact rational coordinate arithmetic
  independently checks the cusp inverse transformation at all eight levels.

Reproduce from the repository root:

```sh
python3 -I -B 05-knowledge/results/joint_metagraph_dyadic_20261113/audit_sporadic/check.py
python3 -I -B -O 05-knowledge/results/joint_metagraph_dyadic_20261113/audit_sporadic/check.py
```

The finite checks support the implementation and equality cases. They do
not replace the all-order arguments or claim any transfer to the unrestricted
Hamilton spectrum. Chan–Pak and a complete tournament-growth recurrence
remain outside this audit.
