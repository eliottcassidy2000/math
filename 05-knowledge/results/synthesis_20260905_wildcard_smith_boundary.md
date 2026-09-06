# The complete-residue two-jet boundary: one lost derivative line costs exactly one p

**Status: PROVED in the stated scope; exact replay and independent proof
audit PASS.** Exact companions verify the consequence object. This
closes the `s=p` boundary explicitly left open by
[THM-4080 / confluent-two-jet-single-scale-smith-partition](../../01-canon/theorems/THM-4080-confluent-two-jet-single-scale-smith-partition.md).
It does not settle larger residue clusters, higher jet orders, or a prize
problem. No external novelty or publication-priority claim is made.

## Inheritance and concept board

The closest proved mechanism is THM-4080's weighted-minor argument. The
canonical hostile is its `p=2,s=2,e=1` profile `(0,0,2,2)`, against the false
continuation `(0,0,1,3)`. The corrected near miss is the same theorem's
restriction `s<p`: this restriction correctly avoids a derivative-rank drop,
but rank modulo p does not measure the depth of that drop. The least-used
sidecar is the **saturation index of the derivative row lattice**. That
index is exactly p here, independently of the geometric scale `p^e`.

The live board considered (i) integral observer lattices and their derivative
filtration; (ii) AMM Rule-A kernel words and the feed/phase sidecar
(THM-4086 / rule-a-transition-clock-and-phase-cocycle); (iii) Mahler's safe
language and terminal integer address (THM-4072 / mahler-safe-terminal-fibre-product-and-finite-state-obstruction);
(iv) rooted tournament response (THM-3729 / rooted-pfaffian-response-and-sign-root-deletion-average); and (v) whole-layer Frobenius
preservation (THM-2022 / gmc2-frobenius-lowest-balanced-face). The first lane
produced a decisive positive signal. The other lanes reinforce the need to
retain the coordinate consumed by the next operation, but this proof supplies
no map from their physical targets to a Smith partition.

## Theorem

Let p be any prime, e>=1, and let

```text
x_i = a + p^e u_i,       0<=i<p,
```

be integer nodes, where the u_i are pairwise distinct modulo p. Thus all
pairwise node differences have exact valuation e. On integer polynomials
of degree less than 2p, take values and first Hasse derivatives:

```text
J(P) = (P(x_i), P'(x_i))_(0<=i<p).
```

The p-primary Smith exponents, in increasing order, are

```text
0, 0,
e, 2e, ..., (p-2)e,
(p-1)e+1,
(p+1)e, (p+2)e, ..., (2p-2)e,
(2p-1)e-1.                                           (1)
```

Empty ranges are omitted. In particular the formula includes p=2:
`(0,0,e+1,3e-1)`. At e=1 the two largest exponents coincide. At p=2,e=1
the two positive exponents also coincide. These are genuine Smith factors,
not merely determinant or rank data.

Relative to the false `s=p` continuation of THM-4080, exactly one unit of
p-adic valuation moves from the last invariant factor to factor p+1.
The total determinant valuation stays `2ep(p-1)`.

## Proof through all determinantal divisors

Translation of the polynomial variable is integral unimodular, so assume
a=0. Let nu_h be the p-adic valuation of the h-th determinantal divisor,
and put nu_0=0. We prove

```text
nu_h = e*(h-1)*(h-2)/2,               1<=h<=p,
nu_h = e*(h*(h-1)/2-p)+1,             p<h<2p,
nu_(2p) = 2ep(p-1).                                  (2)
```

Their consecutive differences are exactly (1).

For an h-minor using degrees q_1<...<q_h and d derivative rows, extract
`p^(e q_j)` from column j after multiplying each derivative row by p^e.
The exact valuation is

```text
e*(sum q_j-d) + v_p(residual minor),                  (3)
```

where the residual entries are `u_i^q` and `q*u_i^(q-1)`.

For h<=p the THM-4080 lower bound and witness still apply: use degrees
`0,...,h-1`, one value row and h-1 derivative rows. Their determinant is a
unit times `(h-1)!`, which remains a p-adic unit even at h=p. Thus the
first line of (2) holds. For h=2p the confluent Vandermonde determinant is
`product_(i<j)(x_j-x_i)^4`, giving the third line directly.

It remains to treat p<h<2p. Write `L=h*(h-1)/2-p`. Since d<=p, every
nonzero h-minor has degree cost `sum q_j-d>=L`. Equality forces precisely
degrees `0,...,h-1` and all p derivative rows. All other choices cost at
least `e(L+1)>=eL+1`.

Consider the p by h derivative matrix R over Z_(p) on those consecutive
degrees. Its reduction modulo p has rank exactly p-1. Indeed, its columns
for q=1,...,p-1 evaluate the independent functions
`1,X,...,X^(p-2)` on every element of F_p. Column q=p is zero. For
`p+1<=q<=h-1<=2p-2`, the function `X^(q-1)` on F_p reduces to
`X^(q-p)` with degree at most p-2, so it adds no direction. Consequently
every minimum-cost residual minor is divisible by p. This proves
`nu_h>=eL+1`.

The p columns q=1,...,p of R have determinant

```text
p! * product_(i<j)(u_j-u_i),                          (4)
```

of exact valuation one. Therefore the rectangular Smith form of R over
Z_(p) has p-1 unit factors and one factor of valuation one. Equivalently,
there is a unimodular operation U on its p rows such that the first p-1
rows of UR remain independent modulo p, while the last row is divisible
by p, and dividing that last row by p gives a matrix R* of full row rank p
modulo p. This is the saturation step. No column operation is needed for
the following extension argument; the row-normal form follows equally from
the rectangular Smith decomposition.

The full matrix of values and derivatives on degree less than h has rank h
modulo p. Its kernel would otherwise contain a nonzero polynomial of degree
less than h<=2p with every element of F_p as a double root; divisibility by
`(X^p-X)^2` of degree 2p makes that impossible. Replacing derivative rows
by R* cannot reduce their span together with the value rows: the span of
the original derivative rows modulo p is contained in that of R*.
Consequently we can choose h-p value rows which, together with all rows of
R*, make an invertible h by h matrix modulo p. Undoing the division of the
last derivative row multiplies its determinant by p; undoing U multiplies
by a unit. Thus those same h-p value rows and the original p derivative
rows give a residual minor of **exact valuation one**. Its scale cost is
eL, so it attains the lower bound. This proves the middle line of (2), and
hence the theorem.

## Consecutive-node corollary: the complete quadratic range

For values and first derivatives on n consecutive integer nodes, write
`n=qp+r`, `0<=r<p`. For every `1<=n<=p^2`, the complete p-primary
partition is the sorted multiset union over the p residue classes of:

- THM-4080's profile `0,0,1,...,s-1,s+1,...,2s-1` for a class of size
  `1<=s<p`;
- formula (1), with e=1, for a class of size s=p;
- the empty profile for an empty class.

There are r classes of size q+1 and p-r of size q. Distinct residue factors
have unit resultants, so the p-local CRT decomposition from
[THM-4010 / confluent-consecutive-hasse-observer-kernel-index-and-smith-firewall](../../01-canon/theorems/THM-4010-confluent-consecutive-hasse-observer-kernel-index-and-smith-firewall.md)
splits the Smith module. Every cluster has size at most p and every pair of
nodes within it differs by p times a p-adic unit. Hence the two single-scale
formulas apply. This extends the previously proved range
`n<=p(p-1)` to `n<=p^2`, including the entire formerly missing last band.

## Connection, destroyed data, and the next boundary

The source object is the integral two-jet observer at one p-adic scale. The
map to the target extracts the deterministic column/row scale weights and
reduces the remaining confluent evaluation matrix modulo p. It preserves
the exact weighted-minor lower bounds and the unit/nonunit distinction. It
forgets the valuation of a nonunit residual minor. The derivative saturation
index p restores precisely that missing coordinate. The cheapest hostile,
`p=2,e=1`, predicts the +1 correction; changing to e=2 distinguishes the
actual p-adic cost +1 from the false geometric cost +e.

The Frobenius kernel `d(X^p)/dX=0` modulo p is the rank-loss mechanism.
The integral lift `d(X^p)/dX=pX^(p-1)` shows that the missing direction
returns after exactly one division by p. The argument turns a modular
vanishing obstruction into a precise integral cost. This is a substantive
transfer from the residual derivative lattice to **every** intermediate
determinantal divisor, not a scalar determinant analogy.

For s>p, pairwise distinct residues u_i are impossible. Thus a larger
single cluster must split at a second p-adic scale; the present proof cannot
be extended by replacing p with s. The next genuine task is a multiscale
cluster law with saturation data, first at n=p^2+1 for consecutive nodes.
For higher jets the Hasse derivative rank and saturation can lose several
directions; the one-line correction established here supplies no general
higher-jet formula.

An independent cross-lane review accepted the row-saturation/minor argument
and the consecutive-node CRT step, and reconstructed the complete p=2
block symbolically. See
[independent audit](synthesis_20260905_lrc_audit_smith_boundary.md).

## Reproduction

Run the standard-library audit and its optimized replay:

```text
python3 04-computation/synthesis_20260905_wildcard_smith_boundary.py
python3 -O 04-computation/synthesis_20260905_wildcard_smith_boundary.py
```

The companion's exact universe, positive controls, hostile controls,
independent determinantal-divisor path, and matching output are recorded in
the source and `synthesis_20260905_wildcard_smith_boundary.out`. Computation
checks the closed formula and witnesses; the proof above establishes the
universal quantifiers. No finite enumeration is promoted to an all-scale
or cross-problem theorem.

The completed audit has 72 boundary matrices (six primes through 13, three
scales, four residue-coordinate systems), 11 inherited positive controls,
70 derivative rank/saturation rows, six exhaustive all-minor cases, and
13 direct consecutive matrices in the final quadratic band. Exact rational
DVR elimination and a separate modular DVR algorithm agree, with the latter
using precision greater than the determinant's known total valuation rather
than assuming any predicted exponent. There are 450,064 exact gates.
Normal and optimized streams agree byte for byte. The script has no Python
assert statements, floating-point literals, or carriage-return bytes.

```text
source_sha256:
5bc8ccd4620d96fc348641ff6e6ed3f700443b7190e1b9b2e53143cbdb02faba
output_sha256:
b5e3ae9d90d6bd116bfcf3d5f9c74e489be37d7f9f9da5302b678a1074aa2285
semantic_sha256:
1503593ff87f7488a820182dd0d3168fd9dac39633bfa931b23ca3aab9272e81
```

On the present Windows host, invoking a `.py` path through the Python
application alias can silently return without executing it. The fully
equivalent command actually used for the normal replay was
`python3 -c 'import runpy; runpy.run_path("04-computation/synthesis_20260905_wildcard_smith_boundary.py",run_name="__main__")'`,
with `-O` added for the optimized replay.
