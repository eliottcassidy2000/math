# The first modular term in the Paley run expansion

2026-10-14. **PROVED exact reduction identities; FINITE-EXACT symbolic
coefficient extraction and independent controls; CITED modular substitution.
The uniform asymptotic remainder is OPEN.**

This extends the incoming S34 packet
[H31, rate and wheel](paley_h31_exact_rate_and_apery_wheel_20261014.md),
commit `294e754ecd`, and HYP-9274
[`paley-rate-three-halves-wheel-w4-apery-form`](../hypotheses/HYP-9274-paley-rate-three-halves-wheel-w4-apery-form.md).
It does not infer a universal theorem from the 27-prime match.

## Inheritance and live concept board

- Closest proved mechanism: run reversal and exact character elimination;
  THM-438 [`paley-cluster-integrals-are-catalan`](../../01-canon/theorems/THM-438-paley-cluster-integrals-are-catalan.md),
  particularly its June 10 joint-run addendum.
- Hostile example: an odd run inside an even-total multi-run sum still
  vanishes; global parity alone misses this.
- Corrected near miss: MISTAKE-061 replaces even cacti by even-series graphs;
  merely classifying cycles loses theta cores.
- Least-used sidecar: keep a formal edge weight `t`, run boundaries and
  orientation signs before specializing to the full Hamiltonian count.
- Board: injectivity / marked runs / terminal cores / modular traces /
  exact reconstruction. The anchor is the explicit Paley coefficient; the
  OEIS recurrence packet is the niche; the source-specific order-budget
  packet is the wildcard.

## 1. Exact statement

Let `q=3 mod4` be an odd prime power, and let chi be the quadratic character
of its finite field, extended by chi(0)=0. Define

`R_q(t) = (1/q!) sum_permutations sigma product_(i=1)^(q-1)
          (1+t chi(sigma(i+1)-sigma(i)))`.

Thus `R_q(1)=2^(q-1)H(P_q)/q!`; `t` retains which subsets of edges contribute.
Let `J_(L1,...,Lr)` sum the corresponding oriented path products over
globally distinct field values on the disjoint marked runs. Let `U_q` be
the **unanchored, injective** wheel sum, with hub-to-rim edges and a directed
four-cycle on the rim. Its orientation agrees with the S34 K2,2 rim gauge:
the two gauges differ by two edge reversals.

For q>=11,

`[t^8]R_q(t) = (P_8(q)+7U_q)/(q)_8`,

where `(q)_8=q(q-1)...(q-7)` and

`P_8(q)=q(q-1)(q-3)
 (q^5-42q^4+766q^3-7452q^2+37561q-75306)/24`.

For an **odd prime** q=3 mod4, the all-prime identity proved by reduction to
the Ahlgren--Ono moment in [the wheel packet](wheel_apery_extension_20261014h.md)
gives `U_q=-q(q-1)(a_q+4)`. Consequently

```text
                 q^6-45q^5+892q^4-9750q^3+59917q^2-187989q+225246-168a_q
[t^8]R_q(t) = ----------------------------------------------------------------.
                               24(q-2)(q-3)(q-4)(q-5)(q-6)(q-7)
```

The exact modular contribution is therefore
`-7a_q/((q-2)...(q-7))`. This is a coefficient of the finite polynomial in t,
not the entire arithmetic contribution to `R_q(1)`. Larger edge sets may
add further terms involving the same or other motives.

## 2. Proof mechanism: terminate at the first irreducible core

For an injective graph sum, parallel copies of an edge cancel in pairs
because chi(x-y)^2=1 when x and y are distinct. Reversing an edge multiplies
the sum by -1. An isolated vertex on a graph of n vertices contributes
the exact factor q-n+1 when removed.

Eliminate a degree-one vertex by summing its character over the whole field
(zero), then subtracting the values already occupied by the other vertices.
Each subtraction identifies the removed vertex with one occupied value,
retaining the orientation of every factor.

For degree two, orient both incident factors outwards, with the appropriate
sign. For distinct neighbours a,b,

`sum_x chi((a-x)(b-x))=-1`.

Subtract occupied values as before. Both operations lower the vertex count;
they never increase the edge count. These are exact identities over every
odd finite field, with the orientation convention specialized to q=3mod4.

After these reductions, a simple terminal core with at most eight edges
has minimum degree at least three, hence at most five vertices. The only
possibilities are K4 and the four-rim wheel. K4 vanishes: a transposition
of two values reverses its Vandermonde character. On five vertices, eight
edges and minimum degree three force degrees (4,3,3,3,3); the complement
is a matching plus an isolated vertex, uniquely identifying the wheel.
Thus all sums in scope reduce to a polynomial plus a wheel multiple.

The [executable recurrence](../../04-computation/wheel_run_interface_20261014h.py)
implements exactly these rules with rational polynomial coefficients. Its
complete reductions give the following wheel coefficients:

| Marked runs | Wheel coefficient |
|---|---:|
| (8) | 0 |
| (6,2) | -14 |
| (4,4) | -6 |
| (4,2,2) | 128 |
| (2,2,2,2) | -960 |

These coefficients are independently recovered by enumerating all covers
of the eight wheel edges by the specified ordered trails, weighting the
induced equality partition by its Mobius value, and dividing by the eight
wheel automorphisms. Canonical restricted-growth partitions independently
check that every orbit has exactly eight labelled covers. Odd-run controls
(7,1) and (5,3) give zero signed coefficient.

A slot subset with E selected edges and r ordered maximal runs has
`binomial(q-E,r)` placements. Its normalized joint sum has denominator
`(q)_(E+r)`, so the combined weight is `1/(r!(q)_E)`.
When E+r>q the injective joint is zero, and the simplified expression still
contributes zero. With E=8 this yields

`S_8=A_8+J_(6,2)+(1/2)J_(4,4)+(1/2)J_(4,2,2)+(1/24)J_(2,2,2,2)`.

The wheel coefficient is exactly `-14-3+64-40=7`.
The same reductions produce P_8 above; full polynomial rows are in the
[machine-readable output](wheel_run_interface_20261014h.json).

In particular, the previously finite-tested single-run patterns satisfy

```text
A_4 = 2q(q-1)(q-3),
A_6 = q(q-1)(q-3)(5q-29),
A_8 = 2q(q-1)(q-3)(q-7)(7q-55),
J_(2,2) = q(q-1)(q-3)(q-9).
```

All four formulas hold over every odd prime-power field q=3mod4. The exact
reducer also gives these smaller cases. Before edge degree eight, odd-run
reversal removes odd degrees and the remaining terminal K4 vanishes, so
the run coefficients in that finite range are rational functions of q.
This is an algebraic decomposition statement, not a claim that a finite
data set disproves every possible polynomial fit to modular traces.

## 3. What the rate calculation still needs

The exact two-cherry formula gives

`J_(2,2)-A_2^2=-11q^3+38q^2-27q`.

This corrects an **asymptotic extrapolation**, not the old exact integer
values: the June addendum's roughly -10 ratio at small q tends to -11.
The incoming S34 two-cherry repulsion factor is consistent with this.
Cherry placement, pair interaction and the quartic run give the formal
first-order candidate c1=3/2. The present coefficient calculation does not
justify exchanging q->infinity with the sum over all run configurations.
That uniform summation bound remains the decisive analytic target.

The finite-field trace is an example of information surviving contraction:
the polynomial part records eliminated vertices, while the wheel retains
a modular statistic. Keeping only low-order polynomial fits erases it.

## 4. H(P31): complete the CRT uniqueness argument

Incoming S34 reports the three residues of the fixed-arc count N1, modulo
`2147483647, 2147483629, 2147483587`. Their product is

`M=9903519940736477367306812281`.

Its reconstructed value is `N1=42745655163116955307561`, so

`H(P31)=465N1=19876729650849384218015865`.

Checking that a reconstructed residue is small does not itself bound the
true count. Here is an independent upper certificate which closes that
gap without assuming the HYP-2371 prediction.

After fixing the arc 0->1, encode each remaining path by its successor
permutation: rows are vertex1 and the 29 unused vertices; columns are the
29 unused vertices and a terminal symbol. Each row has at most sixteen
allowed successors (fifteen out-neighbours plus the terminal). Paths inject
into these permutations; unwanted cycles only enlarge the comparison set.
The Bregman permanent bound gives

`N1 <= (16!)^(30/16) < 10^25 < M`.

The [exact script](../../04-computation/wheel_run_interface_20261014h.py)
checks the integer inequality `(16!)^30 < 10^400 < M^16` and the reported
three residues. Thus the CRT reconstruction is unique. This verifies the
reconstruction interface, **not a second execution of the expensive q31 DP**.

For completeness, the bound follows by exposing a uniform allowed
permutation in a uniformly random row order. Conditional entropy of row i
is at most the log of its unused allowed columns. Given a full permutation,
the number of such columns is uniform on 1,...,d_i as the row order varies.
Thus log(permanent)<=sum_i log(d_i!)/d_i. Monotonicity of (d!)^(1/d) gives
the displayed bound. See also Radhakrishnan's
[entropy proof](https://doi.org/10.1006/jcta.1996.2727) and
[Galvin, Theorem1.1](https://academicweb.nd.edu/~dgalvin1/pdf/bregman.pdf).

## 5. Reproduction and next obligations

Run `python -B -X utf8 04-computation/wheel_run_interface_20261014h.py`.
Independent checks comprise direct injective tuple sums at q=3,7,11,
the earlier A8 values at19,23,31, and a separate Held--Karp evaluation of
the complete weighted permutation polynomial at11. It gives
`[t^8]R_11=2/135` and recovers `H(P11)=95095`. A second agent independently
obtained the same coefficient with a separately written DP.

The next bounded targets are (i) classify exact cores at ten edges while
allowing degree-two elimination to expose smaller cores; (ii) retain trace
symbols in a general run compiler; and (iii) prove a uniform tail bound
strong enough for the proposed c1. A direct Euler-parity obstruction to a
raw wheel quotient does **not** by itself rule out a wheel exposed after
eliminating degree-two vertices at higher edge degrees.
