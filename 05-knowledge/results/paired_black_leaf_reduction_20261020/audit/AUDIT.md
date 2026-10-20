# Independent audit: paired-black leaves and cycle-cut deficit bounds

Audit date: 2026-10-20. Verdict: **APPROVED** for the stated reductions,
all-order arithmetic obstructions, and finite-exact order-nine exclusion.
The strong deficit-two center alphabet and classification of all bridges
remain open.

## 1. Reviewed claim and dependency boundary

Reviewed files in the parent packet directory:

- `BLACK_QUOTIENT_LEAF_CLASSIFICATION.md`
- `CYCLE_CUT_NEAR_PURITY_OBSTRUCTION.md`
- `check_deficit_two.py`
- `verification.json`

The mathematical inputs inherited, rather than reproved from their original
published sources, are the established paired-black degree identity, the
Hamiltonian-bypass theorem used in the rigid-atom classification, and Camion
Hamiltonicity. I read the actual local proofs of the SCC center law,
all-order rigid-atom deduction, and uniform prefix obstruction. The present
audit independently checks their applications and the entire new SC9
computation. It does not replace the prior literature audit.

## 2. Structural reduction

The degree formula uses merged mass and loop-degree two:

    degree_Q = (M-B)/2.

For NS nodes, `M=2f` and `B=0`, so leaves are exactly `f=1`. Hamiltonian
paths and automorphisms multiply over the uniquely ordered strong components.
There are no additional permutations of repeated components because every
arc between two distinct components points in one prescribed direction.
Hence the inherited all-order rigid strong alphabet `{1,C3}` makes the NS
leaf description exact, not merely constructive. Converse reverses the
component word, and each letter is self-converse. Nonpalindromes modulo
reversal are precisely the NS nodes.

For SC nodes, the center law gives `M=t^2 h`, `B=t b`. The half-product `t`
is a positive odd integer. Thus degree one means `t(th-b)=2`, forcing
`t=1` and `h-b=2`. The empty center has `h=b=1` and cannot supply a leaf.
Conversely every strong SC deficit-two center with any rigid half-word
supplies a leaf, uniquely by SCC decomposition.

The bound

    degree_Q >= t(t-1)/2

is correct since `h>=b>=1`. For positive degrees one and two, `t=1`.
If `t>1` and degree is three, equality forces `t=3` and `h=b=1`.
The direct `C5 -> C5^op` control below verifies that the boundary is sharp.
The pure C5 control is retained: nonrigid does not mean nonpure.

The generating function is exact with an unknown series `D`, not an
unproved rational closed formula:

    (1/(1-z-z^3) - (1+z+z^3)/(1-z^2-z^6))/2
      + D(z)/(1-z^2-z^6).

The first numerator counts reversal-nonfixed words and the second arbitrary
rigid half-words around one center. The leaf counts at orders 3 through 9,
`0,2,1,3,3,6,7`, agree. Neither this reduction nor the finite center counts
establishes that every bridge is pendant. At order four the two leaves
share one edge, and the bridge lower bound is only asserted from order five.

## 3. Cycle-cut lemma and arithmetic bounds

For a Hamiltonian successor `c`, let `D=Aut(T) intersect <c>`,
`d=|D|`, and `q=n/d`. A map between two ordered cut paths is uniquely
the corresponding power of `c`; no automorphism outside this cyclic group
can merge more cuts. Thus the cuts meet exactly `q` path orbits.

The reversal for cut `i` is `r_i(j)=2i-1-j`. If one cut `a` is blue, then
the others are blue exactly when `q` divides `2(i-a)`. Thus the number of
blue cut-orbits is exactly `gcd(2,q)`, or zero if there is no blue cut.
Using literal positional reversal avoids confusing fixed paths and fixed
quotient orbits.

Since tournament automorphisms cannot have even order, `d` is odd.
For even `n`, `q>=2^v2(n)` and at most two cut-orbits are blue. This proves

    h-b >= 2^v2(n)-2.

For odd strong tournaments of order at least seven, some Hamiltonian
successor is not an automorphism. If an initially selected successor is
already nonautomorphic there is nothing to show. Otherwise it gives a
circulant presentation, and the uniform prefix proof supplies a closing
nonblue path. Its successor cannot be automorphic, because every cut of
a cycle with automorphic successor is a reflection in a circulant.
For this cycle `q>1` is an odd divisor of `n`, proving

    h-b >= q-1 >= p_min(n)-1.

The hypothesis `n>=7` is essential to this proof; C3 and regular C5 are
positive pure controls. Strongness supplies a Hamiltonian cycle, and SC
supplies the stated `h,b` interpretation. The argument does not assume
every Hamiltonian successor is nonautomorphic.

All prefix witnesses really close: F1 has closing increment 3; F2 has
closing increment `k`; each F3 variant has closing increment 2. The source
hypotheses put these in the connection set. Unit coordinate changes preserve
closure. The four order-seven closing increments are `4,4,3,2`.

Consequences for deficit two are correctly restricted:

- Odd order at least seven must be divisible by three.
- Order divisible by eight is impossible.
- At order `2 mod 4`, a cycle's two cut-orbits may be both blue or both
  black. No invalid automatic-reflection claim is made.
- At order `4 mod 8`, every cycle has four cut-orbits, exactly two black.
  All Hamiltonian cycles belong to one automorphism orbit; the cycle
  stabilizer is `D`, giving `HC=4|Aut|/n`.
- If both global black path orbits occur as cycle cuts, every bypass path
  is blue. This follows because a closing path determines its directed
  cycle, and endpoint orientation is invariant.

The arithmetic obstruction leaves odd multiples of three beginning at
15, or even orders at least 10 with `v2(n)<=2`. No finite cutoff is
proved by this obstruction.

## 4. Independent computation

`independent_audit.py` imports no producer implementation. It uses:

1. A fresh adjacency-matrix reconstruction of normalized anti-involution
   frames.
2. An explicit simultaneous permutation action on the paired indices.
   All supplied representatives have disjoint frame orbits covering
   every frame mask, including all 65,536 order-nine frames.
3. Direct enumeration of symmetric Hamiltonian words and centralizing
   permutations to recompute `B=h_rho/centralizer`. These agree with every
   input row. The anti-involution conjugacy proof in the inherited odd-blue
   audit establishes that normalized frame classes equal unmarked SC classes.
4. A backwards top-down subset recurrence for `H`, not the producer's
   forward extension recurrence.
5. Floyd-Warshall strongness and an independently encoded equitable
   refinement.
6. Direct products of permutations within outdegree cells to compute the
   surviving automorphism groups, deliberately using a coarser restriction
   than the producer's refinement/backtracking implementation.
7. A comparison to the producer's rows only after independent calculation.

Every row matches at orders 3,5,7,9. Order nine gives:

    all classes                  2752
    strong classes               2661
    integral required Aut         184
    odd divisor of 9!              37
    centralizer divisibility       37
    refined-cell divisibility       5
    actual deficit-two classes      0

The five candidates reproduce `(mask,H,B,required Aut,actual Aut)`:

    (32231,15,3,3,5)
    (64526,3,1,1,3)
    (64559,9,1,3,9)
    (64605,27,1,9,27)
    (64607,9,1,3,9)

All are nonstrong and pure. Thus every strong order-nine class is
eliminated by necessary filters before automorphism enumeration.

The filter directions are sound. Equality requires `Aut=H/(B+2)`.
An actual group order is an odd divisor of `n!`, is divisible by the
known subgroup centralizer order, and divides the order of the direct
product of symmetric groups on invariant color cells. These are group
divisibility facts, not numerical upper-bound heuristics.

Direct small controls, enumerating Hamiltonian words and degree-preserving
automorphisms independently:

| Object | H | Aut | B | h | h-B |
|---|---:|---:|---:|---:|---:|
| C3 | 3 | 3 | 1 | 1 | 0 |
| S4 | 5 | 1 | 3 | 5 | 2 |
| Regular C5 | 15 | 5 | 3 | 3 | 0 |
| Paley7, increments {1,2,4} | 189 | 21 | 3 | 9 | 6 |
| C5 -> C5^op | 225 | 25 | 3 | 9 | 6 |

Paley7 attains the prime-order lower bound. It has 24 Hamiltonian cycles:
three with `(q,blue cut-orbits)=(1,1)` and 21 with `(7,1)`.
S4 has one cycle with `(4,2)`. C5 has two, both with `(1,1)`.
The stacked C5 control has paired-black degree three and is nonstrong,
so it is a boundary for the outer-rigidity theorem, not a strong-center
counterexample.

`check_closing_prefixes.py` is a second fresh implementation of the
displayed prefix constructions, not an import of their verifier. It checks
directed cycle closure, a failed reversal arc, and a failed successor arc
for every normalized circulant presentation through order 31, the four
order-seven witnesses, and all three residual patterns for `4<=m<=100`.
These finite checks support the symbolic proof; they are not its basis.

## 5. Reproduction

From the parent packet directory:

    python3 audit/independent_audit.py --source ../odd_blue_spectrum_hollow_cubes_20261020
    python3 audit/check_closing_prefixes.py

Outputs and exact small-control data are in `verification.json`,
`verification.out`, and `check_closing_prefixes.json`.
The original reviewed producer hashes are recorded in
`AUDITED_SOURCE_SHA256SUMS`; the final public packet hashes are in the
parent `SHA256SUMS`. Subsequent changes were presentation-only:
removing work-process language, recording this completed audit, making
reproduction paths portable, and clarifying that this obstruction does
not establish a finite cutoff. No mathematical proof, filter, recurrence
or recorded result was changed.
