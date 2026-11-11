# Endpoint-compatible seeds beyond extension gcd one

**11 November 2026. PROVED IN NOTES and FINITE-EXACT; Chan–Pak OPEN.**
The constructive advance is a finite endpoint certificate for infinite
seed families, with actual adjacency constructors. It gives all counts

    1+(b*3^a-1)*2^s,       b in {30,50,78}, a>=1, s>=1,

and additional binary and mixed binary/ternary coefficient families.
It also realizes the requested 99*2^k+23 for every k>=1, on k+7 vertices.
These are actual Hamiltonian-path counts, not just residue or gcd results.
Seeds for every odd c and a proof of Chan–Pak remain open.

## Results and exact limits

- [E1–E6 and E10 proofs and census](PROOFS.md): D_(n-1)=2 for the extension
  difference lattice (n>=2); the complete one-port endpoint criterion;
  ten explicit inflation presentations; finite-base classification through
  order seven; the remaining source, endpoint and below-ray obligations.
  E10 integrates the incoming strong-core normal form and derives exact
  generating-function transforms that remove all exterior padding.
- [E7–E9 connections](CONNECTIONS.md): the 99 ray and its exact 25-prime
  wheel; source-valid arithmetic recursions; an inflation-only Zeckendorf
  obstruction; concrete controls for the Carmichael fibre analogy; the
  Catalan carry connection and finite path automaton.
- [Exact certificates](verification.json): independent physical-graph
  recovery of every coefficient in the ten named infinite families,
  60 literal tail counts, 30 backtracking controls and ten C_3 substitutions.

For a base B and a selected set O, the inherited closure count is
delta(B,O)+2^t U(B,O). Replacing an unselected vertex v by a transitive
block of size a gives

    delta(a)=sum_j d_j j!S(a,j),    U(a)=sum_j u_j j!S(a,j).

There are at most |B| terms. The seed condition survives for every a
exactly when d=(1,0,...,0); equivalently it survives at a=1,...,|B|.
This is a finite proof certificate for infinitely many actual graphs.
At such a port, substituting an arbitrary tournament M instead changes
the offset to H(M), so an arbitrary seed cannot simply be plugged in.

The complete ordinary-base census through seven vertices finds 1,407
eligible marked presentations, 634 stable ports, and 763 pairs of stable
ports. Every one of those pairs is jointly stable at all positive block
sizes, certified by its two-variable coefficients. This is a statement
about that finite base universe, not all tournaments. The corresponding
marked-presentation sequences are

    P_n: 0,2,9,26,64,152,381,
    Q_n: 0,0,3,19,66,186,489,       n=1,...,7.

Only 95 of those 634 ports and 49 of those 763 pairs lie inside strong
cores. The others follow from three exact padding transforms in E10.
All named infinite constructions grow strong cores and have order at
most 7+floor(log_2(N-1)) for their realized counts N.

Two seven-vertex examples instead preserve the unit seed only at sizes
{1,2} and {1,3}. The second is governed exactly by
3^a-4*2^a+5=0. These pinpoint why finite success without a coefficient
certificate does not establish repeatable recursion.

The 25-prime wheel for 99*2^k+23 has period 428468040. Only eleven of
those primes impose exponent restrictions. Infinitely many terms pass
all 25 tests and are composite with factor 419. The graph realization
holds in both cases; no primality inference is needed.

## Construct the graphs

Run from this directory:

```text
python3 construct.py --family ternary30 --size 3 --tail 1
python3 construct.py --family seven --size 20 --tail 1
python3 construct.py --family ray99 --tail 6
```

The first emits a ten-vertex tournament with H=1619 and seed coefficient
c=809. The second gives c=3670019 and H=14680077 on 26 vertices. The third
gives H=6359 on 13 vertices. Larger parameters require no new search.
Outputs include adjacency rows, selected vertices and the expected count.
`construct.py --help` lists all ten family names; `ternary90` is an audited
duplicate of the later portion of `ternary30`, not a fourth disjoint orbit.

## Reproduce and audit

All scripts use Python's standard library. From this directory:

```text
python3 ray_probe.py
python3 inflation_bank.py
python3 inflation_bank.py --max-n 7
python3 composition_probe.py
python3 joint_ports.py
python3 endpoint_checks.py
python3 core_compression.py
python3 wheel.py
python3 verify.py
```

`inflation_bank.py` exhausts ordinary bases and every selected set/port in
the stated universe. Its quotient-word computation is controlled by 1,276
literal inflated graph counts through order-seven bases. `verify.py`
uses a separate physical-graph DP and reconstructs W from the curvature
of one-vertex counts, rather than importing the endpoint-cover computation.
Triangular inversion at all required sizes checks the entire coefficient
list for each named infinite family. Optimized-mode replay and hashes are
recorded in `AUDIT.md` and `manifest.json`.

The source pins include the two imported counting/canonicalization modules
from the previous [extension-content packet](../modular_centers_collatz_pell_20261111/README.md).
The inherited [all odd c<=4097 seeds](../odd_c_seed_symmetry_20261110/README.md)
remain available and are not claimed as rerun here. Their global below-ray
coverage uses its audited finite range; that argument does not automatically
extend to the new unbounded coefficient families with starts greater than 1.
