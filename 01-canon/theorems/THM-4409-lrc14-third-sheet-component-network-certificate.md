---
id: THM-4409
title: "LRC14 third-sheet component-network certificate"
status: >
  PROVED ANALYTICALLY + VERIFIED-EXACT + INDEPENDENTLY AUDITED. Restoring the
  bipartite contact graph between exact pair components and third-sheet
  components gives a max-flow refinement of THM-4396. At Fejer degree zero it
  is a rational component-length capacity; it closes every primitive distinct
  ternary-unit triple through height 79 and is exact on the equality comb
  (1,5,11). Universal height, entry, synchronization, and LRC(14) remain OPEN.
source: root + lrc_defect3_extra_cleanroom / LRC14 continuation session, 2026-09-03
depends_on:
  - THM-4392-lrc14-raw-carrier-box-spline-fourier-poisson-duality
related:
  - THM-4396-lrc14-finite-dual-exact-pair-hybrid-certificate
  - THM-4386-lrc14-canonical-component-relation-and-zero-defect-incidence
primary_script: 04-computation/lrc14_third_sheet_component_network_thm4409.py
primary_output: 05-knowledge/results/lrc14_third_sheet_component_network_thm4409.out
primary_script_sha256: 03cbe0a2f2da166d53353718179dab93fb174133a4a1b2d276f67c72917266bf
primary_output_sha256: 793049ba29d68d2c5bf43282f531cadf1fb0faf05bf53b4c0eb1d2b6d684577b
independent_audit_script: 04-computation/lrc14_third_sheet_component_network_thm4409_independent_audit.py
independent_audit_output: 05-knowledge/results/lrc14_third_sheet_component_network_thm4409_independent_audit.out
independent_audit_script_sha256: fe7fbe70b7cf7c7fdc278efa07f0b2cf363485d38de92b691e1fbcc8268e1348
independent_audit_output_sha256: 4e121a5d3e3875ecd3e1d28de18cec65cf8ac55f53d18692120c2da71df7d799
hash_basis: raw LF bytes
audit: >
  PASS. The dependency-free referee rederived the signed identity and both
  flow inequalities, supplied exact max-flow/min-cut and conservation checks,
  reconstructed every physical and raw carrier dictionary through height 79,
  and tested nonconstant rational weights on an equality row and a strict-loss
  hostile. Normal, optimized, and fixed-hash-seed outputs byte-match; 3,294,491
  explicit checks remain live under optimization.
---

# THM-4409 -- LRC14 third-sheet component-network certificate

**PROVED ANALYTICALLY + VERIFIED-EXACT + INDEPENDENTLY AUDITED.** This is a
local three-speed theorem and a finite census, not a proof of arbitrary
nonresonance, chart entry, fourteen-runner synchronization, or `LRC(14)`.
The Lonely Runner Conjecture at fourteen runners remains **OPEN**.

## 1. The restored sidecar

Use the notation and radius `lambda=1/14` of THM-4396. Fix a coordinate pair
`i,j`, let `k` be the remaining coordinate, and fix a sheet permutation `pi`.
Put

```text
P_pi = D_(w_i,pi(i)) intersect D_(w_j,pi(j)),
K_pi = D_(w_k,pi(k)),
g_pi = g_(w_i,pi(i);H) g_(w_j,pi(j);K).                (1)
```

Here `0<=g_pi<=1`. Define the exact signed pieces

```text
G_pi = integral_(K_pi) g_pi,
R_pi = integral_(P_pi intersect K_pi) (1-g_pi),
L_pi = integral_(K_pi minus P_pi) g_pi.                 (2)
```

The finite two-factor/exact-third term of THM-4396 is `M=sum_pi G_pi`, and
direct subtraction—not an estimate—gives

```text
mu(F_w)=sum_pi (G_pi+R_pi-L_pi).                        (3)
```

Thus the information missing from the exact-pair positive-part quotient is a
positive deficit on the true intersection and a negative leakage outside it.

Split `P_pi` into rational connected components `I` and `K_pi` into rational
connected components `J`. Join `I` to `J` exactly when their intersection has
positive length. For a nonnegative integrable weight `phi`, let
`Cap_(G_pi)(phi)` be the maximum total nonnegative edge load subject to the
vertex capacities

```text
sum_(e incident to I) x_e <= integral_I phi,
sum_(e incident to J) x_e <= integral_J phi.            (4)
```

The actual loads `x_(I,J)=integral_(I intersect J) phi` are feasible, because
distinct intersections incident to one interval are disjoint up to endpoints.
Therefore

```text
integral_(P_pi intersect K_pi) phi <= Cap_(G_pi)(phi).  (5)
```

Applying `(5)` to `1-g_pi` and `g_pi` gives

```text
R_pi <= Cap_(G_pi)(1-g_pi),
L_pi >= max(0,G_pi-Cap_(G_pi)(g_pi)).                   (6)
```

Substitution in `(3)` proves the finite component-network certificate

```text
mu(F_w) <= U_net(i,j;H,K),                              (7)

U_net=sum_pi [
 G_pi + Cap_(G_pi)(1-g_pi)
      - max(0,G_pi-Cap_(G_pi)(g_pi))].                  (8)
```

This refinement retains only component identities, contact incidence, and
the two vertex-integral capacities. It still forgets overlap endpoints, edge
integrals, and the coupling between the two separately optimized flows.

Source capacities imply `Cap(1-g)<=integral_P(1-g)`, so `(8)` is never weaker
than THM-4396's exact-pair bound. If `T_pi` is the union of touched pair
components and `E_pi` the union of isolated third components, then it is also
never weaker than the contact/disjoint-bit bound

```text
M + sum_pi integral_(T_pi)(1-g_pi)
  - sum_pi integral_(E_pi)g_pi.                         (9)
```

## 2. Exact degree-zero reduction

At `H=K=0`, each smoothed factor is `1/7`, hence `g_pi=1/49` and
`G_pi=1/343`. Let `kappa_pi=Cap_(G_pi)(1)`. Homogeneity gives

```text
Cap_(G_pi)(g_pi)=kappa_pi/49,
Cap_(G_pi)(1-g_pi)=48 kappa_pi/49.                     (10)
```

Since `kappa_pi<=|K_pi|=1/7`, the zero-mode terms in `(8)` cancel exactly:

```text
U_net(i,j;0,0)=sum_(pi in S_3) kappa_pi.               (11)
```

Thus `(11)` is a purely rational maximum-transportation-flow certificate
using only contact graphs and component lengths.

## 3. Why the equality comb is no longer obstructed

For `w=(1,5,11)`, exact rational construction gives, after summing the six
sheet permutations:

| pair | pair pieces / touched | third pieces / isolated | edges / nested / crossing | capacity |
|---|---:|---:|---:|---:|
| `(1,5)` | `6 / 6` | `68 / 62` | `6 / 6 / 0` | `6/77` |
| `(1,11)` | `12 / 6` | `32 / 26` | `6 / 6 / 0` | `6/77` |
| `(5,11)` | `12 / 6` | `8 / 2` | `6 / 6 / 0` | `6/77` |

Every sheet graph is a matching and every edge is nested. For any
nonnegative `phi`, the smaller nested interval is the intersection, so

```text
min(integral_I phi,integral_J phi)
 =integral_(I intersect J) phi.                         (12)
```

There is no competition between edges in a matching. Hence both flow
relaxations in `(8)` are exact for every `0<=g<=1`, in particular for every
coordinate pair and every finite anisotropic pair of Fejer degrees. The
definition-level comb has two raw carriers

```text
(-1,-2,1) and (1,2,-1), each of mass 3/77,             (13)
```

so `(8)` attains `mu(F_(1,5,11))=6/77`. This is precisely the third-sheet
incidence needed to remove THM-4396's all-finite-degree equality obstruction.

## 4. Complete exact census through height 79

The primary and independent implementations exhaust

```text
1<=w_1<w_2<w_3<=79,
each w_i odd and nonzero modulo 3,
gcd(w_1,w_2,w_3)=1.                                    (14)
```

There are `27` eligible speeds and `2,910` triples. Certificate selection
uses only the three degree-zero component networks. After selection, two
definition-level controls independently intersect the physical sheet lists
and enumerate the complete raw-carrier dictionary. They agree on every row.

The exact results are:

```text
physical mass below / on / above 6/77: 2909 / 1 / 0
fixed-pair network successes (01,02,12): 2818 / 2855 / 2859
selected pair counts (01,02,12):          400 /  533 / 1977
minimized network bound exact:                         1747
network bound on 6/77:                        (1,5,11) only
raw carrier-count range:                          0 through 22. (15)
```

The weakest strict selected certificate is

```text
w=(1,11,23), pair=(1,11), U_net=12/161,
6/77-U_net=6/1771.                                     (16)
```

Across all `8,730` pair choices, maximum flow equals the sum of edgewise
minimum capacities. That equality is only **FINITE-EXACT through height 79**;
no universal interval-graph identity is claimed.

The first strict information-loss hostile in lexicographic order is

```text
w=(1,19,79), best pair=(1,79),
mu(F_w)=108/10507 < U_net=8/553.                        (17)
```

Thus the graph sidecar is cheaper than exact overlap geometry but not
consumer-complete for the mass.

## 5. Verification and frontier

The primary uses exact `Fraction` arithmetic and chooses the certificate
before consulting either physical-mass control. The standalone referee does
not import producer or repository mathematics; every rational max flow has a
matching min-cut, capacity, conservation, and terminal-imbalance audit. It
also checks the signed identity and both feasible-flow inequalities for 54
nonconstant rational-weight sheet/pair cases on the strict hostile `(17)`.

Replay from the repository root:

```powershell
python -B 04-computation/lrc14_third_sheet_component_network_thm4409.py
python -B -O 04-computation/lrc14_third_sheet_component_network_thm4409.py
python -B 04-computation/lrc14_third_sheet_component_network_thm4409_independent_audit.py
python -B -O 04-computation/lrc14_third_sheet_component_network_thm4409_independent_audit.py
$env:PYTHONHASHSEED='731'
python -B 04-computation/lrc14_third_sheet_component_network_thm4409_independent_audit.py
```

Normal, optimized, and fixed-hash-seed streams byte-match the frozen LF
outputs. The primary records `43,723` live checks; the referee records
`3,294,491`. Neither relies on `assert`.

The next sharp question is whether the minimum of the three rational
capacities in `(11)` stays at most `6/77` beyond height `79`, or which exact
triple is its first hostile. Even a universal local theorem would still lack
chart entry and synchronization with the other eleven speeds. `LRC(14)`
remains **OPEN**.
