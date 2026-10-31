# Session131: square-sum reflection words and parity seams

**PROVED uniform bounds and exact maps; FINITE-EXACT sharpness certificates.**

If C is any square-sum Hamiltonian cycle on1..N, N>=8, the number of edges
whose square-sum root is even obeys the sharp residuewise uniform bound

```
N mod8:     0  1  2  3  4  5  6  7
minimum:    4  7  8  7  8  7  8  7
```

“Sharp residuewise” means actual cycles at80,81,...,87 attain all eight
bounds. It does **not** assert attainment at every order. The unique G32
cycle has16 such edges. The theorem is proved by cutting even-root edges,
retaining all odd-root path components and their endpoint incidences, then
balancing the eight residue pools. See [PROOFS.md](PROOFS.md), section1.

The packet also proves:

* the marked-endpoint path analogue, retaining both exposed ports;
* H_N, the odd-root subgraph, has exactly four modulo8 components for N>=44;
* an exact compiler for root-labelled partial-reflection cycle words:
  even words have an alternating-shift parameter; odd words are rigid;
* the explicit Zeckendorf odd/even-row interleaving constraints and their
  signed even-Fibonacci value map; numeral reflection requires exact
  decoding/re-encoding, not raw digit reversal;
* the global integer-affine edge-preserver classification: square scaling
  only, distinguished from the phase-dependent Gerbicz block shifts.

## Inheritance boundary

The modulo8 channels and CG/Zeckendorf/forest bridges are recovered from
Session125, with exact paths and hashes in `source_pins.json`. The stronger
endpoint-flow floor is the present refinement. The overnight packet l
correction is retained: full linear-forest states, including singleton
incidences, are required for graph-complete cyclic boundary transfer.
Gerbicz's all-order cycle existence for N>=32 remains inherited and credited;
this packet neither reproves it nor derives it from digit legality.

## Portable reproduction

From this directory, or using absolute script paths from any directory:

```sh
PYTHONDONTWRITEBYTECODE=1 python3 check.py --output /tmp/s131-square-normal.json
PYTHONDONTWRITEBYTECODE=1 python3 -O check.py --output /tmp/s131-square-optimized.json
cmp /tmp/s131-square-normal.json /tmp/s131-square-optimized.json
python3 verify_packet.py
```

The checker uses only Python's standard library and explicit exceptions,
not optimization-sensitive assertions. It checks all eight Hamiltonian
certificates;668 marked-path cases; the long-cycle reflection words;
26,568 complete short root-word cases; all524,287 binary words through
length18;123,504 canonical numeral reflections; and the odd-square
components through order2000. Exact receipts accompany the packet.

The search used a numerical MILP solver only to discover witness cycles.
The portable checker imports no solver. Optimality at80..87 follows from
the elementary theorem plus exact witness checks, not solver messages.

`RELEASE_FILES.json` is the intended release whitelist. Solver scratch
scripts and search logs outside it are not dependencies or release files.
No publication was performed by this lane.
