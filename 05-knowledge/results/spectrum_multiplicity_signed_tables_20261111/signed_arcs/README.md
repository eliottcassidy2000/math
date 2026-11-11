# A native switch and paired rebase from the prime-graph kernel

Three operation-bearing recoveries, with exact scope:

* [KERNEL_SWITCH.md](KERNEL_SWITCH.md): the minimal completion cube of the
  three-port kernel contains the two different H9 cores. Reversing one arc
  after a common marked extension gives23/25 with identical E,U,d. Arbitrary
  flanks give width3A+(A mod2)+2, supplying native gaps6,12,18,… at49/55,
  93/105,155/173. An exact two-coordinate suffix recursion simultaneously
  rebases an infinite family of these unit pairs onto C3, starting823/825.
  Every common future O/Z/J continuation is no more expensive after rebase.

* [PUNCTURED_ATLAS.md](PUNCTURED_ATLAS.md): the HP-supported table atlas turns
  missing moduli21,7 into precisely5→16 on plus and2→1 on minus. Restoring the
  first changes the first-hit channel but no cycle determinant; restoring the
  second inserts the factor1−z². The punctured mixed-sign atlas has a finite
  terminating controller, including a detour whose ten labels all have actual
  tournament witnesses. Equality with the actual HP-supported atlas is the
  still-open spectrum theorem.
* [PRIME_KERNEL.md](PRIME_KERNEL.md): the square-sum graph on the actual first25
  primes has a five-path kernel on ports3,23,47. All arbitrary subdivisions
  compile to eight cycle-parity states. Its six-cycle conflict graph is K6,
  as for the H13 odd cycles, but higher intersections and XOR dependencies
  differ. A directed specialization has cycle-packing value7; completing its
  pair relation to a tournament forces at least9 HPs.

The proofs name exact historical source files and retain their graph
conventions. These are not a new prime wheel for99·2^k+23; that arithmetic
belongs to the already-published S258 packet.

Reproduce from this directory:

    python3 -B check.py
    python3 -O -B check.py
    python3 -B switch_check.py
    python3 -O -B switch_check.py
    shasum -a 256 -c SHA256SUMS

The receipts give literal detour and switch witnesses, matrix-surgery controls, the
25-prime graph and its full kernel, all eight parity states, the directed
carrier's three cycles, and the H13 comparison. The program uses explicit
exceptions instead of assertions, so optimization does not disable checks.
