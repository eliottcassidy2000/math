# Repeated ternary-offset repairs with a field clock and integer carries

**PROVED scoped mechanisms / FINITE-EXACT constructions. Same-agent review
with independent complete algorithms. The claim that 7 and 21 are the only
missing odd Hamilton counts remains OPEN.**

The failed 95/97 outputs can now be repaired in a precise weaker sense:
their successors have the correct integer Hamilton counts and the correct
F49 residue clock, with explicit integer carries retained. The stronger
old requirement, exact commutation with every q2 prefix, remains impossible
through the complete seven-vertex atom bank.

The results are:

- **Four explicit repairs.** Adding just two seven-vertex source charts
  to the old 65 supplies 95->287,289 and 97->293,295. Their full states obey
  u=3v+j(2,-1,0)+7kappa with kappa_H=0. Every source relation is literal.
- **Repeated branching through depth four.** The full bank of strong
  charts through seven realizes all 30 offset branches from H31. Endpoint
  orders are at most 12. All branches are compiled into legal births,
  path complements and named deletions, with peak order 13. This is a
  finite tree, not all-size closure.
- **An exact obstruction.** None of the 74 H95 states or 84 H97 states in
  that full bank supports both original exact affine lifts. Three q2
  observations already force those lifts; arbitrary reencoding cannot
  evade this within the bank. A further complete eight-vertex boundary
  test also fails, with its narrower operation scope stated in THEORY.
- **The missing memory is an integer carry.** Two actual order-eight
  sources have H293 and identical residue data but next q2 counts 3557
  and 3865. The smallest closed integer module forced by the four repair
  differences is 7{(2u,x,y)}. Its exact modular observation quotients are
  proved in CARRIES, with all five Q4 phases retained.
- **A larger infinite offset family.** The previously omitted s78 chart
  preserves the exact correction tau, extending the prefix alphabet to
  {b,q2,s78}. It also sends the initial F49 phase to zero, so the clock
  must retain that singular case rather than assigning it a false phase.

[THEORY.md](THEORY.md) specifies the actual operations and gives the finite
repair table. [CARRIES.md](CARRIES.md) proves the all-size carry module,
Smith quotients, growing-error obstruction and additional infinite family.
[AUDIT.md](AUDIT.md) records the universes, independent checks and remaining
proof obligations. [RECOVERY.md](RECOVERY.md) credits the incoming inverse,
charge, native-compiler and adaptive-port mechanisms.

The finite 67-letter repair bank already fails four of its next required
residue fibres. The larger bank succeeds through four levels; no later
level or complete prime lane is asserted. The next theorem must guarantee
a realized child and return to a closed source family, not just solve the
residue equations.

Reproduce from the repository root (Python standard library and g++17):

```powershell
python -B -X utf8 05-knowledge/results/ternary_repair_clock_20261103/research.py
python -B -O -X utf8 05-knowledge/results/ternary_repair_clock_20261103/research.py
python -B -X utf8 05-knowledge/results/ternary_repair_clock_20261103/seal.py
```

Generated full atlases, state sets and executables remain in a dedicated
temporary directory. The committed receipt retains every selected source,
carry, native command, audit transcript and input/output hash. No A000568
coefficient or unrestricted scalar coverage record is extended here.
