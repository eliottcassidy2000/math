# Joint-entry returns: sharp budgets, a Fibonacci stratum, and phase repair

2026-11-04. **PROVED scoped all-size results / FINITE-EXACT named returns.**
Same-agent review with independent algorithms. The claim that7 and21 are
the only missing positive odd tournament Hamilton counts remains **OPEN**.

Incoming work had already continued the saved H170807 source through
13835447 to41506343/41506345. This packet extends the actual first child:

| Actual parent | Target | Endpoint order | Native complements | Prefix clock |
|---:|---:|---:|---:|---:|
|41506343|124519031=3H+2|29|62|1|
|41506343|124519033=3H+4|29|28|3|

Both programs retain all26 old vertex names and add three permanent
vertices. Each arc correction uses one temporary helper, so peak order
is30. Outer root39, active mark4, previous mark31 and spare6 remain distinct.
The +2 branch moves physical51 from the prefix into the component; the
+4 branch uses prefix[50,51,1]. Every changed old incidence and legal
Hamilton-complement guard is saved in `next_returns.json`. Four joint
port pairs, seventeen continuation queries per pair, and both F49/F343
codecs are retained. Intermediate helper states need not be strong or
preserve H. This is two more actual arrows, not closure of a source grammar.

The main all-size results are:

* **[JE1--JE2, ENTRY](ENTRY.md):** a complete grammar for J(v,w)=1, where
  J counts two-bin covers whose first path runs from v to w. The unique
  path, its transitive complement and a threshold word reconstruct every
  arc. Including coincident ports, marked-class counts have generating
  function x/(1-3x+x^2), hence F_(2n), giving an exact weighted slice of
  A000568 and a scoped minimal marked-class code.
* **[JE3--JE5, ENTRY](ENTRY.md):** a unique-start source reenters through
  the three-birth unit dial iff its marked deletion is transitive. Two
  starts admit only one exception, the marked strong four-vertex atom.
  At every order n>=5 with nontransitive deletion, the minimum starting
  budget permitting unit joint entry is exactly **three**. A uniform
  C3/transitive-prefix adapter attains it, preserving all old incidences.
* **[JR1--JR3, RESPONSE](RESPONSE.md):** all48 declared joint fields of
  that adapter have a three-state recurrence with modes1,2,3. Its integer
  row lattice is Z x6Z x38Z, with Smith factors1,2,114. Realization recovers
  the hidden coordinate modulo19. Modulo7, the F49/F343 charts omit a
  parity bit supplied for free by the retained size clock in this family.
* **[JR4, RESPONSE](RESPONSE.md):** the older bridge--sink star kernel
  applies to every member of this family. Two actual complements expose
  it with amplitude one. The three-state bulk carrier is consequently
  not licensed as a general complement carrier.
* **[RP1--RP3, REPAIR](REPAIR.md):** breaking the newborn triangle gives
  a repeatable entry repair with (h,F,q)->(h+4F,8F-q,q), explicit native
  guards and a guarded inverse. It cannot itself implement either ternary
  branch, even with arbitrary positive input/output prefix clocks.

The twelve-vertex81/27 control is retained literally: equal H462915, but
the same incoming cut256 yields1061775 and1395279. The new field-chart
hostile is separate and has different source sizes; it does not refute
a carrier that already includes size. The full joint matrix, arbitrary
complement closure, and a reusable source selector for both digits remain
open. No global Hamilton-value coverage bound or A000568 record is claimed.

Read [RECOVERY](RECOVERY.md) for the older mechanisms and incoming work,
and [AUDIT](AUDIT.md) for the proof gates and computation domains. Reproduce
after integrating the pinned Session188 skeleton source:

```powershell
python -B -X utf8 05-knowledge/results/joint_entry_return_20261104/verify.py
```

The verifier runs both normal and optimized Python, compares every result
byte for byte, and writes normalized source pins. Python and a C++20
compiler suffice; no third-party mathematical package is required.
