"""
monad-explorer-2026-06-13 — HYP-2461 decisive probe.
Does t=4 (sqrt-15, 18 units, rosette angle 29.0deg, nearest the 30deg bisector,
NOT tested by THM-493) ALSO reach the u(28)=85 crossing, or is the crossing uniquely
the t=3 sqrt-11 Moser resonance?  Contrast t=2 (sqrt-7, 12 units, 41.4deg).
Reuses the exact-integer free annealing of unit_distance_bridge_lattice_family_monad.
"""
import importlib.util, sys
spec = importlib.util.spec_from_file_location(
    "blf", "04-computation/unit_distance_bridge_lattice_family_monad.py")
blf = importlib.util.module_from_spec(spec); spec.loader.exec_module(blf)

ns = list(range(21, 31))
print("t=4 PROBE (sqrt-15, 29.0deg) and t=2 contrast (sqrt-7, 41.4deg)", flush=True)
for t, box in [(2,4), (4,5), (5,5)]:
    blf.run_lattice(t, ns, box=box, iters=110000, restarts=12)
print("\nDONE-PROBE.", flush=True)
