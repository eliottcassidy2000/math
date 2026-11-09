"""S19 audit: Session 206's parameter-capacity necessary condition applied to S19's banks (base 2).
For a seed (h0,R0,D0,hu) with d = hu < R0 the number of parameters m with h(m) <= X is (X/d) log(R0/(R0-d)) + O(sqrt X)
(session206/parameter-capacity/PARAMETER_CAPACITY.md). Seeds with d = 1 have tails inside the one-apex (root) image,
whose size is at most (log 2) X/2 + O(sqrt X) (session206/root-oz-noncofiniteness/ROOT_OBSTRUCTION.md). Sink ports
(R0 = d) have image h0 - d + d * (one-apex image), of size at most (log 2) X/(2d) + O(sqrt X).
Cofinite odd coverage needs the bracket below to be >= 1/2. Usage: python capacity.py CORES_FILE ...
"""
import sys, math
for path in sys.argv[1:]:
    cores = [tuple(map(int, l.split())) for l in open(path)]
    root = any(c[3] == 1 for c in cores)
    br = (math.log(2) / 2 if root else 0.0)
    nsink = 0
    for h0, R0, D0, d in cores:
        if d == 1: continue
        if R0 == d: br += math.log(2) / (2 * d); nsink += 1
        else: br += math.log(R0 / (R0 - d)) / d
    print(f"{path.split('/')[-1]}: {len(cores)} seeds, root-collapsing (d=1) present: {root}, sink ports with d>1: {nsink}, capacity bracket = {br:.4f} (necessary: >= 0.5)")
