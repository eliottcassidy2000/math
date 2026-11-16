# Parents (one-vertex deletions) of the bridge T(a,b) and of T'(a,b) in the near-floor family: is m = phi_a phi_b
# an induced one-vertex parent of T(a,b) (count 2m-1)?  That would make  m -> 2m-1 -> 3m  an induced chain,
# i.e.  k -> 2k -> 3k+1  in the half-count chart n = (H-1)/2.
import importlib.util, sys
spec = importlib.util.spec_from_file_location("nf", "near_floor_family_check.py")
src = open("near_floor_family_check.py").read().split("ok = True")[0]; exec(src)
def sub(T, drop):
    idx = [i for i in range(len(T)) if i != drop]; return [[T[i][j] for j in idx] for i in idx]
for a in range(0, 4):
    for b in range(0, 4):
        if a + b == 0: continue
        pa = 1 if a == 0 else 2 ** a + 1; pb = 1 if b == 0 else 2 ** b + 1; m = pa * pb
        T = build(a, b, False); T1 = build(a, b, True)
        p0 = sorted(set(ham(sub(T, y)) for y in range(len(T)))); p1 = sorted(set(ham(sub(T1, y)) for y in range(len(T1))))
        print("a=%d b=%d m=%3d | T: H=%4d parents %s  (m a parent: %s) | T': H=%4d parents %s" % (a, b, m, 2 * m - 1, p0, m in p0, 3 * m, p1))
