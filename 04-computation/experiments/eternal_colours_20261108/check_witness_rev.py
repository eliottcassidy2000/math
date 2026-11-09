# independent check of mirror-orbit witnesses: "W h n flag rev i-j ..." -> rebuild from TT_n, verify H, strongness, anti-isomorphism
import sys
sys.path.insert(0, sys.argv[1]); from check_witness_codes import H, strong, anti_iso
def build(n, revs):
    A = [[1 if i < j else 0 for j in range(n)] for i in range(n)]
    for (i, j) in revs:
        for (u, v) in {(i, j), (n - 1 - j, n - 1 - i)}:
            A[u][v] ^= 1; A[v][u] ^= 1
    return A
bad = 0; cnt = 0
for line in open(sys.argv[2]):
    if not line.startswith("W "): continue
    t = line.split(); h, n, fl = int(t[1]), int(t[2]), int(t[3])
    if len(sys.argv) > 3 and h not in set(map(int, sys.argv[3].split(","))): continue
    revs = [tuple(map(int, x.split("-"))) for x in t[5:]]
    A = build(n, revs); hh = H(A); st = strong(A); ph = anti_iso(A)
    ok = (hh == h) and (st == bool(fl)) and ph is not None
    cnt += 1; bad += (not ok)
    print(h, n, "H", hh, "strong", st, "anti-iso", ph is not None, "OK" if ok else "FAIL")
print("checked", cnt, "failures", bad)
