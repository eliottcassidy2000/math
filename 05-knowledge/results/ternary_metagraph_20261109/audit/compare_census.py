"""Auditor's comparison of the independent class census (gen_classes.c, ext10.c) with the packet's data.

- order <= 9: strong packing vectors / H_3 values from the auditor's acls_n.txt against S86's alpha_n.txt;
- order 10: the auditor's bitmap my_t10.bin against the union of the session's t10_*.bin and witnesses_order10.txt;
- THM-4950 (5): strong spectrum below 225 and the first value with no strong realiser of order <= 10.
usage: python compare_census.py AUD_DIR S86_DIR S87_DIR PACKET_DIR
"""
import sys
import os
import glob
import numpy as np


def main():
    aud, s86, s87, pk = sys.argv[1:5]
    mine = {}
    for n in range(1, 10):
        vals = set()
        vecs = {}
        for line in open(os.path.join(aud, f"acls_{n}.txt")):
            f = line.split()
            if f[1] == "1":
                vals.add(int(f[7]))
                key = tuple(int(x) for x in f[2:6])
                vecs[key] = vecs.get(key, 0) + 1
        mine[n] = (vals, vecs)
    ok_all = True
    for n in range(3, 10):
        vecs86 = {}
        vals86 = set()
        for line in open(os.path.join(s86, f"alpha_{n}.txt")):
            f = line.split()
            if f[2] == "1":
                a = tuple(int(x) for x in f[4:8])
                vecs86[a] = vecs86.get(a, 0) + 1
                vals86.add(1 + 3 * a[0] + 9 * a[1] + 27 * a[2] + 81 * a[3])
        same = (vecs86 == mine[n][1]) and (vals86 == mine[n][0])
        ok_all &= same
        print(f"order {n}: strong classes {sum(mine[n][1].values())} (S86 {sum(vecs86.values())}); "
              f"multiset of packing vectors identical: {vecs86 == mine[n][1]}; distinct strong H_3 {len(mine[n][0])}")
    le9 = set().union(*(mine[n][0] for n in range(1, 10)))
    print(f"distinct strong H_3 values of orders <= 9 (incl. order 1): {len(le9)}")
    my10 = np.fromfile(os.path.join(aud, "my_t10.bin"), dtype=np.uint8)
    v10 = set(int(x) for x in np.nonzero(my10)[0])
    files = sorted(glob.glob(os.path.join(s87, "t10_*.bin")))
    u = None
    for fn in files:
        b = np.fromfile(fn, dtype=np.uint8)
        u = b if u is None else (u | b)
    s10 = set(int(x) for x in np.nonzero(u)[0])
    w10 = set(int(line.split()[0]) for line in open(os.path.join(pk, "witnesses_order10.txt")))
    print(f"order 10: auditor {len(v10)} values (max {max(v10)}); session bitmaps ({len(files)} slices) {len(s10)}; "
          f"witness file {len(w10)}; identical sets: {v10 == s10 == w10}")
    new10 = v10 - le9
    print(f"order-10 values not attained by strong orders <= 9: {len(new10)}")
    strong = le9 | v10
    exact = [v for v in range(1, 225) if v % 3 == 1 and v not in strong]
    print(f"v = 1 (mod 3), v < 225, with no strong realiser of order <= 10: {exact}")
    miss = [v for v in range(1, 40000) if v % 3 == 1 and v not in strong]
    first = [v for v in miss if v >= 225][0]
    print(f"first v = 1 (mod 3) >= 225 with no strong realiser of order <= 10: {first}")
    # how far does the exact range extend with H_3 - 1 >= (3/2)(H - 1) and H >= f(n) (Busch: f(11) = 225 lower bound
    # 5^((n-1)/3) for every n)?
    print(f"missing below 337 (exactness bound from H_3 >= 1 + 1.5 (H - 1), H >= 225 at n >= 11): "
          f"{[v for v in miss if v < 337]}")
    print(f"order <= 9 strong value sets match S86: {ok_all}")
    with open(os.path.join(aud, "strong_values_le10.txt"), "w") as fo:
        fo.write("\n".join(str(x) for x in sorted(strong)) + "\n")


if __name__ == "__main__":
    main()
