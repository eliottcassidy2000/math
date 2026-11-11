"""a8_census_run.py (AUDIT) - census of two-port kernels at orders 3..9 from the audit's own class lists and DP
(a8_census_dp.c): (1) validation of the DP against the plain walk enumeration (a2_census.c) at orders 3..7;
(2) the table of section 2.4 at orders 8 and 9; (3) all carriers (al = +-1, be = +-2^u, ga = 0, de = d 2^w)."""
import subprocess, os, sys, collections, time
HERE = os.path.dirname(os.path.abspath(__file__))
exe = os.path.join(HERE, "a8_census_dp.exe")
out = []


def P(s):
    print(s)
    out.append(s)
    sys.stdout.flush()


def odd(x):
    w = 0
    while x % 2 == 0:
        x //= 2
        w += 1
    return x, w


# (1) validation, orders 3..7: class codes differ between my list and data/cls_m.txt (other canonical form), so the
# DP is run on the packet's codes (data/cls_m.txt) here and compared line by line with a2_census_m.out (DFS).
for m in range(3, 8):
    r = subprocess.run([exe, str(m), os.path.join(HERE, "..", "data", "cls_%d.txt" % m), "1", "0", "dump"], capture_output=True, text=True)
    dp = {}
    for l in r.stdout.splitlines():
        t = l.split()
        if t[0] == "K":
            dp[(int(t[2]), int(t[3]), int(t[4]))] = tuple(int(z) for z in t[5:10])
    dfs = {}
    for l in open(os.path.join(HERE, "a2_census_%d.out" % m)):
        t = [int(z) for z in l.split()]
        if t[2] < t[3]:
            dfs[(t[1], t[2], t[3])] = (t[6], t[7], t[8], t[9], t[11])
    diff = sum(1 for k in dfs if dp.get(k) != dfs[k])
    P("order %d: DP vs plain walk enumeration on %d kernels (x < y): %d differences" % (m, len(dfs), diff))

# (2), (3)
claimed = {8: (385280, 16962, 602, 6732, 5008, 4620), 9: (13790592, 148596, 3804, 53804, 39306, 51682),
           7: (19152, 2892, 162, 1242, 946, 542), 6: (1680, 630, 60, 290, 218, 62)}
for m in (6, 7, 8, 9):
    t0 = time.time()
    NP = 11 if m == 9 else 1
    procs = [subprocess.Popen([exe, str(m), os.path.join(HERE, "aud_cls_%d.txt" % m), str(NP), str(r)], stdout=subprocess.PIPE, text=True) for r in range(NP)]
    tot = [0] * 6
    carriers = collections.defaultdict(list)
    for p in procs:
        so, _ = p.communicate()
        for l in so.splitlines():
            t = l.split()
            if t[0] == "S":
                for i in range(6):
                    tot[i] += int(t[2 + i])
            elif t[0] == "C":
                code, x, y, al, be, ga, de = int(t[2]), int(t[3]), int(t[4]), int(t[5]), int(t[6]), int(t[7]), int(t[8])
                d, w = odd(de)
                u = abs(be).bit_length() - 1
                carriers[(al, d, u, w)].append((code, x, y))
    P("order %d (own classes, own DP; %.0f s): all %d | <= 2 visits %d | join-like %d | elliptic %d | parabolic %d | hyperbolic %d ; equal to the README row: %s"
      % (m, time.time() - t0, tot[0], tot[1], tot[2], tot[3], tot[4], tot[5], tuple(tot) == claimed[m]))
    if tot[1]:
        P("   share of hyperbolic: %.1f per cent" % (100.0 * tot[5] / tot[1]))
    for key in sorted(carriers):
        al, d, u, w = key
        lst = carriers[key]
        P("   %s carrier d = %d (u = %d, w = %d): %d (class, ordered ports) occurrences; first: code %d ports (%d, %d)"
          % ("plus" if al == 1 else "MINUS", d, u, w, len(lst), lst[0][0], lst[0][1], lst[0][2]))
    P("   distinct d among plus carriers: %s ; among minus carriers: %s"
      % (sorted({k[1] for k in carriers if k[0] == 1}), sorted({k[1] for k in carriers if k[0] == -1})))
open(os.path.join(HERE, "a8_census_run.out"), "w").write("\n".join(out) + "\n")
