"""a7_run_all_five.py (AUDIT) - run a7_enum.exe (fast mode) on ALL 2835 feasible 5-sets of planes, 11 processes,
and summarise: least number of even 4-sets next to each 5-set, by orbit; global histogram.  No symmetry argument
is used here: every feasible 5-set is enumerated completely."""
import subprocess, os, sys, collections, time
HERE = os.path.dirname(os.path.abspath(__file__))
NP = 11
t0 = time.time()
procs = []
for r in range(NP):
    fo = open(os.path.join(HERE, "a7_enum_five_all_%d.part" % r), "w")
    procs.append((subprocess.Popen([os.path.join(HERE, "a7_enum.exe"), os.path.join(HERE, "a7_problems_five.txt"), "0", str(NP), str(r)], stdout=fo), fo))
for p, fo in procs:
    p.wait()
    fo.close()
mins = collections.Counter()
hist = collections.Counter()
nprob = 0
points = 0
for r in range(NP):
    for l in open(os.path.join(HERE, "a7_enum_five_all_%d.part" % r)):
        t = l.split()
        label, k, mn = t[0], int(t[1]), int(t[2])
        orbit = label.split("_")[0]
        mins[(orbit, k, mn)] += 1
        nprob += 1
        for item in t[4:]:
            e, c = item.split(":")
            hist[int(e)] += int(c)
            points += int(c)
out = []
out.append("feasible 5-sets enumerated completely: %d; violation patterns visited: %d; wall time %.0f s" % (nprob, points, time.time() - t0))
out.append("(orbit label, dimension of the pattern space, least number of even 4-sets) -> number of 5-sets: %s" % sorted(mins.items()))
out.append("least number of even 4-sets over all tournaments with at least five odd planes: %d" % min(hist))
out.append("values of E below 30 occurring next to five odd planes (with multiplicity over 5-sets): %s" % sorted((e, c) for e, c in hist.items() if e < 30))
print("\n".join(out))
open(os.path.join(HERE, "a7_enum_five_all.out"), "w").write("\n".join(out) + "\n")
for r in range(NP):
    os.remove(os.path.join(HERE, "a7_enum_five_all_%d.part" % r))
