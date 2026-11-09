"""AUDIT2: direct test of HYP-9610 (A) containment at n = 10, 11 (and the six unwitnessed values), independent of
cover.py's closure: every odd k <= M0(n) other than 3, 7, 13, M0(n) - 2 is searched for with a2_e0hunt.exe (simulated
annealing in the complement), and every hit is re-verified by a2_check.exe (order, d, EVEN, no universal vertex).
usage: python a2_e0full.py n [restarts steps]
"""
import sys, os, subprocess, re
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from a2_lib import *
n = int(sys.argv[1]); restarts = int(sys.argv[2]) if len(sys.argv) > 2 else 4; steps = int(sys.argv[3]) if len(sys.argv) > 3 else 20000
M0 = 2**n - n - 1 if n % 2 == 0 else 2**n - n - 4
targets = [k for k in range(1, M0 + 1, 2) if k not in (3, 7, 13, M0 - 2)]
HUNT = os.path.join(HERE, 'a2_e0hunt.exe')
found = {}
remaining = targets[:]
seed = 1
for rnd_ in range(6):
    if not remaining: break
    chunks = [remaining[i:i + 50] for i in range(0, len(remaining), 50)]
    for ch in chunks:
        out = subprocess.run([HUNT, str(n), str(seed), str(restarts * len(ch)), str(steps)] + list(map(str, ch)),
                             capture_output=True, text=True, check=True).stdout
        seed += 1
        for line in out.splitlines():
            m = re.match(r'FOUND n=(\d+) d=(\d+) naut=(\d+) G-masks: (.*)', line)
            if m:
                found[int(m.group(2))] = (n, list(map(int, m.group(4).split())))
    remaining = [k for k in targets if k not in found]
    print("round %d: found %d of %d; remaining %d %s" % (rnd_, len(found), len(targets), len(remaining), remaining[:20]), flush=True)
    restarts *= 4
ks = sorted(found)
res = run_check([(found[k], -1) for k in ks])
bad = [(k, r) for k, r in zip(ks, res) if not (r[0] == k and r[1] == 'EVEN' and r[3] == 0)]
print("n=%d: M0 = %d; targets %d; realised with re-verified witnesses %d; verification failures %d; not realised: %s"
      % (n, M0, len(targets), len(ks) - len(bad), len(bad), [k for k in targets if k not in found]))
with open(os.path.join(HERE, 'a2_e0full_n%d_witnesses.txt' % n), 'w') as f:
    for k in ks:
        f.write("%d %s\n" % (k, ' '.join(map(str, found[k][1]))))
