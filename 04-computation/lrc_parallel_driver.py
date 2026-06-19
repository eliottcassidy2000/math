#!/usr/bin/env python3
"""
Drives the parallel sigma_2 scan for a given k using subprocess workers.
Splits the search by prefix: for first element f=1 (heaviest) split on (1,g)
for all valid g; for f>=2 split on (f,). Runs a bounded pool of workers,
combines results, writes final answer to results file. Detached-friendly.
Usage: python lrc_parallel_driver.py k
"""
import sys, os, subprocess, tempfile, glob
from fractions import Fraction
from concurrent.futures import ProcessPoolExecutor
import multiprocessing as mp

HERE = os.path.dirname(os.path.abspath(__file__))
WORKER = os.path.join(HERE, "lrc_sigma2_worker.py")

def run_worker(args):
    k, prefix, outdir = args
    tag = "_".join(map(str, prefix))
    out = os.path.join(outdir, f"k{k}_{tag}.txt")
    cmd = [sys.executable, WORKER, str(k), out] + [str(p) for p in prefix]
    subprocess.run(cmd, check=True)
    with open(out) as f:
        return f.read().strip()

def gen_prefixes(k, BOX, plen):
    """All strictly-increasing tuples of length plen that can be the smallest
    elements of a k-subset of {1..BOX}: a1<a2<...<a_plen with room for the rest,
    i.e. a_plen + (k-plen) <= BOX."""
    res = []
    def rec(prefix, start):
        if len(prefix) == plen:
            res.append(tuple(prefix)); return
        rem = k - len(prefix)
        # last fixed element a_plen must leave room: a_plen <= BOX-(k-plen)
        hi = BOX - (k - len(prefix) - 1) - 1  # next element max so that tail fits
        for v in range(start, BOX + 1):
            # ensure room for remaining (plen-len-1)+(k-plen) elements after v
            if v + (k - len(prefix) - 1) > BOX:
                break
            prefix.append(v)
            rec(prefix, v + 1)
            prefix.pop()
    rec([], 1)
    return res

def main():
    k = int(sys.argv[1])
    plen = int(sys.argv[2]) if len(sys.argv) > 2 else 2
    BOX = 3 * k
    floor = Fraction(1, k + 1)
    mediant = Fraction(2, 2 * k + 1)
    outdir = os.path.join(HERE, f"_lrc_k{k}_parts")
    os.makedirs(outdir, exist_ok=True)

    prefixes = gen_prefixes(k, BOX, plen)
    # Adaptive deepening: the heaviest branches are the near-consecutive small
    # prefixes (where the probe filter barely prunes). Expand any prefix that is
    # "tight" (all-consecutive starting at 1, or sum small) one level deeper so
    # no single DFS branch dominates wall-clock.
    def is_tight(p):
        # consecutive starting near 1, e.g. (1,2,3),(1,2,4),(2,3,4),(1,3,4)...
        return p[0] <= 3 and (p[-1] - p[0]) <= len(p) + 1
    expanded = []
    for p in prefixes:
        if is_tight(p) and len(p) < k:
            lo = p[-1] + 1
            hi = BOX - (k - len(p) - 1)
            children = [p + (g,) for g in range(lo, hi + 1)]
            if children:
                expanded.extend(children)
            else:
                expanded.append(p)
        else:
            expanded.append(p)
    prefixes = expanded

    tasks = [(k, p, outdir) for p in prefixes]
    # Shuffle so the heaviest consecutive-prefix tasks don't all run first and
    # clog every core; interleaving light tasks keeps the pool flowing.
    import random
    random.seed(12345)
    random.shuffle(tasks)
    nproc = min(mp.cpu_count(), 16)
    results = []
    with ProcessPoolExecutor(max_workers=nproc) as ex:
        for r in ex.map(run_worker, tasks, chunksize=1):
            results.append(r)

    # Combine
    best = None; bestwit = None; totsurv = 0
    for line in results:
        # prefix=(..) sigma2=X wit=Y surv=N
        # parse sigma2 and surv and wit
        parts = line.split(" sigma2=")
        rest = parts[1]
        sig_str, rest2 = rest.split(" wit=", 1)
        wit_str, surv_str = rest2.rsplit(" surv=", 1)
        totsurv += int(surv_str)
        if sig_str == "None":
            continue
        sig = Fraction(sig_str)
        if best is None or sig < best:
            best = sig; bestwit = wit_str

    resfile = os.path.abspath(os.path.join(HERE, "..", "05-knowledge",
                                           "results", f"lrc_sigma2_k{k}.out"))
    if best is None:
        out_line = (f"k={k} sigma2=None (no sub-mediant) -> sigma_2=mediant={mediant} "
                    f"floor={floor} survivors={totsurv}\n")
    else:
        g = best - floor
        out_line = (f"k={k} sigma2={best}={float(best):.6f} wit={bestwit} "
                    f"floor={floor} mediant={mediant} g={g} survivors={totsurv}\n")
    with open(resfile, "w") as f:
        f.write(out_line); f.flush()
    print(out_line)

if __name__ == "__main__":
    main()
