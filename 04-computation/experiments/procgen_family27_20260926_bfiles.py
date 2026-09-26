#!/usr/bin/env python3
"""procgen_family27_20260926_bfiles.py

OEIS record sequences headed by 27, re-derived from their starting values with exact
big-integer orbits of the shortcut map T(x) = x/2, (3x+1)/2 (trajectory stopped at 1).

  A006877 delay records (standard steps), values A006878
  A006884 path records (standard map maximum), values A006885
  A060412 glide / dropping-time records, values A060413 (T-steps) and A217934 (standard)

The b-files are fetched once from oeis.org with a generic User-Agent into the scratch cache
(no personal data in any request) and their sha256 is reported by the runner.
"""
import hashlib
import json
import math
import os
import subprocess

UA = "Mozilla/5.0 (research; math-repo)"
FILES = {
    "b006877.txt": "https://oeis.org/A006877/b006877.txt",
    "b006878.txt": "https://oeis.org/A006878/b006878.txt",
    "b006884.txt": "https://oeis.org/A006884/b006884.txt",
    "b006885.txt": "https://oeis.org/A006885/b006885.txt",
    "b060412.txt": "https://oeis.org/A060412/b060412.txt",
    "A060413.json": "https://oeis.org/search?q=id:A060413&fmt=json",
    "A217934.json": "https://oeis.org/search?q=id:A217934&fmt=json",
}


def fetch_all(cache):
    os.makedirs(cache, exist_ok=True)
    shas = {}
    for name, url in FILES.items():
        path = os.path.join(cache, name)
        if not os.path.exists(path) or os.path.getsize(path) == 0:
            # curl with a generic User-Agent; no personal data is sent
            subprocess.run(["curl", "-s", "-f", "-L", "-A", UA, "-o", path, url], check=True, timeout=120)
        with open(path, "rb") as f:
            shas[name] = hashlib.sha256(f.read()).hexdigest()
    return shas


def read_bfile(path):
    out = []
    for line in open(path):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        a, b = line.split()[:2]
        out.append((int(a), int(b)))
    return out


def read_json_data(path):
    d = json.load(open(path))
    if isinstance(d, dict):
        d = d.get("results", d)
    return [int(x) for x in d[0]["data"].split(",")]


def T(x):
    return x >> 1 if x % 2 == 0 else (3 * x + 1) >> 1


def orbit(n):
    o = [n]
    while o[-1] != 1:
        o.append(T(o[-1]))
    return o


def stats(n):
    """Exact statistics of n (n >= 1)."""
    o = orbit(n)
    dT = len(o) - 1
    ones = sum(1 for v in o[:-1] if v & 1)
    tT = max(o)
    ipeak = o.index(tT)
    tS = max([n] + [3 * v + 1 for v in o[:-1] if v & 1])
    if n >= 2:
        g = next(j for j in range(1, len(o)) if o[j] < n)
        gC = sum(1 for v in o[:g] if v & 1)
    else:
        g, gC = 0, 0
    return {"n": n, "orbit": o, "dT": dT, "ones": ones, "dS": dT + ones, "tT": tT, "tS": tS,
            "ipeak": ipeak, "glT": g, "glS": g + gC}


def genealogy(records):
    """For each record (in order) find earlier records whose PRE-PEAK orbit segment the record's
    orbit meets (merge point strictly before the earlier record's peak index).  Returns a list of
    (n, parent or None, merge value)."""
    info = []
    pos_maps = []
    for r in records:
        s = stats(r)
        pm = {v: i for i, v in enumerate(s["orbit"])}
        pos_maps.append((s, pm))
    for k, (sk, _) in enumerate(pos_maps):
        parent, mval = None, None
        for j in range(k - 1, -1, -1):
            sj, pmj = pos_maps[j]
            # first point of r_k's orbit lying on r_j's orbit
            for v in sk["orbit"]:
                if v in pmj:
                    if pmj[v] < sj["ipeak"]:
                        parent, mval = sj["n"], v
                    break
            if parent is not None:
                break
        info.append((sk["n"], parent, mval))
    return info


def harmonic(N):
    """H_N for large N (Euler-Maclaurin)."""
    if N < 10 ** 6:
        return sum(1.0 / i for i in range(1, N + 1))
    lnN = math.log(N)
    return lnN + 0.5772156649015329 + 1 / (2 * N) - 1 / (12 * N * N)
