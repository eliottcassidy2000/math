#!/usr/bin/env python3
"""S93 lane seeds, FINISHER check (written after the interruption of the lane).

State of the interrupted re-run of b_census.py (certificates extended from odd c <= 8191 to odd c <= 65535):
  b_census_small.npz and b_census_part_{0,1,2,3}_8.npz are in the NEW format (cert arrays of length 32768),
  b_census_part_{4,5,6,7}_8.npz are still in the OLD format (cert arrays of length 4096), part 4 was being redone.
This script (i) recomputes the merged table minv (least 2-adic valuation of a slope with odd part c, per sign) from
the part files that are present and compares it with b_census_minv.npz (the table behind b_census.out);
(ii) merges the certificates that are present (plus, optionally, files b_census_part_P_32.npz produced later for the
residues P = 4..7 mod 8) and reports for which odd c below the first census gap a certificate is still missing.
Usage: python3 f_cert_coverage.py [write]
   with "write": certs_plus_full.txt / certs_minus_full.txt (format of certs_plus.txt), only certificates present.
"""
import os
import sys

import numpy as np

NC = 1 << 15
small = np.load("b_census_small.npz")
ref = np.load("b_census_minv.npz")

cum = np.full((2, 1 << 20), 255, np.uint8)
cert = dict(cert_v=np.full((2, NC), 255, np.uint8), cert_n=np.full((2, NC), 255, np.uint8),
            cert_c=np.zeros((2, NC), np.int64), cert_O=np.zeros((2, NC), np.int64), cert_M=np.zeros((2, NC), np.int64))


def merge(src, suffix=""):
    L = src["cert_v" + suffix].shape[1]
    for sg in range(2):
        sv = src["cert_v" + suffix][sg]
        sn = src["cert_n" + suffix][sg]
        dv = cert["cert_v"][sg][:L]
        dn = cert["cert_n"][sg][:L]
        better = (sv < dv) | ((sv == dv) & (sn < dn))
        for key in ("cert_v", "cert_n", "cert_c", "cert_O", "cert_M"):
            cert[key][sg][:L][better] = src[key + suffix][sg][better]


for n in range(1, 10):
    cum = np.minimum(cum, small[f"minv_{n}"])
    merge(small, f"_{n}")
print("orders <= 9: merged minv equals b_census_minv.npz minv_9:", bool((cum == ref["minv_9"]).all()))

cnt_p = np.zeros(17, np.int64)
cnt_m = np.zeros(17, np.int64)
ncand = 0
fmt = {}
done32 = {}
for p in range(8):
    fn = f"b_census_part_{p}_8.npz"
    d = np.load(fn)
    L = d["cert_v"].shape[1]
    # prefer a complete set of new-format sub-parts (p, p+8, p+16, p+24 mod 32) if present
    subs = [f"b_census_part_{q}_32.npz" for q in (p, p + 8, p + 16, p + 24)]
    if L < NC and all(os.path.exists(s) for s in subs):
        mv = np.full((2, 1 << 20), 255, np.uint8)
        cp = np.zeros(17, np.int64)
        cm = np.zeros(17, np.int64)
        nc = 0
        for s in subs:
            ds = np.load(s)
            mv = np.minimum(mv, ds["minv"])
            cp += ds["cnt_p"]
            cm += ds["cnt_m"]
            nc += int(ds["ncand"])
            merge(ds)
        same = bool((mv == d["minv"]).all()) and bool((cp == d["cnt_p"]).all()) and bool((cm == d["cnt_m"]).all()) and nc == int(d["ncand"])
        done32[p] = same
        fmt[p] = "new (4 sub-parts mod 32)"
    else:
        merge(d)
        fmt[p] = "new" if L == NC else "OLD (certificates only for odd c <= 8191)"
    cum = np.minimum(cum, d["minv"])
    cnt_p += d["cnt_p"]
    cnt_m += d["cnt_m"]
    ncand += int(d["ncand"])
print("order 10 parts:", fmt)
if done32:
    print("sub-parts mod 32 reproduce minv, the counts by |O| and the number of extensions of the old part files:", done32)
print(f"order 10: {ncand} extensions; seeds by |O| {[(k, int(x)) for k, x in enumerate(cnt_p) if x]}; anti-seeds {[(k, int(x)) for k, x in enumerate(cnt_m) if x]}")
print("orders <= 10: merged minv equals b_census_minv.npz minv_10:", bool((cum == ref["minv_10"]).all()))

for sg, name, gap in ((0, "seeds", 53559), (1, "anti-seeds", 59083)):
    mv = cum[sg][:NC].astype(np.int64)
    cv = cert["cert_v"][sg].astype(np.int64)
    c = 2 * np.arange(NC) + 1
    first_gap = int(c[mv == 255][0]) if (mv == 255).any() else None
    below = c < gap
    have = cv < 255
    print(f"{name}: least odd c without a unit pair of order <= 10 (table): {first_gap} (expected {gap}); largest minv below it: {int(mv[below].max())}")
    print(f"   odd c < {gap}: {int(below.sum())}; with a certificate present: {int((below & have).sum())}; missing: {int((below & ~have).sum())}; "
          f"certificate with v larger than the least v: {int((below & have & (cv > mv)).sum())}; largest certificate v: {int(cv[below & have].max())}")
    print(f"   odd c <= 8191 with a certificate: {int(have[:4096].sum())} of 4096;  odd c in (8191, {gap}) with a certificate: {int((have & below)[4096:].sum())} of {int(below[4096:].sum())}")
    miss = c[below & ~have]
    print(f"   first missing: {miss[:12].tolist()}")
    bad = int(((cert['cert_M'][sg][have] != (c[have] << cv[have]))).sum())
    print(f"   certificates with M != c 2^v: {bad}")
    if len(sys.argv) > 1 and sys.argv[1] == "write":
        fn = "certs_plus_full.txt" if sg == 0 else "certs_minus_full.txt"
        with open(fn, "w") as fh:
            fh.write("# c v order colcode Omask M     (pair (B, O), offset %+d, slope M = c 2^v; B = seedlib.out_from_colcode(order, colcode))\n" % (1 if sg == 0 else -1))
            for i in range(NC):
                if cv[i] < 255:
                    fh.write(f"{2 * i + 1} {cv[i]} {cert['cert_n'][sg, i]} {cert['cert_c'][sg, i]} {cert['cert_O'][sg, i]} {cert['cert_M'][sg, i]}\n")
        print(f"   written {fn}: {int(have.sum())} certificates")
