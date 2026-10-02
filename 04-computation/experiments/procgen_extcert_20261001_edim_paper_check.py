#!/usr/bin/env python3
"""procgen_extcert_20261001_edim_paper_check.py -- the orchestrator's independent audit of the certificate package of
J. Allikvere, "The edge multiset dimension of hypercubes" (arXiv:2608.09983; Zenodo 10.5281/zenodo.21739363,
downloaded with the owner's permission to scratch/lrc_certs/edim/ (not committed); MD5 checked against the Zenodo record).

Only DATA is read from the package (hex masks, orbit representatives, witnesses, the rational bounds JSON, the result
summaries); none of its code is read or executed. Own engines: this script (pure Python / numpy) and
procgen_extcert_20261001_edim_q5.c (Q_5 orbit enumeration).
Checks:
  A. edim_m(Q_d) = infinity for d = 2, 3, 4: no subset of V(Q_d) resolves the edges (exhaustive).
  B. edim_m(Q_5) = infinity: all Aut(Q_5)-orbit representatives with |S| <= 16 (system checked by Burnside and by orbit
     sizes summing to C(32, a)); complements cover |S| >= 16. Plus the paper's intermediate data, recomputed:
     R4 = 14887680 resolving {0,1,2}-weightings of Q_4; 3056640 directional survivors; survivor orbits per size;
     the archived 796 representatives are exactly one per survivor Aut(Q_5)-orbit; all 796 witnesses are genuine collisions.
  C. The archived resolving sets for Q_6 .. Q_10 resolve every edge: the Q_6 descent chain 29 -> 15 (every set), and the
     final sets of sizes 63, 115, 246, 492 for d = 7..10.
  D. The 40 archived rational bounds U_d (11 <= d <= 50) are fractions < 1 agreeing with their decimals (their derivation
     is not redone: THM-4534 proves finiteness for every d >= 6 without them).
"""
import json, os, re, subprocess, sys, tempfile, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.join(HERE, '..', '..', 'scratch', 'lrc_certs', 'edim', 'unpacked')
OKS = []
def ok(c, msg):
    OKS.append(bool(c)); print(('[OK] ' if c else '[FAIL] ') + msg, flush=True)

def edges(d):
    return [(u, u | (1 << i)) for u in range(1 << d) for i in range(d) if not (u >> i) & 1]

def hist(d, e, S):
    h = [0] * (d + 1)
    for s in S:
        h[min(bin(e[0] ^ s).count('1'), bin(e[1] ^ s).count('1'))] += 1
    return tuple(h)

def resolves(d, S):
    return len({hist(d, e, S) for e in edges(d)}) == d * 2 ** (d - 1)

PC = np.array([bin(i).count('1') for i in range(1 << 16)], dtype=np.int16)
def resolves_np(d, S):
    E = edges(d)
    us = np.array([u for (u, v) in E], dtype=np.int64); vs = np.array([v for (u, v) in E], dtype=np.int64)
    H = np.zeros((len(us), d), dtype=np.int16); rows = np.arange(len(us))
    for s in S:
        H[rows, np.minimum(PC[us ^ s], PC[vs ^ s])] += 1
    return len(np.unique(H, axis=0)) == len(us)

def read_text(path):
    b = open(path, 'rb').read()
    if b[:2] in (b'\xff\xfe', b'\xfe\xff'):
        return b.decode('utf-16')
    return b.decode('utf-8', errors='replace')

def members(mask, d):
    return [v for v in range(1 << d) if mask >> v & 1]

# Aut(Q_5) in Python (independent of the C engine) for canonical forms of the archived representatives
import itertools
AUT5 = [(p, t) for p in itertools.permutations(range(5)) for t in range(32)]
def img(x, p, t):
    y = 0
    for i in range(5):
        if x >> i & 1:
            y |= 1 << p[i]
    return y ^ t
IMG = [[img(x, p, t) for x in range(32)] for (p, t) in AUT5]
def canon5(m):
    S = [x for x in range(32) if m >> x & 1]
    return min(sum(1 << g[x] for x in S) for g in IMG)

def check_A():
    for d in (2, 3, 4):
        V = 1 << d
        found = [m for m in range(1, 1 << V) if resolves(d, members(m, d))]
        ok(not found, f'A: Q_{d}: none of the {2 ** V - 1} nonempty vertex subsets resolves the edges (the paper: edim_m(Q_{d}) = infinity)')

def check_B():
    with tempfile.TemporaryDirectory() as td:
        b = os.path.join(td, 'q5')
        subprocess.run(['cc', '-O3', '-o', b, os.path.join(HERE, 'procgen_extcert_20261001_edim_q5.c')], check=True)
        r4 = subprocess.run([b, 'r4'], capture_output=True, text=True, check=True).stdout.strip()
        r = subprocess.run([b], capture_output=True, text=True)
    lines = r.stdout.strip().splitlines()
    for L in lines:
        if not L.startswith('SURV'):
            print('    ' + L)
    summ = [L for L in lines if L.startswith('Q5: ')]
    ok(r.returncode == 0 and 'resolving: 0' in summ[0] and summ[0].endswith(': yes'), 'B: ' + summ[0])
    paper = read_text(os.path.join(PKG, 'results', 'q5', 'q5_result.txt'))
    pr4 = int(re.search(r'R4 built:\s*(\d+)', paper).group(1)); psurv = int(re.search(r'directional survivors=(\d+)', paper).group(1))
    my_r4 = int(re.search(r'R4: (\d+)', r4).group(1)); my_surv = int(re.search(r'among all 2\^32 subsets: (\d+)', summ[1]).group(1))
    ok(my_r4 == pr4 == 14887680, f'B: R4 recomputed = {my_r4} resolving {{0,1,2}}-weightings of Q_4 (paper: {pr4} of 3^16)')
    ok(my_surv == psurv == 3056640, f'B: directional survivors recomputed by orbit sizes = {my_surv} (paper: {psurv}, by a 2^32 scan)')
    # survivor orbits: ours (Aut-orbit canonical reps for |S| <= 16; complements give |S| >= 17)
    surv = {}
    for L in lines:
        if L.startswith('SURV'):
            f = L.split(); surv[int(f[2], 16)] = int(f[1])
    mine = set(surv)
    for m, a in list(surv.items()):
        if a < 16:
            mine.add(canon5(~m & 0xFFFFFFFF))
    summ_txt = read_text(os.path.join(PKG, 'results', 'q5', 'q5_orbit_summary.txt'))
    pdist = {int(a): int(c) for a, c in re.findall(r'\|W\|=(\d+): (\d+) orbits', summ_txt)}
    mdist = {}
    for m in mine:
        a = bin(m).count('1'); mdist[a] = mdist.get(a, 0) + 1
    ok(mdist == pdist and sum(mdist.values()) == 796, f'B: survivor Aut(Q_5)-orbits per size {dict(sorted(mdist.items()))} = the paper\'s table (796 in all)')
    reps = [(int(m, 16), int(a)) for m, a in re.findall(r'^([0-9a-f]{8}) (\d+)\s*$', read_text(os.path.join(PKG, 'results', 'q5', 'q5_orbits.txt')), re.M)]
    theirs = [canon5(m) for m, a in reps]
    ok(len(reps) == 796 and all(bin(m).count('1') == a for m, a in reps) and len(set(theirs)) == 796 and set(theirs) == mine,
       'B: the 796 archived representatives are pairwise inequivalent and hit every survivor Aut(Q_5)-orbit exactly once')
    nwit = 0; bad = 0
    for line in read_text(os.path.join(PKG, 'results', 'q5', 'q5_orbit_witnesses.txt')).splitlines():
        mm = re.match(r'([0-9a-f]{8}) (\d+) e=\((\d+),(\d+)\) f=\((\d+),(\d+)\) hist=\(([\d, ]+)\)', line)
        if not mm:
            continue
        nwit += 1
        W = members(int(mm.group(1), 16), 5); e = (int(mm.group(3)), int(mm.group(4))); f = (int(mm.group(5)), int(mm.group(6)))
        h = tuple(int(x) for x in mm.group(7).split(','))
        good = (len(W) == int(mm.group(2)) and bin(e[0] ^ e[1]).count('1') == 1 and bin(f[0] ^ f[1]).count('1') == 1
                and set(e) != set(f) and hist(5, e, W) == h[:6] == hist(5, f, W) and all(x == 0 for x in h[6:]))
        bad += not good
    ok(nwit == 796 and bad == 0, f'B: all {nwit} archived collision witnesses are genuine (two distinct edges, equal stated histograms)')
    return pdist

def check_C():
    D = os.path.join(PKG, 'results', 'certificates')
    chain = re.findall(r'\|W\|=(\d+): RESOLVING, W=([0-9a-fA-F]+)', read_text(os.path.join(D, 'q6_min_result.txt')))
    sizes = [int(a) for a, w in chain]
    good = sizes == list(range(29, 14, -1)) and all(len(members(int(w, 16), 6)) == int(a) and resolves(6, members(int(w, 16), 6)) for a, w in chain)
    ok(good, f'C: Q_6: all {len(chain)} sets of the archived descent chain (sizes 29 -> 15) resolve the 192 edges; the last is 0x{chain[-1][1]} (= THM-4525\'s set)')
    q = re.search(r'W=([0-9a-fA-F]+) \(\|W\|=(\d+)\)', read_text(os.path.join(D, 'q6_quick.txt')))
    S = members(int(q.group(1), 16), 6)
    ok(len(S) == int(q.group(2)) == 29 and resolves(6, S), 'C: Q_6: q6_quick.txt holds a resolving set of size 29 (the README attributes the size-15 set to this file; it is the last line of q6_min_result.txt)')
    m14 = read_text(os.path.join(D, 'q6_min14_result.txt'))
    print('    q6_min14_result.txt: ' + m14.strip())
    for d, fn, size in ((7, 'q7_result.txt', 63), (8, 'q8_result.txt', 115), (9, 'q9_result.txt', 246), (10, 'q10_result.txt', 492)):
        txt = read_text(os.path.join(D, fn))
        mm = re.search(r'MASK n=%d: ([0-9a-fA-F]+)' % d, txt) or re.search(r'W=([0-9a-fA-F]+)', txt)
        S = members(int(mm.group(1), 16), d)
        ok(len(S) == size and resolves_np(d, S), f'C: Q_{d}: the archived set of size {len(S)} resolves all {d * 2 ** (d - 1)} edges')

def check_D():
    cert = json.load(open(os.path.join(PKG, 'certificates', 'rational_certificates_11_50.json')))
    keys = sorted(int(k) for k in cert)
    good = keys == list(range(11, 51)) and all(int(cert[str(d)]['num']) < int(cert[str(d)]['den']) for d in keys)
    vals = {d: int(cert[str(d)]['num']) / int(cert[str(d)]['den']) for d in keys}
    agree = all(abs(vals[d] - float(cert[str(d)]['decimal_approx'])) < 1e-6 for d in keys)
    dmax = max(vals, key=vals.get)
    ok(good and agree, f'D: all 40 rational bounds U_d (11 <= d <= 50) are < 1 (largest U_{dmax} = {vals[dmax]:.4f}); derivation not redone (non-load-bearing after THM-4534)')

def main():
    t0 = time.time()
    print('==== A ====', flush=True); check_A()
    print('==== B ====', flush=True); check_B()
    print('==== C ====', flush=True); check_C()
    print('==== D ====', flush=True); check_D()
    print(f'elapsed {time.time() - t0:.0f} s')
    print('ALL CHECKS PASSED' if all(OKS) else 'SOME CHECK FAILED')

if __name__ == '__main__':
    main()
