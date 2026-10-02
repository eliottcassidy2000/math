#!/usr/bin/env python3
"""procgen_shave4_20261001_orchestrator_check.py -- the orchestrator's independent audit of the shave4 lane (2026-10-01).

Written from the definitions; the lane's code was not read. Uses the orchestrator's own C embedding engine
procgen_shave4_20261001_orchestrator_check.c (compiled into a temp dir) and nauty's gentourng.
Checks:
  A. D69: path + span-3 has an odd number of embeddings in every tournament for n = 4, 5, 6; at n = 7 it is even for
     174 of the 456 classes and absent from exactly 2; the near-transitive witness has exactly 2 copies.
  B. D_n = P_n + (0, n-2) + (1, n-1) is Redei (odd embedding count in every tournament) for n = 5, 7, 9 (all classes).
  C. u(9) >= 14: the lane's 14-arc graph embeds in all 191536 tournament classes on 9 vertices.
  D. Census of Redei classes (spanning oriented graphs with an odd embedding count in every n-tournament) for
     n = 1..5: 1, 1, 2, 5, 21; the largest Redei graphs have 0, 1, 2, 4, 6 arcs.
"""
import itertools, os, subprocess, sys, tempfile, time

OKS = []
def ok(c, msg):
    OKS.append(bool(c)); print(('[OK] ' if c else '[FAIL] ') + msg, flush=True)

HERE = os.path.dirname(os.path.abspath(__file__))

def path_span3(n):
    return ' '.join([f'{i}>{i+1}' for i in range(n - 1)] + [f'{i}>{i+3}' for i in range(n - 3)])

def Dn(n):
    return ' '.join([f'{i}>{i+1}' for i in range(n - 1)] + [f'0>{n-2}', f'1>{n-1}'])

def run(E, mode, n, spec, stdin_cmd=None, text=None):
    if stdin_cmd:
        p1 = subprocess.Popen(stdin_cmd, stdout=subprocess.PIPE)
        r = subprocess.run([E, mode, str(n), spec], stdin=p1.stdout, capture_output=True, text=True); p1.wait()
    else:
        r = subprocess.run([E, mode, str(n), spec], input=text, capture_output=True, text=True)
    return r.stdout

def census(n):
    """Redei classes among all spanning oriented graphs on n vertices (brute force, n <= 5)"""
    P = [(i, j) for i in range(n) for j in range(i + 1, n)]
    perms = list(itertools.permutations(range(n)))
    tours = subprocess.run(['gentourng', '-q', str(n)], capture_output=True, text=True).stdout.split() if n >= 2 else ['']
    Tadj = []
    for s in tours:
        adj = [[False] * n for _ in range(n)]
        for k, (i, j) in enumerate(P):
            if s[k] == '1': adj[i][j] = True
            else: adj[j][i] = True
        Tadj.append(adj)
    seen = set(); redei = 0; best = -1
    for code in range(3 ** len(P)):
        arcs = []
        x = code
        for (i, j) in P:
            r = x % 3; x //= 3
            if r == 1: arcs.append((i, j))
            elif r == 2: arcs.append((j, i))
        canon = min(tuple(sorted((g[a], g[b]) for (a, b) in arcs)) for g in perms)
        if canon in seen:
            continue
        seen.add(canon)
        good = True
        for adj in Tadj:
            e = sum(1 for g in perms if all(adj[g[a]][g[b]] for (a, b) in arcs))
            if e % 2 == 0:
                good = False; break
        if good:
            redei += 1; best = max(best, len(arcs))
    return redei, best, len(seen)

def main():
    t0 = time.time()
    with tempfile.TemporaryDirectory() as d:
        E = os.path.join(d, 'emb')
        subprocess.run(['cc', '-O2', '-o', E, os.path.join(HERE, 'procgen_shave4_20261001_orchestrator_check.c')], check=True)
        res = {n: run(E, 'parity', n, path_span3(n), ['gentourng', '-q', str(n)]).strip() for n in (4, 5, 6, 7)}
        for n in (4, 5, 6):
            ok('even=0 zero=0' in res[n], f'A: path + span-3, n = {n}: {res[n]}')
        ok('odd=282 even=174 zero=2' in res[7], f'A: path + span-3, n = 7: {res[7]}')
        w = ''.join('1' if (i < j and not (i == 0 and j == 6) and not (i == 3 and j == 6)) else '0'
                    for i in range(7) for j in range(i + 1, 7))
        out = run(E, 'count', 7, path_span3(7), text=w + '\n').split()
        ok(out[-1] == '2', f'A: TT7 with 0->6 and 3->6 reversed has exactly {out[-1]} copies of path + span-3')
        for n in (5, 7, 9):
            r = run(E, 'parity', n, Dn(n), ['gentourng', '-q', str(n)]).strip()
            ok('even=0 zero=0' in r, f'B: D_{n} Redei: {r}')
        r = run(E, 'exist', 9, '0>1 0>3 0>6 1>4 1>5 1>7 2>3 2>5 3>4 3>7 4>7 5>6 5>8 6>8', ['gentourng', '-q', '9']).strip()
        ok('fail=0' in r and 'classes=191536' in r, f'C: the 14-arc graph embeds in every 9-class: {r}')
    exp = {1: (1, 0), 2: (1, 1), 3: (2, 2), 4: (5, 4), 5: (21, 6)}
    for n in range(1, 6):
        rc, best, total = census(n)
        ok((rc, best) == exp[n], f'D: n = {n}: {rc} Redei classes among {total} oriented graphs; largest has {best} arcs')
    print(f'elapsed {time.time() - t0:.0f} s')
    print('ALL CHECKS PASSED' if all(OKS) else 'SOME CHECK FAILED')

if __name__ == '__main__':
    main()
