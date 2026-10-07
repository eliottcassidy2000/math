import subprocess
from start_state import start
T = open('st_1e9.txt').read().split()[0]
# full-histogram comparison of verify2 against verify.c (independent code path) on [1e9, 1e9+1e8)
o1 = subprocess.run(['./verify', str(10**9), str(10**9+10**8), T, open('st_1e9.txt').read().split()[1], 'A'], capture_output=True, text=True).stdout
o2 = subprocess.run(['./verify2', str(10**9), str(10**9+10**8), T, 'B'], capture_output=True, text=True).stdout
h1 = [l for l in o1.splitlines() if l.startswith('A HIST')][0].split()[2:]
h2 = [l for l in o2.splitlines() if l.startswith('B HIST')][0].split()[2:]
print('verify2 histogram == verify histogram on [1e9,1.1e9):', h1 == h2)
r1 = [l.split()[2] + ' ' + l.split()[3] for l in o1.splitlines() if 'RECORD' in l]
r2 = [l.split()[2] + ' ' + l.split()[3] for l in o2.splitlines() if 'RECORD' in l]
print('records equal:', r1 == r2, r2[-3:])
print([l for l in o2.splitlines() if 'DONE' in l])
fin = [l for l in o2.splitlines() if 'FINAL_STATE' in l][0].split()[2]
print('final state ok:', fin == str(pow(2, 10**9 + 10**8, 10**288)).zfill(288))
# spot-check 2000 random exponents' first-zero positions against python by running tiny ranges
import random
random.seed(1)
bad = 0
for _ in range(200):
    n0 = random.randrange(10**9, 10**9 + 10**8)
    t, s = start(n0)
    o = subprocess.run(['./verify2', str(n0), str(n0 + 300), t, 'C'], capture_output=True, text=True).stdout
    h = list(map(int, [l for l in o.splitlines() if l.startswith('C HIST')][0].split()[2:]))
    H = [0]*289
    for n in range(n0, n0 + 300):
        r = str(pow(2, n, 10**288)).zfill(288)[::-1]; pos = r.find('0') + 1
        H[(pos if pos else 289) - 1] += 1
    bad += (H != h)
print('random tiny-range spot checks with mismatches:', bad, 'of 200')
