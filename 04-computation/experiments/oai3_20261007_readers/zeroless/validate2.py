import subprocess
from start_state import start
def run(n0, n1, tag):
    t, s = start(n0)
    out = subprocess.run(['./verify2', str(n0), str(n1), t, tag], capture_output=True, text=True).stdout
    hist = None; fin = None
    for line in out.splitlines():
        if line.startswith(tag + ' HIST'): hist = list(map(int, line.split()[2:]))
        if line.startswith(tag + ' FINAL_STATE'): fin = line.split()[2]
    return out, hist, fin
# 1) exact python comparison on small range
n0, n1 = 957, 30000
out, hist, fin = run(n0, n1, 'V2')
W = 288; H = [0]*(W+1); x = pow(2, n0)
for n in range(n0, n1):
    r = str(x)[::-1]; pos = r.find('0') + 1
    if pos == 0 or pos > W: pos = W + 1
    H[pos-1] += 1; x <<= 1
print('small-range histogram matches python:', H == hist)
print('final state matches pow(2,n1,10^288):', fin == str(pow(2, n1, 10**288)).zfill(288))
# 2) compare with verify.c on [1e9, 1e9+1e8)
out, hist2, fin2 = run(10**9, 10**9 + 10**8, 'V3')
print('final state (1e9+1e8) matches:', fin2 == str(pow(2, 10**9 + 10**8, 10**288)).zfill(288))
print([l for l in out.splitlines() if 'DONE' in l])
