# c5_direct_crosscheck.py -- value-by-value cross-check, p <= PMAX, of the FLINT mod-p^2 values of c5_cw_lists_flint.py against
# two direct O(p^2) recurrences written for this audit: s_3(p) mod p^2 by the T' = S' T^2 recurrence (c1_rigidity_direct.c) and
# cg(p) mod p^2 by the rooted identity b' = T' T^(-1) (c5_connected_direct.c).  Usage: python3 c5_direct_crosscheck.py C1 C5C PMAX
# Second-pass auditor, 2026-10-21.
import sys, subprocess, re, flint
C1, C5C, PMAX = sys.argv[1], sys.argv[2], int(sys.argv[3])
src = open(__file__.replace('c5_direct_crosscheck.py', 'c5_cw_lists_flint.py')).read()
ns = {}; exec('import flint\n' + 'def primes' + src.split('# --- self-test')[0].split('def primes', 1)[1], ns)                      # reuse the FLINT functions only (series, s_mod_p2, cg_mod_p2)
primes = [p for p in ns['primes'](PMAX)]
out = subprocess.run([C5C, "3", str(PMAX)], capture_output=True, text=True, check=True).stdout
cdir = {int(m.group(1)): int(m.group(2)) for m in re.finditer(r'p=(\d+) b\(p\) mod p\^2=(\d+)', out)}
summ = out.strip().splitlines()[-1]
bad_c = sum(1 for p in primes if cdir[p] != ns['cg_mod_p2'](p))
bad_s = 0; hits3 = []
for p in primes:
    if p == 3: continue
    o = subprocess.run([C1, str(p), "1", "explicit", "3"], capture_output=True, text=True, check=True).stdout
    m = re.search(r'eps=([+-]1) base a=3 .*k=1: eps\^k s\(p\^k\)=(\d+), eps\^k D\(p\^k\)=(\d+)', o)
    eps, es, eD = int(m.group(1)), int(m.group(2)), int(m.group(3)); M = p * p
    s = es if eps == 1 else (-es) % M
    if s != ns['s_mod_p2'](3, p): bad_s += 1
    if eD == 0: hits3.append(p)
print(f"p <= {PMAX}: connected cg(p) mod p^2, direct (b' = T' T^-1) vs FLINT: {bad_c} mismatches over {len(primes)} primes; direct summary: {summ}")
print(f"p <= {PMAX}: base-3 s_3(p) mod p^2, direct (T' = S' T^2) vs FLINT: {bad_s} mismatches over {len(primes) - 1} primes; direct base-3 CW primes: {hits3}")
print("direct connected-CW primes:", [int(x) for x in re.findall(r'p=(\d+) .*CONNECTED-CW', out)])
