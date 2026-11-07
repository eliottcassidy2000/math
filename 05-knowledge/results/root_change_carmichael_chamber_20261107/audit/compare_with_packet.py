"""Compare the auditor's raw outputs (this folder) with the packet's outputs (parent folder). Standard library only.
usage (from this folder): python3 -I compare_with_packet.py > COMPARISON.txt
Every line is 'MATCH' or 'MISMATCH' followed by the quantity compared.
"""
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
PKT = HERE.parent
S13 = PKT.parent / 'finite_oz_bank_poisson_root_changes_20261106'


def rd(p):
    return Path(p).read_bytes().decode('utf-8').replace('\r\n', '\n')


out = []


def check(label, a, b):
    out.append('%s  %s' % ('MATCH   ' if a == b else 'MISMATCH', label))
    if a != b:
        out.append('    auditor: %r\n    packet:  %r' % (a, b))


# 1. chamber at 2^20: windows, singletons, counts
mine = rd(HERE / 'chamber_audit_2p20.txt')
theirs = rd(PKT / 'VERIFY_CHAMBER_2p20.txt')
pat = re.compile(r'K=(\d) \[(\d+),(\d+)\): n=(\d+) min=(\d+).*?c0\.\.8=([\d,]+)')
A = {(m[0], m[1]): (m[3], m[4], m[5]) for m in pat.findall(mine)}
B = {(m[0], m[1]): (m[3], m[4], m[5]) for m in pat.findall(theirs)}
check('2^20: all %d window histograms (n, min, c0..8) for K = 0, 1, 2' % len(B), {k: A.get(k) for k in B}, B)
for K in '012':
    a = re.search(r'K=%s singletons (\d+);' % K, mine).group(1)
    b = re.search(r'K=%s: \(C\) parity failures 0; singletons N>1: (\d+)' % K, theirs).group(1)
    check('2^20: K=%s singleton count above 1' % K, a, b)
    a = re.search(r'K=%s singleton values N>1 \(first 60\) with the counted word u = wZ:(.*)' % K, mine).group(1).split()
    b = re.search(r'K=%s singleton values \(first 40\) with word w \(value of wZ\):(.*)' % K, theirs).group(1).split()
    check('2^20: K=%s first %d singleton values and words' % (K, len(b)), a[:len(b)], b)
    a = re.search(r'K=%s: .*?chamber words with value <= X: (\d+)' % K, mine).group(1)
    b = re.search(r'K=%s: .*?chamber words \(dagger-fixed, exactly K J, value <= X\): (\d+)' % K, theirs).group(1)
    check('2^20: K=%s chamber words' % K, a, b)
check('2^20: prefixes with #J <= 2 and F <= 2^20',
      re.search(r'visited prefixes .*?: (\d+)', mine).group(1),
      re.search(r'on (\d+) visited states', theirs).group(1))
a = re.search(r'K=2: .*?odd M_K at N>=3: (\d+) \(above 2\^15: (\d+)\)', mine).groups()
b = re.search(r'odd M_K: (\d+) \(of which above 2\^15: (\d+)\)', rd(PKT / 'explore' / 'odd_multiplicity_counts.txt')).groups()
check('2^20: odd M_2 count (N >= 3) and above 2^15', a, b)
p24 = HERE / 'chamber_audit_K2_2p24.txt'
if p24.exists():
    a = re.search(r'K=2: .*?odd M_K at N>=3: (\d+) \(above 2\^15: (\d+)\)', rd(p24)).groups()
    b = re.findall(r'odd M_K: (\d+) \(of which above 2\^15: (\d+)\)', rd(PKT / 'explore' / 'odd_multiplicity_counts.txt'))[1]
    check('2^24: odd M_2 count (N >= 3) and above 2^15', a, tuple(b))

# 2. K = 0, 1 singleton windows to 2^26 (and full window histograms)
m26 = rd(HERE / 'chamber_audit_K01_2p26.txt')
for K, f in ((0, 'rc_K0_2p26_windows.txt'), (1, 'rc_K1_2p26_windows.txt')):
    t = rd(PKT / 'explore' / f)
    mine_w = {r[0]: (r[1], r[2], r[3]) for r in re.findall(r'K=%d \[(\d+),\d+\): n=(\d+) min=\d+ mean255=([\d.]+) exactmean=[\d.]+ sat255=\d+ c0\.\.8=([\d,]+)' % K, m26)}
    pk_w = {r[0]: (r[1], r[2], r[3]) for r in re.findall(r'\[(\d+),\d+\): n=(\d+) mean=([\d.]+) min=\d+ min_mod3=\([\d,]+\) sat255=\d+ c0\.\.8=([\d,]+)', t)}
    check('2^26: K=%d all %d window rows (n, mean of min(M,255), c0..8)' % (K, len(pk_w)), {k: mine_w.get(k) for k in pk_w}, pk_w)
    s = re.search(r'K=%d singletons (\d+);' % K, m26).group(1)
    tot = sum(int(v[2].split(',')[1]) for v in pk_w.values()) - 1
    check('2^26: K=%d singletons above 1 (auditor count vs sum of the packet windows)' % K, int(s), tot)
check('2^26: K=0 singletons above 1, auditor count vs PROOF section 4 text "3,983"', int(re.search(r'K=0 singletons (\d+);', m26).group(1)), 3983)

# 3. census to 1e8 vs S13 and rc_1e9 (windows below 2^26), the low list
c8 = rd(HERE / 'census_K2_1e8.txt')
s13 = rd(S13 / 'root_changes_1e8.txt')
r9 = rd(PKT / 'rc_1e9.txt')
lev = lambda t: re.findall(r'<=(\d) root changes: odd misses (\d+), largest (\d+)', t)
check('1e8: per-level misses and largest, k = 0, 1, 2 (vs S13 root_changes_1e8.txt)', lev(c8), lev(s13))
kw = lambda t: {m[0]: m[1] for m in re.findall(r'^  \[(\d+),\d+\):((?: k\d:\d+)+)\s*$', t, re.M)}
check('1e8: per-window k-misses (vs S13)', kw(c8), kw(s13))
mm = {m[0]: (m[1], m[2], m[3], m[4], m[5], m[6], m[7], m[8]) for m in re.findall(
    r'\[(\d+),\d+\): n=(\d+) min=(\d+) min_mod3=\((\d+),(\d+),(\d+)\) satcap=\d+ mean255=([\d.]+) sat255=(\d+) c0\.\.8=([\d,]+)', c8)}
pm = {m[0]: (m[1], m[3], m[4], m[5], m[6], m[2], m[7], m[8]) for m in re.findall(
    r'^  \[(\d+),\d+\): n=(\d+) mean=([\d.]+) min=(\d+) min_mod3=\((\d+),(\d+),(\d+)\) sat255=(\d+) c0\.\.8=([\d,]+)', r9, re.M)}
low = [k for k in pm if int(k) < 2 ** 26]
check('1e8: %d multiplicity windows below 2^26 (n, min, min_mod3, mean255, sat255, c0..8) vs rc_1e9.txt' % len(low),
      {k: mm.get(k) for k in low}, {k: pm[k] for k in low})
a = re.search(r'odd N > 10\^5 with M <= 2:(.*)', c8).group(1).split()
b = [x for x in re.search(r'odd N > 10\^5 with multiplicity <= 2:(.*)', r9).group(1).split() if int(x.split(':')[0]) <= 10 ** 8]
check('1e8: odd N > 1e5 with M_2 <= 2', a, b)
a = re.search(r'misses \(<=2\) from 100001:(.*)', c8).group(1).split()
b = re.search(r'misses \(<=2\) above 10\^5:(.*)', r9).group(1).split()
check('F_<=2 misses above 1e5 (auditor to 1e8; packet to 1e9)', a, b)

# 4. (2^26, 1e9] with 4-bit counters: union misses, minima, residue minima (capped at 15), c0..8
big = rd(HERE / 'census_K2_2p26_to_1e9.txt')
check('(2^26,1e9]: odd N with M_2 = 0', re.search(r'union \(M_2 = 0\): odd misses (\d+)', big).group(1), '0')
bw = {m[0]: (m[1], m[2], m[3], m[4], m[5], m[6]) for m in re.findall(
    r'\[(\d+),\d+\): n=(\d+) min=(\d+) min_mod3=\((\d+),(\d+),(\d+)\) satcap=\d+ c0\.\.8=([\d,]+)', big)}
cap = lambda x: str(min(int(x), 15))
pw = {k: (pm[k][0], pm[k][1], cap(pm[k][2]), cap(pm[k][3]), cap(pm[k][4]), pm[k][7]) for k in pm if int(k) >= 2 ** 26}
check('(2^26,1e9]: 4 windows (n, min, min_mod3 capped at 15, c0..8) vs rc_1e9.txt', bw, pw)
nm = sorted(set(re.search(r'odd N > 2\^22 with M <= 12 \(\d+ values\):(.*)', c8).group(1).split())
            | set(re.search(r'odd N > 2\^22 with M <= 12 \(\d+ values\):(.*)', big).group(1).split()),
            key=lambda x: int(x.split(':')[0]))
pk = re.search(r'odd N > 2\^22 with multiplicity <= 12:(.*)', rd(PKT / 'explore' / 'rc_low_2p28.txt')).group(1).split()
check('(2^22,2^28): values with M_2 <= 12 (packet explore/rc_low_2p28.txt)', [x for x in nm if int(x.split(':')[0]) < 2 ** 28], pk)
out.append('INFO      values with M_2 <= 12 in [2^28, 1e9] (beyond the packet): %s' % [x for x in nm if int(x.split(':')[0]) >= 2 ** 28])

# 5. depth tables
d7 = rd(PKT / 'depth_K7_1e6.txt')
for k in range(3, 8):
    t = rd(HERE / ('census_K%d_1e6.txt' % k))
    check('1e6: <=%d root changes misses/largest' % k, lev(t)[k], lev(d7)[k])
check('1e6: <=7 miss list', re.search(r'misses \(<=7\) from 1:(.*)', rd(HERE / 'census_K7_1e6.txt')).group(1).split(),
      re.search(r'misses \(<=7\) from 1:(.*)', d7).group(1).split())
c1p = [(m[0], m[1]) for m in re.findall(r'^  \[(\d+),\d+\): n=\d+ mean=[\d.]+ min=\d+ min_mod3=\([\d,]+\) sat255=\d+ c0\.\.8=\d+,(\d+)', d7, re.M) if m[0] != '1' and m[1] != '0']
c1a = [(m[0], m[1]) for m in re.findall(r'^  \[(\d+),\d+\): n=\d+ min=\d+ min_mod3=\([\d,]+\) satcap=\d+ mean255=[\d.]+ sat255=\d+ c0\.\.8=\d+,(\d+)', rd(HERE / 'census_K7_1e6.txt'), re.M) if m[0] != '1' and m[1] != '0']
check('1e6: windows above [1,2) containing an N with M_7(N) = 1 (HYP-9494 (C) says none for K >= 2)', c1a, c1p)
out.append('INFO      that window holds N = 1897, M_7(1897) = 1 (k_singletons_2048.out): HYP-9494 (C) is false for K = 7')
d3 = rd(PKT / 'depth_K3_4e6.txt')
check('4e6: per-level misses k = 0..3', lev(rd(HERE / 'census_K3_4e6.txt')), lev(d3))
check('4e6: <=3 miss list', re.search(r'misses \(<=3\) from 1:(.*)', rd(HERE / 'census_K3_4e6.txt')).group(1).split(),
      re.search(r'misses \(<=3\) from 1:(.*)', d3).group(1).split())

# 6. literal witnesses
lw = rd(PKT / 'LITERAL_WITNESSES.txt')
dl = rd(HERE / 'depth_literal_30001.txt')
for m in re.finditer(r'^(\d+) (\d+) ([OZJ]+) (\d+) \((\d+), (\d+)\)$', lw, re.M):
    N, k, w, n, H, C = m.groups()
    mm2 = re.search(r'^%s depth (\d+) my witness \S+ \((\d+) vertices\).*?packet witness %s: #J (\d+), F (\d+), (\d+) vertices, literal \(H, C\) = \((\d+), (\d+)\)' % (N, w), dl, re.M)
    check('literal witness %s = %s (depth, #J, F, vertices, H, C)' % (N, w),
          (mm2.group(1), mm2.group(3), mm2.group(4), mm2.group(5), mm2.group(6), mm2.group(7)), (k, k, N, n, H, C))
c1 = '7 13 21 41 43 55 67 77 83 163 173 231 233 247 251 277 299 421 449 503 533 569 583 637 821 1067 1339 1861 2647'.split()
check('odd N <= 30001 with no word at all (any number of J) = Session197 list',
      re.search(r'with no word at all:(.*)', dl).group(1).split(), c1)

# 7. exploratory means at 2^24
ex = rd(HERE / 'census_K2_2p24_exact32.txt')
v3m = {m[0]: m[1].split() for m in re.findall(r'\[(\d+),\d+\) mean by v3\(N\)=0\.\.5,>=6:([\d. ]+)\|', ex)}
v3p = {m[0]: m[1].split() for m in re.findall(r'\[(\d+),\d+\) mean by v3\(N\)=0\.\.6:([\d. ]+)', rd(PKT / 'explore' / 'v3_means_2p24.txt'))}
check('2^24: mean M_2 by v3(N) (exact 32-bit) vs explore/v3_means_2p24.txt', {k: v3m.get(k) for k in v3p}, v3p)
m3m = {m[0]: (m[1].split(), m[2].split()) for m in re.findall(r'\[(\d+),\d+\) mean by v3.*?\| mod 3: ([\d. ]+)\| mod 9:([\d. ]+)', ex)}
m3p = {m[0]: (m[1].split(), m[2].split()) for m in re.findall(r'\[(\d+),\d+\) mean by N mod 3: ([\d. ]+)\| mod 9:([\d. ]+)', rd(PKT / 'explore' / 'mod3_means_2p24.txt'))}
check('2^24: mean M_2 by N mod 3 and mod 9 (windows >= 2^20)', {k: m3m.get(k) for k in m3p if k in m3m}, {k: m3p[k] for k in m3p if k in m3m})
out.append('INFO      max M_2 below 2^24: %s' % re.search(r'max seen (\d+)', ex).group(1))
print('\n'.join(out))
print('SUMMARY: %d MATCH, %d MISMATCH' % (sum(l.startswith('MATCH') for l in out), sum(l.startswith('MISMATCH') for l in out)))
