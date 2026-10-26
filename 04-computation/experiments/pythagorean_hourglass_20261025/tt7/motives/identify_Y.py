# Identification on the cheap triple products Y02, Y05, Y09, Y18 (all primes q = 3 mod 4 in cheapY_*.txt).
# Base: 1, q, q^2, q^3, a84, q a84, a86, a46, |b|^2.  Candidates: products |b|^2 * Tr_g and q^2 * Tr_g for weight-2
# orbits g (2-power and odd levels), and every single orbit as before.  Report the new dimension after each.
import re, os, subprocess, tempfile, sys
from sympy import Matrix
rows = []
for fn in sorted(f for f in os.listdir('.') if f.startswith('cheapY_') and f.endswith('.txt')):
    rows += [list(map(int, l.split())) for l in open(fn)]
rows = sorted({r[0]: r for r in rows}.values())
qs = [r[0] for r in rows]
qmax = int(sys.argv[1]) if len(sys.argv) > 1 else 820
qs = [q for q in qs if q <= qmax]; rows = [r for r in rows if r[0] <= qmax]
Ycols = [[r[k] for r in rows] for k in (1, 3, 5, 6)]      # Y02, Y05, Y09, Y18
P = max(qs) + 2
gpf = os.path.join(tempfile.gettempdir(), 'forms_idY.gp')
open(gpf, 'w').write(f'default(parisize,800000000);\nf=mfeigenbasis(mfinit([8,4,1],0))[1];g=mfeigenbasis(mfinit([8,6,1],0))[1];h=mfeigenbasis(mfinit([4,6,1],0))[1];b=mfeigenbasis(mfinit([32,3,-4],0))[1];\n'
                     f'print(mfcoefs(f,{P}));print(mfcoefs(g,{P}));print(mfcoefs(h,{P}));c=mfcoefs(b,{P});print(vector({P}+1,n,norml2(subst(lift(c[n]),y,I))));\nquit;\n')
out = subprocess.run(['gp', '-q', gpf], capture_output=True, text=True).stdout.strip().splitlines()
f84, f86, f46, B2 = (eval(l.replace('I', '1j')) for l in out[-4:])
B2 = [int(round(abs(x))) for x in B2]
assert B2[3] == 16 and B2[7] == 64 and f86[3] == 20
base = [[1] * len(qs), qs, [q * q for q in qs], [q ** 3 for q in qs], [f84[q] for q in qs], [q * f84[q] for q in qs],
        [f86[q] for q in qs], [f46[q] for q in qs], [B2[q] for q in qs]]
def rank(cols): return Matrix(cols).T.rank()
rB = rank(base); rBY = rank(base + Ycols)
print(f"{len(qs)} primes up to {max(qs)}; new dimension of the Y's: {rBY - rB}")
forms = {}
for fn in ('../candidate_forms.out', 'candidate_forms_odd.out'):
    for l in open(fn):
        lab, p, tr, nr = l.split(); forms.setdefault(lab, {})[int(p)] = int(tr)
def k_of(lab): return int(re.match(r'N\d+k(\d+)', lab).group(1))
hits = []
for lab, d in forms.items():
    if any(q not in d for q in qs): continue
    k = k_of(lab)
    cands = []
    if k == 2:
        cands.append(('|b|^2*' + lab, [B2[q] * d[q] for q in qs]))
        cands.append(('q^2*' + lab, [q * q * d[q] for q in qs]))
    for name, col in cands:
        if not any(col): continue
        rc = rank(base + [col]); dd = rank(base + [col] + Ycols) - rc
        if dd < rBY - rB: hits.append((dd, name))
print("products that reduce the new dimension:", sorted(hits) if hits else 'none')
