# Seven-vertex census: add the weight-6 forms of levels 4 and 8 (4.6.a.a = eta(2z)^12, 8.6.a.a), together and with
# the weight-6 orbits of levels 16..128 (Frechette-Ono-Papanikolas: 6F5(1) is a Hecke trace on S_6(Gamma0(8))).
# mac-mini-2026-10-25-pyth.
import sys
sys.argv = [sys.argv[0]]
exec(open('tt7_candidate_fit.py').read().split("print(f\"{len(ps)} values")[0])
e12 = eta_product([(2, 12)])
import subprocess, tempfile, os
gpf = os.path.join(tempfile.gettempdir(), 'c6_8_6_a_a.gp')
open(gpf, 'w').write('default(parisize,200000000);\nmf=mfinit([8,6,1],0);L=mfeigenbasis(mf);c=mfcoefs(L[1],%d);print(c);\nquit;\n' % max(ps))
out = subprocess.run(['gp', '-q', gpf], capture_output=True, text=True).stdout.strip().splitlines()[-1]
c6 = eval(out)
assert c6[3] == 20 and e12[3] == -12
cols = dict(base)
cols['q(q-1)c6'] = {p: p * (p - 1) * c6[p] for p in ps}
cols['q(q-1)e12'] = {p: p * (p - 1) * e12[p] for p in ps}
sol = in_span(cols)
print("base + 8.6.a.a + 4.6.a.a:", 'EXACT ' + str({k: str(v) for k, v in sol.items() if v}) if sol else 'NOT in span')
for lab in forms:
    if lab_k(lab) == 6: cols.update(orbit_cols(lab))
sol = in_span(cols)
print(f"base + all weight-6 newforms of levels 4..128 ({len(cols)} columns):", 'EXACT' if sol else 'NOT in span')
for lab in forms:
    if lab_k(lab) == 4 and 'c1' in lab: cols.update(orbit_cols(lab))
print(f"  + weight-4 trivial-character orbits: {len(cols)} columns for {len(ps)} values")
if len(cols) < len(ps) - 3:
    sol = in_span(cols); print("  ->", 'EXACT' if sol else 'NOT in span')
