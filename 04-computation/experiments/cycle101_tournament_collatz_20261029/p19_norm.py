"""Partial factorisation of the quadratic norm N_{Q(sqrt-19)}(det B_1) of the two-run observer of P19, read from
out/paley_fourier2.json (written by paley_fourier2.py 19), with the splitting test kronecker(-19, p) for each prime
found below 10^7.  Writes out/p19_norm.out."""
import json, sys
sys.path.insert(0, '../growth_gf_runs_fermat_20261029')
from runs import gp_batch
raw = json.load(open('out/paley_fourier2.json'))['19']
fields = raw.strip('[]').split(', ')
qn = int(fields[2])
out = gp_batch(['my(f=factor(%d, 10^7)); [f, vector(#f[,1], i, kronecker(-19, f[i,1])), vector(#f[,1], i, ispseudoprime(f[i,1]))]' % qn])
with open('out/p19_norm.out', 'w') as fo:
    fo.write('N_K(det B_1) for P19 has %d digits\n' % len(str(qn)))
    fo.write('factor(N, 10^7), kronecker(-19, p), ispseudoprime:\n' + out[0] + '\n')
print(open('out/p19_norm.out').read()[:1500])
