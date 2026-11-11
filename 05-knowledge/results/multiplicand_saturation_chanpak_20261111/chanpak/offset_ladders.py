"""Exact run-basis certificates for physical generalized-offset ladders.

Default: certify the three positive carriers through the proved run cutoff.
--search7: the narrowly scoped order-seven discovery menu, with compact output.
All masks are outgoing adjacency rows; marks stay outside the cloned module.
"""
from pathlib import Path
import argparse
import json
import math
import subprocess
from head_lift_probe import inflate, ends, converse
from realize35_s2 import seed_metrics

HERE = Path(__file__).resolve().parent
CASES = [
    {'name': 'three', 'rows': (40,33,3,6,15,28), 'selected': (3,), 'port': 4},
    {'name': 'twentyone', 'rows': (42,28,57,48,33,2), 'selected': (3,4), 'port': 2},
    {'name': 'thirtyfive', 'rows': (86,124,40,113,100,65,4), 'selected': (6,), 'port': 1},
]

def run_basis(n):
    """F_j(m)=j! S(m,j), indexed by m,j including the empty boundary."""
    f = [[0]*(n+1) for _ in range(n+1)]
    f[0][0] = 1
    for m in range(1,n+1):
        for j in range(1,m+1):
            f[m][j] = j*(f[m-1][j]+f[m-1][j-1])
    return f

def invert(values):
    """Positive chain lengths determine coefficients with >=1 runs."""
    f = run_basis(len(values))
    a = [0]
    for m,q in enumerate(values,1):
        numerator = q-sum(a[j]*f[m][j] for j in range(1,m))
        denominator = math.factorial(m)
        if numerator % denominator:
            raise ValueError('nonintegral skeleton coefficient')
        a.append(numerator//denominator)
    return a[1:]

def certify(case):
    b,o,v = case['rows'],case['selected'],case['port']
    if v in o:
        raise ValueError('selected port is outside this theorem')
    cutoff = len(b)+1 # k=n-1 exterior vertices; at most two paths.
    values = []
    for m in range(1,cutoff+1):
        t = inflate(b,v,m)
        marks = [i+(m-1)*(i>v) for i in o]
        h,u,w = seed_metrics(t,marks)
        values.append([h,u,w,h-w])
    coefficients = {key: invert([row[j] for row in values])
                    for j,key in enumerate(('H','U','W','D'))}
    d,u = values[0][3],values[0][1]
    expected_d = [d,(d-1)//2]+[0]*(cutoff-2)
    expected_u = [u,u//2]+[0]*(cutoff-2)
    if (d-1)%2 or u%2 or coefficients['D'] != expected_d or coefficients['U'] != expected_u:
        raise ValueError(f'not an offset-doubling ladder: {case}')
    return {**case, 'cutoff': cutoff, 'positive_length_H_U_W_D': values,
            'run_coefficients_from_r1': coefficients,
            'offset': d, 'capacity': u,
            'D_formula': f'1+{d-1}*2^(m-1)',
            'U_formula': f'{u}*2^(m-1)'}

def all_metrics(b):
    n=len(b); full=(1<<n)-1
    e=ends(b); s=ends(converse(b)); h=[sum(x) for x in e]; h[0]=1
    u=[sum(s[m][i]*h[full^m] for m in range(1,full+1)) for i in range(n)]
    k=[[sum(e[m][i]*s[full^m][j] for m in range(1,full)) for j in range(n)] for i in range(n)]
    return h[full],u,k

def selected(data,indices):
    h,u,k=data
    return h-sum(k[i][j] for i in indices for j in indices),sum(u[i] for i in indices)

def search7():
    path='05-knowledge/results/actual_hamilton_multiplicities_20261021/session24-multiplicities/extended_census.json'
    bank=json.loads(subprocess.check_output(['git','show','f92e070445:'+path]))['records']
    hits=[]; count=0
    for record in bank:
        if record['n'] != 7:
            continue
        row=tuple(record['adjacency_outmasks'])
        rows=[row]+([tuple(converse(row))] if record['ordinary_multiplicity']==2 else [])
        for b in rows:
            count+=1; old=all_metrics(b)
            for v in range(len(b)):
                new=all_metrics(inflate(b,v,2)); third=None
                for mask in range(1,1<<len(b)):
                    if mask>>v&1:
                        continue
                    o=[i for i in range(len(b)) if mask>>i&1]
                    d,u=selected(old,o)
                    if d<=1 or (d-1)&(d-2):
                        continue
                    if selected(new,[i+(i>v) for i in o]) != (2*d-1,2*u):
                        continue
                    if third is None:
                        third=all_metrics(inflate(b,v,3))
                    if selected(third,[i+2*(i>v) for i in o]) != (4*d-3,4*u):
                        continue
                    odd=u
                    while odd and odd%2==0:
                        odd//=2
                    hits.append({'rows': b, 'selected': o, 'port': v,
                                 'D': d, 'U': u, 'odd_part_U': odd})
    by_odd={}
    for h in hits:
        by_odd.setdefault(str(h['odd_part_U']),[]).append(h)
    output={'status':'FINITE-EXACT three-length candidates only',
            'ordinary_order7_representatives': count, 'hits':hits,
            'by_odd_part': {k:len(v) for k,v in by_odd.items()}}
    (HERE/'offset_ladders_search7.json').write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps({k:v for k,v in output.items() if k!='hits'},indent=2))
    print('First representative by odd part:')
    for key,group in by_odd.items():
        print(key, json.dumps(group[0]))

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--search7',action='store_true')
    args=ap.parse_args()
    if args.search7:
        search7(); return
    data={'status':'PROVED by finite complete ordered-run coefficient certificate',
          'cases':[certify(c) for c in CASES]}
    (HERE/'offset_ladders_receipt.json').write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps(data,indent=2))

if __name__=='__main__':
    main()
