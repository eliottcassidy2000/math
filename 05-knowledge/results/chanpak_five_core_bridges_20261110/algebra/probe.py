"""Independent endpoint and seven-context probes for the fixed-core compiler."""
import importlib.util
import itertools
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
SRC = ROOT / '05-knowledge/results/session251_creative_completion_20261111/ternary/compiler.py'
spec = importlib.util.spec_from_file_location('inherited_compiler', SRC)
c = importlib.util.module_from_spec(spec)
spec.loader.exec_module(c)


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def literal_endpoints(rows):
    n = len(rows)
    starts, ends = [0]*n, [0]*n
    for p in itertools.permutations(range(n)):
        if all(rows[p[i]] >> p[i+1] & 1 for i in range(n-1)):
            starts[p[0]] += 1
            ends[p[-1]] += 1
    return starts, ends


def reduced(sig):
    d, kp, km, ip, im, op, om, jp, jm = sig
    return (d+jp+jm, d, d+kp, ip, op, -im, -om)


def evaluate(z, a, b):
    h, d, L, i, o, s, e = z
    return h+d*(c.g(a)+c.g(b)-2)+L*a*b+i*a+o*b+s*(a&1)+e*(b&1)


def recover(values):
    h, p10, p30, p20, p01, p03, p11 = values
    i = (p30-p10)//2
    s = p10-h-i
    d = (p20-h-2*i)//2
    o = (p03-p01)//2
    e = p01-h-o
    L = p11-h-i-o-s-e
    return (h,d,L,i,o,s,e)


def rank(rows):
    from fractions import Fraction
    basis = {}
    chosen = []
    for idx, row in enumerate(rows):
        v = list(map(Fraction,row))
        for pivot, old in sorted(basis.items()):
            f=v[pivot]
            if f:
                v=[x-f*y for x,y in zip(v,old)]
        pivot=next((j for j,x in enumerate(v) if x),None)
        if pivot is not None:
            f=v[pivot]
            basis[pivot]=[x/f for x in v]
            chosen.append(idx)
    return len(basis), chosen


def main():
    cores=c.bank()
    signatures=[]
    contexts=[(0,0),(1,0),(3,0),(2,0),(0,1),(0,3),(1,1)]
    checks=0
    coredata=[]
    for rows in cores:
        starts,ends=literal_endpoints(rows)
        full=(1<<len(rows))-1
        d=sum(starts)
        coredata.append({'rows':rows,'H':d,'starts':starts,'ends':ends})
        for O in range(full+1):
            sig=c.signature(rows,O)
            require(sig[2]==-d,'total alternating identity')
            require(sig[4]==-sum(starts[v] for v in range(len(rows)) if not O>>v&1),'incoming alternating identity')
            require(sig[6]==-sum(ends[v] for v in range(len(rows)) if O>>v&1),'outgoing alternating identity')
            z=reduced(sig)
            require(all(x>=0 for x in z),'nonnegative coefficients')
            require(z[2]%2==0 and (z[3]+z[5])%2==0 and (z[4]+z[6])%2==0,'parity conditions')
            vals=[c.hp_subsets(c.construct(rows,O,a,b))[-1] for a,b in contexts]
            require(recover(vals)==z,'seven physical samples')
            for a,b in itertools.product(range(5),repeat=2):
                require(evaluate(z,a,b)==c.count(sig,a,b),'reduced formula')
                if a<3 and b<3:
                    require(evaluate(z,a,b)==c.hp_subsets(c.construct(rows,O,a,b))[-1],'direct DP')
                checks+=1
            signatures.append({'rows':rows,'O':O,'reduced':z,'samples':vals})
    r,chosen=rank([q['reduced'] for q in signatures])
    report={'status':'FINITE-EXACT controls; identities are proved in companion note',
            'core_classes':len(cores),'coloured_cores':len(signatures),
            'formula_comparisons':checks,'signature_rank':r,
            'rank_witnesses':[signatures[j] for j in chosen],
            'core_endpoint_profiles':coredata}
    out=Path(__file__).with_name('probe.json')
    out.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in report.items() if k not in ('core_endpoint_profiles','rank_witnesses')},sort_keys=True))


if __name__=='__main__':
    main()
