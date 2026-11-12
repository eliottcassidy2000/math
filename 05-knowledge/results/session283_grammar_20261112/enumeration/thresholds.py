"""All-backbone first rational-ambiguity boundary for odd b>=27."""
from math import gcd,lcm
from collections import Counter
import itertools,json

def patterns():
    return sorted((x,y,z,w) for x,y,z,w in itertools.combinations(range(0,24,2),4) if x+y+z+w<=28)
P=patterns()
def inspect(a):
    best=None
    for ds in P:
        p=tuple(a+d for d in ds)
        defects=[min(x,y)-gcd(x,y) for x,y in itertools.combinations(p,2)]
        r=sum(d>min(defects) for d in defects)
        if r==5:
            candidate=(sum(ds),ds,defects)
            if best is None or candidate<best:best=candidate
    if best is None:raise AssertionError('constructive bound failed')
    return best

def schedule(n):
    if n<0:raise ValueError('n must be nonnegative')
    b=max(27,(n-36)//4+1)
    if b%2==0:b+=1
    candidates=0
    while True:
        candidates+=1
        N=4*(b+2)+inspect(b+2)[0]
        if n<N:return {'n':n,'backbone_top':b,'first_rational_ambiguity':N,'candidates':candidates}
        b+=2

def run():
    period=2*lcm(*((y-x)//((y-x)&-(y-x)) for ds in P for x,y in itertools.combinations(ds,2)))
    hist=Counter();rows=[]
    for a in range(29,29+period,2):
        t,ds,defects=inspect(a)
        if inspect(a+period)[:2]!=(t,ds):raise AssertionError('period')
        hist[t]+=1
        rows.append(P.index(ds))
    return {'status':'PASS','a_domain':'odd a>=29; backbone top b=a-2',
            'period':period,'finite_patterns':len(P),'residue_classes':len(rows),
            'offset_histogram':dict(sorted(hist.items())),
            'examples':[{'backbone_top':b,'first_rational_ambiguity':4*(b+2)+inspect(b+2)[0],'witness':[b+2+d for d in inspect(b+2)[1]]} for b in [27,29,31,33,35,37,39,41,99,999]],
            'schedule_examples':[schedule(n) for n in [150,200,250,1000,10000,1000000]],
            'classification_order':'a=29+2i for i=0,...,3464; reduce a modulo6930',
            'patterns':P,'witness_pattern_indices':rows}
if __name__=='__main__':print(json.dumps(run(),indent=2))
