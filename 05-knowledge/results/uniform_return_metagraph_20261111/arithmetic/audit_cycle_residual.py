"""Independent first-vertex DP for the residual-C3 minimum branching sector.

This audits the literal count, strongness, small complete rooted-class census,
and the exact rational series read from the producer's productive language.
The proof of the pruning filters is reviewed separately in the receipt.
"""
from collections import Counter
from hashlib import sha256
from itertools import product
from pathlib import Path
import json,random,sys
import audit_defect as old

HERE=Path(__file__).resolve().parent
DEFECT=HERE.parent/'defect'
sys.path.insert(0,str(DEFECT))
import cycle_residual as producer


def graph(left,contact,right):
    rows=old.graph(left+contact+right)
    first=len(left)
    for v in range(3):
        rows[first+v]&=~(7<<first)
        rows[first+v]|=1<<(first+(v+1)%3)
    return rows


def literal(rows):
    n=len(rows);full=(1<<n)-1;root=n-1
    forward=old.starts(rows)
    backward=old.starts([sum(((rows[w]>>v)&1)<<w for w in range(n)) for v in range(n)])
    hp=[1]+[sum(x) for x in forward[1:]]
    F=sum(backward[mask][root]*hp[full^mask] for mask in range(1,full+1) if mask>>root&1)
    reach=[row|(1<<v) for v,row in enumerate(rows)]
    for v in range(n):
        for w in range(n):
            if reach[w]>>v&1:reach[w]|=reach[v]
    return F-hp[full],all(x==full for x in reach)


def check_case(left,contact,right):
    defect,strong=literal(graph(left,contact,right))
    weights=producer.two.START
    for x in left:weights=producer.two.advance(weights,x)
    weights=producer.apply_block(weights,contact)
    for x in right:weights=producer.two.advance(weights,x)
    old.check(defect==producer.two.value(weights),'independent C3 first-vertex DP')
    return defect,strong


def main():
    contacts=sorted({min(c,c[1:]+c[:1],c[2:]+c[:2]) for c in product(range(4),repeat=3)})
    counts=Counter();tested=0
    for extra in range(3):
        for ell in range(extra+1):
            for left in product(range(4),repeat=ell):
                for right in product(range(4),repeat=extra-ell):
                    for contact in contacts:
                        defect,strong=check_case(left,contact,right)
                        tested+=1
                        if defect==6 and strong:counts[extra+6]+=1
    old.check(dict(counts)=={6:1,7:4,8:9},'complete C3-residual census through order8')
    rng=random.Random(1414256)
    long_checks=0
    for extra in range(3,6):
        for _ in range(8):
            ell=rng.randrange(extra+1)
            left=tuple(rng.randrange(4) for j in range(ell))
            right=tuple(rng.randrange(4) for j in range(extra-ell))
            contact=contacts[rng.randrange(len(contacts))]
            check_case(left,contact,right);long_checks+=1
    receipt=DEFECT/'cycle_residual.json'
    data=json.loads(receipt.read_text())
    old.check(data['acyclic_after_self_loops'],'productive DAG certificate')
    coeff=[0]*41
    for regex in data['exact_disjoint_path_regex']:
        stars=sum(x.endswith('*') for x in regex)
        # Root and two inports contribute3; unique C marker contributes3.
        base=3+sum(3 if x.startswith('C') else 1 for x in regex if not x.endswith('*'))
        old.check(stars<=2,'at most two independent gaps')
        for n in range(base,41):
            coeff[n]+=(int(n==base) if stars==0 else 1 if stars==1 else n-base+1)
    expected=[0]*41
    for n in range(6,41):
        t=n-6
        expected[n]=t+1+2*max(t,0)+2*max(t-1,0)
    old.check(coeff==expected==data['counts_by_total_order'],'exact all-length rational series')
    report={'status':'PASS; scoped independent audit',
        'complete_literal_presentations_through_order8':tested,
        'C3_contact_orbits':len(contacts),'b6_strong_rooted_classes_by_order':dict(counts),
        'additional_random_presentations_orders9_through11':long_checks,
        'series':'x^6(1+2x+2x^2)/(1-x)^2',
        'exact_productive_regex_paths':len(data['exact_disjoint_path_regex']),
        'scope':'Own graph labels, first-vertex/reverse DP, reachability strongness and complete small orbit census; literal controls are independent, productive DFA discovery is inherited.',
        'proof_review':'q-s=0 fixes an initial residual SCC: first left contact11 or empty-left C3 contacts all11. Such prefixes keep q-s=0. Tagged capacity is induced-extension monotone; a disjoint C3 packing multiplies prefix capacity by3. Thus left bound2 and post-core bound6 are necessary. Last nonzero contact is exact strongness under the fixed initial SCC. The C3 rotation quotient is intrinsic after removing root and its ordered inports.',
        'reviewed_sha256':{p.name:sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in [DEFECT/'cycle_residual.py',receipt]}}
    output=json.dumps(report,indent=2,sort_keys=True)+'\n'
    (HERE/'cycle_residual_audit.json').write_text(output,encoding='utf-8')
    print(output)


if __name__=='__main__':main()
