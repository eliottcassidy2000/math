"""Direct exhaustive four-slot catalogue and exact separating invariant.

Discovery used a rational nullspace. Verification needs only integer sums:
all distinct-source pair labels are enumerated and all compatible disjoint
pair trades checked. This module does not depend on a symbolic package.
"""
import json,sys
from collections import defaultdict
from itertools import combinations
from pathlib import Path
HERE=Path(__file__).resolve().parent
COEFFICIENTS=[-1,-1,0,-1,-2,-2,-1,-1,0,0,-1,-1,1,-1,1,-1,0,-2,0,0,0,1,-2,1,1,0,-1,0]
def main():
    d=json.loads((HERE/'domains.json').read_text());rows=sorted(d['sources'],key=lambda r:r['source']);sources=[r['source'] for r in rows]
    f=json.loads((HERE/'fibre_receipt.json').read_text())
    code=[tuple(r['code']) for r in rows]
    if len(rows)!=28 or len(set(code))!=28:raise RuntimeError('distinct full paired labels required')
    bank=defaultdict(list)
    for i,j in combinations(range(28),2):
        bank[((code[i][0]+code[j][0])%128,)+tuple(a^b for a,b in zip(code[i][1:],code[j][1:]))].append((i,j))
    trades=[]
    for pairs in bank.values():
        for p,q in combinations(pairs,2):
            if set(p)&set(q):continue
            if sum(COEFFICIENTS[i] for i in p)!=sum(COEFFICIENTS[i] for i in q):raise RuntimeError('invariant failed on actual four-trade')
            trades.append(dict(positive=[sources[i] for i in p],negative=[sources[i] for i in q]))
    trades.sort(key=lambda t:(t['positive'],t['negative']))
    if len(trades)!=27:raise RuntimeError('complete four-trade count')
    def value(mask):return sum(c for i,c in enumerate(COEFFICIENTS) if mask>>i&1)
    initial,terminal=value(f['initial']),value(f['terminal'])
    if initial==terminal:raise RuntimeError('endpoints must be separated')
    # Incoming Session46's simpler invariant concerns a different16source pair.
    critical=[sources.index(s) for s in (1638,1974)]
    if any(sum(s in t['positive'] for s in (1638,1974))!=sum(s in t['negative'] for s in (1638,1974)) for t in trades):raise RuntimeError('incoming critical invariant')
    result=dict(status='PROVED exact invariant after finite authenticated label catalogue',sources=sources,coefficients=COEFFICIENTS,
                complete_distinct_pair_universe=378,complete_four_slot_trades=trades,
                endpoint_values=[initial,terminal],incoming_critical_endpoint_values=[sum(mask>>i&1 for i in critical) for mask in (f['initial'],f['terminal'])],
                scope='All fixed-margin original-dyadic moves of support at most4 in this complete28source binary supply, across all current signatures; arbitrary sequences preserve the displayed integer statistic. Budgets4/5 cannot connect these endpoints.')
    out=HERE/'four_move_receipt.json'
    if '--write' in sys.argv:out.write_text(json.dumps(result,indent=2))
    elif json.loads(out.read_text())!=result:raise RuntimeError('saved invariant receipt')
    print('PASS378pairs;27trades;invariant endpoints',initial,terminal,';incoming critical endpoints',result['incoming_critical_endpoint_values'])
if __name__=='__main__':main()
