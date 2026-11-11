"""Small exact controls; the all-round proof is SATURATION_OBSTRUCTION.md."""
from pathlib import Path
import itertools,json

HERE=Path(__file__).resolve().parent


def need(x,label):
    if not x:raise RuntimeError(label)


def closure_below(a,b,limit):
    values={x for x in (a,b) if x<=limit};x=a
    while x*b<=limit:
        y=b
        while x*y<=limit:values.add(x*y);y*=b
        x*=a
    return values


def run(a,b,rounds=4):
    live={a,b};rows=[]
    for t in range(rounds+1):
        hole=a**(t+1)*b;maximum=max(live)
        full=closure_below(a,b,max(hole,maximum))
        actual=min(full-live)
        need(actual==hole,'exact first missing closure element')
        need(all(x in live for x in full if x<hole),'complete lower frontier')
        if t>=3:need(hole<maximum,'permanent frontier hostile control')
        rows.append(dict(t=t,live=sorted(live),maximum=maximum,first_missing=actual,
                         saturated_to_max=not any(x<=maximum for x in full-live)))
        if t<rounds:live|={x*y for x,y in itertools.combinations(sorted(live),2)}
    return dict(seeds=[a,b],rounds=rows)


def main():
    cases=[run(a,b) for a,b in ((2,3),(2,5),(3,4),(4,9))]
    dependent={2,4}
    for _ in range(2):dependent|={x*y for x,y in itertools.combinations(sorted(dependent),2)}
    need(32 in dependent,'dependent seeds erase one-b obstruction')
    receipt=dict(status='PASS; tiny initial-state controls only',
                 independent_cases=cases,
                 dependent_hostile=dict(seeds=[2,4],round2=sorted(dependent),
                                        collision='2^3*4 = 2*4^2 = 32'))
    (HERE/'saturation_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print('PASS: exact first holes in20 initial states; independent rank2 saturation fails, dependent hostile retained')


if __name__=='__main__':main()
