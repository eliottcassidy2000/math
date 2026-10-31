"""Independent full labelled enumeration for the capacity recurrence through six."""
from itertools import combinations, product
from pathlib import Path
import argparse
import json

HERE = Path(__file__).resolve().parent

def run():
    receipt = json.loads((HERE/'capacity_receipt.json').read_text())
    rows_by_name = {
        'join2': [2, 0], 'C3': [2, 4, 1],
        'P7': [sum(1 << j for j in range(7) if (j-i)%7 in (1,2,4)) for i in range(7)]}
    sets = {}
    for name, rows in rows_by_name.items():
        n = len(rows)
        cycles = [sum(1 << i for i in t) for t in combinations(range(n),3)
                  if all(sum(rows[i]>>j&1 for j in t)==1 for i in t)]
        sets[name] = [tuple(i for i in range(n) if s>>i&1) for s in range(1<<n)
                      if not any(s&c==c for c in cycles)]
    values = [0,1]; accepted = []; checks = 0
    for k in range(2,7):
        best = values[-1]; count = 0
        for name, rows in rows_by_name.items():
            for weights in product(range(k), repeat=len(rows)):
                if sum(a>0 for a in weights)<2:
                    continue
                if any(sum(weights[i] for i in s)>k for s in sets[name]):
                    continue
                count += 1
                best = max(best, sum(values[a] for a in weights))
        values.append(best); accepted.append(count)
        checks += 1
        if best != receipt['values'][k]:
            raise RuntimeError('Full labelled enumeration disagrees')
    # Independently replay every stored witness, including budgets beyond six.
    for k, witness in enumerate(receipt['witnesses'][2:], 2):
        if 'reuse_budget' in witness:
            valid = witness['reuse_budget']==k-1 and receipt['values'][k]==receipt['values'][k-1]
        else:
            a = witness['demands']; name = witness['template']
            valid = (len(a)==len(rows_by_name[name]) and all(0<=x<k for x in a)
                     and sum(x>0 for x in a)>=2
                     and all(sum(a[i] for i in s)<=k for s in sets[name])
                     and sum(receipt['values'][x] for x in a)==receipt['values'][k])
        checks += 1
        if not valid:
            raise RuntimeError('Invalid construction witness')
    return {'checks': checks, 'independent_values': values, 'accepted_labelled_candidates': accepted,
            'witness_budgets': len(receipt['values'])-1,
            'scope': 'Independent exhaustive upper/lower comparison through six; witness validation through twelve.'}

if __name__=='__main__':
    ap=argparse.ArgumentParser(); ap.add_argument('--output',type=Path); args=ap.parse_args()
    text=json.dumps(run(),indent=2,sort_keys=True)+'\n'
    if args.output: args.output.write_text(text,encoding='utf-8',newline='\n')
    else: print(text,end='')
