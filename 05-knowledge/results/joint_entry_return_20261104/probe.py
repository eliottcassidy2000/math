"""Four named-source next-target bank probes, not a closure assertion."""
from pathlib import Path
import json,sys
HERE=Path(__file__).resolve().parent
OLD=HERE.parent/'recursive_seed_switching_20261103/selector'
sys.path.insert(0,str(OLD))
from check import bank_accept

def main():
    bank=json.loads((OLD/'bank.json').read_text())['seeds']
    records=[]
    for parent in (41506343,41506345):
        for digit in (2,4):
            target=3*parent+digit
            hits,stats=bank_accept(bank,target)
            records.append(dict(parent=parent,digit=digit,H=target,hits=hits,statistics=stats))
            print(parent,digit,target,'hits',len(hits),flush=True)
    (HERE/'probe.json').write_text(json.dumps(dict(status='FINITE_EXACT_BANK_INVERSION',records=records),indent=2)+'\n')

if __name__=='__main__':main()
