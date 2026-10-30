"""Independent receipt checks: score-based supports and literal parity banks.

Does not import the producer. This is same-agent review with a separate path.
"""
from pathlib import Path
from functools import lru_cache
from itertools import combinations
import json

HERE=Path(__file__).resolve().parent
GATES=0

def need(ok,why):
    global GATES
    GATES+=1
    if not ok:raise RuntimeError(why)

def transitive_masks(rows):
    good=[]
    for mask in range(1<<len(rows)):
        scores=sorted((row&mask).bit_count() for i,row in enumerate(rows) if mask>>i&1)
        if scores==list(range(mask.bit_count())):good.append(mask)
    return good

def run():
    receipt=json.loads((HERE/'receipt.json').read_text(encoding='utf8'))
    total=0
    for n in (5,7,9):
        rows=[sum(1<<j for j in range(n) if 0<(j-i)%n<=(n-1)//2) for i in range(n)]
        good=transitive_masks(rows)
        maximal=[s for s in good if not any(s!=t and s&t==s for t in good)]
        if n==5:
            records=receipt['rounding_charts']
            need([(r['epsilon'],r['odd_mask']) for r in records]==
                 [(e,m) for e in (0,1) for m in range(1<<n)],'T5 chart universe')
            bank=[r['colors'] for r in records]
        else:bank=receipt['carousel_rounding_charts'][str(n)]
        need(len(bank)==1<<(n+1),'complete larger chart universe')
        for index,colors in enumerate(bank):
            epsilon,odd=divmod(index,1<<n)
            need(len(colors)==n and all((colors[i] in (0,1,2))==bool(odd>>i&1)
                                       and (colors[i] in (-1,0,1,2)) for i in range(n)),
                 'chart domain exactly the odd support')
            classes=[sum(1<<i for i,c in enumerate(colors) if c==j) for j in range(3)]
            for support in maximal:
                r=(support&odd).bit_count()
                for j,part in enumerate(classes):
                    need((support&part).bit_count()<=(r-epsilon+1)//2+epsilon*int(j==2),
                         'joint support parity inequality')
            total+=1
    # Closed form from most significant first-difference counting.
    def prefix_formula(n):
        return sum(2**(n>>(r+1)).bit_count()*3**r for r in range(n.bit_length()) if n>>r&1)
    @lru_cache(None)
    def halves(n):
        if n<=1:return n
        return 2*halves(n//2)+halves((n+1)//2)
    for n in list(range(65))+[(1<<r)+j for r in range(1,513,17) for j in (-1,0,1,3)]:
        need(prefix_formula(n)==halves(n),'separate exact arithmetic clocks')
        if n<=64:need(prefix_formula(n)==receipt['capacities'][n],'recorded capacity values')
    # Same Ramsey observer on both literal twelve-vertex hostile relations.
    profiles=[]
    for converse in (False,True):
        rows=[0]*12
        for i,j in combinations(range(12),2):
            block,u=divmod(i,4);other,v=divmod(j,4)
            if block!=other:forward=(other-block)%3==1
            else:
                forward=(u==0 or (v-u)%3==1)
                if converse and block==2:forward=not forward
            rows[i if forward else j]|=1<<(j if forward else i)
        profiles.append(transitive_masks(rows))
    need(profiles[0]==profiles[1],'all4096 named Ramsey subset predicates agree')
    need(max(s.bit_count() for s in profiles[0])==6,'hostile transitive budget')
    result={'independent_gates':GATES,'complete_parity_charts':total,
            'hostile_transitive':6,'common_transitive_subsets':len(profiles[0]),
            'status':'Same-agent independent algorithm; no external referee'}
    print(json.dumps(result,sort_keys=True))

if __name__=='__main__':run()
