"""Fixed-path n=7 Boolean chart: a literal test of the 17/71 pairing."""
from collections import Counter
from itertools import combinations
from pathlib import Path
import json
from scout import rows, h_count

HERE=Path(__file__).resolve().parent

def main():
    n=7
    edges=list(combinations(range(n),2))
    free=[i for i,(a,b) in enumerate(edges) if b>a+1]
    fixed=sum(1<<i for i,(a,b) in enumerate(edges) if b==a+1)
    def expand(x):
        return fixed+sum((x>>i&1)<<e for i,e in enumerate(free))
    counts=[h_count(rows(n,expand(x))) for x in range(1<<len(free))]
    full=(1<<len(free))-1
    hits=[x for x,h in enumerate(counts) if h==17 and counts[x^full]==71]
    result={'n':n,'fixed_mask':fixed,'free_edge_indices':free,
            'histogram':dict(sorted(Counter(counts).items())),
            '17_to_71_tile_complement_hits':hits,
            'H17_complement_values':sorted(set(counts[x^full] for x,h in enumerate(counts) if h==17)),
            'H71_tile':next(x for x,h in enumerate(counts) if h==71),
            'H71_mask':next(expand(x) for x,h in enumerate(counts) if h==71)}
    if hits:
        x=hits[0]
        result['witness']={'tiles':[x,x^full],'masks':[expand(x),expand(x^full)],
                           'out_rows':[rows(n,expand(x)),rows(n,expand(x^full))]}
    (HERE/'chart_scout.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in result.items() if k not in ('histogram','17_to_71_tile_complement_hits')}))
    print('pair_hits',len(hits),flush=True)

if __name__=='__main__':main()
