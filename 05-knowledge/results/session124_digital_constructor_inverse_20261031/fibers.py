"""Finite exact boundary data for infinite actual-seed singleton fibers."""
from collections import defaultdict
import argparse,json
from pathlib import Path
import digital
BASE=Path(__file__).resolve().parent
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output',type=Path,required=True)
args=parser.parse_args()
raw=json.loads((BASE/'actual_seed_roots.json').read_text())
states=[(r['state'][0],r['state'][1],r['state'][2],r['state'][0]-r['state'][3]) for r in raw]
def consecutive(values):
    best=run=0;previous=None
    for x in sorted(values):
        run=run+1 if previous is not None and x==previous+2 else 1
        best=max(best,run);previous=x
    return best
result=[]
for h in sorted({s[0] for s in states}):
    sources=defaultdict(set)
    for state in states:
        if state[0]!=h:continue
        for q,kappa,b in digital.kappa_sources(state):sources[q].add(kappa)
    L=max(sources)
    offsets={}
    for q in range(2*L+h+1):
        own=set(sources.get(q,()))
        previous=digital.predecessor(h,q)
        if previous:own|=offsets[previous[0]]
        offsets[q]=own
    boundary={q:sorted(offsets[q]) for q in range(L+1,2*L+h+1) if offsets[q] and digital.predecessor(h,q)[0]<=L}
    signatures={tuple(K) for K in boundary.values()}
    allK=set().union(*sources.values())
    record={'h':h,'maximum_source_q':L,'global_offsets':sorted(allK),'eventual_boundary_nodes':len(boundary),'eventual_distinct_fibers':len(signatures),'maximum_eventual_fiber_size':max(map(len,signatures)),'maximum_eventual_consecutive_odd_width':max(map(consecutive,signatures)),'boundary_nodes':[{'q':q,'offsets':K} for q,K in boundary.items()]}
    result.append(record)
    print(h,'offsets',len(allK),'max_q',L,'boundary',len(boundary),'signatures',len(signatures),'max_fiber',record['maximum_eventual_fiber_size'],'max_odd_interval',record['maximum_eventual_consecutive_odd_width'])
args.output.write_text(json.dumps(result,indent=2)+'\n')
print('EXACT_INFINITE_FIBERS_FROM_FINITE_ACTUAL_SEED_ANCESTRY_PASS')
