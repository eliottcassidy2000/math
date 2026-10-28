"""Check concurrent implementations without duplicating their production kernels."""
import contextlib,io,json,sys
from pathlib import Path
from collections import defaultdict
from time import perf_counter
import native

HERE=Path(__file__).resolve().parent
def run(bound):
    fast=native.source('incoming_fast','05-knowledge/results/session87_four_vertices_enumeration_20261028/fast_counts.py')
    start=perf_counter();a,stats=fast.degree_bucket_dp(bound);elapsed=perf_counter()-start
    reference=json.loads((HERE/'ordinary-300.json').read_text())['values']
    native.need(a==reference[:bound+1],'incoming degree-bucket complete prefix')
    selected=native.source('incoming_selector','05-knowledge/results/kernel125_closed_carriers_20261028/root_minimum.py')
    from root_switch_skip import generate
    parent_a=[()];parent_b=[()];roots=[]
    for n in range(1,8):
        parent_a,_=generate(parent_a,n)
        with contextlib.redirect_stdout(io.StringIO()):
            parent_b,row=selected.generate(parent_b,n,65536,True,True)
        buckets=defaultdict(list)
        for i,t in enumerate(parent_a):buckets[native.oracle.colours(t)[1]].append(i)
        seen=set()
        for t in parent_b:
            matches=[i for i in buckets[native.oracle.colours(t)[1]] if native.oracle.isomorphisms(t,parent_a[i],False)]
            native.need(len(matches)==1 and matches[0] not in seen,'incoming selected class bijection')
            seen.add(matches[0])
        native.need(len(seen)==len(parent_a),'incoming selected coverage')
        roots.append(dict(n=n,classes=len(seen),reroots_examined=row['reroots_examined'],bijective=True))
    return dict(degree_bucket_bound=bound,degree_bucket_complete_agreement=True,
        degree_bucket_single_run_seconds=elapsed,degree_bucket_statistics=stats,
        selector_class_bijections=roots,note='Single run time is not a paired benchmark; incoming algorithms retain their original credit.')

if __name__=='__main__':
    if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(0)
    result=run(300)
    (HERE/'incoming.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in result.items() if k!='degree_bucket_statistics'}))
