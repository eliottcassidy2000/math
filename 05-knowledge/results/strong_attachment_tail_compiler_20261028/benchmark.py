"""Same-process alternating timings; all coefficient equality checks are exact."""
import argparse,json,platform,statistics,sys
from pathlib import Path
from time import perf_counter
from counting import contracted,contracted13,inherited,need

def main():
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);args=p.parse_args()
    rows=[]
    for n,core in [(100,False),(200,False),(250,False),(50,True)]:
        timings={'baseline':[],'unit_tail':[],'tail13':[]};counts={};target=None
        for repeat in range(3):
            order=['baseline','unit_tail'] if repeat%2==0 else ['unit_tail','baseline']
            if not core and n==100:order.append('tail13')
            for mode in order:
                start=perf_counter()
                if mode=='baseline':values,stats=inherited.local_dp(n,core=core)
                elif mode=='unit_tail':values,stats=contracted(n,core=core)
                else:values,stats=contracted13(n)
                elapsed=perf_counter()-start
                if target is None:target=values
                need(values==target,'alternating benchmark coefficient agreement')
                timings[mode].append(elapsed)
                counts[mode]={k:v for k,v in stats.items() if k!='rows'}
        median={k:statistics.median(v) for k,v in timings.items() if v}
        row=dict(bound=n,core=core,timings=timings,median_seconds=median,
                 baseline_over_unit_tail=median['baseline']/median['unit_tail'],counts=counts)
        rows.append(row);print(json.dumps({k:v for k,v in row.items() if k not in ('counts','timings')}),flush=True)
    Path(args.output).write_text(json.dumps(dict(python=sys.version,platform=platform.platform(),
        scope='same process; three alternating runs; wall clock; no historical cross-machine comparison',rows=rows),indent=2)+'\n',encoding='utf-8')

if __name__=='__main__':main()
