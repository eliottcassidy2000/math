"""Matched root-switch loops, including the same parent setup, excluding H audits."""
import json,statistics,sys
from pathlib import Path
from root_switch_skip import generate,need

def main():
    parents=[()]
    for n in range(1,8):parents,_=generate(parents,n)
    samples={'all_rows':[],'one_per_block':[]};stats={};target=None
    for repeat in range(3):
        order=[False,True] if repeat%2==0 else [True,False]
        for skip in order:
            children,r=generate(parents,8,skip_blocks=skip)
            if target is None:target=children
            need(children==target,'same complete ordered order-eight bank')
            name='one_per_block' if skip else 'all_rows';samples[name].append(r['seconds'])
            stats[name]={k:v for k,v in r.items() if k!='seconds'}
    medians={k:statistics.median(v) for k,v in samples.items()}
    result=dict(order=8,samples=samples,median_seconds=medians,
        all_rows_over_one_per_block=medians['all_rows']/medians['one_per_block'],stats=stats,
        scope='three alternating runs, same setup and checks, no extra H/prime audit; all representatives equal')
    Path(sys.argv[1]).write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in result.items() if k!='stats'}))

if __name__=='__main__':main()
