"""One local seed step after one selected fixed-bank failure; no target-layer census."""
from pathlib import Path
import sys,json,time
HERE=Path(__file__).resolve().parent;RESULTS=HERE.parents[1]
sys.path.insert(0,str(HERE));import inverse
sys.path.insert(0,str(RESULTS/'zero_ready_component_returns_20261103/inverse'));import reader
PARENT_ROWS=(38,60,24,49,33,4);PARENT_ROOT=5;TARGET=13835447

def options():
    n=len(PARENT_ROWS)
    for a in range(n):
        for b in range(a+1,n):
            rows=list(PARENT_ROWS);rows[a]^=1<<b;rows[b]^=1<<a
            yield dict(kind='one_arc_reversal',pair=[a,b],rows=tuple(rows))
    for mask in range(1<<n):
        yield dict(kind='one_birth',newborn_out_mask=mask,rows=reader.birth(PARENT_ROWS,mask))

def main():
    all_hits=[];tried=0;by_kind={}
    for operation in options():
        rows=operation['rows'];fields=reader.all_fields(rows);is_strong=reader.strong(rows)
        for w,f in enumerate(fields):
            found,clocks=inverse.inverse(f,TARGET);tried+=1
            for z in found:
                output=reader.run(f,z['word']);component=reader.graph(rows,w,z['word'])
                n=len(component)+z['prefix']+1
                all_hits.append(dict(operation=operation,root=w,seed_strong=is_strong,seed_fields=f,component_strong=reader.strong(component),fields=output,output_order=n,**z))
                by_kind[operation['kind']]=by_kind.get(operation['kind'],0)+1
    all_hits.sort(key=lambda z:(not z['seed_strong'],z['root']!=PARENT_ROOT,z['prefix']<2,z['output_order']<23,z['output_order'],len(z['word']),z['word']))
    report=dict(target=TARGET,parent_seed_rows=PARENT_ROWS,parent_root=PARENT_ROOT,rooted_candidates=tried,hits=len(all_hits),by_kind=by_kind,witnesses=all_hits)
    (HERE/'enlargement.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(rooted_candidates=tried,hits=len(all_hits),by_kind=by_kind,first=all_hits[:3])),flush=True)
if __name__=='__main__':main()
