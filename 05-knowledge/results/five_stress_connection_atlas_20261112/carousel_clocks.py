"""All positive lengths in one-clock carousel families, by exact dominance."""
from pathlib import Path
from itertools import combinations
import json,sys
from tournament_tools import switched_tt
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'seed_recursion_births_20261111'))
import port_transport as pt
ib=pt.ib

def need(x,msg):
    if not x:raise RuntimeError(msg)
def chart(rows,port,mandatory,optional,offset,kind):
    h,u,w=ib.coefficients(rows,port);found=[];maxbound=0;tests=0
    for mask in range(1<<len(optional)):
        O=mandatory+[v for i,v in enumerate(optional) if mask>>i&1]
        if not O:continue
        tests+=1
        d={j:h[j]-sum(w[j][a][b] for a,b in combinations(O,2)) for j in range(1,len(h))}
        poly=pt.exponential(d);poly[1]=poly.get(1,0)-1;poly={b:c for b,c in poly.items() if c}
        U={j:sum(u[j][v] for v in O) for j in range(1,len(h))};ue=pt.exponential(U)
        if not poly:roots='all';bound=None
        else:
            bound=pt.dominance(poly);maxbound=max(maxbound,bound)
            roots=[a for a in range(1,bound) if pt.evaluate(poly,a)==0]
        if roots:
            found.append({'selected':O,'clock_sizes':roots,'big_block_sizes':roots if roots=='all' else [a+offset for a in roots],'no_further_roots_from_clock':bound,'defect_minus1_exponentials':sorted(poly.items()),'U_exponentials':sorted(ue.items())})
    return {'kind':kind,'base_order':len(rows),'port':port,'mark_patterns':tests,'largest_dominance_bound':maxbound,'unit_solutions':found}

def run():
    out=[]
    for k in [3,5,7,9,11]:
        rows=switched_tt(k,sum((i%2)<<i for i in range(k)))
        charts=[chart(rows,0,[],list(range(1,k)),0,'big block unselected')]
        labels,child=ib.inflate(rows,0,2)
        clones=[i for i,v in enumerate(labels) if v==0]
        optional=[i for i,v in enumerate(labels) if v!=0]
        for chosen in clones:
            port=next(i for i in clones if i!=chosen)
            charts.append(chart(child,port,[chosen],optional,1,'only first selected' if chosen==clones[0] else 'only last selected'))
        result={'carousel_order':k,'scope':'All positive big-block lengths; the block is unselected or contains exactly one selected endpoint. Every subset of the other singleton blocks. Interior marks and multiple selected clone vertices excluded.','charts':charts}
        out.append(result)
        print(json.dumps({'k':k,'charts':[(v['kind'],v['mark_patterns'],v['largest_dominance_bound'],len(v['unit_solutions'])) for v in charts]}),flush=True)
    dest=HERE/('carousel_clocks.optimized.json' if not __debug__ else 'carousel_clocks.json')
    dest.write_text(json.dumps(out,indent=2)+'\n',encoding='utf8')

if __name__=='__main__':run()
