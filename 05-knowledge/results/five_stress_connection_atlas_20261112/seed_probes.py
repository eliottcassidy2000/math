"""Full coefficient probes of the four C5-core seeds in the round census."""
from pathlib import Path
import json,sys
from itertools import combinations
from collections import defaultdict
from math import comb
from tournament_tools import switched_tt
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'seed_recursion_births_20261111'))
import cones
ib=cones.pt.ib

def need(x,msg):
    if not x:raise RuntimeError(msg)
def exp2(poly):
    out=defaultdict(int)
    for (i,j),c in poly.items():
        for a in range(1,i+1):
            for b in range(1,j+1):out[a,b]+=c*(-1)**(i+j-a-b)*comb(i,a)*comb(j,b)
    return {k:v for k,v in out.items() if v}
def ev(poly,a,b):return sum(c*p**a*q**b for (p,q),c in poly.items())
def encode(poly):return [[*k,c] for k,c in sorted(poly.items())]
def stabled(rows,mask,v):
    h,u,w=ib.coefficients(rows,v)
    O=[x for x in range(len(rows)) if mask>>x&1]
    defect=[h[j]-sum(w[j][a][b] for a,b in combinations(O,2)) for j in range(len(h))]
    slope=[sum(u[j][a] for a in O) for j in range(len(h))]
    return defect,slope

def run():
    data=json.loads((HERE/'tournament_receipt.json').read_text())
    selected=[r for r in data['seed_records'] if not r['parabola_core']]
    records=[]
    for item in selected:
        rows=item['rows'];mask=item['selected'];n=len(rows)
        ports=[]
        for v in range(n):
            if mask>>v&1:continue
            d,u=stabled(rows,mask,v)
            ports.append({'port':v,'stable_all_positive':d==[0,1]+[0]*(len(d)-2),'defect_coefficients':d,'slope_coefficients':u})
        gates=[]
        for v in range(n):
            if not mask>>v&1:continue
            labels,child=ib.inflate(rows,v,2)
            positions=[i for i,x in enumerate(labels) if x==v]
            old=[i for i,x in enumerate(labels) if x!=v and mask>>x&1]
            for chosen in positions:
                selected_child=old+[chosen];cm=sum(1<<i for i in selected_child)
                port=next(i for i in positions if i!=chosen)
                d,u=stabled(child,cm,port)
                gates.append({'old_selected':v,'selected_end':'first' if chosen==positions[0] else 'last','rows':child,'selected':cm,'port':port,'stable_all_positive':d==[0,1]+[0]*(len(d)-2),'defect_coefficients':d,'slope_coefficients':u})
        records.append({'source':item,'unselected_ports':ports,'selected_clone_gates':gates})
        print(json.dumps({'n':n,'selected':mask,'stable_ports':[p['port'] for p in ports if p['stable_all_positive']],'stable_selected_gates':[(g['old_selected'],g['selected_end']) for g in gates if g['stable_all_positive']]}),flush=True)
    # Split a C5 port as left chain, selected middle, right chain.
    word='0001010';rows=switched_tt(7,sum(int(c)<<i for i,c in enumerate(word)))
    grid=[];literal=0
    for pair in combinations(range(3,7),2):
        O=[1,*pair];d,u=cones.kernel(rows,[0,2],O);de=exp2(d);ue=exp2(u)
        roots=[(a,b) for a in range(1,41) for b in range(1,41) if ev(de,a,b)==1]
        for a in range(1,4):
            for b in range(1,4):
                child,sel,_=cones.expand(rows,O,[0,2],[a,b]);H,U,W,*_=ib.tc.path_tables(child)
                dd=H-sum(W[x][y] for x,y in combinations(sel,2));uu=sum(U[x] for x in sel)
                need(ev(de,a,b)==dd and ev(ue,a,b)==uu,'two-clock literal check');literal+=1
        grid.append({'selected':O,'defect_exponentials':encode(de),'U_exponentials':encode(ue),'unit_roots_positive_40_square':roots,'U_at_roots':[[a,b,ev(ue,a,b)] for a,b in roots]})
        print(json.dumps(grid[-1]),flush=True)
    out={'scope':'Every unselected port and first/last one-selected-clone gate of the four non-parabola unit seeds through order12; full finite Stirling coefficients. Six two-clock charts on one C5 port; roots searched only in [1,40]^2 unless separately proved.','records':records,'two_clock_charts':grid,'two_clock_literal_controls':literal}
    dest=HERE/('seed_probes.optimized.json' if not __debug__ else 'seed_probes.json')
    dest.write_text(json.dumps(out,indent=2)+'\n',encoding='utf8')

if __name__=='__main__':run()
