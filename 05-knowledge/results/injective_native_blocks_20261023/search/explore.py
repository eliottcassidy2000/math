"""Bounded actual native transition-menu generator; discovery only."""
import csv,gzip,json,sys,time
from array import array
from collections import defaultdict,Counter
from itertools import combinations
from pathlib import Path
HERE=Path(__file__).resolve().parent
RESULTS=HERE.parents[1]
NATIVE=RESULTS/'metagraph_order8_20261023'
D=json.loads((RESULTS/'metagraph_order8_electrical_20261023/receipt.json').read_text())['order8']
X=D['higher_layers'][-1]['lifts'][0];Y=D['radical_lifts_mod4'];Z=D['eighth_radical_lifts_mod8'][0]
SIG=[tuple([Z[i]%4]+[Y[j][i]%2 for j in range(1,5)]) for i in range(3528)]
NODES=list(csv.DictReader((NATIVE/'merged_nodes.tsv').open(),delimiter='\t'))
MERGED={int(r[k]):int(r['merged_id']) for r in NODES for k in ('class_a','class_b')}
OWNER=array('H');OWNER.frombytes(gzip.decompress((NATIVE/'owner_map.u16.gz').read_bytes()))
if sys.byteorder!='little':OWNER.byteswap()
CHORDS=[(i,j) for i in range(8) for j in range(i+2,8)]
FULL=(1<<21)-1

def paths(source):
    canon=int(NODES[source]['canon_a']);A=[[0]*8 for _ in range(8)]
    for bit,(i,j) in enumerate(combinations(range(8),2)):
        A[i][j]=(canon>>(27-bit))&1;A[j][i]=1-A[i][j]
    rows=[]
    def visit(p,unused):
        if not unused:
            mask=sum(A[p[i]][p[j]]<<b for b,(i,j) in enumerate(CHORDS))
            if MERGED[OWNER[mask]]!=source:raise RuntimeError('actual source owner mismatch')
            t=MERGED[OWNER[FULL^mask]]
            rows.append({'path':p,'target':t,'mask':mask,'phase':(X[source]-X[t])%8})
            return
        for v in range(8):
            if (unused>>v)&1 and (not p or A[p[-1]][v]):visit(p+[v],unused^(1<<v))
    visit([],255)
    return rows

def menu(rows):
    banks=defaultdict(list)
    for i,r in enumerate(rows):
        p=r['path']
        for k in range(1,6):banks[(p[0],p[-1],r['phase'],k,p[k],p[k+1])].append(i)
    out={}
    for pins,ids in banks.items():
        for i,j in combinations(ids,2):
            a,b=rows[i],rows[j];u,v=a['target'],b['target']
            if u==v or SIG[u]!=SIG[v]:continue
            if u>v:u,v,i,j=v,u,j,i
            cost=sum(x!=y for x,y in zip(a['path'],b['path']))
            rec=(cost,pins,i,j)
            if (u,v) not in out or rec<out[u,v]:out[u,v]=rec
    return out

def main():
    mode=sys.argv[1] if len(sys.argv)>1 else 'two_codes'
    sources=list(range(3528)) if mode=='all' else [s for s in range(3528) if SIG[s] in (SIG[1238],SIG[2246])]
    results=[];pairs=defaultdict(list);start=time.time()
    for s in sources:
        if NODES[s]['class_a']==NODES[s]['class_b']:continue
        rows=paths(s);m=menu(rows)
        results.append({'source':s,'canon':int(NODES[s]['canon_a']),'sigma':SIG[s],'x':X[s],'rows':rows,
                        'menu':[[u,v,list(w)] for (u,v),w in sorted(m.items())]})
        for pair,w in m.items():pairs[pair].append(s)
    summary={'sources':len(results),'paths':sum(len(r['rows']) for r in results),'menus':sum(len(r['menu']) for r in results),
             'paired_targets':len(pairs),'shared_pairs':sum(len(v)>1 for v in pairs.values()),
             'target_zero_256':[[*p,v] for p,v in pairs.items() if (X[p[0]]-X[p[1]])%128==0 and len(v)>1],
             'elapsed':time.time()-start}
    (HERE/('cache_'+mode+'.json')).write_text(json.dumps(results,separators=(',',':')))
    (HERE/('summary_'+mode+'.json')).write_text(json.dumps(summary,indent=2))
    print(json.dumps(summary,indent=2))

if __name__=='__main__':main()
