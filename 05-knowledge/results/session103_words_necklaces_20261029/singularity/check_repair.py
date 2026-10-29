"""Exact repair control: singularity is not hereditary under arbitrary unique-HC extension."""
from native_kernel import *
from check_insertion_update import insertion_update
from pathlib import Path
import json,hashlib

def main():
    path=Path(__file__).resolve().parents[1]/'necklaces/results.json';raw=path.read_bytes()
    A=json.loads(raw)['hostile']['singular_child'];inc=[0,0,1,1,1,0]
    B=[row+[inc[u]] for u,row in enumerate(A)]+[[1-v for v in inc]+[0]]
    ports,M=observer(B);need((ports,M)==direct_observer(B),'independent observer mismatch')
    need((ports,M)==insertion_update(A,inc),'boundary insertion mismatch')
    H=hp_table(B)[127];hc=sum(H[0][t]*B[t][0] for t in range(7))
    cyc=[0,1,2,3,4,6,5]
    need(hc==1 and all(B[cyc[i]][cyc[(i+1)%7]] for i in range(7)),'unique cycle check')
    need(det(M)==-1 and det(observer(A)[1])==0,'repair determinant check')
    relay_rows=[]
    for k in range(4):
        n=6+k
        T=[row+[0]*k for row in A]+[[0]*n for _ in range(k)]
        for v in range(6,n):
            for u in range(6):
                T[u][v]=int(u in (2,3,4));T[v][u]=1-T[u][v]
            for w in range(6,v):T[w][v]=1
        pt,MT=observer(T)
        need((pt,MT)==direct_observer(T),'relay observer disagreement')
        HT=hp_table(T)[(1<<n)-1]
        hct=sum(HT[0][t]*T[t][0] for t in range(n))
        d=det(MT)
        need(hct==1 and d==[0,-1,-1,1][k],'bounded relay discrepancy')
        relay_rows.append({'transitive_relay_vertices':k,'order':n,'hamilton_cycles':hct,'determinant':d})
    out={'status':'FINITE-EXACT','input_sha256':hashlib.sha256(raw).hexdigest(),'bounded_relay_family':relay_rows,'old_order':6,'old_determinant':0,'new_order':7,'new_vertex_outneighbours':[0,1,5],'new_vertex_inneighbours':[2,3,4],'new_adjacency':B,'new_unique_hamilton_cycle':cyc,'new_determinant':-1,'new_nullity':0,'new_smith_form':{'1':42},'scope':'Singularity is not hereditary under arbitrary vertex extensions retaining unique Hamilton cycle. No identification or refutation of missing b/e append grammar.'}
    Path(__file__).with_name('repair_receipt.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
