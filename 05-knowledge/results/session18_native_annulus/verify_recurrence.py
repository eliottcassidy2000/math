#!/usr/bin/env python3
"""Independent grouped-edit audit of exact coherence recurrence (G)--(H)."""
import argparse,hashlib,json
from collections import defaultdict
from pathlib import Path
import verify_annulus as a
HERE=Path(__file__).resolve().parent

def independent_g(sw):
    M=sw.N;delta=defaultdict(lambda:[0,0])
    m=2*M+1
    for m,v in [(m,(a.L[m%3][0]-a.U[m%3][0],a.L[m%3][1]-a.U[m%3][1])),(4*M+1,a.U[(4*M+1)%3]),(4*M+3,a.U[(4*M+3)%3])]:
        while m%8==5:m=(m-1)//4
        delta[m][0]+=v[0];delta[m][1]+=v[1]
    sources=[]
    for s,(da,db) in delta.items():
        u,t,d=sw.routes.get(s);sg=1-2*(t%2);v=a.charge(M,s)
        sources.append((u,d,(sg*v[0],sg*v[1]),(sg*da,sg*db)))
    out=[]
    for h in range(a.DEPTH_CAP+1):
        D=defaultdict(lambda:[0,0]);sub=0
        for u,d,v,dv in sources:
            if d<=h:
                D[u][0]+=dv[0];D[u][1]+=dv[1];sub+=a.Q(dv)+a.beta(v,dv)
        g=sum(a.Q(v)+a.beta((sw.W[u][0][h],sw.W[u][1][h]),v) for u,v in D.items())-sub
        out.append(g)
    return out

def audited_step(sw):
    g=independent_g(sw);before=[e-m for e,m in zip(sw.E,sw.M)];sw.step()
    a.check(g==[e-m-old for e,m,old in zip(sw.E,sw.M,before)],'grouped coherence derivative')
    return g

def main():
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');args=p.parse_args()
    data=json.loads((HERE/'frontier.json').read_text());routes=a.Routes(data['frontier'],dict(data['frontier_ranks']))
    old=a.Sweep(4096,routes);large=a.Sweep(64*4096,routes)
    G0=[e-m-f+n for e,m,f,n in zip(large.E,large.M,old.E,old.M)]
    oldg=audited_step(old);sumg=[0]*(a.DEPTH_CAP+1)
    for j in range(64):
        g=audited_step(large);sumg=[x+y for x,y in zip(sumg,g)]
    G1=[e-m-f+n for e,m,f,n in zip(large.E,large.M,old.E,old.M)]
    a.check([y-x for x,y in zip(G0,G1)]==[x-y for x,y in zip(sumg,oldg)],'64-scale gap recurrence')
    for sw in (old,large):
        e,m,info=a.direct(sw.N,routes)
        a.check(e==sw.E and m==sw.M,'final recurrence source rebuild')
    a.check((G0[-1],G1[-1])==(-201674,-201616),'nonmonotone scale coherence hostile')
    result={'status':'PROVED elementary recurrence audited in fixed finite control','N_pair':[4096,4097],'G0':G0,'G1':G1,'old_g':oldg,'large_64step_g':sumg,'checks':a.CHECKS}
    blob=a.encode(result);out=HERE/'recurrence_receipt.json'
    if args.write:out.write_bytes(blob)
    else:a.check(out.read_bytes()==blob,'recurrence saved receipt equality')
    print(json.dumps({'status':'PASS','checks':result['checks'],'receipt_sha256':hashlib.sha256(blob).hexdigest(),'final_G0_G1':[G0[-1],G1[-1]]},sort_keys=True))
if __name__=='__main__':main()
