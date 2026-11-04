#!/usr/bin/env python3
from pathlib import Path
import json,subprocess,hashlib
import audit as A
import physical_audit as B
P=Path(__file__).resolve().parent

def main():
    q=(2,4,1);T,blocks=B.inflate(q,(4,4,4));U,_=B.inflate(A.opposite(q),(4,4,4));V=list(U);child=(2,4,9,3)
    for i in range(4):V[i]=(V[i]&~15)|child[i]
    V=tuple(V);core_queries=[(3,1),(6,4)];graphs=[T,U,V];req=[(T,5),((14,12,8,0),4),(child,4)]
    for g in graphs:
        for a,c in core_queries:req.append((B.ports(g,B.lift(a,blocks),B.lift(c,blocks)),5))
    for g in graphs:
        n=len(g);amb=list(g)+[1<<(n+1),(1<<n)-1]
        for i in range(n):amb[i]|=1<<n
        for a,c in core_queries:req.append((B.ports(amb,B.lift(a,blocks)|(1<<(n+1)),B.lift(c,blocks)|(1<<(n+1))),1))
    inp='\n'.join(' '.join(map(str,[len(g),h,*g])) for g,h in req)+'\n';out=subprocess.check_output([str(P/'physical')],input=inp,text=True);values=[list(map(int,z.split())) for z in out.splitlines()]
    source_profile=values[0];A.require([x*pow(source_profile[0],-1,16)%16 for x in source_profile]==[1,14,4,8,0],'recursive actual profile closure');A.require(all(x>=y for x,y in zip(values[2],values[1])) and values[2]!=values[1],'actual strict child profile dominance')
    core=values[3:9];ambient=values[9:15];contrasts=[[core[z][h]-core[2+z][h] for h in range(5)] for z in range(2)];A.require(contrasts[0]==[-x for x in contrasts[1]],'complete core duality');A.require(contrasts[0][0]<0,'negative baseline chosen');A.require(all(core[4+z][h]>=core[2+z][h] for z in range(2) for h in range(5)),'all tested target profiles monotone');A.require(core[0][0]<core[4][0],'negative separator persists for changed target')
    ad=[ambient[z][0]-ambient[2+z][0] for z in range(2)];A.require(ad[0]==-ad[1] and ad[0]<0 and ambient[0][0]<ambient[4][0],'ambient duality and dominating target')
    data={'source_rows':T,'baseline_rows':U,'dominating_target_rows':V,'child_sets':blocks,'core_queries':core_queries,'closed_source_profile':source_profile,'normalized_mod16':[x*pow(source_profile[0],-1,16)%16 for x in source_profile],'source_child_profile':values[1],'target_child_profile':values[2],'core_augmented_profiles':core,'ambient_H':ambient,'baseline_core_contrasts':contrasts,'baseline_ambient_contrasts':ad}
    (P/'extensions_controls.json').write_text(json.dumps(data,sort_keys=True,indent=2)+'\n');(P/'extensions_input.txt').write_text(inp);(P/'extensions_output.txt').write_text(out)
    receipt={'status':'PASS','literal_graph_queries':len(req),'closure_normalized_profile':[1,14,4,8,0],'baseline_ambient_contrasts':ad,'sources':{f:hashlib.sha256((P/f).read_bytes()).hexdigest() for f in ['extensions_audit.py','audit.py','physical_audit.py','physical.cpp']},'artifacts':{f:hashlib.sha256((P/f).read_bytes()).hexdigest() for f in ['extensions_controls.json','extensions_input.txt','extensions_output.txt']}}
    (P/'extensions_receipt.json').write_text(json.dumps(receipt,sort_keys=True,indent=2)+'\n');print(json.dumps(receipt,sort_keys=True))
if __name__=='__main__':main()
