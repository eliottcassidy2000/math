#!/usr/bin/env python3
from pathlib import Path
from collections import Counter
import json,sys,re
from audit_core import *
D=Path(__file__).resolve().parent
PACKET=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else D.parent

def printed_words(text,n):
    formulas=re.findall(r'Q[123]=\(([^)]+)\)',text)
    need(len(formulas)==6,'six printed A/B words')
    symbols={'a':[0],'b':[1],'x':[2],'W':list(range(3,n-3)),'y':[n-3],'c':[n-2],'z':[n-1]}
    return [sum((symbols[item.strip()] for item in f.split(',')),[]) for f in formulas]

def proof_controls(n,text):
    A,B,C=printed_words(text,n)[:3];r1,r2,r3=map(reversal,[A,B,C]);z=n-1;c=n-2;y=n-3
    need((r1[z],r1[3])==(2,c),'Q1 printed anti witness')
    need((r2[5],r2[0])==(2,3),'Q2 ordinary printed witness')
    need((r3[3],r3[0])==(2,z),'Q3 ordinary printed witness')
    need((r3[z],r3[n-4])==(0,5),'Q3 h3 printed witness')
    chain=[(0,4),(1,3),(n-4,c),(7,z),(0,n-8)];E=path_edges(A);parities=[]
    for k,((u,v),target) in enumerate(zip(chain,chain[1:])):
        if k%2:need(target==(n-1-v,n-1-u),'P printed chain link');parities.append(0)
        else:
            p,q=r2[v],r2[u];need(pair(p,q)==target,'Q2 printed chain link')
            parities.append(int(p>q)^int((u,v) in E)^int(target in E))
    need(parities==[0,0,1,0],'h5 chain exact individual parities')
    need(all(v not in {5,y,c,z} for v in [4,n-8]),'h5 actual port contradiction')

def main():
    text=(PACKET/'README.md').read_text();counts=Counter()
    for n in range(15,81):
        qs=printed_words(text,n)[:3]
        for h in range(3,n-3):
            if h==n-4:continue
            I={h,n-3,n-2};A=symbolic_source(n,I,[(n-2,1,1),(n-2,2,1),(2,n-3,1)])
            replay_family(A,qs);counts['all_h_affine_cells']+=1
    for n in range(15,301):proof_controls(n,text)
    hostiles=[]
    for hostile in json.loads((PACKET/'witnesses.json').read_text()):
        n=hostile['n'];h=hostile['h'];A=literal_rows(hostile['rows']);qs=printed_words(text,n)[3:]
        need(qs==hostile['words'],'printed B words differ from witness words')
        need(A[n-2][1]==A[n-2][2]==A[n-3][2]==1,'opposite-chord guards')
        flags=check_literal(A,{h,n-3,n-2},qs)
        need(flags[-1] and not any(blue_equations(k_move(A,list(range(n))),qs[-1])),'literal final word not blue at both endpoints')
        hostiles.append({'n':n,'h':h,'actual_blue_flags':flags,'printed_words_match':True})
    print(json.dumps({'status':'PASS','method':'independent GF(2) affine source families, printed-formula parsing, and exact labelled rows','all_h_orders':[15,80],'counts':dict(counts),'individual_witness_and_chain_orders':[15,300],'literal_hostiles':hostiles,'scope':'Guarded c->x and x->y only. Opposite-chord examples refute only the printed final-word template.'},sort_keys=True,indent=2))
if __name__=='__main__':main()
