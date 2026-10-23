#!/usr/bin/env python3
"""Correct lower bits to lift the unique n8 q2 radical from mod4 to mod8."""
import gzip,hashlib,json,argparse
from pathlib import Path
from check_electrical import projection,nullbasis2,rank2,bits,need
from check_n8_lift import Av,solve2
HERE=Path(__file__).resolve().parent

def main(base,out):
    source=json.loads((HERE/'n8_lift_receipt.json').read_text());n=source['vertices']
    with gzip.open(base/'edges.tsv.gz','rt') as f:
        next(f);tab=[list(map(int,line.split())) for line in f]
    blue=[(a,b,w) for a,b,w,z in tab if w];black=[(a,b,z//2) for a,b,w,z in tab if z];edges=[(a,b,w+z) for a,b,w,z in tab]
    p=projection(n,blue,black);Q=p['_q'];basis=p['_basis'];keep=p['_keep'];ix={v:i for i,v in enumerate(keep)};r=len(keep)
    A2=[0]*r
    for a,b,w in blue:
        if w%2:
            d=(1<<ix[a] if a in ix else 0)^(1<<ix[b] if b in ix else 0)
            if a in ix:A2[ix[a]]^=d
            if b in ix:A2[ix[b]]^=d
    prev=[]
    for rec in source['witness_lifts']:
        v=[0]*n
        for i in rec['binary_one_nodes']:v[i]+=1
        for i in rec['correction_one_nodes']:v[i]+=2
        need(Av(n,edges,v)==[4*z for z in rec['A_lift_divided_by_four']],'previous lift replay');prev.append(v)
    q2=[sum(v<<j for j,v in enumerate(row)) for row in source['q2_matrix']];rad=nullbasis2(q2,len(q2));need(len(rad)==1,'one final radical')
    v=[sum(prev[j][i] for j in bits(rad[0])) for i in range(n)];av=Av(n,edges,v);need(all(z%4==0 for z in av),'mod4 radical lift');b=[z//4 for z in av]
    rhs=[]
    for z in basis:rhs.append(sum(b[keep[i]] for i in bits(z))%2)
    coeff=solve2(Q,rhs,len(Q))
    z=[0]*n
    for j in bits(coeff):
        for i in bits(basis[j]):z[keep[i]]+=1
    az=Av(n,edges,z);need(all(a%2==0 for a in az),'kernel combination')
    w=solve2(A2,[(-b[i]-az[i]//2)%2 for i in keep],r)
    t=z[:]
    for i in bits(w):t[keep[i]]+=2
    at=Av(n,edges,t);need(all((at[i]+2*b[i])%4==0 for i in range(n)),'correct mod-four correction equation')
    lift=[(v[i]+2*t[i])%8 for i in range(n)];alift=Av(n,edges,lift);need(all(z%8==0 for z in alift),'certified actual mod-eight lift')
    energy=sum(x*y for x,y in zip(lift,alift));need(energy%8==0,'q3 integral');q3=(energy//8)%2
    result={'n':8,'q3_matrix':[[q3]],'factors_divisible_by_sixteen':1-q3,'complete_two_primary':q3==1,'certified_modulus':8,'q2_radical_coefficients':list(bits(rad[0])),'lift_mod8':lift,'A_lift_divided_by_eight':[z//8 for z in alift],'source_sha256':{str(file):hashlib.sha256(file.read_bytes()).hexdigest() for file in [base/'edges.tsv.gz',HERE/'n8_lift_receipt.json',Path(__file__)]}}
    if q3:
        result['two_primary_exponent_counts']={'1':3358,'2':4,'3':1};result['spanning_tree_valuation_two']=3369
    out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k not in ['lift_mod8','A_lift_divided_by_eight']},indent=2));print('ALL CHECKS PASSED')
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--base',type=Path,default=HERE.parent/'session38-wiring');ap.add_argument('--out',type=Path,default=HERE/'n8_final_radical_receipt.json');a=ap.parse_args();main(a.base,a.out)
