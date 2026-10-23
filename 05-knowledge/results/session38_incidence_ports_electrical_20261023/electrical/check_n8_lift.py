#!/usr/bin/env python3
"""Certify mod-four lifts and q2 directly on the full actual n8 ledger.
Imports only this session's independent graph arithmetic, no inherited producer.
"""
import gzip,json,hashlib,time,argparse
from pathlib import Path
from check_electrical import projection,bits,nullbasis2,rank2,need,public
HERE=Path(__file__).resolve().parent

def solve2(rows,b,n):
    piv={}
    for row,rhs in zip(rows,b):
        z=row|(int(rhs)<<n)
        while z&((1<<n)-1):
            p=((z&((1<<n)-1))&-(z&((1<<n)-1))).bit_length()-1
            if p in piv:z^=piv[p]
            else:piv[p]=z;break
        else:need(not(z>>n&1),'inconsistent binary lift equation')
    x=0
    for p,z in sorted(piv.items(),reverse=True):
        if ((z>>n)&1)^((z&x).bit_count()&1):x|=1<<p
    need(all(((row&x).bit_count()&1)==int(rhs) for row,rhs in zip(rows,b)),'binary solve replay')
    return x

def Av(n,edges,v):
    out=[0]*n
    for a,b,w in edges:
        z=w*(v[a]-v[b]);out[a]+=z;out[b]-=z
    return out

def main(base,out):
    with gzip.open(base/'nodes.tsv.gz','rt') as f:n=sum(1 for _ in f)-1
    with gzip.open(base/'edges.tsv.gz','rt') as f:
        next(f);tab=[list(map(int,line.split())) for line in f]
    blue=[(a,b,w) for a,b,w,z in tab if w];black=[(a,b,z//2) for a,b,w,z in tab if z];edges=[(a,b,w+z) for a,b,w,z in tab]
    p=projection(n,blue,black);Q=p['_q'];basis=p['_basis'];keep=p['_keep'];ix={v:i for i,v in enumerate(keep)};r=len(keep)
    rad=nullbasis2(Q,len(Q));need(len(rad)==p['factors_divisible_by_four'],'q1 nullity')
    A2=[0]*r
    for a,b,w in blue:
        if w%2:
            d=(1<<ix[a] if a in ix else 0)^(1<<ix[b] if b in ix else 0)
            if a in ix:A2[ix[a]]^=d
            if b in ix:A2[ix[b]]^=d
    lifts=[];witness=[]
    for c in rad:
        raw=0
        for i in bits(c):raw^=basis[i]
        x=[0]*n
        for i in bits(raw):x[keep[i]]=1
        ax=Av(n,edges,x);need(all(z%2==0 for z in ax),'V1 kernel')
        y=solve2(A2,[(-ax[v]//2)%2 for v in keep],r)
        lift=x[:]
        for i in bits(y):lift[keep[i]]+=2
        alift=Av(n,edges,lift);need(all(z%4==0 for z in alift),'certified actual mod-four lift')
        lifts.append(lift);witness.append({'binary_one_nodes':[i for i,v in enumerate(x) if v],'correction_one_nodes':[keep[i] for i in bits(y)],'A_lift_divided_by_four':alift and [z//4 for z in alift]})
    q2=[]
    for x in lifts:
        row=0
        for j,y in enumerate(lifts):
            ay=Av(n,edges,y);num=sum(a*b for a,b in zip(x,ay));need(num%4==0,'q2 integral')
            if num//4%2:row|=1<<j
        q2.append(row)
    q2mat=[[(row>>j)&1 for j in range(len(q2))] for row in q2]
    need(all(q2mat[i][j]==q2mat[j][i] for i in range(len(q2)) for j in range(len(q2))),'q2 symmetry')
    rank=rank2(q2);remaining=len(q2)-rank
    result={'status':'FINITE-EXACT certified mod-four lift and second intrinsic half-energy','n':8,**public(p),'q2_matrix':q2mat,'q2_rank':rank,'factors_divisible_by_eight':remaining,'complete_two_primary':remaining==0,'certified_modulus':4,'witness_lifts':witness,'source_sha256':{str(file):hashlib.sha256(file.read_bytes()).hexdigest() for file in [base/'nodes.tsv.gz',base/'edges.tsv.gz',Path(__file__),HERE/'check_electrical.py']}}
    if remaining==0:
        result['two_primary_exponent_counts']={'1':p['half_form_rank'],'2':rank}
        result['spanning_tree_valuation_two']=p['half_form_rank']+2*rank
    out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='witness_lifts'},indent=2));print('ALL CHECKS PASSED')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--base',type=Path,default=HERE.parent/'session38-wiring');ap.add_argument('--out',type=Path,default=HERE/'n8_lift_receipt.json');args=ap.parse_args();main(args.base,args.out)
