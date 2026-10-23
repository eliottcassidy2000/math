#!/usr/bin/env python3
"""Finish n8's unique dyadic tail with certified integer lattice corrections.
M_k = 2 M_(k-1) + span_Z(chosen lifts of V_k). This gives exact lower-bit
corrections using images A U_i / 2^i mod2, never uncorrected divisions.
"""
import gzip,json,hashlib,argparse
from pathlib import Path
from check_electrical import projection,bits,rank2,need
from check_n8_lift import Av
HERE=Path(__file__).resolve().parent

def main(base,out,maxlevel):
    rec2=json.loads((HERE/'n8_lift_receipt.json').read_text());rec3=json.loads((HERE/'n8_final_radical_receipt.json').read_text());n=rec2['vertices']
    with gzip.open(base/'edges.tsv.gz','rt') as f:
        next(f);tab=[list(map(int,line.split())) for line in f]
    blue=[(a,b,w) for a,b,w,z in tab if w];black=[(a,b,z//2) for a,b,w,z in tab if z];edges=[(a,b,w+z) for a,b,w,z in tab]
    p=projection(n,blue,black);keep=p['_keep'];ix={v:i for i,v in enumerate(keep)};basis1=p['_basis'];r=len(keep)
    adj=[{} for _ in range(n)]
    for a,b,w in edges:
        if a!=b:adj[a][a]=adj[a].get(a,0)+w;adj[b][b]=adj[b].get(b,0)+w;adj[a][b]=adj[a].get(b,0)-w;adj[b][a]=adj[b].get(a,0)-w
    lifts2=[]
    for rec in rec2['witness_lifts']:
        v=[0]*n
        for i in rec['binary_one_nodes']:v[i]+=1
        for i in rec['correction_one_nodes']:v[i]+=2
        need(Av(n,edges,v)==[4*z for z in rec['A_lift_divided_by_four']],'replay U2');lifts2.append(v)
    # Generators retain (level, sparse integer vector, image bitmask).
    generators=[]
    for v in keep:
        image=sum((w%2)<<ix[u] for u,w in adj[v].items() if u in ix)
        generators.append((0,[(v,1)],image))
    for z in basis1:
        supp=[(keep[i],1) for i in bits(z)];value={}
        for v,coeff in supp:
            for u,w in adj[v].items():value[u]=value.get(u,0)+coeff*w
        need(all(v%2==0 for v in value.values()),'U1 image divisible2')
        image=sum(((w//2)%2)<<ix[u] for u,w in value.items() if u in ix)
        generators.append((1,supp,image))
    for v in lifts2:
        av=Av(n,edges,v);need(all(w%4==0 for w in av),'U2 image divisible4')
        generators.append((2,[(i,w) for i,w in enumerate(v) if w],sum(((av[u]//4)%2)<<ix[u] for u in keep)))
    # Binary image elimination also records exactly which lattice generators were used.
    piv={}
    for j,(_,_,image) in enumerate(generators):
        witness=1<<j;z=image
        while z:
            top=z.bit_length()-1
            if top in piv:z^=piv[top][0];witness^=piv[top][1]
            else:piv[top]=(z,witness);break
    need(len(piv)==r-1,'one missing higher-valuation image')
    def solve(image):
        witness=0;z=image
        while z:
            top=z.bit_length()-1
            need(top in piv,'higher correction not in certified image span')
            z^=piv[top][0];witness^=piv[top][1]
        need(image==__import__('functools').reduce(int.__xor__,(generators[j][2] for j in bits(witness)),0),'image solve replay')
        return witness
    lift=rec3['lift_mod8'];need(Av(n,edges,lift)==[8*z for z in rec3['A_lift_divided_by_eight']],'replay U3');k=3;layers=[{'k':3,'q':0}];certs=[];complete=False
    while k<maxlevel:
        av=Av(n,edges,lift);mod=1<<k;need(all(z%mod==0 for z in av),'previous qualifying lift')
        image=sum(((av[u]//mod)%2)<<ix[u] for u in keep);selected=solve(image)
        correction=[0]*n
        for j in bits(selected):
            level,support,_=generators[j]
            for v,w in support:correction[v]+=(1<<(k-level))*w
        need(all(z%2==0 for z in correction),'binary potential unchanged')
        k+=1;mod=1<<k;lift=[(v-c)%mod for v,c in zip(lift,correction)];av=Av(n,edges,lift)
        need(all(z%mod==0 for z in av),'certified next dyadic lift')
        energy=sum(a*b for a,b in zip(lift,av));need(energy%mod==0,'next half form integral');q=energy//mod%2
        layers.append({'k':k,'q':q});certs.append({'k':k,'lift':lift,'A_lift_divided_by_modulus':[z//mod for z in av]})
        if q:complete=True;break
    result={'n':8,'complete_two_primary':complete,'certified_through_k':k,'tail_layers':layers,'generator_image_rank':len(piv),'generator_counts':{str(level):sum(i==level for i,_,_ in generators) for level in range(3)},'witness_lifts':certs,'source_sha256':{file.name:hashlib.sha256(file.read_bytes()).hexdigest() for file in [base/'edges.tsv.gz',HERE/'n8_lift_receipt.json',HERE/'n8_final_radical_receipt.json',Path(__file__)]}}
    if complete:result['two_primary_exponent_counts']={'1':3358,'2':4,str(k):1};result['spanning_tree_valuation_two']=3358+8+k
    else:result['remaining_factor_valuation_at_least']=k+1
    out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='witness_lifts'},indent=2));need(complete,'declared precision exhausted without completion');print('ALL CHECKS PASSED')
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--base',type=Path,default=HERE.parent/'session38-wiring');ap.add_argument('--out',type=Path,default=HERE/'n8_complete_tower_receipt.json');ap.add_argument('--max-level',type=int,default=16);a=ap.parse_args();main(a.base,a.out,a.max_level)
