#!/usr/bin/env python3
"""Replay native c3 port current independently against the full n8 Laplacian."""
import argparse,csv,gzip,hashlib,json
from pathlib import Path
HERE=Path(__file__).resolve().parent

def need(p,msg):
    if not p:raise RuntimeError(msg)

def main(base,out):
    with gzip.open(base/'nodes.tsv.gz','rt') as f:rows=[{k:int(v) for k,v in row.items()} for row in csv.DictReader(f,delimiter='\t')]
    n=len(rows);image=[0]*n;degree=[0]*n;bluedegree=[0]*n;curr=[r['current_c3'] for r in rows];phi=[r['c3'] for r in rows]
    with gzip.open(base/'edges.tsv.gz','rt') as f:
        next(f)
        for line in f:
            u,v,b,k=map(int,line.split());need(k%2==0,'literal black twins')
            degree[u]+=b+k;degree[v]+=b+k;bluedegree[u]+=b;bluedegree[v]+=b
            flow=(b+k)*(phi[u]-phi[v]);image[u]+=flow;image[v]-=flow
            if b:need((phi[u]-phi[v])%2==0,'even-order blue c3 parity')
    need(all(image[i]==-curr[i] for i in range(n)),'intrinsic port current matches negative Laplacian')
    need(sum(curr)==0,'net source current zero')
    need(all(degree[i]==rows[i]['tilings'] for i in range(n)),'full fibre degree')
    need(all(bluedegree[i]==rows[i]['blue_tilings'] for i in range(n)),'blue fibre degree')
    need(all(bluedegree[i]%2==rows[i]['SC'] for i in range(n)),'SC parity boundary')
    need(all(r['c3']%2==0 for r in rows if r['SC']),'even-order SC cyclic triangles even')
    result={'n':8,'nodes':n,'native_c3_currents_verified':n,'net_current':sum(curr),'fibre_and_blue_degrees_verified':n,'self_converse_parity_verified':n,'SC_c3_even_count':sum(r['SC'] for r in rows),'source_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [base/'nodes.tsv.gz',base/'edges.tsv.gz',Path(__file__)]}}
    out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2));print('ALL CHECKS PASSED')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--base',type=Path,default=HERE.parent/'session38-wiring');p.add_argument('--out',type=Path,default=HERE/'n8_currents_receipt.json');a=p.parse_args();main(a.base,a.out)
