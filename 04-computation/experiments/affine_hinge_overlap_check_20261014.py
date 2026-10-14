#!/usr/bin/env python3
"""Verify all exact edges and retained colour certificates. Standard library only.
Requires the earlier independent_hinge_seam_four_colouring_20261014.py beside this file.
"""
import json,sys
from itertools import combinations
from pathlib import Path
from independent_hinge_seam_four_colouring_20261014 import vertices,unit,colour

ROOT=Path(__file__).resolve().parents[2]
CERT=ROOT/'05-knowledge/results/affine_hinge_overlap_certificates_20261014.json'

def require(test,message):
    if not test:raise RuntimeError(message)

def build(m):
    p,labs,D=vertices(m);points={}
    for side in [0,1]:
        for v,lab in zip(p,labs):
            x=list(v[0]);x[0]+=side*D;w=(tuple(x),tuple(v[1]))
            points.setdefault(w,[]).append([side,lab[0]])
    vs=list(points);labels=list(points.values())
    edges=[(i,j) for i,j in combinations(range(len(vs)),2) if unit(vs[i],vs[j],D*D)]
    overlap=[i for i,lab in enumerate(labels) if len(lab)>1]
    return labels,edges,overlap

def verify_colouring(cols,n,edges):
    require(len(cols)==n and all(c in range(4) for c in cols),'invalid colour vector')
    require(all(cols[u]!=cols[v] for u,v in edges),'improper colour certificate')

def base_colour(labs):
    for side,(a,b,i,j) in labs:
        if a==0 and b==0:return (i+side-j)%3
    return None

def glued_conflicts(labels,edges):
    cols=[]
    for labs in labels:
        cs=[]
        for side,lab in labs:
            c=colour(*lab)
            if side:c={0:1,1:2,2:0,3:3}[c]
            cs.append(c)
        require(len(set(cs))==1,'glued overlap ambiguity')
        cols.append(cs[0])
    bad=[(u,v) for u,v in edges if cols[u]==cols[v]]
    require(all(cols[u]==3 for u,v in bad),'non-fourth-colour formula conflict')
    return len(bad)

def bare_patterns(overlap,edges,triangle):
    adj={v:set() for v in overlap}
    for u,v in edges:
        if u in adj and v in adj:adj[u].add(v);adj[v].add(u)
    require(all(v in adj[u] for u,v in combinations(triangle,2)),'normalizing triangle missing')
    col={v:c for c,v in enumerate(triangle)};out=[]
    def visit():
        if len(col)==len(overlap):out.append(tuple(col[v] for v in overlap));return
        v=max((v for v in overlap if v not in col),key=lambda v:(len({col[u] for u in adj[v] if u in col}),len(adj[v])))
        bad={col[u] for u in adj[v] if u in col}
        for c in range(4):
            if c not in bad:col[v]=c;visit();del col[v]
    visit();return set(out)

def main():
    cert=json.loads(CERT.read_text());out=[]
    for m in [1,2]:
        labs,edges,ov=build(m);n=len(labs)
        expected=(140,367,6) if m==1 else (1182,3653,20)
        require((n,len(edges),len(ov))==expected,'graph count mismatch')
        row={'m':m,'vertices':n,'edges':len(edges),'overlap_vertices':len(ov),'exact_pairs':n*(n-1)//2,'glued_formula_conflicts':glued_conflicts(labs,edges)}
        if m==1:
            c=cert['Y1'];require(ov==c['overlap_vertices'],'overlap order mismatch')
            patterns=bare_patterns(ov,edges,c['normalizing_triangle'])
            recorded={tuple(r['overlap_colouring']) for r in c['patterns']}
            require(patterns==recorded and len(patterns)==8,'incomplete overlap relation')
            for r in c['patterns']:
                require(r['extends'],'unexpected nonextension')
                cols=r['extension'];verify_colouring(cols,n,edges)
                require([cols[v] for v in ov]==r['overlap_colouring'],'wrong overlap extension')
            row['all_normalized_overlap_patterns']=8;row['all_extend']=True
        else:
            c=cert['Y2_all_base_canonical_three'];cols=c['colouring'];verify_colouring(cols,n,edges)
            fixed=0
            for v,lab in enumerate(labs):
                desired=base_colour(lab)
                if desired is not None:require(cols[v]==desired,'canonical base condition failed');fixed+=1
            require(fixed==30,'base count mismatch')
            row['canonical_base_vertices_fixed']=fixed
        row['verified']=True;out.append(row)
    print(json.dumps({'scope':'FINITE-EXACT; Y1 full overlap relation, Y2 one stronger fixed-base extension; no global affine colouring theorem claimed','tests':out,'no_float_filters':True,'no_sat_dependency':True},indent=2,sort_keys=True))

if __name__=='__main__':main()
