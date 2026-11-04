#!/usr/bin/env python3
"""Finite truncation diagnostic; full cubic and quartic data are never discarded from the bank."""
from pathlib import Path
from collections import defaultdict,Counter
from itertools import combinations
import json,gzip
from check_six import count,need
from verify_ordered import tournament,modules,maximal
P=Path(__file__).resolve().parent;N=32768

def main():
    f=P/'mod4_six_vectors.bin';raw=f.read_bytes() if f.exists() else gzip.decompress((P/'mod4_six_vectors.bin.gz').read_bytes())
    need(len(raw)==15859712,'bank size')
    B=[m for m in range(4096) if (m&63).bit_count()<=2 and (m>>6).bit_count()<=2]
    classes=[];pairs=[]
    for degree in range(1,5):
        ix=[i for i,m in enumerate(B) if m.bit_count()<=degree];groups=defaultdict(list)
        for code in range(N):groups[bytes(raw[i*N+code] for i in ix)].append(code)
        classes.append(dict(maximum_degree=degree,coordinates=len(ix),classes=len(groups),class_size_histogram=dict(sorted(Counter(map(len,groups.values())).items()))))
        if degree==2:
            pairs=[v for v in groups.values() if len(v)>1];need(all(len(v)==2 for v in pairs),'quadratic collision pair size')
    need(len(pairs)==240 and classes[2]['classes']==N and classes[3]['classes']==N,'finite degree counts')
    categories=Counter();module_records=[]
    for a,b in pairs:
        T=tournament(6,a);mods=modules(T);basis={}
        for i,m in enumerate(mods):
            x=sum(1<<j for j,(u,v) in enumerate(combinations(range(6),2)) if m>>u&1 and m>>v&1);z=1<<i
            while x:
                h=x.bit_length()-1
                if h in basis:x^=basis[h][0];z^=basis[h][1]
                else:basis[h]=(x,z);break
            need(bool(x),'module edge-incidence independence')
        x=a^b;z=0
        while x:
            h=x.bit_length()-1
            if h not in basis:break
            x^=basis[h][0];z^=basis[h][1]
        if x:categories['outside_original_module_span']+=1;continue
        fs=[m for i,m in enumerate(mods) if z>>i&1]
        lam=not any(u&v and u&v not in (u,v) for u,v in combinations(fs,2))
        covered=lam and all(sum(maximal([a for a in fs if a!=D and a&D==a]))==D for D in maximal(fs))
        need(covered,'module-span pair is fully covered laminar')
        categories['fully_covered_laminar']+=1;module_records.append([a,b,fs])
    need(categories=={'fully_covered_laminar':120,'outside_original_module_span':120},'module family classification')
    a,b=pairs[0];difference=[[m&63,m>>6,(raw[i*N+a]-raw[i*N+b])%4] for i,m in enumerate(B) if raw[i*N+a]!=raw[i*N+b]]
    need((a,b)==(408,15462),'canonical first pair')
    fa,fb=count(a,6,1),count(b,6,1);need((fa,fb)==(455,477),'literal cubic witness')
    return dict(status='PASS',scope='Finite modulo-four truncation diagnostic only; the complete response bank is retained.',degree_classes=classes,quadratic_collision_pairs=pairs,original_module_classification=dict(categories),module_family_records=module_records,first_pair=dict(source_code=a,target_code=b,source_rows=tournament(6,a),target_rows=tournament(6,b),modules=[3,12,48,63],complete_nonzero_mod4_difference=difference,literal_query=[6,1],literal_counts=[fa,fb]))
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,separators=(',',':')))
