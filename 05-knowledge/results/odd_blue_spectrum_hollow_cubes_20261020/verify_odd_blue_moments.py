#!/usr/bin/env python3
"""Exact auxiliary-cube triangle polynomial and moment checks.

Reads only the frozen pairing_n3/5/7/9.json files. Imports no producer code.
Reconstructs every canonical tournament's adjacency, c3 and matching slopes.
"""
import argparse
from collections import Counter, defaultdict
from itertools import combinations
from math import comb, factorial
from pathlib import Path
import hashlib, json


def frame_data(n, mask):
    k=(n-1)//2
    edges=list(combinations(range(k),2));e=len(edges)
    a=[0]*n
    def arc(i,j,bit):
        if bit:a[i]|=1<<j
        else:a[j]|=1<<i
    for i in range(k):
        arc(i,k,1);arc(k,k+1+i,1)
        arc(i,k+1+i,(mask>>i)&1)
    for z,(i,j) in enumerate(edges):
        bit=(mask>>(k+z))&1
        arc(i,j,bit);arc(k+1+j,k+1+i,bit)
        bit=(mask>>(k+e+z))&1
        arc(i,k+1+j,bit);arc(j,k+1+i,bit)
    c3=comb(n,3)-sum(comb(row.bit_count(),2) for row in a)
    slopes=[2*k-1-2*(a[i].bit_count()-((mask>>i)&1)) for i in range(k)]
    return c3,slopes


def run(n,data_dir):
    path=data_dir/f'pairing_n{n}.json'
    data=json.loads(path.read_text());rows=data['rows'];k=(n-1)//2
    packets=defaultdict(list)
    for row in rows:
        direct,coefficients=frame_data(n,row['mask'])
        assert direct==row['c3']
        assert coefficients==row['triangle_coefficients']
        packets[row['cube_packet']].append(row)
    reports=[]
    for root,vertices in sorted(packets.items()):
        anchor=next(row for row in vertices if row['id']==root)
        coefficients=anchor['triangle_coefficients'];orbits=anchor['hollow_pair_orbits']
        steps=[]
        for orbit_mask in orbits:
            indices=[i for i in range(k) if (orbit_mask>>i)&1]
            assert len(indices)%2==1
            assert len({coefficients[i] for i in indices})==1
            signed_step=sum(coefficients[i]*(1-2*((anchor['mask']>>i)&1)) for i in indices)
            assert signed_step%2==1
            steps.append(signed_step)
        dim=len(steps);positive_steps=sorted(abs(step) for step in steps)
        cmin=anchor['c3']+sum(min(0,step) for step in steps)
        polynomial=Counter({cmin:1})
        for step in positive_steps:
            shifted=Counter({c+step:count for c,count in polynomial.items()})
            polynomial.update(shifted)
        actual=Counter(row['c3'] for row in vertices)
        assert polynomial==actual
        assert cmin==min(actual) and all(step>0 and step%2==1 for step in positive_steps)
        assert len(vertices)==1<<dim
        assert len({row['blue'] for row in vertices})==1
        signed_moments=[sum((-1)**row['c3']*row['c3']**degree for row in vertices)
                        for degree in range(dim+1)]
        assert all(v==0 for v in signed_moments[:dim])
        sharp=(-1)**(cmin+dim)*factorial(dim)
        for step in positive_steps:sharp*=step
        assert signed_moments[dim]==sharp and sharp!=0
        # Exact integer polynomial division by q+1, r times, plus nonzero residue.
        dense=[polynomial.get(c,0) for c in range(max(polynomial)+1)]
        for _ in range(dim):
            quot=[0]*(len(dense)-1)
            quot[0]=dense[0]
            for i in range(1,len(quot)):quot[i]=dense[i]-quot[i-1]
            assert dense[-1]==quot[-1]
            dense=quot
        assert sum(coef*(-1)**i for i,coef in enumerate(dense))!=0
        reports.append({'packet_id':root,'B':anchor['blue'],'dimension':dim,
                        'class_count':len(vertices),'c_min':cmin,
                        'positive_odd_steps':positive_steps,'signed_moments_0_through_r':signed_moments,
                        'first_nonzero_signed_moment':sharp})
    by_B=[]
    for B in sorted({row['blue'] for row in rows}):
        selected=[row for row in rows if row['blue']==B]
        min_dim=min(row['cube_dimension'] for row in selected)
        shores=[[row for row in selected if row['c3']%2==parity] for parity in [0,1]]
        counts=[len(shore) for shore in shores]
        totals=[sum(row['c3'] for row in shore) for shore in shores]
        moments=[sum((-1)**row['c3']*row['c3']**degree for row in selected)
                 for degree in range(min_dim+1)]
        assert counts[0]==counts[1]
        assert all(value==0 for value in moments[:min_dim])
        if n%4==1:
            assert min_dim>=2 and totals[0]==totals[1]
        by_B.append({'B':B,'minimum_packet_dimension':min_dim,
                     'shore_class_counts':counts,'shore_c3_totals':totals,
                     'signed_moments_0_through_min_dimension':moments})
    if n==3:
        assert by_B[0]['shore_c3_totals']==[0,1]
        assert reports[0]['signed_moments_0_through_r']==[0,-1]
    if n==5:
        assert by_B[0]['shore_c3_totals']==[4,4]
        assert by_B[0]['signed_moments_0_through_min_dimension']==[0,0,6]
    return {'n':n,'input_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
            'classes_checked':len(rows),'packets_checked':len(reports),
            'packet_dimensions':dict(sorted(Counter(r['dimension'] for r in reports).items())),
            'all_packet_factorizations_exact':True,'all_vanishing_moments_exact':True,
            'all_packet_root_multiplicities_exact':True,
            'by_blue_degree':by_B,'packets':reports}


if __name__=='__main__':
    here=Path(__file__).resolve().parent
    default_data=here if (here/'pairing_n3.json').exists() else here.parent/'session14-odd-blue-pairing'
    ap=argparse.ArgumentParser();ap.add_argument('--data-dir',type=Path,default=default_data)
    ap.add_argument('--output',type=Path,default=here/'odd_blue_moment_verification.json')
    args=ap.parse_args()
    reports=[run(n,args.data_dir) for n in [3,5,7,9]]
    result={'status':'PASS','scope':'exact checks of frozen SC-class cube data at n3,5,7,9',
            'reports':reports,'boundaries':{'n3_first_signed_moment':-1,
                                          'n5_B1_second_signed_moment':6}}
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    for row in reports:
        print(f"n={row['n']}: {row['classes_checked']} classes, {row['packets_checked']} cube polynomials; all exact checks PASS")
