#!/usr/bin/env python3
"""Scoped literal controls for the graded boundary theorems."""
import argparse
import hashlib
import itertools
import json
from pathlib import Path


def require(p,msg):
    if not p:
        raise RuntimeError(msg)


def constants(K):
    c=[0,0]
    splits={}
    for k in range(2,K+1):
        best=min(c[i]+c[k-i]+(2*i==k) for i in range(1,k//2+1))
        c.append(k+best)
        splits[k]=[i for i in range(1,k//2+1) if c[i]+c[k-i]+(2*i==k)==best]
    return c,splits


def step(points,K):
    result=set(points)
    for u,v in itertools.combinations(sorted(points),2):
        if u[1]+v[1] <= K:
            result.add((u[0]+v[0],u[1]+v[1]))
    return result


def run():
    K,T=12,24
    c,splits=constants(64)
    points={(1,0),(0,1)}
    row_checks=0
    boundary=[]
    for t in range(T+1):
        if t:
            points=step(points,K)
        for u,k in points:
            if k:
                require(u <= k*t-c[k],'all-round height bound')
        for k in range(1,min(t,K)+1):
            row={u for u,j in points if j==k and u>0}
            require(row == set(range(1,k*t-c[k]+1)),'complete affine row')
            row_checks += 1
        if t in [0,1,2,3,4,8,12,24]:
            boundary.append({'t':t,'rows':{str(k):sorted(u for u,j in points if j==k)
                                        for k in range(1,min(K,4)+1)}})
    # Independent integer multiplication, without using valuation transitions.
    integers={2,3}
    vectors={(1,0),(0,1)}
    for t in range(5):
        require(integers == {2**u*3**v for u,v in vectors},'literal integer conjugacy')
        if t<4:
            integers |= {a*b for a,b in itertools.combinations(sorted(integers),2)}
            vectors=step(vectors,100)

    controls=[]
    for seed in [{(1,0),(0,1),(2,1)},
                 {(2,0),(0,1),(1,1)},
                 {(2,0),(0,1),(3,2),(1,3)}]:
        a=next(x for x,y in seed if y==0)
        points=set(seed)
        late={}
        for t in range(17):
            if t:
                points=step(points,4)
            if t>=12:
                rows={}
                for x,k in points:
                    if k:
                        r=x%a
                        rows.setdefault((k,r),set()).add((x-r)//a)
                for key,row in rows.items():
                    require(row == set(range(min(row),max(row)+1)),'late complete coset row')
                late[t]={key:(min(row),max(row)) for key,row in rows.items()}
        require(all(set(late[t])==set(late[12]) for t in late),'stable finite coset set')
        forecasts=[]
        for key,(lo,hi) in sorted(late[12].items()):
            slope=late[13][key][1]-hi
            intercept=hi-12*slope
            require(1 <= slope <= key[0],'grade-bounded affine slope')
            for t in [14,15,16]:
                require(late[t][key] == (lo,slope*t+intercept),'finite affine forecast control')
            forecasts.append({'grade':key[0],'residue':key[1],'minimum':lo,
                              'slope':slope,'intercept':intercept})
        controls.append({'seeds':sorted(seed),'zero_face_seed':a,'forecasts':forecasts})

    return {'status':'FINITE-EXACT controls; all-time scope is proved in BOUNDARY_COMPILER.md',
            'universe':{'free_strip_grade':K,'free_generations':T,
                        'general_strip_grade':4,'general_generations':16,
                        'forecast_fit_times':[12,13],'forecast_check_times':[14,15,16]},
            'complete_row_checks':row_checks,'constants_through64':c[1:],
            'optimal_splits_through16':{str(k):splits[k] for k in range(2,17)},
            'selected_free_rows':boundary,'general_coset_controls':controls,
            'dropped_diagonal_tax_hostile':{'target':[2,2],'predicted_without_tax':2,'actual_first_time':3}}


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--write',action='store_true')
    args=p.parse_args()
    data=json.dumps(run(),indent=2,sort_keys=True)+'\n'
    target=Path(__file__).with_name('receipt.json')
    if args.write:
        target.write_text(data)
    else:
        require(target.read_text()==data,'receipt drift')
    print(json.dumps({'passed':True,'receipt_sha256':hashlib.sha256(data.encode()).hexdigest()}))
