#!/usr/bin/env python3
"""Independent exact polynomial interpolation and blue integer Smith controls."""
from check_metagraph_independent import laplacian,determinant,check
from collections import Counter
from fractions import Fraction as F
from pathlib import Path
import json,math,time
import sympy as S
from sympy.matrices.normalforms import smith_normal_form
from symmetric_dyadic import valuations
HERE=Path(__file__).resolve().parent

def parts(n,edges):
    a=[set([i]) for i in range(n)]
    for u,v in edges:
        if a[u] is a[v]:continue
        merged=a[u]|a[v]
        for i in merged:a[i]=merged
    return sorted({tuple(sorted(c)) for c in a})

def tree(n,edges):return determinant([r[:-1] for r in laplacian(n,edges)[:-1]])

def component_product(n,edges,cc):
    result=1
    for c in cc:
        where={v:i for i,v in enumerate(c)}
        sub=Counter({(where[u],where[v]):w for (u,v),w in edges.items() if u in where and v in where})
        result*=tree(len(c),sub)
    return result

def quotient(edges,cc):
    where={v:i for i,c in enumerate(cc) for v in c};out=Counter()
    for (u,v),w in edges.items():out[tuple(sorted((where[u],where[v])))]+=w
    return out

def interpolate_integer_values(values):
    diffs=list(values);leading=[]
    while diffs:
        leading.append(diffs[0]);diffs=[b-a for a,b in zip(diffs,diffs[1:])]
    result=[F(0)]*len(values);falling=[F(1)]
    for k,delta in enumerate(leading):
        for j,x in enumerate(falling):result[j]+=F(delta,math.factorial(k))*x
        new=[F(0)]*(len(falling)+1)
        for j,x in enumerate(falling):new[j]-=k*x;new[j+1]+=x
        falling=new
    check(all(x.denominator==1 for x in result),'nonintegral interpolation')
    return list(map(int,result))

def run():
    source=json.loads((HERE.parent/'metagraph/dyadic_network_receipt.json').read_text())
    expected=json.loads((HERE.parent/'metagraph/tree_scales_receipt.json').read_text())
    polys={r['n']:r for r in expected['tree_polynomials']};bluevals={r['n']:r['blue_two_primary_exponents'] for r in expected['blue_critical_groups']}
    out=[]
    for data in source:
        start=time.time();n=data['n'];N=data['nodes'];L=Counter({(u,v):w for u,v,w in data['blue_edges']});B=Counter({(u,v):w for u,v,w in data['paired_black_edges']})
        CL=parts(N,L);CB=parts(N,B)
        vals=[];blue_taus=[]
        for c in CL:
            if len(c)==1:continue
            M=laplacian(N,L);A=S.Matrix([[M[i][j] for j in c[:-1]] for i in c[:-1]])
            if n<=6:
                D=smith_normal_form(A,domain=S.ZZ)
                diag=[abs(int(D[i,i])) for i in range(D.rows)]
                check(all(x>0 for x in diag),'blue Smith zero')
                v=[(x&-x).bit_length()-1 for x in diag]
                tau=math.prod(diag)
            else:
                v=valuations(A.tolist());tau=int(determinant(A.tolist()))
                check(sum(v)==(tau&-tau).bit_length()-1,'symmetric dyadic vs exact determinant mismatch')
            vals.extend(x for x in v if x)
            blue_taus.append(str(tau))
        vals.sort();check(vals==bluevals[n],'blue Smith values disagree')
        row=dict(n=n,blue_two_primary_exponents=vals,blue_component_tree_counts=blue_taus)
        if n<=6:
            lo=len(CB)-1;hi=N-len(CL)
            M=laplacian(N,L);K=laplacian(N,B)
            values=[determinant([[t*M[i][j]+2*K[i][j] for j in range(N-1)] for i in range(N-1)]) if t else (0 if len(CB)>1 else tree(N,Counter({e:2*w for e,w in B.items()}))) for t in range(hi+1)]
            # At t=0, disconnected paired-black graph can make Bareiss pivot zero.
            coef=interpolate_integer_values(values)
            check(coef==polys[n]['coefficients_by_blue_count'],'tree polynomial disagreement')
            low=2**(N-len(CB))*component_product(N,B,CB)*tree(len(CB),quotient(L,CB))
            high=2**(len(CL)-1)*component_product(N,L,CL)*tree(len(CL),quotient(B,CL))
            check(coef[lo]==low and coef[hi]==high,'endpoint mismatch')
            row.update(coefficients=coef,lower_endpoint=low,upper_endpoint=high)
        row['elapsed_seconds']=round(time.time()-start,3);out.append(row);print(row,flush=True)
        (HERE/'tree_controls_receipt.json').write_text(json.dumps(dict(status='PASS',rows=out),indent=2)+'\n')
    print('ALL INDEPENDENT TREE CONTROLS PASSED')
if __name__=='__main__':run()
