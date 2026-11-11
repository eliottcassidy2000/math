"""Exact all-H9-family extension content: independent subset-path controls."""
from itertools import combinations
from math import gcd
from pathlib import Path
import json
from h9_fibres import rows as apex_rows, need

C3=(2,4,1)
S4=(2,4,9,3)

def tt(n):
    return tuple(sum(1<<j for j in range(i+1,n)) for i in range(n))

def join(blocks):
    n=sum(map(len,blocks)); out=[]; offset=0
    for block in blocks:
        later=((1<<n)-1)^((1<<(offset+len(block)))-1)
        out.extend((row<<offset)|later for row in block)
        offset+=len(block)
    return tuple(out)

def path_table(adj):
    n=len(adj); size=1<<n; dp=[[0]*n for _ in range(size)]
    for i in range(n): dp[1<<i][i]=1
    for mask in range(1,size):
        for i,value in enumerate(dp[mask]):
            if value:
                bits=adj[i] & (size-1-mask)
                while bits:
                    bit=bits & -bits; bits-=bit
                    dp[mask|bit][bit.bit_length()-1]+=value
    return dp

def h(adj): return sum(path_table(adj)[-1])

def extend(adj,mask):
    n=len(adj)
    return tuple(row | ((1<<n) if not(mask>>i&1) else 0)
                 for i,row in enumerate(adj))+(mask,)

def coefficients(adj):
    n=len(adj); full=(1<<n)-1
    end=path_table(adj)
    reverse=tuple(sum(1<<j for j in range(n) if adj[j]>>i&1) for i in range(n))
    start=path_table(reverse); k=[[0]*n for _ in range(n)]
    for mask in range(1,full):
        e=end[mask]; s=start[full^mask]
        for i,ei in enumerate(e):
            if ei:
                for j,sj in enumerate(s):
                    if sj: k[i][j]+=ei*sj
    base=sum(end[-1])
    linear=[start[-1][j]-end[-1][j]+sum(k[i][j] for i in range(n)) for j in range(n)]
    quadratic={(i,j):-k[i][j]-k[j][i] for i in range(n) for j in range(i+1,n)}
    content=base
    for value in linear+list(quadratic.values()): content=gcd(content,value)
    return base,linear,quadratic,content,k

def evaluate(coeff,mask):
    base,lin,quad,*_=coeff
    return base+sum(v for i,v in enumerate(lin) if mask>>i&1)+sum(
        v for (i,j),v in quad.items() if mask>>i&1 and mask>>j&1)

def main():
    census=[]; matrices=[]; exact_graph_checks=0; family_checks=0
    for n in range(5,13):
        records=[]
        for v in ('OIII','OIOI'):
            for a in range(n-4):
                b=n-5-a; adj=join((tt(a),apex_rows(v),tt(b)))
                cf=coefficients(adj)
                need(cf[0]==9 and cf[3]==1,'one-apex source content')
                records.append({'kind':v,'padding':[a,b],'C1':1})
                family_checks+=1
        for a in range(n-5):
            for b in range(n-5-a):
                c=n-6-a-b
                adj=join((tt(a),C3,tt(b),C3,tt(c)))
                cf=coefficients(adj); expected=3 if b else 1
                need(cf[0]==9 and cf[3]==expected,'two-triangle source content')
                i=a; j=a+3+b
                need(cf[4][i][j]==(1 if b==0 else 0),'forward cross cover')
                need(cf[4][j][i]==9*2**b,'reverse cross cover')
                need(cf[2][(a,a+1)]==-6,'internal triangle coefficient')
                # Independent actual graph counts for every cut at n<=7;
                # selected low-weight and hostile cuts at n8,9.
                masks=range(1<<n) if n<=7 else (
                    [0,(1<<n)-1,1<<i,1<<j,(1<<i)|(1<<j)] if n<=9 else [])
                for mask in masks:
                    need(evaluate(cf,mask)==h(extend(adj,mask)),'literal extension polynomial')
                    exact_graph_checks+=1
                if b>=1:
                    escaped=join((tt(a),S4,tt(b),S4,tt(c)))
                    need(h(escaped)==25,'actual two-birth escape')
                    # Delete the two fourth vertices to recover the exact source.
                    keep=[x for x in range(n+2) if x not in (a+3,a+4+b+3)]
                    induced=tuple(sum((escaped[u]>>v&1)<<q for q,v in enumerate(keep))
                                  for u in keep)
                    need(induced==adj,'old graph retained by escape')
                    exact_graph_checks+=1
                records.append({'kind':'C3_C3','padding':[a,b,c],'C1':expected})
                if n<=8: matrices.append({'padding':[a,b,c],'rows':adj,'C1':expected})
                family_checks+=1
        expected_total=2*(n-4)+(n-4)*(n-5)//2
        expected_nonunit=(n-5)*(n-6)//2 if n>=7 else 0
        need(len(records)==expected_total,'ordinary H9 count')
        need(sum(x['C1']==3 for x in records)==expected_nonunit,'nonunit source count')
        census.append({'n':n,'classes':len(records),'C1_3':expected_nonunit,
                       'C1_1':len(records)-expected_nonunit})
    # Explicit one-apex extensions and the minimal hostile.
    need(h(apex_rows('OIIII'))==17 and h(apex_rows('OIOII'))==19,'core escape controls')
    triangle_controls=[]
    gap_cases={1:[(0,0),(1,2)],2:[(0,0,0),(0,1,0),(1,2,1)],
               3:[(0,0,0,0),(0,1,0,0),(0,1,1,0),(1,1,1,1)],
               4:[(0,0,0,0,0),(0,1,0,0,0),(0,1,1,1,0)]}
    for r, cases in gap_cases.items():
        for gaps in cases:
            blocks=[tt(gaps[0])]
            for j in range(r): blocks.extend((C3,tt(gaps[j+1])))
            adj=join(blocks); cf=coefficients(adj)
            adjacent=any(g==0 for g in gaps[1:-1])
            expected=3**(r-2 if r>=2 and adjacent else r-1)
            need(cf[0]==3**r and cf[3]==expected,'arbitrary triangle-chain content')
            triangle_controls.append({'triangles':r,'gaps':gaps,'vertices':len(adj),
                                      'H':cf[0],'C1':cf[3]})
    out={'status':'FINITE-EXACT controls for all-order structural proof',
      'family_checks':family_checks,'independent_literal_extension_checks':exact_graph_checks,
      'census':census,'small_source_matrices':matrices,
      'nonunit_H9_generating_function':'z^7/(1-z)^3',
      'all_H9_generating_function':'2*z^5/(1-z)^2+z^6/(1-z)^3',
      'two_birth_escape_counts':[9,25],'terminal_content':1}
    out['arbitrary_triangle_chain_controls']=triangle_controls
    Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'result':'PASS','families':family_checks,'literal_checks':exact_graph_checks}))

if __name__=='__main__': main()
