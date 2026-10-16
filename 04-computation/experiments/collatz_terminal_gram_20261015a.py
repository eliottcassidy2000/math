#!/usr/bin/env python3
"""Joint terminal Gram certificates, signed-graph controls and exact cube maximum."""
import itertools
import json
from fractions import Fraction as F
from pathlib import Path
from collatz_terminal_mesh_20261015a import certificate,T
from collatz_clock_mesh_audit_20261015a import envelope,need

def gram(ledger,cols):
    rows=[[row.get(str(e),0) for e in cols] for row in ledger]
    return [[sum(row[i]*row[j] for row in rows) for j in range(len(cols))] for i in range(len(cols))]

def balance(G):
    sig={}
    for root in range(len(G)):
        if root in sig:continue
        sig[root]=1;stack=[root]
        while stack:
            i=stack.pop()
            for j in range(len(G)):
                if i==j or not G[i][j]:continue
                sign=sig[i]*(1 if G[i][j]>0 else -1)
                if j in sig:
                    if sig[j]!=sign:return None
                else:sig[j]=sign;stack.append(j)
    # Column zero is ROOT; flip its component if needed. Root starts at +1 here.
    return [sig[i] for i in range(len(G))]

def cube_max(G):
    n=len(G);need(n<=20,'cube enumeration budget')
    x=[1]*n;h=[sum(G[i][j] for j in range(n) if j!=i) for i in range(n)]
    value=sum(sum(row) for row in G);best=value;bestx=x[:]
    # ROOT is fixed to one, the remaining signs run through Gray code.
    for step in range(1,2**(n-1)):
        i=(step&-step).bit_length();old=x[i]
        value-=4*old*h[i];x[i]=-old
        for j in range(n):
            if j!=i:h[j]-=2*old*G[j][i]
        if value>best:best=value;bestx=x[:]
    need(best==sum(bestx[i]*G[i][j]*bestx[j] for i in range(n) for j in range(n)),'cube witness')
    # Independent full quadratic evaluation: no Gray-code updates or fields.
    diagonal=sum(G[i][i] for i in range(n))
    edges=[(i,j,2*G[i][j]) for i in range(n) for j in range(i+1,n) if G[i][j]]
    direct=None
    for tail in itertools.product((-1,1),repeat=n-1):
        y=(1,)+tail
        value=diagonal+sum(w*y[i]*y[j] for i,j,w in edges)
        direct=value if direct is None else max(direct,value)
    need(best==direct,'independent exhaustive cube maximum')
    return best,bestx

def audit_case(k,depth,enumerate_cube=False):
    row=certificate(k,depth,keep_ledger=True)
    cols=sorted({1}|{int(e) for r in row['terminal_weight_ledger'] for e in r})
    G=gram(row['terminal_weight_ledger'],cols);upper=sum(abs(v) for r in G for v in r)
    b=row['b'];X=2**k;R=sum(envelope(a,b) for a in range(X,2*X,b))
    Z=F(upper,4*b)+R;A=2*X*X-4*Z-X
    sig=balance(G)
    if sig is not None:
        need(sig[0]==1,'marked root orientation')
        need(sum(sig[i]*G[i][j]*sig[j] for i in range(len(G)) for j in range(len(G)))==upper,'balanced exact maximum')
    triangles=[]
    for i,j,h in itertools.combinations(range(len(G)),3):
        if G[i][j]*G[i][h]*G[j][h]<0:
            triangles.append((min(abs(G[i][j]),abs(G[i][h]),abs(G[j][h])),i,j,h))
    triangles.sort(reverse=True);used=set();paid=0;chosen=[]
    for w,i,j,h in triangles:
        edges={(i,j),(i,h),(j,h)}
        if not used&edges:used|=edges;paid+=w;chosen.append([i,j,h,w])
    refined=upper-4*paid
    exact=None;witness=None
    if enumerate_cube:
        exact,witness=cube_max(G);need(exact<=refined<=upper,'frustration upper bound')
    # Check the actual shared terminal push, including absorption at ROOT.
    nextrow=certificate(k,depth+1,keep_ledger=True)
    transported=[]
    for r in row['terminal_weight_ledger']:
        new={}
        for es,v in r.items():
            e=int(es);target=1 if e==1 else T(e);sign=1 if e==1 else -1
            key=str(target);new[key]=new.get(key,0)+sign*v
        transported.append({e:v for e,v in new.items() if v})
    need(transported==nextrow['terminal_weight_ledger'],'sourcewise common terminal extension')
    nextcols=sorted({1}|{int(e) for r in transported for e in r})
    nextG=gram(transported,nextcols)
    need(sum(abs(v) for r in nextG for v in r)<=upper,'entrywise Gram norm monotonicity')
    result={'k':k,'depth':depth,'b':b,'sources_unresolved':row['sources_not_grounded_by_this_controller'],
            'terminal_columns':cols,'gram_integer_numerator':G,'gram_denominator':4*b,
            'rowwise_S_upper_over_X2':row['current_energy_upper_over_X2'],
            'rowwise_paid':row['certifies_E_le_2X2'],
            'joint_S_upper_over_X2':float(Z/X**2),'joint_S_upper_exact':[Z.numerator,Z.denominator],
            'joint_paid':A>=0 and A*A>=16*X*Z,'balanced_sign_witness':sig,
            'edge_disjoint_negative_triangles':chosen,'frustration_bound_numerator':refined,
            'cube_exact_numerator':exact,'cube_maximizer':witness}
    print({key:result[key] for key in ('k','depth','sources_unresolved','rowwise_paid','joint_paid','joint_S_upper_over_X2')},flush=True)
    return result

def main():
    rows=[audit_case(k,d,e) for k,d,e in ((10,89,False),(12,110,False),(14,134,False),(12,112,False),(12,128,True))]
    result={'status':'FINITE-EXACT joint terminal certificates; general monotonicity proved in note','cases':rows}
    target=Path(__file__).resolve().parents[2]/'05-knowledge/results/collatz_terminal_gram_20261015a.json'
    target.write_text(json.dumps(result,indent=2)+'\n')

if __name__=='__main__':main()
