#!/usr/bin/env python3
"""Independent subset-path polynomial and literal mixed-quartic audit."""
from pathlib import Path
from itertools import combinations
from math import comb,factorial
import json,hashlib,subprocess,sys
P=Path(__file__).resolve().parent

def need(x,m):
    if not x:raise RuntimeError(m)
def trim(p):
    while p and not p[-1]:p.pop()
    return p
def add(p,q,scale=1):
    r=p+[0]*max(0,len(q)-len(p))
    for i,x in enumerate(q):r[i]+=scale*x
    return trim(r)
def mul(p,q):
    r=[0]*(len(p)+len(q)-1) if p and q else []
    for i,x in enumerate(p):
        for j,y in enumerate(q):r[i+j]+=x*y
    return trim(r)
def fa(F,G,scale=1):
    out={s:p[:] for s,p in F.items()}
    for s,q in G.items():out[s]=add(out.get(s,[]),q,scale)
    return {s:p for s,p in out.items() if p}
def fm(F,G):
    out={}
    for s,p in F.items():
        for t,q in G.items():
            if not s&t:out[s|t]=add(out.get(s|t,[]),mul(p,q))
    return {s:p for s,p in out.items() if p}
def scale(F,z):return {s:[z*x for x in p] for s,p in F.items() if z}
def paths(rows):
    n=len(rows);M=[[{} for _ in range(n)] for _ in range(n)]
    for st in range(n):
        d={(1<<st,st):[1]}
        for mask in range(1,1<<n):
            if not mask>>st&1:continue
            for j in range(n):
                pol=d.get((mask,j))
                if not pol:continue
                M[st][j][mask]=pol
                for k in range(n):
                    if mask>>k&1:continue
                    key=(mask|1<<k,k);d[key]=add(d.get(key,[]),mul(pol,[int(rows[j]>>k&1),1]))
    C=[{}for _ in rows];R=[{}for _ in rows]
    for i in range(n):
        for j in range(n):R[i]=fa(R[i],M[i][j]);C[j]=fa(C[j],M[i][j])
    H={}
    for F in C:H=fa(H,F)
    return M,C,R,H

def kernel(rows,p,q,r,cache=None):
    M,C,R,H=cache if cache is not None else paths(rows)
    X=fa(fa(fm(C[p],M[q][r]),fm(C[q],M[p][r])),fa(fm(M[r][p],R[q]),fm(M[r][q],R[p])),-1)
    Y=fa(fm(H,X),fm(fa(fm(C[p],R[q]),fm(C[q],R[p])),fa(C[r],R[r],-1)))
    full=(1<<len(rows))-1;x=X.get(full,[]);y=Y.get(full,[])
    need(trim([(i+1)*x[i+1] for i in range(len(x)-1)])==y,'independent derivative channel')
    return x,y

def tt(n):return tuple(sum(1<<j for j in range(i+1,n)) for i in range(n))
def op(rows):return tuple(((1<<len(rows))-1)^r^(1<<i) for i,r in enumerate(rows))
def graph(n,mask):
    out=[0]*n
    for z,(i,j) in enumerate(combinations(range(n),2)):
        if mask>>z&1:out[i]|=1<<j
        else:out[j]|=1<<i
    return tuple(out)
def euler(m,r):
    L=m-2;u=r-1;poly=[0]*(L+1)
    for j in range(u+1):poly[j+L-u]+=comb(u,j)
    for j in range(L-u+1):poly[j+u]-=comb(L-u,j)
    S=[[0]*(L+1) for _ in range(L+1)];S[0][0]=1
    for a in range(1,L+1):
        for b in range(1,a+1):S[a][b]=S[a-1][b-1]+b*S[a-1][b]
    out=[2*factorial(j)*sum(poly[a]*S[a][j] for a in range(j,L+1)) for j in range(L+1)]
    return trim(out),trim(poly)
def weights(m):
    out={m-1:1}
    if m%2:
        for r in range(2,(m-1)//2+1):out[r]=(-1)**r*comb(m-2,r-1)
    else:
        for r in range(2,(m-2)//2+1):out[r]=(-1)**r*(comb(m-3,r-1)-comb(m-3,r-2))
    return out

def external(erows,into):
    M,C,S,H=paths(erows)
    def constant(F):return {s:[p[0]] for s,p in F.items() if p and p[0]}
    M=[[constant(F)for F in row] for row in M];C=list(map(constant,C));S=list(map(constant,S));alpha={0:[1]};beta={0:[1]};tau={}
    for i in range(len(erows)):
        if into>>i&1:alpha=fa(alpha,C[i])
        else:beta=fa(beta,S[i])
    for i in range(len(erows)):
        for j in range(len(erows)):
            if not(into>>i&1) and into>>j&1:tau=fa(tau,M[i][j])
    out=[]
    for i in range(len(erows)):
        ell={};rho={}
        for j in range(len(erows)):
            if into>>j&1:ell=fa(ell,M[i][j])
            else:rho=fa(rho,M[j][i])
        out.append((fa(fm(alpha,S[i]),fm(beta,C[i]),-1),fm(fm(alpha,beta),fa(ell,rho,-1))))
    return tau,out

def substitute(A,E,into):
    m=len(A);out=list(A)+[r<<m for r in E]
    for i in range(m):
        for j in range(len(E)):
            if into>>j&1:out[m+j]|=1<<i
            else:out[i]|=1<<(m+j)
    return tuple(out)
def ports(rows,a,c):
    n=len(rows);out=list(rows)+[1<<(n+1),0]
    for i in range(n):
        if a>>i&1:out[i]|=1<<n
        else:out[n]|=1<<i
        if c>>i&1:out[i]|=1<<(n+1)
        else:out[n+1]|=1<<i
    return out

def main():
    inner=[];checks=0
    qs=[graph(3,z) for z in range(8)]+[graph(4,z) for z in range(64)]+[(6,12,24,17,3),(2,4,9,3)]+[tt(n)for n in range(5,9)]
    for rows in qs:
        cache=paths(rows)
        for r in range(2,len(rows)):
            x,y=kernel(rows,0,1,r,cache);checks+=1
            if rows==tt(len(rows)):need(x==euler(len(rows),r)[0],'Eulerian kernel via Stirling monomials')
            inner.append({'rows':rows,'marks':[0,1,r],'X':x,'derivative_channel':y})
    combos=[]
    for m in range(3,62):
        x=[];ordinary=[]
        for r,w in weights(m).items():
            e,p=euler(m,r);x=add(x,e,w);ordinary=add(ordinary,p,w)
        need(ordinary==([1] if m%2 else [1,2]),'ordinary binomial identity');need(x==([2] if m%2 else [2,4]),'transformed binomial identity')
        combos.append({'m':m,'weights':weights(m),'ordinary_combination':ordinary,'transformed_combination':x})
    # Generic, nontransitive internal children: actual full named quartics.
    models=[];requests=[];spec=[]
    for A0 in [(2,4,1),(2,4,9,3),(6,12,24,17,3),tt(5),op((2,4,9,3))]:
        cache=paths(A0)
        for E,into in [((2,4,1),3),((2,4,9,3),5),((2,4,9,3),9)]:
            tau,channels=external(E,into);T=substitute(A0,E,into);U=substitute(op(A0),E,into);m=len(A0);full=(1<<len(E))-1;x,y=kernel(A0,0,1,2,cache);powers=[{0:[1]}]
            for j in range(max(len(x),len(y))):powers.append(fm(powers[-1],tau))
            for i,(K,G) in enumerate(channels):
                pred=0
                for coeff,F in [(x,K),(y,G)]:
                    for j,z in enumerate(coeff):pred+=z*fm(F,powers[j]).get(full,[0])[0]
                record={'internal_rows':A0,'exterior_rows':E,'into_mask':into,'source_rows':T,'target_rows':U,'marks':[0,1,2,m+i],'X':x,'derivative_channel':y,'predicted_contrast':pred,'query_start':len(requests)}
                for graph0 in [T,U]:
                    for a0 in range(4):
                        for c0 in range(4):
                            a=(a0&1)|((a0>>1)<<1);c=((c0&1)<<2)|((c0>>1)<<(m+i));g=ports(graph0,a,c);requests.append(' '.join(map(str,[len(g),*g])))
                models.append(record)
    subprocess.run(['g++','-O3','-std=c++17',str(P/'count.cpp'),'-o',str(P/'count')],check=True)
    out=subprocess.check_output([str(P/'count')],input='\n'.join(requests)+'\n',text=True);vals=list(map(int,out.split()));need(len(vals)==len(requests),'literal query count')
    for rec in models:
        start=rec['query_start'];values=[vals[start:start+16],vals[start+16:start+32]];coeff=[]
        for vector in values:coeff.append(sum((-1)**(4-a.bit_count()-c.bit_count())*vector[4*a+c] for a in range(4)for c in range(4)))
        need(coeff[0]-coeff[1]==rec['predicted_contrast'],'generic actual mixed quartic identity');rec['actual_values']=values;rec['actual_coefficients']=coeff
    for file,data in [('independent_kernels.json',inner),('independent_binomials.json',combos),('generic_mixed_controls.json',models)]: (P/file).write_text(json.dumps(data,separators=(',',':'))+'\n')
    receipt={'status':'PASS','independent_subset_path_kernels':checks,'binomial_orders':list(range(3,62)),'generic_mixed_coordinates':len(models),'literal_physical_queries':len(requests),'sources':{f:hashlib.sha256((P/f).read_bytes()).hexdigest() for f in ['kernel_audit.py','count.cpp']},'artifacts':{f:hashlib.sha256((P/f).read_bytes()).hexdigest() for f in ['independent_kernels.json','independent_binomials.json','generic_mixed_controls.json']}}
    (P/'kernel_receipt.json').write_text(json.dumps(receipt,sort_keys=True,indent=2)+'\n');print(json.dumps(receipt,sort_keys=True))
if __name__=='__main__':main()
