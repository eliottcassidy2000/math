"""Exact Eisenstein-ring and prime-power controls for the eight-step affine clock."""
from pathlib import Path
from itertools import product
from collections import Counter
import argparse,json
HERE=Path(__file__).resolve().parent
CHECKS=0
def need(ok,msg):
    global CHECKS
    CHECKS+=1
    if not ok:raise ValueError(msg)
Z=(0,0);I=(1,0);W=(0,1);W2=(-1,-1)
def add(x,y):return(x[0]+y[0],x[1]+y[1])
def mul(x,y):return(x[0]*y[0]-x[1]*y[1],x[0]*y[1]+x[1]*y[0]-x[1]*y[1])
def affine(x,y):return(mul(x[0],y[0]),add(x[1],mul(x[0],y[1])))
def formal():
    records=[]
    for a,c in product((I,W,W2),repeat=2):
        if(a,c)==(I,I):continue
        cols=[]
        for b,d in((I,Z),(Z,I)):
            x,y=(a,b),(c,d)
            for _ in range(8):x,y=y,affine(x,y)
            need((x[0],y[0])==(a,c),'formal slope return');cols.append((x[1],y[1]))
        M=[[cols[j][i] for j in range(2)] for i in range(2)]
        N=[[add(M[i][j],(-int(i==j),0)) for j in range(2)] for i in range(2)]
        need(all(x%3==0 for row in N for v in row for x in v),'missing factor3')
        Q=[[(u//3,v//3) for u,v in row] for row in N]
        need(all(add(mul(Q[i][0],Q[0][j]),mul(Q[i][1],Q[1][j]))==Z for i,j in product(range(2),repeat=2)),'formal square-zero failure')
        need(any(x in(I,W,W2,(-1,0),(0,-1),(1,1)) for row in Q for x in row),'no cyclotomic unit entry')
        det=add(mul(Q[0][0],Q[1][1]),tuple(-t for t in mul(Q[0][1],Q[1][0])))
        need(det==Z,'formal determinant not zero')
        records.append({'slopes':[a,c],'Q_with_M8_equal_I_plus_3Q':Q})
    return records

def root(p,k):
    w=next(a for a in range(2,p) if (a*a+a+1)%p==0);m=p
    for _ in range(1,k):w=next(w+t*m for t in range(p) if ((w+t*m)**2+w+t*m+1)%(m*p)==0);m*=p
    need((w*w+w+1)%m==0 and (w-1)%p!=0,'root lift')
    return m,w

def step(x,y,m):return y,(x[0]*y[0]%m,(x[1]+x[0]*y[1])%m)
def scan(p,k,exhaustive=True):
    m,w=root(p,k);counts=Counter();samples=0
    vectors=list(product(range(m),repeat=2)) if exhaustive else [(i,(i*i+17*i+3)%m) for i in range(64)]
    for a,c in product((1,w,w*w%m),repeat=2):
        if(a,c)==(1,1):continue
        for b,d in vectors:
            x,y=(a,b),(c,d)
            for _ in range(8):x,y=step(x,y,m)
            need((x[0],y[0])==(a,c),'lifted slopes')
            nb,nd=(x[1]-b)%m,(y[1]-d)%m
            order=1
            while (order*nb)%m or (order*nd)%m:order*=p
            predicted=8*order
            x,y=(a,b),(c,d);first=None
            for t in range(1,predicted+1):
                x,y=step(x,y,m)
                if (x,y)==((a,b),(c,d)):first=t;break
            need(first==predicted,'predicted exact period failed');counts[first]+=1;samples+=1
    if exhaustive:
        expected={8:8*m}
        for r in range(1,k+1):expected[8*p**r]=8*m*(p**r-p**(r-1))
        need(dict(counts)==expected,'valuation-stratum count')
    return {'p':p,'k':k,'modulus':m,'primitive_cube_root':w,'exhaustive_translations':exhaustive,'pairs_checked':samples,'period_histogram':sorted(counts.items())}
def matrix_clock_controls():
    def mm(A,B,m):return [[sum(A[i][t]*B[t][j] for t in range(2))%m for j in range(2)] for i in range(2)]
    def power(A,n,m):
        B=[[1,0],[0,1]]
        while n:
            if n%2:B=mm(B,A,m)
            A=mm(A,A,m);n//=2
        return B
    out=[];M=[[0,1],[1,1]];I2=[[1,0],[0,1]]
    for k in range(1,13):
        m=7**k;order=16*7**(k-1)
        need(power(M,order,m)==I2,'matrix clock upper bound')
        need(power(M,order//2,m)!=I2,'matrix clock factor2')
        if k>1:need(power(M,order//7,m)!=I2,'matrix clock factor7')
        out.append({'depth':k,'Fibonacci_matrix_order':order,'whole_affine_pair_order':16*7**k})
    return out

if __name__=='__main__':
    pa=argparse.ArgumentParser();pa.add_argument('--write',action='store_true');args=pa.parse_args()
    out={'status':'FINITE-EXACT controls; formal ring proof and scope in ZERO-SLOT.md','formal_matrices':formal(),
         'seventh_power_clock_controls':matrix_clock_controls(),'prime_power_controls':[scan(7,1),scan(13,1),scan(19,1),scan(7,2),scan(7,3,False)],'checks':CHECKS}
    out=json.loads(json.dumps(out));f=HERE/'clock_lifts.json'
    if args.write:f.write_text(json.dumps(out,indent=2)+'\n')
    else:need(json.loads(f.read_text())==out,'receipt mismatch')
    print(json.dumps({k:v for k,v in out.items() if k!='formal_matrices'},indent=2))
