"""Exact modular spectral exclusions for bounded strong two-port periods."""
from endpoint_profiles import OBS,save,need,np
from fractions import Fraction
from numba import njit
import json,sys,time
from itertools import product

LETTERS=[3,4,5,6,7,8]
PRIME=101

@njit(cache=True)
def detzero(A,lam,p):
    B=A.copy();n=len(B)
    for i in range(n):B[i,i]=(B[i,i]-lam)%p
    for k in range(n):
        j=k
        while j<n and B[j,k]==0:j+=1
        if j==n:return True
        if j!=k:
            for l in range(k,n):B[k,l],B[j,l]=B[j,l],B[k,l]
        inv=1
        # Fermat inversion with exact modular exponentiation.
        base=B[k,k];power=p-2
        while power:
            if power&1:inv=inv*base%p
            base=base*base%p;power//=2
        for j in range(k+1,n):
            q=B[j,k]*inv%p
            for l in range(k+1,n):B[j,l]=(B[j,l]-q*B[k,l])%p
            B[j,k]=0
    return False

@njit(cache=True)
def scan(M,length,stride,p):
    N=6**length;hits=[];tested=0
    stack=np.zeros((length+1,17,17),dtype=np.int64)
    for j in range(17):stack[0,j,j]=1
    digits=np.zeros(length,dtype=np.int64)
    depth=0;scale=pow(2,length)%p
    eigs=np.array([scale,scale*(2**stride)%p,scale*(3**stride)%p],dtype=np.int64)
    while depth>=0:
        if digits[depth]>=6:
            digits[depth]=0;depth-=1
            if depth>=0:digits[depth]+=1
            continue
        a=digits[depth]
        for i in range(17):
            for j in range(17):
                total=0
                for k in range(17):total+=stack[depth,i,k]*M[a,k,j]
                stack[depth+1,i,j]=total%p
        if depth+1<length:depth+=1;continue
        tested+=1
        survives=True
        for x in eigs:
            if not detzero(stack[length],x,p):survives=False;break
        if survives:
            code=0
            for x in digits:code=code*6+x
            hits.append(code)
        digits[depth]+=1
    return tested,hits

def run(maxstride=3):
    M=np.array([[[int(2*Fraction(s))%PRIME for s in row] for row in OBS['rank']['query_matrices'][a]] for a in LETTERS],dtype=np.int64)
    start=time.time();out=[]
    for stride in range(1,maxstride+1):
        # Avoid floating log in the mathematical cutoff: (3/2)^L <=3^r.
        bound=0
        while 3**(bound+1)<=3**stride*2**(bound+1):bound+=1
        for length in range(1,bound+1):
            tested,hits=scan(M,length,stride,PRIME);words=[]
            for code in hits:
                w=[0]*length
                for i in range(length-1,-1,-1):w[i]=LETTERS[code%6];code//=6
                words.append(w)
            out.append({'stride':stride,'length':length,'tested':tested,'modular_survivors':words})
            print('stride',stride,'length',length,'tested',tested,'survivors',len(words),'seconds',round(time.time()-start,1),flush=True)
    save('endpoint_periods.json.gz',{'scope':'Six strong-preserving one-exit two-port commands; target with nonzero1,2^a,3^a modes; arithmetic progression strides tested explicitly. Growth Hnew>=3H/2 bounds period length. Modular failures exclude exact eigenvalues; surviving words still need physical-source analysis.','prime':PRIME,'records':out})
if __name__=='__main__':run(int(sys.argv[1]) if len(sys.argv)>1 else 3)
