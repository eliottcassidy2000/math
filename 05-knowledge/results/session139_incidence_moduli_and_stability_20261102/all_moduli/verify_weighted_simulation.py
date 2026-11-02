#!/usr/bin/env python3
from fractions import Fraction as Q
from itertools import permutations,combinations
from math import factorial
from collections import deque
from array import array
import json

def check(test,label):
    if not test:raise RuntimeError(label)

def add(a,b):
    c=[0]*max(len(a),len(b))
    for i,v in enumerate(a):c[i]+=v
    for i,v in enumerate(b):c[i]+=v
    return c

def mul(a,b):
    c=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):c[i+j]+=x*y
    return c

def shift(a,t):
    from math import comb
    c=[0]*len(a)
    for i,v in enumerate(a):
        for j in range(i+1):c[j]+=v*comb(i,j)*t**(i-j)
    return c

def evalpoly(a,x):
    z=0
    for c in reversed(a):z=z*x+c
    return z

def op(a,g):
    if g=='P':return [0]+a
    plus,minus=shift(a,1),shift(a,-1)
    if g=='C':return add([0,0]+a,[0]+add(plus,[-x for x in minus]))
    return add(add(mul([0,-2,0,1],a),mul([0,1,2],plus)),mul([0,1,-2],minus))

def weighted_polynomial(A):
    n=len(A);pc=[Q(0)]*(n+1)
    for perm in permutations(range(n)):
        for cut in range(1<<(n-1)):
            w=Q(1)
            for i in range(n-1):
                if not (cut>>i&1):w*=A[perm[i]][perm[i+1]]
            pc[cut.bit_count()+1]+=w
    pc=[v/factorial(c) for c,v in enumerate(pc)]
    F=[Q(0)]*(n+1);fall=[1]
    for c in range(1,n+1):
        fall=mul(fall,[-(c-1),1]);F=add(F,[pc[c]*v for v in fall])
    return F,pc

def residue(x,p):
    x=Q(x);return x.numerator*pow(x.denominator,-1,p)%p

def find_words(p,targets):
    d=(p-1)//2;N=p**d
    def encode(v,eps):return sum(x*p**i for i,x in enumerate(v))+(N if eps==-1 else 0)
    def decode(code):
        eps=-1 if code>=N else 1;code%=N;v=[]
        for _ in range(d):v.append(code%p);code//=p
        return v,eps
    start=encode(list(range(1,d+1)),-1);par=array('i',[-2])*(2*N);kind=bytearray(2*N);par[start]=-1;queue=array('i',[start]);i=0
    wanted={encode(t,-1) for t in targets}
    while not all(par[t]!=-2 for t in wanted):
        code=queue[i];i+=1;v,eps=decode(code)
        def f(x):
            x%=p
            return 0 if not x else v[x-1] if x<=d else eps*v[p-x-1]
        outs=[([-eps,[(x*f(x))%p for x in range(1,d+1)]]),([eps,[(x*x*f(x)+x*(f(x+1)-f(x-1)))%p for x in range(1,d+1)]]),([-eps,[((x**3-2*x)*f(x)+(2*x*x+x)*f(x+1)+(x-2*x*x)*f(x-1))%p for x in range(1,d+1)]])]
        for j,(e,w) in enumerate(outs):
            z=encode(w,e)
            if par[z]==-2:par[z]=code;kind[z]=j;queue.append(z)
    words=[]
    for target in targets:
        z=encode(target,-1);w=[]
        while par[z]>=0:w.append('PCD'[kind[z]]);z=par[z]
        words.append(''.join(reversed(w)))
    return words

def matmul(A,B):
    return [[sum(A[i][k]*B[k][j] for k in range(len(B))) for j in range(len(B[0]))] for i in range(len(A))]

def transpose(A):return [list(x) for x in zip(*A)]

def substitute(outer,A):
    n=len(A);q=len(outer);N=n+q-1
    B=[[Q(0)]*N for _ in range(N)]
    for i in range(n):
        for j in range(n):B[i][j]=A[i][j]
    for i in range(1,q):
        for v in range(n):B[v][n+i-1]=outer[0][i];B[n+i-1][v]=outer[i][0]
        for j in range(1,q):B[n+i-1][n+j-1]=outer[i][j]
    return B

def hamilton(A):
    n=len(A);dp=[[Q(0)]*n for _ in range(1<<n)]
    for i in range(n):dp[1<<i][i]=Q(1)
    for S in range(1,1<<n):
        for i in range(n):
            if not dp[S][i]:continue
            for j in range(n):
                if not(S>>j&1) and A[i][j]:dp[S|1<<j][j]+=dp[S][i]*A[i][j]
    return sum(dp[-1])

def contextual_checks(A,B):
    w=[1,1,-2,0,0];U=[[Q(int(i==j))-Q(w[i]*w[j],3) for j in range(5)] for i in range(5)]
    contexts=[]
    for q in range(1,5):
        pairs=list(combinations(range(q),2))
        for bits in range(1<<len(pairs)):
            out=[[Q(0)]*q for _ in range(q)]
            for k,(i,j) in enumerate(pairs):out[i][j]=Q(bits>>k&1);out[j][i]=1-out[i][j]
            contexts.append(out)
    for rows in [(16,17,3,7,12),(16,25,3,5,12)]:contexts.append([[Q(rows[i]>>j&1) for j in range(5)] for i in range(5)])
    records=[]
    for out in contexts:
        X,Y=substitute(out,A),substitute(out,B);N=len(X)
        V=[[Q(int(i==j)) for j in range(N)] for i in range(N)]
        for i in range(5):
            for j in range(5):V[i][j]=U[i][j]
        check(matmul(matmul(V,X),transpose(V))==Y,'contextual similarity')
        check(all(sum(row)==1 for row in V),'all-one fixed')
        hx,hy=hamilton(X),hamilton(Y)
        check(hy-hx==hamilton(out)/3456,'contextual first-channel response')
        records.append((hx,hy))
    rectangle=(records[-1][1]-records[-2][1])-(records[-1][0]-records[-2][0])
    check(rectangle==Q(1,1728),'mixed reflector rectangle')
    return {'contexts':len(contexts),'max_substituted_order':9,'mixed_rectangle':str(rectangle)}


def main():
    A=[[Q(0) if i==j else Q(1,2) for j in range(5)] for i in range(5)]
    A[0][3],A[3][0],A[1][4],A[4][1]=Q(5,8),Q(3,8),Q(3,4),Q(1,4)
    B=[[Q(x) for x in row] for row in [['0','1/2','1/2','7/12','5/12'],['1/2','0','1/2','11/24','2/3'],['1/2','1/2','0','7/12','2/3'],['5/12','13/24','5/12','0','1/2'],['7/12','1/3','1/3','1/2','0']]]
    FA,pa=weighted_polynomial(A);FB,pb=weighted_polynomial(B)
    check(FA==[0,Q(3,2),0,5,0,1],'weighted A F')
    check(FB==[0,Q(5185,3456),0,5,0,1],'weighted B F')
    check(pa[1]==Q(15,2) and pb[1]==Q(25921,3456),'weighted H')
    Z5=mul([0,0,4,0,1],[0,2,0,8,0,1])
    check(all(residue(evalpoly(FA,x),5)==evalpoly(Z5,x)%5 for x in range(5)),'A mod5 actual simulation')
    check(all(residue(evalpoly(FB,x),5)==x%5 for x in range(5)),'B mod5 actual simulation')
    targets=[tuple(residue(evalpoly(F,x),11) for x in range(1,6)) for F in [FA,FB]]
    words=find_words(11,targets);records=[]
    for name,F,target,word in zip(['A','B'],[FA,FB],targets,words):
        poly=[0,1]
        for g in word:poly=op(poly,g)
        check(all(evalpoly(poly,x)%11==residue(evalpoly(F,x),11) for x in range(11)),('all function values',name))
        records.append({'weighted_input':name,'word':word,'order':1+sum('PCD'.index(g)+1 for g in word),'actual_H':evalpoly(poly,1),'positive_half_values':list(target)})
    print(json.dumps({'status':'PASS','contextual_checks':contextual_checks(A,B),'weighted_profiles':{'A':[str(x) for x in pa[1:]],'B':[str(x) for x in pb[1:]]},'weighted_polynomials':{'A':[str(x) for x in FA],'B':[str(x) for x in FB]},'mod5_simulation':{'A':'J(W5,W11)','B':'K1'},'mod11_simulation':records},sort_keys=True,indent=2))
if __name__=='__main__':main()
