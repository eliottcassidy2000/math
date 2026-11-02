#!/usr/bin/env python3
"""Independent content + weighted/context audit. No research imports.
Binary profiles via disjoint odd-cycle collections, weighted profiles via
subset endpoint/partition DP, literal contextual block similarity checks.
"""
from fractions import Fraction as F
from itertools import combinations,permutations
from math import factorial,gcd,comb
import json

def ck(ok,msg):
    if not ok:raise RuntimeError(msg)

def cyc(n,steps=None):
    if steps is None:steps=range(1,(n+1)//2)
    return [[int((j-i)%n in steps) for j in range(n)] for i in range(n)]

def block(r,s):
    n=r+s+1;a=[[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i<r and j<r:a[i][j]=int(i<j)
            elif r<=i<r+s and r<=j<r+s:a[i][j]=int(i<j)
            elif i<r and r<=j<r+s:a[i][j]=1
            elif r<=i<r+s and j==n-1:a[i][j]=1
            elif i==n-1 and j<r:a[i][j]=1
    return a

def join(a,b):
    n=len(a);m=len(b)
    return [row+[1]*m for row in a]+[[0]*n+row for row in b]

def oddcycle_profile(a):
    n=len(a);N=1<<n;cycles=[0]*N
    # Every subset's directed cycles have their least vertex fixed first.
    for root in range(n):
        dp=[[0]*n for _ in range(N)];dp[1<<root][root]=1
        for S in range(1<<root,N):
            if not S>>root&1 or S&((1<<root)-1):continue
            if S.bit_count()>=3 and S.bit_count()%2:
                cycles[S]=sum(dp[S][v]*a[v][root] for v in range(n))
            for v in range(n):
                if not dp[S][v]:continue
                for w in range(root+1,n):
                    if not S>>w&1 and a[v][w]:dp[S|1<<w][w]+=dp[S][v]
    polys=[[0]*(n+1) for _ in range(N)];polys[0][0]=1
    for S in range(1,N):
        root=S&-S
        for k,z in enumerate(polys[S^root][:-1]):polys[S][k+1]+=z
        C=S
        while C:
            if C&root and cycles[C]:
                for k,z in enumerate(polys[S^C][:-1]):polys[S][k+1]+=2*cycles[C]*z
            C=(C-1)&S
    # x^k = sum_c S(k,c)(x)_c.
    stir=[0]*(n+1);stir[0]=1;prof=[0]*(n+1)
    for k in range(1,n+1):
        stir=[0]+[c*stir[c]+stir[c-1] for c in range(1,n+1)]
        for c in range(1,n+1):prof[c]+=polys[-1][k]*stir[c]
    return prof[1:]

def content(p):return gcd(*(factorial(c)*factorial(c-1)*z for c,z in enumerate(p,1)))

def weighted_profile(a):
    n=len(a);N=1<<n;dp=[[F(0)]*n for _ in range(N)];h=[F(0)]*N;h[0]=F(1)
    for v in range(n):dp[1<<v][v]=F(1)
    for S in range(1,N):
        for v in range(n):
            if not S>>v&1:continue
            R=S^(1<<v)
            if R:dp[S][v]=sum(dp[R][w]*a[w][v] for w in range(n) if R>>w&1)
        h[S]=sum(dp[S])
    ps=[[F(0)]*(n+1) for _ in range(N)];ps[0][0]=F(1)
    for S in range(1,N):
        root=S&-S;C=S
        while C:
            if C&root:
                for k in range(1,n+1):ps[S][k]+=h[C]*ps[S^C][k-1]
            C=(C-1)&S
    return ps[-1][1:]

def ham(a):
    n=len(a);dp=[[F(0)]*n for _ in range(1<<n)]
    for v in range(n):dp[1<<v][v]=F(1)
    for S in range(1,1<<n):
        for v in range(n):
            if not S>>v&1:continue
            R=S^(1<<v)
            if R:dp[S][v]=sum(dp[R][w]*a[w][v] for w in range(n) if R>>w&1)
    return sum(dp[-1])

def substitute(q,port,a):
    n=len(a);ext=[j for j in range(len(q)) if j!=port];N=n+len(ext)
    out=[[F(0)]*N for _ in range(N)]
    for i in range(n):
        for j in range(n):out[i][j]=a[i][j]
        for k,j in enumerate(ext):out[i][n+k]=q[port][j];out[n+k][i]=q[j][port]
    for k,i in enumerate(ext):
        for l,j in enumerate(ext):out[n+k][n+l]=q[i][j]
    return out

def mm(a,b):return [[sum(x*y for x,y in zip(row,col)) for col in zip(*b)] for row in a]
def tr(a):return [list(x) for x in zip(*a)]

def literal_gate(poly,op):
    if op=='P':return [0]+poly
    s=1 if op=='C' else 2;out=[0]*(len(poly)+s+1)
    # On x^k, enumerate odd cycles in C3[L_k,L_s,K1].
    for k,z in enumerate(poly):
        if not z:continue
        out[k+s+1]+=z
        for u in range(1,k+1):
            for v in range(1,s+1):
                if (u+v)%2==0:out[k+s+1-u-v]+=2*z*comb(k,u)*comb(s,v)
    return out

def val(poly,x):return sum(c*x**k for k,c in enumerate(poly))
def mod(x,p):
    x=F(x);return x.numerator*pow(x.denominator,-1,p)%p

def main():
    w5,w11=block(2,1),block(2,2)
    graphs={'W5':w5,'W11':w11,'Z5':join(w5,w11),'R5':cyc(5),'R11':cyc(11)}
    wants={'W5':1,'W11':1,'Z5':5,'R5':15,'R11':33};records={}
    for name,a in graphs.items():
        p=oddcycle_profile(a);g=content(p);ck(g==wants[name],('content',name,g))
        ck(gcd(p[0],factorial(len(a))*factorial(len(a)-1))%g==0,('factorial restriction',name))
        records[name]={'profile':p,'content':g}
    ck(records['R11']['profile']==[93027,787479,2034813,2376110,1481546,537075,117975,15840,1265,55,1],'R11 literal profile')
    a=graphs['R5'];aut=sum(all(a[i][j]==a[g[i]][g[j]] for i in range(5) for j in range(5)) for g in permutations(range(5)))
    ck(aut==5,'R5 automorphisms')
    S={1,3,4,5,9};p11=cyc(11,S);group={tuple((a*x+b)%11 for x in range(11)) for a in S for b in range(11)}
    ck(len(group)==55 and all(p11[i][j]==p11[g[i]][g[j]] for g in group for i in range(11) for j in range(11)),'affine action')
    ck(all(tuple(g[h[x]] for x in range(11)) in group for g in group for h in group),'affine closure')
    A=[[F(0) if i==j else F(1,2) for j in range(5)] for i in range(5)]
    A[0][3],A[3][0],A[1][4],A[4][1]=F(5,8),F(3,8),F(3,4),F(1,4)
    w=[1,1,-2,0,0];U=[[F(i==j)-F(w[i]*w[j],3) for j in range(5)] for i in range(5)]
    B=mm(mm(U,A),tr(U));pa,pb=weighted_profile(A),weighted_profile(B);delta=F(1,3456)
    ck(pa==[F(15,2),30,30,10,1] and pb==[F(25921,3456),30,30,10,1],'weighted profiles')
    ck(all(sum(row)==1 for row in U) and all(sum(col)==1 for col in zip(*U)),'two-sided ones')
    # Row-ones preservation alone does not preserve an incoming homogeneous block.
    bad=[[F(1),F(0)],[F(1,2),F(1,2)]];bad_inv=[[F(1),F(0)],[F(-1),F(2)]]
    ck(all(sum(row)==1 for row in bad) and mm([[F(1),F(1)]],bad_inv)!=[[1,1]],'one-sided similarity hostile')
    contexts=[]
    for n in range(1,4):
        edges=list(combinations(range(n),2))
        for bits in range(1<<len(edges)):
            q=[[0]*n for _ in range(n)]
            for k,(i,j) in enumerate(edges):q[i][j]=bits>>k&1;q[j][i]=1-q[i][j]
            for port in range(n):contexts.append((q,port))
    reflect=[]
    for rows in [(16,17,3,7,12),(16,25,3,5,12)]:
        q=[[row>>j&1 for j in range(5)] for row in rows];contexts.append((q,0));reflect.append(q)
    differences=[]
    for q,port in contexts:
        qa,qb=substitute(q,port,A),substitute(q,port,B);N=len(qa)
        D=[[U[i][j] if i<5 and j<5 else F(i==j) for j in range(N)] for i in range(N)]
        ck(mm(mm(D,qa),tr(D))==qb,('context block similarity',len(q),port))
        d=ham(qb)-ham(qa);ck(d==delta*ham(q),('context first-channel',len(q),port));differences.append(d)
    ck(differences[-1]-differences[-2]==2*delta,'mixed reflector rectangle')
    FA=[0,F(3,2),0,5,0,1];FB=[0,F(5185,3456),0,5,0,1]
    wordrecords=[]
    for name,Fpoly,word,wantH,wantn in [('A',FA,'PPDDCDDCDPCDP',100671287595046269,29),('B',FB,'PPCDCCPCDDDD',13291517610922779,27)]:
        poly=[0,1]
        for op in word:poly=literal_gate(poly,op)
        ck(len(poly)-1==wantn and val(poly,1)==wantH,('positive word integer',name))
        ck(all(val(poly,x)%11==mod(val(Fpoly,x),11) for x in range(11)),('positive word function',name))
        wordrecords.append({'input':name,'word':word,'H':wantH,'order':wantn})
    ck(all(mod(val(FA,x),5)==0 and mod(val(FB,x),5)==x for x in range(5)),'mod5 functions')
    print(json.dumps({'status':'PASS','binary_profile_method':'odd-cycle collections then Stirling transform','content':records,'R5_automorphism_order':aut,'Paley11_affine_subgroup':len(group),'weighted_profile_method':'endpoint and partition subset DP','weighted_profile_delta':[str(b-a) for a,b in zip(pa,pb)],'literal_contexts':len(contexts),'mixed_reflector_rectangle':str(differences[-1]-differences[-2]),'mod11_words':wordrecords,'orthogonal_two_sided_ones_checked':True,'one_sided_ones_similarity_hostile':'rejected'},sort_keys=True,indent=2))
if __name__=='__main__':main()
