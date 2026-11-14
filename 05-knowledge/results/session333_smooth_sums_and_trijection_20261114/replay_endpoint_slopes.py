from math import gcd
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parent.parent /
                      "session328_fermat_gluing_and_phase_repair_20261114"))
from fermat_profiles import ends, strong, rooted, hp_all

def chk(x):
    if not x:raise RuntimeError("endpoint and slope replay")

def sub(A,S):
    return tuple(sum((A[v]>>w&1)<<j for j,w in enumerate(S)) for v in S)

def inflate(A,v,k):
    seq=[i for i in range(len(A)) if i!=v]+[v]*k
    return tuple(sum(1<<j for j,b in enumerate(seq)
                    if (a==b and i<j) or (a!=b and A[a]>>b&1))
                 for i,a in enumerate(seq))

def endpaths(A,s,t):
    F=[[0]*len(A) for _ in range(1<<len(A))];F[1<<s][s]=1
    for S,fs in enumerate(F):
        for i,c in enumerate(fs):
            if c:
                C=A[i]&~S
                while C:
                    v=C&-C;C^=v
                    F[S|v][v.bit_length()-1]+=c
    return F[-1][t]

def binary(N):
    r=N.bit_length()
    m=r+2
    y,x=m,m+1
    A=[sum(1<<j for j in range(i+1,m)) for i in range(m)]+[(1<<m)-1,1<<y|1]
    A[m-1]|=1<<x
    for i in range(1,r+1):
        if N>>(i-1)&1:A[i]|=1<<x
        else:A[x]|=1<<i
    return tuple(A)

def main():
    for N in range(1,130):
        A=binary(N)
        chk(strong(A) and endpaths(A,0,len(A)-3)==N)
    print("all fixed-endpoint counts 1–129 verified directly")
    B=[
    ((258,276,9,291,301,327,415,319,4),8,2,135,(115,41,20,176)),
    ((262,28,424,337,293,331,279,379,2),8,1,405,(379,83,56,488)),
    ((760,1013,769,518,44,652,956,540,697,16),9,4,1215,(1095,309,188,1524)),
    ]
    for A,v,z,D,t in B:
        V=[i for i in range(len(A)) if i!=v]
        U=sub(A,V);s=V.index(z);P=hp_all(U)
        h=P[((1<<len(U))-1)^(1<<s)]
        H,L,E=P[-1],list(rooted(U))[s][1],ends(U)[-1][s]
        chk((H,h,E,L)==t and D==L-h and H+2*h-E-L==1)
        for j in range(1,5):
            X=inflate(A,v,j)
            chk(sum(ends(X)[-1])==1+D*2**j and strong(X))
        print("literal offset-one slope",D)
    A=(2,8,3,21,103,79,143,63)
    chk(A[6]>>7&1 and A[-1].bit_count()==6)
    for j in range(1,5):
        X=inflate(A,7,j)
        chk(sum(ends(X)[-1])==1+34*2**j and strong(X))
    print("literal sole-entry coefficient 34")

if __name__=="__main__":
    main()
