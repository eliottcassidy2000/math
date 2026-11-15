"""Exact dyadic input cylinders; all conclusions inside the indicated table/menu."""
from collections import defaultdict
from pathlib import Path
from sys import argv

def suffixes(k):
    old=[(1,0,0,0,0)]
    yield old
    for d in range(k):
        new=[]
        for D,R,A,B,C in old:
            new.append((2*D,2*R,A+D,B+(d==0),C+R+(R%2)))
            new.append((2*D,2*R+1,A,B,C))
        yield new
        old=new

def rows():
    fn=Path(__file__).with_name('forms_orders3_5.txt')
    out=defaultdict(list)
    for line in fn.read_text().splitlines():
        if line.startswith('#'):
            continue
        words,loc,_=line.split('|')
        a,q,s,e,t,f,n=map(int,words.split())
        out[a,q].append(((a,q,s,e,t,f,n),tuple(map(int,loc.split()))))
    return out

def signature(coeff,X,Y):
    a,q,s,e,t,f,n=coeff
    d,u,l,j,k=X
    D,v,L,J,K=Y
    x=q*d*v+t*d+a*l
    y=q*D*u+s*D+a*L
    i=a*j+f*(d==1)
    I=a*J+e*(D==1)
    c=q*u*v+t*u+s*v+a*(k+K)+n+f*(u%2)+e*(v%2)
    return (a,q*d*D,x,y,i,I),c

def main(K):
    SUF=list(suffixes(K))
    best=0
    for k in range(K+1):
        counts=0
        for (a,q),group in sorted(rows().items()):
            index={}
            for coeff,loc in group:
                for w in range(k+1):
                    for X in SUF[w]:
                        for Y in SUF[k-w]:
                            sig,c=signature(coeff,X,Y)
                            flip=(*sig[:2],sig[3],sig[2],sig[5],sig[4])
                            descriptor=(loc,*X[:2],*Y[:2],int(sig>flip))
                            K0=min(sig,flip)
                            if K0 not in index:
                                index[K0]={}
                            if c not in index[K0]:
                                index[K0][c]=descriptor
            counts+=len(index)
            for sig, val in index.items():
                for c in sorted(val):
                    if c-2 in val:
                        continue
                    j=c
                    while j in val:
                        j+=2
                    L=(j-c)//2
                    if L>best:
                        best=L
                        print('record bits',k,'length',L,'signature',sig,flush=True)
                        for d in range(c,j,2):
                            print(d,val[d],flush=True)
        print('level',k,'distinct signatures',counts,'record',best,flush=True)
if __name__=='__main__':
    main(int(argv[1]) if len(argv)>1 else 9)
