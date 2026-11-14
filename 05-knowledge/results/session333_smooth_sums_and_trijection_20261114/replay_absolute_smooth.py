"""Exact finite identity seeds and low cases of the full smooth rail."""
from itertools import permutations
from pathlib import Path
import sys

sys.path.insert(0,str(Path(__file__).resolve().parent.parent /
                      "session328_fermat_gluing_and_phase_repair_20261114"))
from fermat_profiles import ends,strong

def verify(x,msg):
    if not x:raise RuntimeError(msg)

def graph(w):
    n=len(w); x=n;y=n+1
    A=[((1<<n)-1)^((1<<(i+1))-1) for i in range(n)]+[1<<y,0]
    for i,d in enumerate(map(int,w)):
        for v,b in [(x,(d>>1)&1),(y,d&1)]:
            if b:A[v]|=1<<i
            else:A[i]|=1<<v
    return tuple(A)

def visits(A,v):
    C=[0]*len(A)
    W=[i for i in range(len(A)) if i!=v]
    for P in permutations(W):
        for mask in range(1,1<<(1+len(W))):
            paths=[]
            for j in range(len(W)+1):
                if (mask>>j)&1:paths.append(v)
                if j<len(W):paths.append(P[j])
            if all(A[i]>>j&1 for i,j in zip(paths,paths[1:])):
                C[mask.bit_count()]+=1
    return C

def main():
    vals=[]
    for k in range(1,4):
        line=[]
        for l in range(1,4):
            T=graph("11"+"0"*k+"2"*l+"311")
            n=sum(ends(T)[-1])
            verify(n==1+36*2**k*3**l and strong(T),"literal bivariate seed")
            line.append(n)
        vals.append(line)
    for A,v,h,d in [((2,12,9,17,7),0,13,6),((2,4,41,19,7,27),4,19,9)]:
        V=visits(A,v)
        verify(V[1:3]==[h,d] and not any(V[3:]) and strong(A),"low ternary branch")
    print("9 direct Hamilton counts:",vals)
    print("literal two-visit coefficients for B=1,2 passed")

if __name__=="__main__":main()
