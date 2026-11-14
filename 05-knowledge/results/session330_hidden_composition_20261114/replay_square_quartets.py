from itertools import combinations
from math import isqrt

W=(1,3,78,66,55,89,80,64,57,87,82,62,38,43,21,60,84,85,36,13,68,32,49,15,34,30,51,93,28,72,9,91,53,47,17,83,61,39,42,79,65,16,20,44,37,63,81,19,45,76,5,95,26,74,70,99,97,24,40,41,59,22,27,94,50,71,10,6,58,86,35,14,11,25,56,88,33,67,54,90,31,18,46,75,69,12,4,96,73,8,92,77,23,98,2,7,29,52,48)
def square(x): return isqrt(x)**2==x
def canon(v):return min(w[i:]+w[:i] for w in (v,v[::-1]) for i in range(len(w)))
def algebra(n):
    S=[i*i for i in range(2,isqrt(2*n-1)+1)]
    ans={}
    for A,B,D in combinations(S,3):
        C=A+D-B
        if not B<=C<D or C not in S:continue
        M=min(A-1,n-B+A)
        if B==C:M=min(M,(A-1)//2)
        bad={A,2*A-B,C,A+C-B}
        for t in range(max(1,C-n),M+1):
            if 2*t not in bad:
                v=(t,A-t,t+B-A,C-t)
                c=canon(v)
                assert c not in ans
                ans[c]=(A,B,C,D)
    return ans
def graph(n):
    N={x:{y for y in range(1,n+1) if y!=x and square(y+x)} for x in range(1,n+1)}
    c=set()
    for a,b in combinations(N,2):
        for x,y in combinations(N[a]&N[b],2):c.add(canon((a,x,b,y)))
    return c
def edges(v):
    return {frozenset(e) for e in zip(v,v[1:]+v[:1])}
def hc(v):
    return set(v)==set(range(1,len(v)+1)) and all(square(sum(e)) for e in edges(v))
def req(ok,msg):
    if not ok:raise RuntimeError(msg)

def main():
    for n in range(2,181):
        v=algebra(n)
        req(set(v)==graph(n),f'graph/formula mismatch at {n}')
        rep={w for w,(_,B,C,_) in v.items() if B==C}
        if n<99:req(not rep,'earlier repeated label')
        if n==99:req(rep=={(1,3,97,99)},'wrong first repeated label')
        if n in (46,98,99,100,145,146,180):print(n,len(v),len(rep))
    i=W.index(99)
    V=W[:1]+W[1:i+1][::-1]+W[i+1:]
    req(hc(W) and hc(V),"invalid Hamiltonian switch")
    req(set(W[::2])==set(V[::2]),"vertex parity changed")
    req(edges(W)-edges(V)=={frozenset((1,3)),frozenset((97,99))},"wrong removed edges")
    req(edges(V)-edges(W)=={frozenset((1,99)),frozenset((3,97))},"wrong new edges")
    print("99: minimum repeated-label Hamiltonian 2-switch, same position parity")
    # Two newly activated repeated templates (not a census of every cycle at this N).
    c=245
    roots=[(a,isqrt((a*a+c*c)//2),c) for a in range(1,c)
           if (a*a+c*c)%2==0 and square((a*a+c*c)//2)]
    req(roots==[(35,175,245),(155,205,245)],"wrong c=245 progression templates")
    V=(612,613,12012,12013,30012,30013)
    E={frozenset((x,y)) for x,y in combinations(V,2) if square(x+y)}
    T=edges((612,613,30012,30013))|edges((12012,12013,30012,30013))
    req(E==T and len(E)==7,"wrong induced theta")
    print("30013: two newly activated repeated templates form an induced theta")


if __name__=="__main__":main()
