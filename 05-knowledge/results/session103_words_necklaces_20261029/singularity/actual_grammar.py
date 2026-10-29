"""Actual abcde insertion semantics, pinned source d3132eafbfa9c43c6d9f1110566b1e9c3a3286e1.
Keeps chronological vertex names and an explicit unique-cycle mark; no rotation quotient.
"""
from native_kernel import need

def build(word):
    A=[[0]*4 for _ in range(4)]
    for u,v in [(0,1),(1,2),(2,3),(3,0),(0,2),(1,3)]:A[u][v]=1
    cycle=list(range(4));history=[]
    for letter in word:
        n=len(A);top=[u for u in range(n) if sum(A[u])==n-2]
        need(len(top) in (1,2),'not grammar state')
        pred=lambda u:cycle[(cycle.index(u)-1)%n]
        if len(top)==1:
            a=top[0];ports={'a':pred(a),'b':a}
        else:
            a,b=top
            if not A[a][b]:a,b=b,a
            need(pred(b)==a,'maxima not adjacent')
            ports={'c':pred(a),'d':a,'e':b}
        need(letter in ports,'illegal letter')
        p=ports[letter];history.append({'letter':letter,'new_vertex':n,'parent':p,'cycle_before':cycle[:]})
        A=[row+[int(u==p)] for u,row in enumerate(A)]+[[int(v!=p) for v in range(n)]+[0]]
        cycle.insert(cycle.index(p)+1,n)
    return A,cycle,history
