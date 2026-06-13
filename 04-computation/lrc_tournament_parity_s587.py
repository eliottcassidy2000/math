#!/usr/bin/env python3
"""Tournaments encode even/odd & add/mult. Verify: (1) Redei -- H=#Ham-paths is ALWAYS
ODD; (2) one arc flip changes H by an EVEN amount; (3) self-converse classes vs chiral
(converse) pairs -- the converse-orbit parity. Connect to the LRC worry-set (=self-
converse). opus-2026-06-03-S587."""
from itertools import permutations, combinations
def tournaments(n):
    E=list(combinations(range(n),2)); m=len(E)
    for mask in range(2**m):
        adj=[[0]*n for _ in range(n)]
        for b,(i,j) in enumerate(E):
            if mask>>b&1: adj[i][j]=1
            else: adj[j][i]=1
        yield adj,E,mask
def ham_paths(adj,n):
    cnt=0
    for p in permutations(range(n)):
        if all(adj[p[i]][p[i+1]] for i in range(n-1)): cnt+=1
    return cnt
def canon(adj,n):
    best=None
    for p in permutations(range(n)):
        k=tuple(adj[p[i]][p[j]] for i in range(n) for j in range(n))
        if best is None or k<best: best=k
    return best
def converse(adj,n): return [[adj[j][i] for j in range(n)] for i in range(n)]
def main():
    for n in [4,5,6]:
        Hs=[]; classes={}; flip_changes=set()
        seen={}
        for adj,E,mask in tournaments(n):
            H=ham_paths(adj,n); Hs.append(H)
            c=canon(adj,n)
            if c not in classes: classes[c]=(H,canon(converse(adj,n),n))
        allodd=all(h%2==1 for h in Hs)
        # flip parity: take one tournament, flip each arc, record H-change parity
        adj0=next(tournaments(n))[0]; H0=ham_paths(adj0,n)
        even_flip=True
        for i,j in combinations(range(n),2):
            a=[r[:] for r in adj0]; a[i][j],a[j][i]=a[j][i],a[i][j]
            if (ham_paths(a,n)-H0)%2!=0: even_flip=False
        # self-converse vs chiral
        sc=sum(1 for c,(H,cv) in classes.items() if c==cv)
        chiral=len(classes)-sc
        sc_H_odd=all(H%2==1 for c,(H,cv) in classes.items() if c==cv)
        print(f"n={n}: #classes={len(classes)} (self-converse={sc}, chiral pairs->{chiral//2}); "
              f"H always odd (Redei): {allodd}; one-flip changes H by EVEN: {even_flip}; "
              f"self-converse H all odd: {sc_H_odd}")
if __name__=='__main__': main()
