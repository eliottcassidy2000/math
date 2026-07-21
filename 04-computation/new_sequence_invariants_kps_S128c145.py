#!/usr/bin/env python3
"""
MORE sequence-like invariants — hunt for NEW integer sequences from our objects.
kind-pasteur-2026-07-21-S128c145.  Owner: look for more sequence-like invariants.

Differentiated from the census/reciprocal-sum sweep (opus THM-1985/1990): here the sequences are
(a) INVARIANT-VALUE-COUNTS |X|(n) = # distinct values invariant X takes over n-tournaments
    (a sequence ABOUT an invariant -- likely uncataloged), and
(b) EXTREMAL sequences maxX(n), and
(c) STRUCTURAL COUNTS (#strong, #self-complementary, #regular), and
(d) METAGRAPH sequences (edges, 1-WL colors).
Exhaustive over iso classes n=3..6 (exact integer arithmetic). First terms -> OEIS lookup.
"""
import itertools, math
from fractions import Fraction
from collections import defaultdict

def setup(n):
    pairs=list(itertools.combinations(range(n),2)); pos={pr:k for k,pr in enumerate(pairs)}
    maps=[[ (pos[(p[i],p[j])],0) if p[i]<p[j] else (pos[(p[j],p[i])],1) for (i,j) in pairs]
          for p in itertools.permutations(range(n))]
    return pairs,maps
def canon_orbit(bits,maps,npairs):
    seen=set()
    for m in maps:
        v=0
        for k in range(npairs):
            b=(bits>>k)&1; tb,inv=m[k]
            if inv:b^=1
            v|=b<<tb
        seen.add(v)
    return min(seen),len(seen)
def adj(bits,pairs,n):
    A=[[0]*n for _ in range(n)]
    for k,(i,j) in enumerate(pairs):
        if (bits>>k)&1:A[i][j]=1
        else:A[j][i]=1
    return A
def matmul(A,B,n): return [[sum(A[i][k]*B[k][j] for k in range(n)) for j in range(n)] for i in range(n)]
def charpoly_int(M,n):
    c=[1]; Mk=[[1 if i==j else 0 for j in range(n)] for i in range(n)]
    for k in range(1,n+1):
        Mk=matmul(M,Mk,n); tr=sum(Mk[i][i] for i in range(n)); ck=-tr//k; c.append(ck)
        for i in range(n): Mk[i][i]+=ck
    return tuple(c)
def bareiss(M,n):
    if n==0: return 1
    M=[r[:] for r in M]; sign=1; prev=1
    for k in range(n-1):
        if M[k][k]==0:
            sw=next((r for r in range(k+1,n) if M[r][k]!=0),None)
            if sw is None: return 0
            M[k],M[sw]=M[sw],M[k]; sign=-sign
        for i in range(k+1,n):
            for j in range(k+1,n):
                M[i][j]=(M[i][j]*M[k][k]-M[i][k]*M[k][j])//prev
        prev=M[k][k]
    return sign*M[n-1][n-1]
def score(A,n): return tuple(sorted(sum(A[i]) for i in range(n)))
def Kmat(A,n): return [[A[i][j]-A[j][i] for j in range(n)] for i in range(n)]
def disc(A,n):
    cK=charpoly_int(Kmat(A,n),n); val=sum(cK[k]*((-1)**(n-k)) for k in range(n+1))
    return Fraction(abs(((-1)**n)*val),2**(n-1))
def var_lam2(A,n):
    # tr(S^4) etc; use Sum lambda^4 - (mean)^2 * n ; lambda^2 are -eigs of S^2. Use tr(K^2),tr(K^4).
    K=Kmat(A,n); K2=matmul(K,K,n); K4=matmul(K2,K2,n)
    s2=-sum(K2[i][i] for i in range(n)); s4=sum(K4[i][i] for i in range(n))  # Sum lam^2, Sum lam^4
    mean=Fraction(s2,n); return Fraction(s4,n)-mean*mean   # variance of lambda^2
def arb_inv(A,n):
    indeg=[sum(A[x][c] for x in range(n)) for c in range(n)]
    Lm=[[ (indeg[i] if i==j else 0)-A[j][i] for j in range(n)] for i in range(n)]
    def root(r):
        minor=[[Lm[i][j] for j in range(n) if j!=r] for i in range(n) if i!=r]
        return bareiss(minor,n-1)
    return tuple(sorted(root(r) for r in range(n)))
def total_arb(A,n):
    ai=arb_inv(A,n); return sum(ai)
def Hcount(A,n):
    dp=[[0]*n for _ in range(1<<n)]
    for v in range(n): dp[1<<v][v]=1
    for mask in range(1<<n):
        for v in range(n):
            if not (mask>>v)&1 or dp[mask][v]==0: continue
            for w in range(n):
                if (mask>>w)&1 or not A[v][w]: continue
                dp[mask|(1<<w)][w]+=dp[mask][v]
    return sum(dp[(1<<n)-1][v] for v in range(n))
def Rsigned(A,n):
    dp=[defaultdict(int) for _ in range(1<<n)]
    for v in range(n): dp[1<<v][v]=1
    for mask in range(1<<n):
        if not dp[mask]: continue
        for v,val in list(dp[mask].items()):
            if val==0: continue
            for w in range(n):
                if (mask>>w)&1 or not A[v][w]: continue
                inv=bin(mask>>(w+1)).count("1"); sgn=-1 if inv&1 else 1
                dp[mask|(1<<w)][w]+=val*sgn
    return abs(sum(dp[(1<<n)-1].values()))
def cyc(A,n):
    out=[]
    for k in range(3,n+1):
        cnt=0
        for S in itertools.combinations(range(n),k):
            s0=S[0]
            for p in itertools.permutations(S[1:]):
                seq=(s0,)+p
                if all(A[seq[i]][seq[(i+1)%k]] for i in range(k)): cnt+=1
        out.append(cnt)
    return tuple(out)
def is_strong(A,n):
    # strongly connected: from 0 reach all AND all reach 0
    def reach(start,mat):
        seen={start}; st=[start]
        while st:
            u=st.pop()
            for v in range(n):
                if mat[u][v] and v not in seen: seen.add(v); st.append(v)
        return len(seen)==n
    AT=[[A[j][i] for j in range(n)] for i in range(n)]
    return reach(0,A) and reach(0,AT)
def is_regular(A,n):
    if n%2==0: return False
    return len(set(sum(A[i]) for i in range(n)))==1

SEQS=defaultdict(list)
META_E=[]; META_WL=[]
for n in range(3,7):
    pairs,maps=setup(n); npairs=len(pairs); reps={}
    for bits in range(2**npairs):
        c,osz=canon_orbit(bits,maps,npairs)
        if c not in reps: reps[c]=(bits,osz)
    Hs=[]; Rs=[]; specs=set(); scs=set(); cycs=set(); discs=set(); arbs=set(); auts=set()
    vars_=[]; tarbs=[]; nstrong=0; nSC=0; nreg=0
    canon_of={}
    for c,(bits,osz) in reps.items():
        A=adj(bits,pairs,n)
        H=Hcount(A,n); R=Rsigned(A,n)
        Hs.append(H); Rs.append(R)
        specs.add(charpoly_int(A,n)); scs.add(score(A,n)); cycs.add(cyc(A,n))
        discs.add(disc(A,n)); arbs.add(arb_inv(A,n)); auts.add(math.factorial(n)//osz)
        vars_.append(var_lam2(A,n)); tarbs.append(total_arb(A,n))
        if is_strong(A,n): nstrong+=1
        # self-complementary: canon(T) == canon(T^op)
        AT=[[A[j][i] for j in range(n)] for i in range(n)]
        bitsT=0
        for k,(i,j) in enumerate(pairs):
            if AT[i][j]: bitsT|=1<<k
        cT,_=canon_orbit(bitsT,maps,npairs)
        if cT==c: nSC+=1
        if is_regular(A,n): nreg+=1
    SEQS["iso |V(G_n)| A000568"].append(len(reps))
    SEQS["|score| distinct"].append(len(scs))
    SEQS["|specA| distinct (spectral cls)"].append(len(specs))
    SEQS["|cyc| distinct"].append(len(cycs))
    SEQS["|H| distinct (Redei spectrum)"].append(len(set(Hs)))
    SEQS["|R| distinct"].append(len(set(Rs)))
    SEQS["|disc| distinct"].append(len(discs))
    SEQS["|arb_inv| distinct"].append(len(arbs))
    SEQS["|var(lam2)| distinct"].append(len(set(vars_)))
    SEQS["|Aut| distinct values"].append(len(auts))
    SEQS["max H (Szele)"].append(max(Hs))
    SEQS["max |R|"].append(max(Rs))
    SEQS["max total-arb"].append(max(tarbs))
    SEQS["#strong iso classes"].append(nstrong)
    SEQS["#self-complementary iso"].append(nSC)
    SEQS["#regular iso classes"].append(nreg)
    # metagraph edges + 1-WL colors
    canon_bits={}
    for bits in range(2**npairs):
        c,_=canon_orbit(bits,maps,npairs); canon_bits[bits]=c
    verts=sorted(set(canon_bits.values())); gadj=defaultdict(set)
    for bits in range(2**npairs):
        cu=canon_bits[bits]
        for k in range(npairs):
            cv=canon_bits[bits^(1<<k)]
            if cv!=cu: gadj[cu].add(cv)
    E=sum(len(gadj[v]) for v in verts)//2
    META_E.append(E)
    # 1-WL
    color={v:0 for v in verts}
    for _ in range(len(verts)+2):
        sig={v:(color[v],tuple(sorted(color[w] for w in gadj[v]))) for v in verts}
        nc={}; nxt={}
        for v in verts:
            if sig[v] not in nc: nc[sig[v]]=len(nc)
            nxt[v]=nc[sig[v]]
        if len(set(nxt.values()))==len(set(color.values())): color=nxt; break
        color=nxt
    META_WL.append(len(set(color.values())))
SEQS["metagraph |E(G_n)| (wiggly)"]=META_E
SEQS["metagraph 1-WL color classes"]=META_WL

print("="*84)
print("SEQUENCE-LIKE INVARIANTS (n=3,4,5,6) — first terms for OEIS lookup")
print("="*84)
for name,vals in SEQS.items():
    print(f"  {name:36s}: {vals}")
print()
print("Candidates likely UNCATALOGED (invariant-value-counts / metagraph): flagged for OEIS search.")
print("DONE.")
