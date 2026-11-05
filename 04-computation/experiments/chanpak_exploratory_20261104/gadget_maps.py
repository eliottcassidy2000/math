# Near-join gadget maps: T = T1 => G except back arc u -> v (u in G, v = port of T1); new port z in G.
# Claim: (H(T), e_z(T), B_z(T)) = M_{G,u,z} (h1, e_v, B_v). Compute M from 3 probes, verify on more probes.
import itertools
from fractions import Fraction
def stats(n,out):
    full=(1<<n)-1
    end=[[0]*n for _ in range(1<<n)]
    for v in range(n): end[1<<v][v]=1
    for S in range(1,1<<n):
        for v in range(n):
            c=end[S][v]
            if not c: continue
            nb=out[v]&~S
            while nb:
                w=(nb&-nb).bit_length()-1; nb&=nb-1
                end[S|(1<<w)][w]+=c
    hs=[sum(end[S]) if S else 1 for S in range(1<<n)]
    inn=[0]*n
    for a in range(n):
        for b in range(n):
            if out[a]>>b&1: inn[b]|=1<<a
    st=[[0]*n for _ in range(1<<n)]
    for v in range(n): st[1<<v][v]=1
    for S in range(1,1<<n):
        for v in range(n):
            c=st[S][v]
            if not c: continue
            nb=inn[v]&~S
            while nb:
                w=(nb&-nb).bit_length()-1; nb&=nb-1
                st[S|(1<<w)][w]+=c
    h=hs[full]
    e=[end[full][v] for v in range(n)]
    B=[sum(hs[S]*st[full&~S][t] for S in range(1<<n) if not S>>t&1) for t in range(n)]
    return h,e,B
def nearjoin(n1,out1,v,n2,out2,u):
    n=n1+n2; out=[0]*n
    for a in range(n1): out[a]=out1[a]
    for a in range(n2): out[n1+a]=out2[a]<<n1
    for a in range(n1):
        for b in range(n2):
            if a==v and b==u: out[n1+b]|=1<<a
            else: out[a]|=1<<(n1+b)
    return n,out
def tour_from_code(n,code):
    out=[0]*n; k=0
    for i in range(n):
        for j in range(i+1,n):
            if code>>k&1: out[i]|=1<<j
            else: out[j]|=1<<i
            k+=1
    return out
# probes: (n, out, port)
probes=[(1,[0],0),(2,[2,0],0),(3,[2,4,1],0),(3,[6,4,0],1),(4,tour_from_code(4,0b101100),2),(4,tour_from_code(4,0b011010),1),(5,tour_from_code(5,0b1011001101),3)]
ptrip=[]
for (n,out,v) in probes:
    h,e,B=stats(n,out); ptrip.append((h,e[v],B[v]))
print('probe triples',ptrip)
def solve3(A,b):
    # solve A x = b, A 3x3 Fractions
    M=[list(map(Fraction,A[i]))+[Fraction(b[i])] for i in range(3)]
    for c in range(3):
        p=next(r for r in range(c,3) if M[r][c]!=0); M[c],M[p]=M[p],M[c]
        for r in range(3):
            if r!=c and M[r][c]!=0:
                f=M[r][c]/M[c][c]; M[r]=[M[r][k]-f*M[c][k] for k in range(4)]
    return [M[i][3]/M[i][i] for i in range(3)]
base=[0,1,2]  # probes with independent triples
A=[list(ptrip[i]) for i in base]
import sys
nbad=0; ngad=0; mats={}
for n2 in range(1,5):
    for code in range(1<<(n2*(n2-1)//2)):
        G=tour_from_code(n2,code)
        for u in range(n2):
            for z in range(n2):
                outs=[]
                for (n,out,v) in probes:
                    N,O=nearjoin(n,out,v,n2,G,u)
                    h,e,B=stats(N,O); outs.append((h,e[n+z],B[n+z]))
                # rows of M: M[r] . triple = outs[.][r]
                Mrows=[solve3(A,[outs[i][r] for i in base]) for r in range(3)]
                ok=all(sum(Mrows[r][k]*ptrip[i][k] for k in range(3))==outs[i][r] for i in range(len(probes)) for r in range(3))
                integral=all(x.denominator==1 for row in Mrows for x in row)
                ngad+=1
                if not ok or not integral: nbad+=1
                mats[(n2,code,u,z)]=[[int(x) for x in row] for row in Mrows]
print('gadgets (order<=4, all marks):',ngad,' linear-and-integral failures:',nbad)
for key in [(1,0,0,0),(2,0,0,0),(2,0,0,1),(2,0,1,0),(2,0,1,1)]:
    print(key,mats[key])
# compare with C formulas
cm={}
for line in open('C:/tmp/s78/gm4.txt'):
    a=list(map(int,line.split())); cm[(a[0],a[1],a[2],a[3])]=[a[4:7],a[7:10],a[10:13]]
mism=sum(1 for k in mats if mats[k]!=cm[k])
print('C formula vs probe matrices: compared',len(mats),'mismatches',mism)
