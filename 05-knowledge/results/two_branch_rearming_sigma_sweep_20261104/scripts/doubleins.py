"""Two-cycle tournaments as double insertions into unique-cycle tournaments.
For each digit-2 class: Y = Phi(T,a) with cycles C1 = P S1 S2 R, C2 = P S2 S1 R.  Let S be the short segment
(|S2| <= 2 -> S = S2, else |S1| = 1 -> S = S1), D = Y - S (unique-cycle).  Record (|S|, order of D) and, per canonical
form of D (labelg), how many distinct two-cycle Y arise from it.  Also test: is the block of D's cycle that beats S
contiguous, and are the two insertion positions of S into D's cycle exactly the two positions of C1/C2?"""
import sys, subprocess
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from douglas_chords import d6_to_rows, ham_paths_ending, ham_cycles, phi
from collections import Counter, defaultdict
LABELG='/opt/homebrew/bin/labelg'
def to_d6(rows):
    n=len(rows); bits=[]
    for i in range(n):
        for j in range(n): bits.append(rows[i]>>j&1)
    while len(bits)%6: bits.append(0)
    s='&'+chr(n+63)
    for k in range(0,len(bits),6):
        v=0
        for b in bits[k:k+6]: v=v*2+b
        s+=chr(v+63)
    return s
path=sys.argv[1]; order=int(sys.argv[2])
sizes=Counter(); perD=defaultdict(set); Ys=0; Dlist=[]; recs=[]
for line in open(path):
    parts=line.split()
    if len(parts)==4:
        if int(parts[1])!=1: continue
        H,d6=int(parts[2]),parts[3]
    else:
        H,d6=int(parts[0]),parts[1]
    rows=d6_to_rows(d6); n=len(rows)
    if n!=order: continue
    for a in range(n):
        outs=[j for j in range(n) if rows[a]>>j&1]
        if len(outs)!=1: continue
        r0=outs[0]
        if ham_paths_ending(rows,[v for v in range(n) if v!=a])[r0]!=2: continue
        Y=phi(rows,a)
        if a!=0:
            perm=[a]+[v for v in range(n) if v!=a]; inv={v:k for k,v in enumerate(perm)}
            Y2=[0]*n
            for i in range(n):
                for j in range(n):
                    if Y[i]>>j&1: Y2[inv[i]]|=1<<inv[j]
            Y=Y2
        c1,c2=ham_cycles(Y); pos={v:i for i,v in enumerate(c1)}
        perm2=[pos[v] for v in c2]
        p=next(i for i in range(n) if perm2[i]!=i)-1; q=perm2[p+1]-1
        L2=1
        while p+1+L2<n and perm2[p+1+L2]==q+1+L2: L2+=1
        r=q+L2
        S1=[c1[i] for i in range(p+1,q+1)]; S2=[c1[i] for i in range(q+1,r+1)]
        S = S2 if len(S2)<=2 else S1
        if len(S)>2: print("unexpected", H); continue
        keep=[v for v in range(n) if v not in S]; idx={v:k for k,v in enumerate(keep)}
        D=[0]*len(keep)
        for v in keep:
            for u in keep:
                if Y[v]>>u&1: D[idx[v]]|=1<<idx[u]
        sizes[(len(S), len(S1), len(S2))]+=1
        recs.append((to_d6(D), to_d6(Y)))
        Ys+=1
can=subprocess.run([LABELG,'-q'],input='\n'.join(r[0] for r in recs)+'\n',capture_output=True,text=True).stdout.split()
canY=subprocess.run([LABELG,'-q'],input='\n'.join(r[1] for r in recs)+'\n',capture_output=True,text=True).stdout.split()
for (d,y),cd,cy in zip(recs,can,canY): perD[cd].add(cy)
mult=Counter(len(v) for v in perD.values())
print(f"order {order}: {Ys} two-cycle classes; short segment sizes (|S|,|S1|,|S2|): {sorted(sizes.items())}")
print(f"   distinct unique-cycle D's: {len(perD)} (orders {sorted(set(len(d6_to_rows(k)) for k in perD))}); number of two-cycle Y per D: {sorted(mult.items())}")
