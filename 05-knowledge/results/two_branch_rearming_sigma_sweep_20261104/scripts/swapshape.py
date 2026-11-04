"""For each digit-2 class: Y = Phi(T, a) has exactly two Hamiltonian cycles C1, C2 (both through a = position 0 of C1).
Test: C2 is C1 with two consecutive segments swapped: C1 = P [p+1..q] [q+1..r] S  ->  C2 = P [q+1..r] [p+1..q] S.
Tabulate (len1, len2) = (q-p, r-q) and the position p."""
import sys
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from douglas_chords import d6_to_rows, ham_paths_ending, ham_cycles, phi
from collections import Counter
path=sys.argv[1]; order=int(sys.argv[2])
lens=Counter(); pstat=Counter(); bad=0; seen=0; both=Counter()
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
        cyc=ham_cycles(Y); assert len(cyc)==2
        c1,c2=cyc; pos={v:i for i,v in enumerate(c1)}
        perm2=[pos[v] for v in c2]   # c2 as positions of c1, starts at 0
        # find p: first index where perm2 deviates from identity
        p=next(i for i in range(n) if perm2[i]!=i) - 1   # perm2[p] == p is the last agreeing position
        # after p, c2 runs q+1..r then p+1..q then r+1.. : so perm2[p+1] = q+1
        q=perm2[p+1]-1
        L2=1
        while p+1+L2<n and perm2[p+1+L2]==q+1+L2: L2+=1
        r=q+L2
        expect=list(range(0,p+1))+list(range(q+1,r+1))+list(range(p+1,q+1))+list(range(r+1,n))
        if expect!=perm2:
            bad+=1
            if bad<=3: print("not a segment swap:", H, perm2)
            continue
        lens[(q-p, r-q)]+=1; pstat[p]+=1; both[(p,q-p,r-q)]+=1; seen+=1
print(f"order {order}: {seen} digit-2 pairs are segment swaps, {bad} are not")
print("  (len1,len2) distribution:", sorted(lens.items()))
print("  p distribution:", sorted(pstat.items()))
print("  (p,len1,len2) most common:", both.most_common(8))
