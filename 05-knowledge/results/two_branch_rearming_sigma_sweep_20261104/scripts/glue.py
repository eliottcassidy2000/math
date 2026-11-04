"""For each digit-2 class (Y = Phi(T,a) with two HCs C1 = P S1 S2 R, C2 = P S2 S1 R): count the Hamiltonian cycles of
Y - S1, Y - S2, Y - (S1 u S2), and of the sub-tournaments induced on P u S1 u S2 ("left", closing with the arc c_r -> c_0?)"""
import sys
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from douglas_chords import d6_to_rows, ham_paths_ending, ham_cycles, phi
from collections import Counter
path=sys.argv[1]; order=int(sys.argv[2])
def cycles_sub(Y, verts):
    verts=list(verts); idx={v:k for k,v in enumerate(verts)}
    sub=[0]*len(verts)
    for v in verts:
        for u in verts:
            if Y[v]>>u&1: sub[idx[v]]|=1<<idx[u]
    return len(ham_cycles(sub))
stats=Counter(); n_pairs=0
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
        V=set(range(n))
        k1=cycles_sub(Y, V-set(S1)); k2=cycles_sub(Y, V-set(S2)); k12=cycles_sub(Y, V-set(S1)-set(S2))
        stats[(k1,k2,k12)]+=1; n_pairs+=1
print(f"order {order}: {n_pairs} pairs; (HCs of Y-S1, Y-S2, Y-S1-S2) distribution:", sorted(stats.items()))
