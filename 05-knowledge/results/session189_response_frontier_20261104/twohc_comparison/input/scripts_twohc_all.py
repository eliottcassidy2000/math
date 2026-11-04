"""All tournaments (gentourng) with exactly two Hamiltonian cycles, orders 4..8: count them, and test whether the two
cycles always differ in exactly three arcs (a segment swap).  Also count unique-HC tournaments (control: F(2n-6))."""
import subprocess, sys
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from douglas_chords import ham_cycles
from collections import Counter
for n in range(int(sys.argv[1]),int(sys.argv[2])+1):
    out=subprocess.run(['gentourng','-q',str(n)],capture_output=True,text=True).stdout.split()
    one=two=0; diffs=Counter(); sole_exit_two=0; outdeg1_in_two=Counter()
    for line in out:
        if len(line)!=n*(n-1)//2: continue
        rows=[0]*n; idx=0
        for i in range(n):
            for j in range(i+1,n):
                if line[idx]=='1': rows[i]|=1<<j
                else: rows[j]|=1<<i
                idx+=1
        cyc=ham_cycles(rows)
        if len(cyc)==1: one+=1
        elif len(cyc)==2:
            two+=1
            c1,c2=cyc
            a1={(c1[i],c1[(i+1)%n]) for i in range(n)}; a2={(c2[i],c2[(i+1)%n]) for i in range(n)}
            diffs[len(a1-a2)]+=1
            k=sum(1 for v in range(n) if bin(rows[v]).count('1')==1); outdeg1_in_two[k]+=1
    print(f"order {n}: unique-HC tournaments {one}; exactly-two-HC tournaments {two}; arc differences {sorted(diffs.items())}; #vertices of out-degree 1 among two-HC tournaments {sorted(outdeg1_in_two.items())}")
