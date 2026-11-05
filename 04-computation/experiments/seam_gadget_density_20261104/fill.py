import sys
from collections import defaultdict
P=defaultdict(set)
for line in open(sys.argv[1]):
    h,r=map(int,line.split()); P[h].add(r)
def fill(lo,hi,a,b):
    tot=have=0
    for h in range(lo,hi+1,2):
        for r in range(2*((int(a*h)+1)//2),int(b*h)+1,2):
            tot+=1; have+= r in P[h]
    return have/tot, tot
f1=fill(801,999,2.0,2.75); f2=fill(1601,1799,2.0,2.75)
print(f"fill[801,999]x[2,2.75]={f1[0]:.4f} (n={f1[1]})  fill[1601,1799]x[2,2.75]={f2[0]:.4f} (n={f2[1]})  ratio={f2[0]/f1[0] if f1[0] else float('nan'):.4f}")
