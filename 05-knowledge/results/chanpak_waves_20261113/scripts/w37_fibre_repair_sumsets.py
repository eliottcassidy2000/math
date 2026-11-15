from collections import defaultdict
F=defaultdict(set)
for line in open('pairs10.txt'):
    h,p=map(int,line.split()); F[h].add(p)
big=sorted(F.items(),key=lambda t:-len(t[1]))[:9]
def offs(ps): m=min(ps); return sorted(set((p-m)//24 for p in ps))
def fill(S,lo,hi): return sum(1 for x in range(lo,hi) if x in S)/max(hi-lo,1)
for i in range(0,9,3):
    (h1,p1),(h2,p2),(h3,p3)=big[i],big[i+1],big[i+2]
    D1,D2,D3=offs(p1),offs(p2),offs(p3)
    one=len(D1)/(D1[-1]+1)
    S2=set(a+b for a in D1 for b in D2); top2=D1[-1]+D2[-1]
    S3=set(s+c for s in S2 for c in D3); top3=top2+D3[-1]
    print(f'H={h1},{h2},{h3}: one-fibre fill {one:.3f} (range {D1[-1]} units of 24); two-fibre sumset fill of central half {fill(S2,top2//4,3*top2//4):.4f}; three-fibre {fill(S3,top3//4,3*top3//4):.4f}')
# check the mod-24 law globally: is P2 mod 24 a function of H mod 24 at order 10?
from collections import Counter
law=defaultdict(set)
for h,ps in F.items():
    for p in ps: law[h%24].add(p%24)
print('P2 mod 24 values for each H mod 24 (order 10):',{k:sorted(v) for k,v in sorted(law.items())})
