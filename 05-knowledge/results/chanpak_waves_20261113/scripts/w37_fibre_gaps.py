from collections import defaultdict, Counter
import random
F=defaultdict(set)
for line in open('pairs10.txt'):
    h,p=map(int,line.split()); F[h].add(p)
g=Counter(); mods=Counter()
for h,ps in F.items():
    s=sorted(ps)
    for a,b in zip(s,s[1:]): g[b-a]+=1
    for p in s: mods[(h%4,p%8)]+=1
print('smallest P2 gaps inside H-fibres (gap: count):',sorted(g.items())[:8])
print('P2 mod 8 by H mod 4:',sorted(mods.items())[:12])
# two-fibre repair sumsets: take the 6 largest fibres, offsets D = {p - min}; test sumset coverage of multiples of 4
big=sorted(F.items(),key=lambda t:-len(t[1]))[:6]
for (h1,p1),(h2,p2) in [(big[0],big[1]),(big[2],big[3]),(big[4],big[5])]:
    D1=sorted(set(p-min(p1) for p in p1)); D2=sorted(set(p-min(p2) for p in p2))
    step=4
    S=set(a+b for a in D1 for b in D2)
    top=min(D1[-1],D2[-1])
    lo,hi=top//4,top     # central window of the sumset
    want=[x for x in range(lo - lo%step, hi, step)]
    hit=sum(1 for x in want if x in S)
    print(f'fibres H={h1} ({len(p1)} values) and H={h2} ({len(p2)}): one-fibre fill (step 4) {len(D1)/(D1[-1]//4+1):.3f}; sumset fill of [{lo},{hi}) in steps of 4: {hit/len(want):.4f}')
