import itertools, random, sys
sys.stdout.reconfigure(line_buffering=True)
def maxgap_after(P, killers):
    ivs=[]
    for v in list(P)+list(killers):
        w=1.0/(14*v); inv=1.0/v
        for j in range(v):
            c=j*inv; lo=(c-w)%1.0; hi=(c+w)%1.0
            if lo<hi: ivs.append((lo,hi))
            else: ivs.append((lo,1.0)); ivs.append((0.0,hi))
    ivs.sort()
    mlo,mhi=ivs[0]; merged=[]
    for lo,hi in ivs[1:]:
        if lo<=mhi+1e-15: mhi=hi if hi>mhi else mhi
        else: merged.append((mlo,mhi)); mlo,mhi=lo,hi
    merged.append((mlo,mhi))
    if len(merged)==1 and merged[0][0]<=1e-15 and merged[0][1]>=1-1e-15: return 0.0
    mg=0.0
    for i in range(len(merged)):
        a,b=merged[i]; nlo=merged[(i+1)%len(merged)][0]+(1.0 if i+1==len(merged) else 0.0)
        g=nlo-b
        if g>mg: mg=g
    return mg
def Rsharp(P,rem):
    L=maxgap_after(P,rem)
    return (1.0/(7*L*max(rem)),L) if L>0 else (999.0,0.0)
random.seed(7)
ALL=[list(c) for c in itertools.combinations(range(1,13),7)]
worst=0.0; wc=None
def upd(P,rem):
    global worst,wc
    R,L=Rsharp(P,rem)
    if R>worst and L>0: worst=R; wc=(list(P),sorted(rem),round(L,7))
for P in ALL:
    lo=13*max(P)+1
    for off in range(0,24):
        for step in (1,2,3,4,5):
            upd(P,[lo+off+step*i for i in range(5)])
print("(a) consec small all cores: maxR=%.4f %s"%(worst,wc))
for P in random.sample(ALL,200):
    for base in (300,500,900,1500,3000):
        for step in (1,2,3,5,8):
            upd(P,[base+step*i for i in range(5)])
print("(b) consec moderate scale: maxR=%.4f %s"%(worst,wc))
for _ in range(40000):
    P=random.choice(ALL); base=random.choice([333,500,1000,3000])
    rem=sorted(random.sample(range(base,base+random.choice([5,8,12,20,40])),5))
    upd(P,rem)
print("(c) clustered large: maxR=%.4f %s"%(worst,wc))
for _ in range(40000):
    P=random.choice(ALL); lo=13*max(P)+1
    rem=sorted(random.sample(range(lo,lo+random.choice([40,120,400,1500])),5))
    upd(P,rem)
print("(d) random spread: maxR=%.4f %s"%(worst,wc))
print("FINAL sampled max R_sharp=%.4f -> sampled bank %s; UNIFORM r=6 OPEN"%(
    worst, "below 1" if worst <= 1 else "reaches/exceeds 1"))
