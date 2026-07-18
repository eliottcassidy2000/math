# float scan of R_sharp over the worst shape, wide b-range
P=[1,2,4,7,9,11,12]
def maxgap(P,killers):
    ivs=[]
    for v in list(P)+list(killers):
        w=1.0/(14*v); inv=1.0/v
        for j in range(v):
            c=j*inv; lo=(c-w)%1.0; hi=(c+w)%1.0
            if lo<hi: ivs.append((lo,hi))
            else: ivs.append((lo,1.0)); ivs.append((0.0,hi))
    ivs.sort(); mlo,mhi=ivs[0]; merged=[]
    for lo,hi in ivs[1:]:
        if lo<=mhi+1e-15: mhi=hi if hi>mhi else mhi
        else: merged.append((mlo,mhi)); mlo,mhi=lo,hi
    merged.append((mlo,mhi))
    mg=0.0
    for i in range(len(merged)):
        a=merged[i][1]; b=merged[(i+1)%len(merged)][0]+(1.0 if i+1==len(merged) else 0.0)
        if b-a>mg: mg=b-a
    return mg
worst=0.0; wb=None
for b in range(157,4001):
    K=[b,b+2,b+4,b+6,b+8]
    L=maxgap(P,K)
    R=1.0/(7*L*(b+8))
    if R>worst: worst=R; wb=b
print("worst core [1,2,4,7,9,11,12], consecutive step-2 killers, b in [157,4000]:")
print("  max R_sharp = %.6f at b=%d ; ALL < 1: %s"%(worst,wb,worst<1))
# also step 1,3,4 and a few other cores for cross-check
for step in (1,3,4):
    w2=0.0; wb2=None
    for b in range(157,2001):
        K=[b+step*i for i in range(5)]; L=maxgap(P,K); R=1.0/(7*L*(K[-1]))
        if R>w2: w2=R; wb2=b
    print("  step=%d: max R_sharp=%.5f at b=%d (<1: %s)"%(step,w2,wb2,w2<1))
