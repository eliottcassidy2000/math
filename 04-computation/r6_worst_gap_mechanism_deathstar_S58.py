from fractions import Fraction as F
P=[1,2,4,7,9,11,12]; K=[171,173,175,177,179]
# build all danger arcs with their owner speed; find complement gaps and their bounding owners
arcs=[]  # (lo,hi,owner)
for v in P+K:
    w=F(1,14*v)
    for j in range(v):
        c=F(j,v); lo=c-w; hi=c+w
        lo%=1; hi%=1
        if lo<hi: arcs.append((lo,hi,v))
        else: arcs.append((lo,F(1),v)); arcs.append((F(0),hi,v))
# sweep: sort by lo, merge to covered intervals but track which owner ends at a gap's left and starts at right
arcs.sort()
# compute union intervals with boundary owners
# event points
covered=[]
cur_lo,cur_hi,left_owner,right_owner_of_hi=None,None,None,None
merged=[]  # (lo,hi, owner_at_lo, owner_at_hi)
for (lo,hi,v) in arcs:
    if cur_lo is None:
        cur_lo,cur_hi,lo_owner,hi_owner=lo,hi,v,v
    elif lo<=cur_hi:
        if hi>cur_hi: cur_hi=hi; hi_owner=v
    else:
        merged.append((cur_lo,cur_hi,lo_owner,hi_owner))
        cur_lo,cur_hi,lo_owner,hi_owner=lo,hi,v,v
merged.append((cur_lo,cur_hi,lo_owner,hi_owner))
# gaps between merged covered intervals (circular). gap i: from merged[i].hi to merged[i+1].lo
best=None
for i in range(len(merged)):
    a=merged[i]; b=merged[(i+1)%len(merged)]
    gap_lo=a[1]; gap_hi=b[0]+(1 if i+1==len(merged) else 0)
    g=gap_hi-gap_lo
    if best is None or g>best[0]:
        best=(g, gap_lo, gap_hi, a[3], b[2])  # width, left end (arc of owner a.hi_owner), right end (owner b.lo_owner)
g,glo,ghi,leftowner,rightowner=best
print("largest gap L =",g,"=",float(g))
print("  spans (",float(glo),",",float(ghi),")")
print("  LEFT end = right edge of a danger arc owned by speed", leftowner)
print("  RIGHT end = left edge of a danger arc owned by speed", rightowner)
print("  threshold 1/(7*179) =", F(1,7*179), "=", float(F(1,7*179)))
print("  L / threshold =", float(g/F(1,7*179)), " => R_sharp =", float(F(1,7*179)/g))
# which core-safe arc does it sit in? compute S(P) arcs
def danger(v):
    w=F(1,14*v); out=[]
    for j in range(v):
        c=F(j,v); lo=(c-w)%1; hi=(c+w)%1
        if lo<hi: out.append((lo,hi))
        else: out.append((lo,F(1))); out.append((F(0),hi))
    return out
def sub(safe,arcs):
    for clo,chi in sorted(arcs):
        new=[]
        for lo,hi in safe:
            if chi<=lo or clo>=hi: new.append((lo,hi)); continue
            if clo>lo: new.append((lo,clo))
            if chi<hi: new.append((chi,hi))
        safe=new
    return safe
S=[(F(0),F(1))]
for v in P: S=sub(S,danger(v))
S.sort(key=lambda x:x[1]-x[0])
print("core-safe arcs: count",len(S),"; widest =",float(S[-1][1]-S[-1][0]),"at",(float(S[-1][0]),float(S[-1][1])))
# which core arc contains the gap center?
gc=(glo+ghi)/2 % 1
for lo,hi in S:
    if lo<=gc<=hi: print("  gap center",float(gc),"is inside core-safe arc",(float(lo),float(hi)),"width",float(hi-lo))
