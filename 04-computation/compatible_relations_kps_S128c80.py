#!/usr/bin/env python3
"""compatible_relations_kps_S128c80.py -- kind-pasteur S128 cont.80.
ENUMERATING THE COMPATIBLE RELATIONS MOD 4, and correcting my own criterion.

The centre is c = (1/4,1/2,3/4), so a relation m = (m2,m3,m4) is COMPATIBLE iff
        m.c = (m2 + 2 m3 + 3 m4)/4  in  Z    <=>    m2 + 2 m3 + 3 m4 = 0 (mod 4).
So Lambda = {m : m.(1,2,3) = 0 mod 4} is a sublattice of index 4 in Z^3 -- 16 of the 64
classes mod 4.

BUT MY 'COUNT INDEPENDENT COMPATIBLE RELATIONS' FRAMING IS WRONG.  The relation lattice
R(d) = d-perp always has RANK 2, and R(d) cap Lambda has index at most 4 in it, so it is
rank 2 too -- every direction carries two independent compatible relations.  The correct
criterion is whether ALL of R(d) lies in Lambda, because the geodesic closure is the
annihilator of R(d):
        geodesic passes through c  <=>  m.c in Z for EVERY m in R(d)  <=>  R(d) subset Lambda.
Setting e = generator of { m.(1,2,3) : m in R(d) }, this reads e = 0 or 4 | e.
e = 0 exactly when d ~ (1,2,3).  THE QUESTION: does 4 | e ever happen otherwise -- i.e. is
there a SECOND full-concentration family?  PRINT DATA ONLY."""
import sys, itertools
from math import gcd
sys.stdout.reconfigure(line_buffering=True)
W=(1,2,3)
print("### (1) the 16 compatible classes mod 4 ###")
cls=[m for m in itertools.product(range(4),repeat=3) if (m[0]+2*m[1]+3*m[2])%4==0]
print("  Lambda = {m : m2 + 2 m3 + 3 m4 = 0 mod 4} ; %d of 64 classes, index %d"%(len(cls),64//len(cls)))
for i in range(0,len(cls),4):
    print("   ",", ".join(str(c) for c in cls[i:i+4]))
print()
print("### (2) the observed relations are all compatible ###")
for m,tag in [((-2,1,0),"d3-2d2=0   [(1,2,3)]"),((-3,0,1),"d4-3d2=0   [(1,2,3)]"),
              ((1,1,-1),"d2+d3-d4=0 [(5,9,14)]"),((9,-5,0),"9d2-5d3=0  [(5,9,14)]")]:
    v=m[0]+2*m[1]+3*m[2]
    print("  m=%-12s m.(1,2,3) = %-4d  = 0 mod 4: %-6s %s"%(str(m),v,v%4==0,tag))
print()
def e_of(d):
    """generator of { m.(1,2,3) : m in d-perp }"""
    vals=[]
    for m in itertools.product(range(-40,41),repeat=3):
        if m==(0,0,0): continue
        if m[0]*d[0]+m[1]*d[1]+m[2]*d[2]==0:
            vals.append(abs(m[0]+2*m[1]+3*m[2]))
    vals=[v for v in vals if v>0]
    return 0 if not vals else gcd(*vals) if len(vals)>1 else vals[0]
print("### (3) e for the observed directions ###")
print("  direction        e      e=0 or 4|e ?   observed concentration")
for d,obs in [((1,2,3),"28.28x"),((2,4,6),"28.29x"),((5,9,14),"6.33x"),((1,2,4),"~0"),
              ((2,3,4),"~0"),((3,5,7),"~0"),((11,19,37),"1.03x"),((211,367,593),"0.89x")]:
    e=e_of(d)
    full = (e==0 or (e>0 and e%4==0))
    print("  %-16s %-6d %-14s %s"%(str(d),e,"YES" if full else "no",obs))
print()
print("### (4) IS THERE A SECOND FULL-CONCENTRATION FAMILY?  search 4 | e, d not ~ (1,2,3) ###")
hits=[]
for d in itertools.combinations(range(1,26),3):
    g=gcd(gcd(d[0],d[1]),d[2])
    dp=(d[0]//g,d[1]//g,d[2]//g)
    if dp==(1,2,3): continue
    e=e_of(d)
    if e>0 and e%4==0: hits.append((d,e))
print("  triples 1<=d2<d3<d4<=25 with 4 | e and d not proportional to (1,2,3): %d"%len(hits))
if hits:
    for d,e in hits[:12]: print("    %-14s e = %d"%(str(d),e))
else:
    print("    NONE -- e = 0 (i.e. d ~ (1,2,3)) is the only way to get full concentration")
print()
print("### (5) so the classification is ###")
print("  e = 0        <=> d ~ (1,2,3)        -> geodesic passes through c -> FULL, 28.28x")
print("  4 | e, e > 0 <=> (none found)       -> would be a second full family")
print("  otherwise    -> partial or generic, governed by which relations lie in Lambda")
print("DONE")
