import itertools, math, cmath
# triangular lattice point -> complex coordinate (Eisenstein: 1 and w=e^{i pi/3})
w = cmath.exp(1j*math.pi/3)
def Z(i,j): return i + j*w
def disk(R): return [(i,j) for i in range(-R,R+1) for j in range(-R,R+1) if i*i+i*j+j*j<=R*R]
def count_unit(P, tol=1e-7):
    e=0
    for a,b in itertools.combinations(P,2):
        if abs(abs(a-b)-1)<tol: e+=1
    return e

# unit-modulus elements u = a/conj(a) of Z[w'] (Eisenstein, w'=e^{2 pi i/3}) -> rotations
# that map the lattice's directions; small-height ones add the most unit distances.
we = cmath.exp(2j*math.pi/3)
def eisen(a,b): return a + b*we
units=[]
for a in range(-4,5):
    for b in range(-4,5):
        z=eisen(a,b)
        if abs(z)<1e-9: continue
        u=z/z.conjugate()           # |u|=1, an algebraic number of modulus 1, bounded height
        ang=math.degrees(cmath.phase(u))%60
        if 0.5<ang<59.5: units.append((round(abs(z)**2),a,b,u,ang))
units=sorted(set((h,round(ang,4)) for h,a,b,u,ang in units))
# build base patch (the 49-edge-ish compact region), overlay a rotated copy, pick best 22-subset by edges
L=disk(4)
base=[Z(i,j) for (i,j) in L]
best=(0,None)
for a in range(-3,4):
  for b in range(-3,4):
    z=eisen(a,b)
    if abs(z)<1e-9: continue
    u=z/z.conjugate()
    # overlay: base ∪ u*base (rotated about origin); merge near-coincident points
    pts=list(base)
    for p in base:
        q=u*p
        if all(abs(q-r)>1e-6 for r in pts): pts.append(q)
    # take the 22 points of highest local unit-degree
    deg={p: sum(abs(abs(p-r)-1)<1e-7 for r in pts) for p in pts}
    sub=sorted(pts,key=lambda p:-deg[p])[:22]
    e=count_unit(sub)
    if e>best[0]: best=(e,(a,b,round(math.degrees(cmath.phase(u))%360,3)))
print("distinct small modulus-1 rotation angles (deg, mod 60) from Z[w]:",[a for _,a in units[:12]])
print("best 22-pt UNION (lattice + one CM rotation) unit distances:",best[0],"via rotation (a,b,deg)=",best[1])
print("vs pure-lattice baseline 49; optimum u(22) in {60,61}")
