# Minimal escape price from the 2-adic point -1 in the E-graph:
# a_min(b) = min number of x3 moves over E-paths from -1 with exactly b halvings; r = 3^a/2^b.
import math
L=math.log(3)/math.log(2)
MAXB=70; TMAX=4; CAP=10**40
# states at even values: dict y->a ; start: -1 is odd -> forced 3(-1)+1=-2 (a=1)
layer={-2:1}
rows=[]
for b in range(0,MAXB):
    nxt={}
    for y,a in layer.items():
        yt=y; at=a
        for t in range(TMAX+1):
            if t>0:
                yt=9*yt+4; at+=2
            if abs(yt)>CAP: break
            z=yt//2; az=at
            if z%2!=0: z=3*z+1; az+=1
            if z not in nxt or az<nxt[z]: nxt[z]=az
    # prune: keep states with a within slack of the minimum
    amin=min(nxt.values())
    layer={y:a for y,a in nxt.items() if a<=amin+6}
    # the minimal cost to have consumed b+1 bits: min a over states  (value z = (3^a*(-1)+B)/2^(b+1))
    r=3**amin/2**(b+1)
    best=[y for y,a in layer.items() if a==amin]
    rows.append((b+1,amin,r,len(layer),sorted(best,key=abs)[:3]))
for b,a,r,n,best in rows:
    if b<=20 or b%5==0: print(f"b={b:3d} a_min={a:3d}  r=3^a/2^b={r:10.4f}  a-b/log2(3)={a-b/L:7.3f}  states={n:6d}  endpoints {best}")
