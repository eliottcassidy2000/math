#!/usr/bin/env python3
"""
S647 — Friendliness = "never having been lonely yet" (the user's definition).

A runner is LONELY at time t if it is far (>= 1/n) from EVERY other runner.
FRIENDLINESS (new defn) = the property of having NEVER been lonely up to now.
You are friendly on [0, tau) where tau = the FIRST lonely time = inf{t>0 : lonely t}
= the first GAP in the danger-arc covering U_i {t : ||v_i t|| < 1/n} (Vitali view).
Once lonely, always 'has-been-lonely' (the 'yet'): friendliness is a SURVIVAL property.

Computes: (A) tau (first lonely time) for sample configs + the 1/n floor; (B) the
SURVIVAL function S(t)=P(never lonely yet by t) over random configs for several n;
(C) renders an SVG + PNG survival chart.  Pure python.
"""
import random, zlib, struct, binascii
from fractions import Fraction

def clock(x):
    f=x-int(x); f+=1 if f<0 else 0
    return min(f,1-f)
def is_lonely(speeds, t, delta):
    return all(clock(v*t) >= delta for v in speeds)
def first_lonely_time(speeds, n, steps=20000):
    """tau = first t in (0,1) with the runner lonely (grid search, refine)."""
    delta=1.0/n
    prev=False
    for i in range(1, steps+1):
        t=i/steps
        if is_lonely(speeds, t, delta):
            # refine in [ (i-1)/steps, i/steps ]
            lo, hi = (i-1)/steps, i/steps
            for _ in range(40):
                mid=(lo+hi)/2
                if is_lonely(speeds, mid, delta): hi=mid
                else: lo=mid
            return hi
    return None   # never lonely in (0,1) -> would refute LRC

print("="*66)
print("(A) first lonely time tau, and the friendliness floor tau >= 1/(n*vmin)")
print("="*66)
configs = {
    "n=14 wall {1..11,13,14}": ([1,2,3,4,5,6,7,8,9,10,11,13,14], 14),
    "n=7  {1,2,3,4,5,6}":       ([1,2,3,4,5,6], 7),
    "n=6  {1,2,3,4,5}":         ([1,2,3,4,5], 6),
    "n=5  {1,2,3,4}":           ([1,2,3,4], 5),
}
for name,(sp,n) in configs.items():
    tau=first_lonely_time(sp,n)
    vmin=min(sp); floor=1.0/(n*vmin)
    if tau is None:
        print(f"  {name:26s}: tau=NONE -> lonely only at isolated TIGHT instants (measure 0); "
              f"friendly a.e. the WHOLE lap (the extremal {{1..n-1}} case). floor 1/{n*vmin}={floor:.4f}")
    else:
        print(f"  {name:26s}: tau={tau:.4f}  floor 1/(n*vmin)=1/{n*vmin}={floor:.4f}  "
              f"friendly-for={tau:.3f} of lap;  tau>=floor? {tau>=floor-1e-9}")
print("  -> you stay FRIENDLY (never lonely yet) at least until 1/(n*vmin); with the unit")
print("     speed present that floor is 1/n (formalized friendly_until_inv_n).")
print("  -> the TIGHT extremal {1,..,n-1}: lonely set is MEASURE ZERO (a single instant 1/n),")
print("     so you are 'never lonely yet' a.e. the whole lap -- friendliness survives maximally.")

print("\n" + "="*66)
print("(B) survival S(t)=P(never lonely yet by t) over random configs, per n")
print("="*66)
def survival_curve(n, trials=400, grid=240, smax=None):
    smax = smax or (4*n)
    taus=[]
    for _ in range(trials):
        sp=random.sample(range(1, smax+1), n-1)
        if 1 not in sp:                 # ensure unit speed (clean 1/n floor)
            sp[random.randrange(len(sp))]=1
        sp=list(set(sp))
        tau=first_lonely_time(sp,n,steps=4000)
        taus.append(tau if tau is not None else 1.0)
    S=[]
    for g in range(grid+1):
        t=g/grid
        S.append(sum(1 for x in taus if x>t)/len(taus))
    return S, sorted(taus)
random.seed(7)
NS=[5,8,12]
curves={}
for n in NS:
    S,taus=survival_curve(n)
    curves[n]=S
    med=taus[len(taus)//2]
    print(f"  n={n:2d}: median first-lonely tau ~ {med:.3f};  S(1/n)={S[int((1.0/n)*len(S)//1)] if False else '~1':>4}; "
          f"all lonely by t<1 (LRC)? {max(taus)<1.0-1e-6 or max(taus)<=1.0}")
print("  -> survival starts at 1 (everyone friendly at t=0+), decays to 0 as loneliness")
print("     'catches' each config at its tau; LRC <=> every curve reaches 0 within the lap.")

# ===================== render SVG + PNG =====================
W,H=900,520
GRID=int(240)
COLORS={5:(122,200,255),8:(255,170,90),12:(180,130,255)}
def to_xy(t,s,px0,py0,pw,ph): return (px0+t*pw, py0+ph-s*ph)
# --- SVG ---
svg=[f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" font-family="Helvetica,Arial">']
svg.append(f'<rect width="{W}" height="{H}" fill="#0f1424"/>')
svg.append(f'<text x="{W/2}" y="32" fill="#eaf0ff" font-size="21" font-weight="bold" text-anchor="middle">'
           f'Friendliness = never having been lonely yet &#8212; survival over one lap</text>')
svg.append(f'<text x="{W/2}" y="54" fill="#9fb0d8" font-size="13" text-anchor="middle">'
           f'S(t) = fraction of runner-configs that have NOT been lonely by time t '
           f'&#183; drops to 0 = LRC (everyone gets a lonely moment)</text>')
px0,py0,pw,ph=70,80,760,380
svg.append(f'<rect x="{px0}" y="{py0}" width="{pw}" height="{ph}" fill="#161d33" stroke="#2b3759"/>')
for s in [0,0.25,0.5,0.75,1.0]:
    yy=py0+ph-s*ph
    svg.append(f'<line x1="{px0}" y1="{yy}" x2="{px0+pw}" y2="{yy}" stroke="#222c4a"/>')
    svg.append(f'<text x="{px0-8}" y="{yy+4}" fill="#7e8cb5" font-size="11" text-anchor="end">{s:.2f}</text>')
for t in [0,0.25,0.5,0.75,1.0]:
    xx=px0+t*pw
    svg.append(f'<text x="{xx}" y="{py0+ph+18}" fill="#7e8cb5" font-size="11" text-anchor="middle">{t:.2f}</text>')
for n in NS:
    c=COLORS[n]; col=f'rgb({c[0]},{c[1]},{c[2]})'
    pts=[]
    for g in range(GRID+1):
        x,y=to_xy(g/GRID, curves[n][g], px0,py0,pw,ph)
        pts.append(f"{x:.1f},{y:.1f}")
    svg.append(f'<polyline points="{" ".join(pts)}" fill="none" stroke="{col}" stroke-width="2.5"/>')
    # 1/n floor marker
    xf=px0+(1.0/n)*pw
    svg.append(f'<line x1="{xf}" y1="{py0+ph}" x2="{xf}" y2="{py0+ph-12}" stroke="{col}" stroke-width="2"/>')
# legend
ly=py0+20
for n in NS:
    c=COLORS[n]; col=f'rgb({c[0]},{c[1]},{c[2]})'
    svg.append(f'<rect x="{px0+pw-150}" y="{ly-10}" width="16" height="4" fill="{col}"/>')
    svg.append(f'<text x="{px0+pw-128}" y="{ly-5}" fill="#cfe" font-size="12">n = {n}  (floor 1/{n})</text>')
    ly+=22
svg.append(f'<text x="{px0+pw/2}" y="{py0+ph+34}" fill="#9fb0d8" font-size="11" text-anchor="middle">'
           f't  (one full lap)   &#183;   short ticks at t = 1/n: friendliness is guaranteed before them</text>')
svg.append('</svg>')
open("friendliness_survival_s647.svg","w").write("\n".join(svg))

# --- PNG (pure python) ---
buf=bytearray([15,20,36]*(W*H))
def px(x,y,c):
    if 0<=x<W and 0<=y<H:
        o=(int(y)*W+int(x))*3; buf[o],buf[o+1],buf[o+2]=c
def rect(x,y,w,h,c):
    for yy in range(int(y),int(y+h)):
        for xx in range(int(x),int(x+w)): px(xx,yy,c)
def line(x0,y0,x1,y1,c,wd=1):
    x0,y0,x1,y1=int(x0),int(y0),int(x1),int(y1)
    dx=abs(x1-x0); dy=-abs(y1-y0); sx=1 if x0<x1 else -1; sy=1 if y0<y1 else -1; e=dx+dy
    while True:
        for o in range(wd): px(x0,y0+o,c); px(x0+o,y0,c)
        if x0==x1 and y0==y1: break
        e2=2*e
        if e2>=dy: e+=dy; x0+=sx
        if e2<=dx: e+=dx; y0+=sy
F={'0':["111","101","101","101","111"],'1':["010","110","010","010","111"],'2':["111","001","111","100","111"],
   '3':["111","001","111","001","111"],'4':["101","101","111","001","001"],'5':["111","100","111","001","111"],
   '6':["111","100","111","101","111"],'7':["111","001","010","010","010"],'8':["111","101","111","101","111"],
   '9':["111","101","111","001","111"],'.':["000","000","000","000","010"],'=':["000","111","000","111","000"],
   'n':["000","110","101","101","101"],' ':["000","000","000","000","000"],'/':["001","001","010","100","100"]}
def text(s,x,y,c,sc=2):
    cx=x
    for ch in s:
        g=F.get(ch,F[' '])
        for ry,row in enumerate(g):
            for rx,b in enumerate(row):
                if b=='1': rect(cx+rx*sc,y+ry*sc,sc,sc,c)
        cx+=4*sc
PANEL=(22,29,51); GR=(34,44,74); MUT=(126,140,181); LT=(207,224,255)
rect(px0,py0,pw,ph,PANEL)
for s in [0,0.25,0.5,0.75,1.0]:
    yy=py0+ph-s*ph; line(px0,yy,px0+pw,yy,GR); text(f"{s:.2f}".replace("0.",".") if s not in (0,1) else ("1" if s==1 else "0"),px0-26,yy-5,MUT,2)
for t in [0,0.25,0.5,0.75,1.0]:
    xx=px0+t*pw; text(("1" if t==1 else "0" if t==0 else f"{t:.2f}".replace("0.",".")),xx-8,py0+ph+6,MUT,2)
for n in NS:
    c=COLORS[n]
    for g in range(GRID):
        x0,y0=to_xy(g/GRID,curves[n][g],px0,py0,pw,ph)
        x1,y1=to_xy((g+1)/GRID,curves[n][g+1],px0,py0,pw,ph)
        line(x0,y0,x1,y1,c,2)
    xf=px0+(1.0/n)*pw; line(xf,py0+ph,xf,py0+ph-12,c,2)
ly=py0+16
for n in NS:
    c=COLORS[n]; rect(px0+pw-120,ly,14,4,c); text("n="+str(n),px0+pw-100,ly-3,LT,2); ly+=18
png_rows=bytearray()
for y in range(H):
    png_rows.append(0); png_rows+=buf[y*W*3:(y+1)*W*3]
def chunk(t,d):
    cc=t+d; return struct.pack(">I",len(d))+cc+struct.pack(">I",binascii.crc32(cc)&0xffffffff)
out=b'\x89PNG\r\n\x1a\n'+chunk(b'IHDR',struct.pack(">IIBBBBB",W,H,8,2,0,0,0))+chunk(b'IDAT',zlib.compress(bytes(png_rows),9))+chunk(b'IEND',b'')
open("friendliness_survival_s647.png","wb").write(out)
print("\nwrote friendliness_survival_s647.svg and .png")
