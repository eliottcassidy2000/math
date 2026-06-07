#!/usr/bin/env python3
"""S647 — friendliness chart as a PNG, pure python (zlib only, no matplotlib)."""
import zlib, struct, binascii

# ---------------- the math: friendliness = covering depth ----------------
SPEEDS=[1,2,3,4,5,6,7,8,9,10,11,13,14]; N=14; DELTA=1.0/N
def clock(x):
    f=x-int(x);  f+=1 if f<0 else 0
    return min(f,1-f)
def friend(t): return sum(1 for v in SPEEDS if clock(v*t)<DELTA)
M=200000; hist=[0]*(len(SPEEDS)+1); series=[]
for i in range(M):
    t=i/M; f=friend(t); hist[f]+=1
    if i%(M//1400)==0: series.append((t,f))
p=[c/M for c in hist]; p0=p[0]; maxk=max(k for k in range(len(p)) if p[k]>0)
meanf=sum(k*p[k] for k in range(len(p)))

# ---------------- tiny framebuffer ----------------
W,H=960,520
buf=bytearray([15,20,36]*(W*H))     # dark navy bg
def px(x,y,c):
    if 0<=x<W and 0<=y<H:
        o=(y*W+x)*3; buf[o]=c[0]; buf[o+1]=c[1]; buf[o+2]=c[2]
def rect(x,y,w,h,c):
    for yy in range(int(y),int(y+h)):
        for xx in range(int(x),int(x+w)): px(xx,yy,c)
def hline(x0,x1,y,c):
    for x in range(int(x0),int(x1)): px(x,int(y),c)
def line(x0,y0,x1,y1,c):
    x0,y0,x1,y1=int(x0),int(y0),int(x1),int(y1)
    dx=abs(x1-x0); dy=-abs(y1-y0); sx=1 if x0<x1 else -1; sy=1 if y0<y1 else -1; e=dx+dy
    while True:
        px(x0,y0,c)
        if x0==x1 and y0==y1: break
        e2=2*e
        if e2>=dy: e+=dy; x0+=sx
        if e2<=dx: e+=dx; y0+=sy

# 3x5 font (digits + few symbols), scaled
F={'0':["111","101","101","101","111"],'1':["010","110","010","010","111"],
   '2':["111","001","111","100","111"],'3':["111","001","111","001","111"],
   '4':["101","101","111","001","001"],'5':["111","100","111","001","111"],
   '6':["111","100","111","101","111"],'7':["111","001","010","010","010"],
   '8':["111","101","111","101","111"],'9':["111","101","111","001","111"],
   '%':["101","001","010","100","101"],'.':["000","000","000","000","010"],
   ' ':["000","000","000","000","000"]}
def text(s,x,y,c,sc=2):
    cx=x
    for ch in s:
        g=F.get(ch,F[' '])
        for ry,row in enumerate(g):
            for rx,bit in enumerate(row):
                if bit=='1': rect(cx+rx*sc,y+ry*sc,sc,sc,c)
        cx+=4*sc

BLUE=(90,169,255); RED=(255,90,122); GRID=(34,44,74); PANEL=(22,29,51)
LT=(207,224,255); MUT=(126,140,181)

# ---------------- Panel A: friendliness over one lap ----------------
ax,ay,aw,ah=64,86,520,370
rect(ax,ay,aw,ah,PANEL)
ymax=maxk+1
def AX(t): return ax+t*aw
def AY(f): return ay+ah-(f/ymax)*ah
for k in range(0,ymax+1):
    yy=AY(k); hline(ax,ax+aw,yy,GRID)
    if k%2==0: text(str(k),ax-22,yy-5,MUT,2)
# filled step area
last=0
for idx,(t,f) in enumerate(series):
    x0=AX(series[idx-1][0]) if idx>0 else AX(0)
    x1=AX(t)
    yv=AY(f)
    # vertical fill from baseline to yv across [x0,x1]
    for xx in range(int(AX(0) if idx==0 else x0),int(x1)+1):
        for yy in range(int(yv),int(AY(0))):
            o=(yy*W+xx)*3
            if 0<=xx<W and 0<=yy<H:
                buf[o]=min(255,(buf[o]+BLUE[0])//2+30); buf[o+1]=(buf[o+1]+BLUE[1])//2; buf[o+2]=(buf[o+2]+BLUE[2])//2
    line(x0,yv,x1,yv,BLUE)
hline(ax,ax+aw,AY(0),RED)               # the lonely baseline (friendliness 0)
rect(ax,ay,aw,ah,PANEL) if False else None

# ---------------- Panel B: the distribution ----------------
bx,by,bw,bh=636,86,288,370
rect(bx,by,bw,bh,PANEL)
nb=maxk+1; pmax=max(p[:nb]); gap=4
bwid=(bw-gap*(nb+1))/nb
for k in range(nb):
    hh=(p[k]/pmax)*(bh-26)
    xx=bx+gap+k*(bwid+gap); yy=by+bh-hh
    rect(xx,yy,bwid,hh,RED if k==0 else BLUE)
    if k%1==0 and bwid>=10: text(str(k),xx+bwid/2-3,by+bh+4,MUT,2)
    elif k in (0,5,10,13): text(str(k),xx,by+bh+4,MUT,2)
# percent labels on the two tall friendly bars + the lonely bar
for k in (0,1,2,13):
    hh=(p[k]/pmax)*(bh-26); xx=bx+gap+k*(bwid+gap); yy=by+bh-hh
    text(str(int(round(p[k]*100)))+"%",xx-2,yy-20,LT,2)

png_rows=bytearray()
for y in range(H):
    png_rows.append(0); png_rows+=buf[y*W*3:(y+1)*W*3]
def chunk(t,d):
    c=t+d; return struct.pack(">I",len(d))+c+struct.pack(">I",binascii.crc32(c)&0xffffffff)
out=b'\x89PNG\r\n\x1a\n'+chunk(b'IHDR',struct.pack(">IIBBBBB",W,H,8,2,0,0,0))+\
    chunk(b'IDAT',zlib.compress(bytes(png_rows),9))+chunk(b'IEND',b'')
open("friendliness_chart_s647.png","wb").write(out)
print(f"p0(lonely)={p0:.4f}  friendly={1-p0:.4f}  mean={meanf:.2f}  maxk={maxk}")
print("wrote friendliness_chart_s647.png", len(out),"bytes")
