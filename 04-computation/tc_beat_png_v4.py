#!/usr/bin/env python3
"""
Tournament Codec v4: BEAT PNG AT ALL SIZES.

Key fix: use zlib.compressobj with Z_FILTERED strategy (same as PNG).
PNG uses this strategy to optimize for filtered byte streams.
Also: try Z_DEFAULT_STRATEGY and Z_HUFFMAN_ONLY, pick smallest.

kind-pasteur-2026-03-25-S2
"""
import sys, io, struct, zlib, time
import numpy as np
from PIL import Image

sys.stdout.reconfigure(line_buffering=True)

# ============================================================
# PREDICTORS
# ============================================================
def pix(a,x,d=0):
    if a is not None and 0<=x<len(a): return int(a[x])
    return d

def _pa(a,b,c):
    p=a+b-c; pa,pb,pc=abs(p-a),abs(p-b),abs(p-c)
    if pa<=pb and pa<=pc: return a
    if pb<=pc: return b
    return c
def _md(a,b,c):
    if c>=max(a,b): return min(a,b)
    if c<=min(a,b): return max(a,b)
    return a+b-c
def _gp(r,p,x,w):
    a,b,c=pix(r,x-1),pix(p,x),pix(p,x-1)
    d,e=pix(p,x+1,b),pix(r,x-2)
    dh,dv=abs(a-c)+abs(b-d),abs(b-c)+abs(a-e)
    if dh>2*dv: return b
    if dv>2*dh: return a
    return (a+b)>>1

FILTERS=[
    lambda r,p,x,w:0,
    lambda r,p,x,w:pix(r,x-1),
    lambda r,p,x,w:pix(p,x),
    lambda r,p,x,w:(pix(r,x-1)+pix(p,x))>>1,
    lambda r,p,x,w:_pa(pix(r,x-1),pix(p,x),pix(p,x-1)),
    lambda r,p,x,w:_md(pix(r,x-1),pix(p,x),pix(p,x-1)),
    lambda r,p,x,w:_gp(r,p,x,w),
    lambda r,p,x,w:pix(p,x-1),
]
NF=len(FILTERS)

def frow(row,prev,fid,w):
    pf=FILTERS[fid]; o=np.empty(w,dtype=np.uint8)
    for x in range(w): o[x]=(int(row[x])-pf(row,prev,x,w))&0xFF
    return o

def zc(data, level=9):
    """Try multiple zlib strategies, return smallest."""
    results = []
    for strategy in [zlib.Z_DEFAULT_STRATEGY, zlib.Z_FILTERED, zlib.Z_HUFFMAN_ONLY]:
        try:
            obj = zlib.compressobj(level, zlib.DEFLATED, 15, 9, strategy)
            c = obj.compress(data) + obj.flush()
            results.append(c)
        except:
            pass
    # Also try standard compress
    results.append(zlib.compress(data, level))
    return min(results, key=len)

def compress_plane(plane):
    h,w=plane.shape; cands=[]
    # A: per-row optimal (all filters)
    buf=bytearray()
    for y in range(h):
        row,prev=plane[y],(plane[y-1] if y>0 else None)
        bs=1<<30; bf=0; br=None
        for f in range(NF):
            r2=frow(row,prev,f,w); s=r2.astype(np.int16); s[s>128]-=256
            a2=int(np.sum(np.abs(s)))
            if a2<bs: bs=a2; bf=f; br=r2
        buf.append(bf); buf.extend(br)
    cands.append(zc(bytes(buf)))
    # B: whole-frame per filter
    for f in range(NF):
        b2=bytearray()
        for y in range(h): b2.extend(frow(plane[y],plane[y-1] if y>0 else None,f,w))
        cands.append(zc(bytes(b2)))
    # C: raw
    cands.append(zc(plane.tobytes()))
    # D: transpose raw
    cands.append(zc(plane.T.copy().tobytes()))
    # E: delta-row
    be=bytearray(plane[0])
    for y in range(1,h): be.extend(((plane[y].astype(int)-plane[y-1].astype(int))&0xFF).astype(np.uint8))
    cands.append(zc(bytes(be)))
    # F: PNG-style interleaved, per-row optimal (5 PNG filters only)
    bf2=bytearray()
    for y in range(h):
        row,prev=plane[y],(plane[y-1] if y>0 else None)
        bs2=1<<30; bi2=0; br2=None
        for f in range(5):
            r2=frow(row,prev,f,w); s=r2.astype(np.int16); s[s>128]-=256
            a2=int(np.sum(np.abs(s)))
            if a2<bs2: bs2=a2; bi2=f; br2=r2
        bf2.append(bi2); bf2.extend(br2)
    cands.append(zc(bytes(bf2)))
    # G: transpose + per-row
    pt=plane.T.copy(); ht2,wt2=pt.shape
    bg=bytearray()
    for y in range(ht2):
        row,prev=pt[y],(pt[y-1] if y>0 else None)
        bs3=1<<30; bi3=0; br3=None
        for f in range(NF):
            r2=frow(row,prev,f,wt2); s=r2.astype(np.int16); s[s>128]-=256
            a2=int(np.sum(np.abs(s)))
            if a2<bs3: bs3=a2; bi3=f; br3=r2
        bg.append(bi3); bg.extend(br3)
    cands.append(zc(bytes(bg)))

    bi=min(range(len(cands)),key=lambda i:len(cands[i]))
    return bytes([bi])+cands[bi]

def ct_ycocg(img):
    r,g,b=img[:,:,0].astype(int),img[:,:,1].astype(int),img[:,:,2].astype(int)
    Co=r-b; tmp=b+(Co>>1); Cg=g-tmp; Y=tmp+(Cg>>1)
    return np.stack([Y&0xFF,Co&0xFF,Cg&0xFF],axis=-1).astype(np.uint8)
def ct_rct(img):
    r,g,b=img[:,:,0].astype(int),img[:,:,1].astype(int),img[:,:,2].astype(int)
    return np.stack([(r+2*g+b)>>2&0xFF,(b-g)&0xFF,(r-g)&0xFF],axis=-1).astype(np.uint8)
def ct_grd(img):
    r,g,b=img[:,:,0].astype(int),img[:,:,1].astype(int),img[:,:,2].astype(int)
    return np.stack([g&0xFF,(r-g)&0xFF,(b-g)&0xFF],axis=-1).astype(np.uint8)
CTS=[lambda i:i, ct_ycocg, ct_rct, ct_grd]

def compress_image(img):
    h,w=img.shape[:2]; rgb=img.ndim==3 and img.shape[2]==3
    if rgb:
        bt=1<<60; bd=None
        for ci,ct in enumerate(CTS):
            t=ct(img); pd=[compress_plane(t[:,:,c]) for c in range(3)]
            tt2=sum(len(d) for d in pd)
            if tt2<bt:
                bt=tt2; hdr=struct.pack('<HHB',w,h,0x30|ci)
                bd=hdr+b''.join(struct.pack('<I',len(d))+d for d in pd)
        return bd
    else:
        p=img if img.ndim==2 else img[:,:,0]
        pd=compress_plane(p); return struct.pack('<HH',w,h)+pd

def png_size(img):
    pil=Image.fromarray(img,'L' if img.ndim==2 else 'RGB')
    buf=io.BytesIO(); pil.save(buf,format='PNG',optimize=True,compress_level=9)
    return buf.tell()

# ============================================================
# TESTS
# ============================================================
print("="*76)
print("  TC v4 vs PNG: MULTI-SIZE STRESS TEST")
print("  kind-pasteur-2026-03-25-S2")
print("="*76)

results=[]
W,T,L=0,0,0

def test(name,img):
    global W,T,L
    ps=png_size(img); ts=len(compress_image(img))
    r=ps/ts if ts>0 else 0
    if r>1.001: res="WIN"; W+=1
    elif r<0.999: res="LOSS"; L+=1
    else: res="TIE"; T+=1
    print(f"  {name:<28} PNG={ps:>8} TC={ts:>8} {r:>7.4f} {res}")
    results.append((name,ps,ts,r,res))

for SZ in [128, 256]:
    print(f"\n--- {SZ}x{SZ} ---")
    np.random.seed(42)
    test(f"{SZ}_solid", np.full((SZ,SZ),128,dtype=np.uint8))
    test(f"{SZ}_grad_h", np.tile(np.linspace(0,255,SZ,dtype=np.uint8),(SZ,1)))
    x,y=np.meshgrid(np.arange(SZ),np.arange(SZ))
    test(f"{SZ}_checker8", ((x//8+y//8)%2*255).astype(np.uint8))
    test(f"{SZ}_circles", (np.sin(np.sqrt((x-SZ//2)**2+(y-SZ//2)**2)/5)*127+128).astype(np.uint8))
    test(f"{SZ}_random", np.random.randint(0,256,(SZ,SZ),dtype=np.uint8))
    test(f"{SZ}_dithered_bw", np.random.choice([0,255],(SZ,SZ)).astype(np.uint8))

    for d in [0.01, 0.05, 0.10, 0.20]:
        im=np.zeros((SZ,SZ),dtype=np.uint8); n=int(SZ*SZ*d)
        im[np.random.randint(0,SZ,n),np.random.randint(0,SZ,n)]=np.random.randint(1,256,n).astype(np.uint8)
        test(f"{SZ}_sparse_{d:.0%}", im)

    sm=np.random.randint(0,256,(max(SZ//16,2),max(SZ//16,2)),dtype=np.uint8)
    base=np.array(Image.fromarray(sm).resize((SZ,SZ),Image.BILINEAR)).astype(float)
    test(f"{SZ}_natural_n10", np.clip(base+np.random.normal(0,10,(SZ,SZ)),0,255).astype(np.uint8))
    test(f"{SZ}_natural_n25", np.clip(base+np.random.normal(0,25,(SZ,SZ)),0,255).astype(np.uint8))

    test(f"{SZ}_rgb_nat", np.stack([np.clip(np.array(Image.fromarray(np.random.randint(0,256,(max(SZ//16,2),max(SZ//16,2)),dtype=np.uint8)).resize((SZ,SZ),Image.BILINEAR)).astype(float)+np.random.normal(0,10,(SZ,SZ)),0,255).astype(np.uint8) for _ in range(3)],axis=-1))

print(f"\n{'='*76}")
tp=sum(r[1] for r in results); tt=sum(r[2] for r in results)
print(f"  {W} WINS / {T} TIES / {L} LOSSES out of {len(results)}")
print(f"  Aggregate: {tp/tt:.3f}x")

if L>0:
    print(f"\n  LOSSES:"); [print(f"    {n}: PNG={p} TC={t} ({p/t:.4f})") for n,p,t,r,res in results if res=="LOSS"]
if T>0:
    print(f"\n  TIES:"); [print(f"    {n}: PNG={p} TC={t}") for n,p,t,r,res in results if res=="TIE"]
