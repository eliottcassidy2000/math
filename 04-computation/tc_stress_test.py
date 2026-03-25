#!/usr/bin/env python3
"""
TC v3 Stress Test: verify 100% win rate at 256x256 and with adversarial images.

kind-pasteur-2026-03-25-S2
"""
import sys, io, struct, zlib, time
import numpy as np
from PIL import Image

sys.stdout.reconfigure(line_buffering=True)

# Import the codec from v3 (inline the key functions for speed)
def pix(arr, x, default=0):
    if arr is not None and 0 <= x < len(arr): return int(arr[x])
    return default

def _paeth(a,b,c):
    p=a+b-c; pa,pb,pc=abs(p-a),abs(p-b),abs(p-c)
    if pa<=pb and pa<=pc: return a
    if pb<=pc: return b
    return c
def _med(a,b,c):
    if c>=max(a,b): return min(a,b)
    if c<=min(a,b): return max(a,b)
    return a+b-c
def _gap(r,p,x,w):
    a,b,c=pix(r,x-1),pix(p,x),pix(p,x-1)
    d=pix(p,x+1,b); e=pix(r,x-2)
    dh=abs(a-c)+abs(b-d); dv=abs(b-c)+abs(a-e)
    if dh>2*dv: return b
    if dv>2*dh: return a
    return (a+b)>>1
def _avg4(r,p,x,w):
    v=[pix(r,x-1)]
    if p is not None:
        v.append(pix(p,x)); v.append(pix(p,x-1))
        if x+1<w: v.append(pix(p,x+1))
    return sum(v)//len(v)

FILTERS = [
    lambda r,p,x,w: 0, lambda r,p,x,w: pix(r,x-1), lambda r,p,x,w: pix(p,x),
    lambda r,p,x,w: (pix(r,x-1)+pix(p,x))>>1, lambda r,p,x,w: _paeth(pix(r,x-1),pix(p,x),pix(p,x-1)),
    lambda r,p,x,w: _med(pix(r,x-1),pix(p,x),pix(p,x-1)), lambda r,p,x,w: _gap(r,p,x,w),
    lambda r,p,x,w: pix(p,x-1), lambda r,p,x,w: _avg4(r,p,x,w),
]
NF = len(FILTERS)

def filter_row(row, prev, fid, w):
    pf=FILTERS[fid]; out=np.empty(w,dtype=np.uint8)
    for x in range(w): out[x]=(int(row[x])-pf(row,prev,x,w))&0xFF
    return out

def compress_plane(plane, level=9):
    h,w=plane.shape; candidates=[]
    # A: per-row optimal
    buf=bytearray()
    for y in range(h):
        row,prev=plane[y],(plane[y-1] if y>0 else None)
        bs=float('inf'); bf=0; br=None
        for fid in range(NF):
            res=filter_row(row,prev,fid,w); s=res.astype(np.int16); s[s>128]-=256
            a2=np.sum(np.abs(s))
            if a2<bs: bs=a2; bf=fid; br=res
        buf.append(bf); buf.extend(br)
    candidates.append(zlib.compress(bytes(buf),level))
    # B: whole-frame each filter
    for fid in range(NF):
        b2=bytearray()
        for y in range(h):
            b2.extend(filter_row(plane[y],plane[y-1] if y>0 else None,fid,w))
        candidates.append(zlib.compress(bytes(b2),level))
    # C: transpose + per-row
    pt=plane.T.copy(); ht2,wt2=pt.shape
    buf3=bytearray()
    for y in range(ht2):
        row,prev=pt[y],(pt[y-1] if y>0 else None)
        bs=float('inf'); bf=0; br=None
        for fid in range(NF):
            res=filter_row(row,prev,fid,wt2); s=res.astype(np.int16); s[s>128]-=256
            a2=np.sum(np.abs(s))
            if a2<bs: bs=a2; bf=fid; br=res
        buf3.append(bf); buf3.extend(br)
    candidates.append(zlib.compress(bytes(buf3),level))
    # D: raw
    candidates.append(zlib.compress(plane.tobytes(),level))
    # E: delta-row
    be=bytearray(plane[0])
    for y in range(1,h): be.extend(((plane[y].astype(int)-plane[y-1].astype(int))&0xFF).astype(np.uint8))
    candidates.append(zlib.compress(bytes(be),level))
    # F: PNG-style none filter
    bf2=bytearray()
    for y in range(h): bf2.append(0); bf2.extend(bytes(plane[y]))
    candidates.append(zlib.compress(bytes(bf2),level))
    # G: PNG-style per-row (5 filters)
    bg=bytearray()
    for y in range(h):
        row,prev=plane[y],(plane[y-1] if y>0 else None)
        bs=float('inf'); bfi=0; br2=None
        for fid in [0,1,2,3,4]:
            res=filter_row(row,prev,fid,w); s=res.astype(np.int16); s[s>128]-=256
            a2=np.sum(np.abs(s))
            if a2<bs: bs=a2; bfi=fid; br2=res
        bg.append(bfi); bg.extend(br2)
    candidates.append(zlib.compress(bytes(bg),level))
    # H: transpose raw
    candidates.append(zlib.compress(plane.T.copy().tobytes(),level))

    bi=min(range(len(candidates)),key=lambda i:len(candidates[i]))
    return bytes([bi])+candidates[bi]

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
        bt=float('inf'); bd=None
        for ci,ct in enumerate(CTS):
            t=ct(img); pd=[compress_plane(t[:,:,c]) for c in range(3)]
            tt2=sum(len(d) for d in pd)
            if tt2<bt:
                bt=tt2; hdr=struct.pack('<HHB',w,h,0x30|ci)
                bd=hdr+b''.join(struct.pack('<H',len(d))+d for d in pd)
        return bd
    else:
        p=img if img.ndim==2 else img[:,:,0]
        pd=compress_plane(p); return struct.pack('<HH',w,h)+pd

def png_size(img):
    pil=Image.fromarray(img,'L' if img.ndim==2 else 'RGB')
    buf=io.BytesIO(); pil.save(buf,format='PNG',optimize=True,compress_level=9)
    return buf.tell()

# ============================================================
# STRESS TESTS at 256x256
# ============================================================

print("=" * 76)
print("  STRESS TEST: TC v3 vs PNG at 256x256 + adversarial images")
print("  kind-pasteur-2026-03-25-S2")
print("=" * 76)

SZ = 256
results = []
W,T,L = 0,0,0

def test(name, img):
    global W,T,L
    ps=png_size(img); ts=len(compress_image(img))
    r=ps/ts if ts>0 else 0
    if r>1.001: res="WIN"; W+=1
    elif r<0.999: res="LOSS"; L+=1
    else: res="TIE"; T+=1
    print(f"  {name:<25} PNG={ps:>8} TC={ts:>8} ratio={r:.4f} {res}")
    results.append((name,ps,ts,r,res))

print(f"\n--- Standard tests at {SZ}x{SZ} ---")
np.random.seed(42)
test("solid", np.full((SZ,SZ),128,dtype=np.uint8))
test("grad_h", np.tile(np.linspace(0,255,SZ,dtype=np.uint8),(SZ,1)))
test("grad_v", np.tile(np.linspace(0,255,SZ,dtype=np.uint8).reshape(-1,1),(1,SZ)))
x,y=np.meshgrid(np.arange(SZ),np.arange(SZ))
test("grad_diag", ((x+y)*255//(2*SZ-2)).astype(np.uint8))
test("checker8", ((x//8+y//8)%2*255).astype(np.uint8))
test("checker1", ((x+y)%2*255).astype(np.uint8))
test("circles", (np.sin(np.sqrt((x-SZ//2)**2+(y-SZ//2)**2)/5)*127+128).astype(np.uint8))

print(f"\n--- Adversarial tests ---")
np.random.seed(123)
# Nearly incompressible
test("random_uniform", np.random.randint(0,256,(SZ,SZ),dtype=np.uint8))
# Very low entropy
test("binary_noise", np.random.choice([0,255],(SZ,SZ)).astype(np.uint8))
# Sparse at various densities
for d in [0.001, 0.01, 0.05, 0.10, 0.20, 0.50]:
    im=np.zeros((SZ,SZ),dtype=np.uint8); n=int(SZ*SZ*d)
    im[np.random.randint(0,SZ,n),np.random.randint(0,SZ,n)]=np.random.randint(1,256,n).astype(np.uint8)
    test(f"sparse_{d:.1%}", im)

# Natural photo simulations
np.random.seed(555)
sm=np.random.randint(0,256,(SZ//8,SZ//8),dtype=np.uint8)
base=np.array(Image.fromarray(sm).resize((SZ,SZ),Image.BILINEAR)).astype(float)
test("smooth_16x", np.clip(base,0,255).astype(np.uint8))
test("natural_noise5", np.clip(base+np.random.normal(0,5,(SZ,SZ)),0,255).astype(np.uint8))
test("natural_noise15", np.clip(base+np.random.normal(0,15,(SZ,SZ)),0,255).astype(np.uint8))
test("natural_noise30", np.clip(base+np.random.normal(0,30,(SZ,SZ)),0,255).astype(np.uint8))

# Dithered
g=np.linspace(0,255,SZ).reshape(1,-1).repeat(SZ,axis=0).astype(float)
out=np.zeros((SZ,SZ),dtype=np.uint8)
for y2 in range(SZ):
    for x2 in range(SZ):
        o=g[y2,x2]; n=255 if o>128 else 0; out[y2,x2]=n; e=o-n
        if x2+1<SZ: g[y2,x2+1]+=e*7/16
        if y2+1<SZ:
            if x2>0: g[y2+1,x2-1]+=e*3/16
            g[y2+1,x2]+=e*5/16
            if x2+1<SZ: g[y2+1,x2+1]+=e/16
test("dithered", out)

# Text
im=np.full((SZ,SZ),250,dtype=np.uint8)
for _ in range(SZ//4):
    yy=np.random.randint(0,SZ); xs2=np.random.randint(0,SZ//2)
    im[yy,xs2:min(xs2+np.random.randint(5,SZ//2),SZ)]=np.random.randint(0,20)
test("text", im)

# Blocks
im2=np.zeros((SZ,SZ),dtype=np.uint8)
for _ in range(30):
    x0,y0=np.random.randint(0,SZ-16),np.random.randint(0,SZ-16)
    bw,bh=np.random.randint(8,64),np.random.randint(8,64)
    im2[y0:min(y0+bh,SZ),x0:min(x0+bw,SZ)]=np.random.randint(0,256)
test("blocks", im2)

print(f"\n--- RGB tests ---")
np.random.seed(777)
test("rgb_natural", np.stack([np.clip(np.array(Image.fromarray(np.random.randint(0,256,(SZ//16,SZ//16),dtype=np.uint8)).resize((SZ,SZ),Image.BILINEAR)).astype(float)+np.random.normal(0,10,(SZ,SZ)),0,255).astype(np.uint8) for _ in range(3)],axis=-1))

im3=np.zeros((SZ,SZ,3),dtype=np.uint8)
for _ in range(20):
    x0,y0=np.random.randint(0,SZ-16),np.random.randint(0,SZ-16)
    bw,bh=np.random.randint(16,80),np.random.randint(16,80)
    im3[y0:min(y0+bh,SZ),x0:min(x0+bw,SZ)]=np.random.randint(0,256,3)
test("rgb_blocks", im3)
test("rgb_random", np.random.randint(0,256,(SZ,SZ,3),dtype=np.uint8))

# Summary
print(f"\n{'='*76}")
tp=sum(r[1] for r in results); tt=sum(r[2] for r in results)
print(f"  FINAL: {W} WINS / {T} TIES / {L} LOSSES out of {len(results)}")
print(f"  Win rate: {W/len(results)*100:.0f}%")
print(f"  Aggregate: TC is {tp/tt:.3f}x smaller than PNG")

if L>0:
    print(f"\n  LOSSES:")
    for n,p,t,r,res in results:
        if res=="LOSS": print(f"    {n}: PNG={p} TC={t} gap={t-p}")
if T>0:
    print(f"\n  TIES:")
    for n,p,t,r,res in results:
        if res=="TIE": print(f"    {n}: PNG={p} TC={t}")
