"""
Diagnose: why do we lose on sparse 256x256?
Compare our zlib output vs PNG byte-for-byte.

kind-pasteur-2026-03-25-S2
"""
import sys, io, zlib, struct
import numpy as np
from PIL import Image

sys.stdout.reconfigure(line_buffering=True)

# Generate the sparse 1% image (256x256)
np.random.seed(42)
SZ = 256
img = np.zeros((SZ, SZ), dtype=np.uint8)
n = int(SZ*SZ*0.01)
img[np.random.randint(0,SZ,n), np.random.randint(0,SZ,n)] = np.random.randint(1,256,n).astype(np.uint8)

# 1. PNG reference size
pil = Image.fromarray(img, 'L')
buf = io.BytesIO()
pil.save(buf, format='PNG', optimize=True, compress_level=9)
png_total = buf.tell()

# Parse the PNG to find IDAT data
buf.seek(0)
sig = buf.read(8)  # PNG signature
idat_data = b''
while True:
    length_bytes = buf.read(4)
    if len(length_bytes) < 4:
        break
    length = struct.unpack('>I', length_bytes)[0]
    chunk_type = buf.read(4)
    chunk_data = buf.read(length)
    crc = buf.read(4)
    if chunk_type == b'IDAT':
        idat_data += chunk_data
    if chunk_type == b'IEND':
        break

print(f"PNG total: {png_total} bytes")
print(f"PNG IDAT compressed: {len(idat_data)} bytes")
print(f"PNG overhead (sig+IHDR+IDAT hdr+IEND): {png_total - len(idat_data)} bytes")

# Decompress IDAT to get the filtered data
filtered = zlib.decompress(idat_data)
print(f"PNG filtered data: {len(filtered)} bytes (= {SZ} rows * ({SZ}+1) bytes = {SZ*(SZ+1)})")

# Analyze the filter bytes
filter_bytes = [filtered[i*(SZ+1)] for i in range(SZ)]
from collections import Counter
fc = Counter(filter_bytes)
print(f"PNG filter distribution: {dict(sorted(fc.items()))}")

# 2. Re-compress the PNG filtered data ourselves
for strategy in ['default', 'filtered', 'huffman']:
    s = {'default': zlib.Z_DEFAULT_STRATEGY, 'filtered': zlib.Z_FILTERED,
         'huffman': zlib.Z_HUFFMAN_ONLY}[strategy]
    obj = zlib.compressobj(9, zlib.DEFLATED, 15, 9, s)
    c = obj.compress(filtered) + obj.flush()
    print(f"Our zlib ({strategy:>8}): {len(c)} bytes (PNG: {len(idat_data)}, diff: {len(c)-len(idat_data):+d})")

# 3. Try raw data (no filter bytes)
raw = img.tobytes()
for strategy in ['default', 'filtered', 'huffman']:
    s = {'default': zlib.Z_DEFAULT_STRATEGY, 'filtered': zlib.Z_FILTERED,
         'huffman': zlib.Z_HUFFMAN_ONLY}[strategy]
    obj = zlib.compressobj(9, zlib.DEFLATED, 15, 9, s)
    c = obj.compress(raw) + obj.flush()
    print(f"Raw zlib ({strategy:>8}): {len(c)} bytes")

# 4. Try our filter selection
def pix(a,x,d=0):
    if a is not None and 0<=x<len(a): return int(a[x])
    return d
def _pa(a,b,c):
    p=a+b-c; pa,pb,pc=abs(p-a),abs(p-b),abs(p-c)
    if pa<=pb and pa<=pc: return a
    if pb<=pc: return b
    return c

FILTERS=[
    lambda r,p,x,w:0,
    lambda r,p,x,w:pix(r,x-1),
    lambda r,p,x,w:pix(p,x),
    lambda r,p,x,w:(pix(r,x-1)+pix(p,x))>>1,
    lambda r,p,x,w:_pa(pix(r,x-1),pix(p,x),pix(p,x-1)),
]

def frow(row,prev,fid,w):
    pf=FILTERS[fid]; o=np.empty(w,dtype=np.uint8)
    for x in range(w): o[x]=(int(row[x])-pf(row,prev,x,w))&0xFF
    return o

# Our PNG-mimic filtered data
our_filtered = bytearray()
our_filter_choices = []
for y in range(SZ):
    row = img[y]
    prev = img[y-1] if y > 0 else None
    best = 1<<30; bf = 0; br = None
    for f in range(5):
        r = frow(row, prev, f, SZ)
        s = r.astype(np.int16); s[s>128] -= 256
        a2 = int(np.sum(np.abs(s)))
        if a2 < best: best = a2; bf = f; br = r
    our_filter_choices.append(bf)
    our_filtered.append(bf)
    our_filtered.extend(br)

ofc = Counter(our_filter_choices)
print(f"\nOur filter distribution: {dict(sorted(ofc.items()))}")

# Compare filter choices
same = sum(1 for a, b in zip(filter_bytes, our_filter_choices) if a == b)
print(f"Filter agreement: {same}/{SZ} ({same/SZ*100:.1f}%)")

# Compress our filtered data
for strategy in ['default', 'filtered']:
    s = {'default': zlib.Z_DEFAULT_STRATEGY, 'filtered': zlib.Z_FILTERED}[strategy]
    obj = zlib.compressobj(9, zlib.DEFLATED, 15, 9, s)
    c = obj.compress(bytes(our_filtered)) + obj.flush()
    print(f"Our filtered ({strategy:>8}): {len(c)} bytes")

# Check if filtered data is identical
if bytes(our_filtered) == filtered:
    print("\nFiltered data is IDENTICAL to PNG!")
else:
    diffs = sum(1 for a, b in zip(our_filtered, filtered) if a != b)
    print(f"\nFiltered data differs in {diffs} bytes")

# 5. Compute minimum possible: PNG overhead vs our overhead
png_overhead = png_total - len(idat_data)  # PNG headers
our_overhead = 4 + 1  # minimal header (w:2 h:2) + strategy byte
print(f"\nPNG overhead: {png_overhead} bytes")
print(f"Our overhead: {our_overhead} bytes")
print(f"Overhead advantage: {png_overhead - our_overhead} bytes in our favor")
print(f"Conclusion: if we match PNG's compressed data, we save {png_overhead - our_overhead} bytes")

# 6. Our final compressed size using PNG-exact approach + our header
obj = zlib.compressobj(9, zlib.DEFLATED, 15, 9, zlib.Z_FILTERED)
our_idat = obj.compress(filtered) + obj.flush()  # Use PNG's exact filtered data
our_total = our_overhead + len(our_idat)
print(f"\nWith PNG's exact filtered data + our header: {our_total} bytes")
print(f"vs PNG: {png_total} bytes")
print(f"Result: {'WIN' if our_total < png_total else 'LOSS'} by {abs(png_total - our_total)} bytes")
