#!/usr/bin/env python3
"""Literal integer subset-DP controls for the lossless six-name modulo-four bank."""
from pathlib import Path
from itertools import combinations
import hashlib,json,gzip
ROOT=Path(__file__).resolve().parent
N=32768
B=[s for s in range(4096) if (s&63).bit_count()<=2 and (s>>6).bit_count()<=2]

def need(ok,why):
    if not ok:raise RuntimeError(why)

def count(code,a,c):
    rows=[0]*8
    for k,(i,j) in enumerate(combinations(range(6),2)):
        if code>>k&1:rows[i]|=1<<j
        else:rows[j]|=1<<i
    rows[6]|=1<<7
    for i in range(6):
        if a>>i&1:rows[i]|=1<<6
        else:rows[6]|=1<<i
        if c>>i&1:rows[i]|=1<<7
        else:rows[7]|=1<<i
    dp=[[0]*8 for _ in range(256)]
    for i in range(8):dp[1<<i][i]=1
    for mask in range(1,256):
        for u,z in enumerate(dp[mask]):
            if not z:continue
            rem=rows[u]&~mask
            while rem:
                b=rem&-rem;rem^=b;dp[mask|b][b.bit_length()-1]+=z
    return sum(dp[255])

def main():
    p=ROOT/'mod4_six_vectors.bin'
    raw=p.read_bytes() if p.exists() else gzip.decompress((ROOT/'mod4_six_vectors.bin.gz').read_bytes())
    need(len(raw)==N*len(B),'bank length')
    need(all(c<4 for c in raw),'canonical residues')
    records=[]
    codes=sorted({0,32767,21845,10922}|{(7919*j+17)&32767 for j in range(64)})
    for code in codes:
        for a,c in [(0,0),(63,63),((code*23+7)&63,(code*53+19)&63),(21,42)]:
            value=count(code,a,c);m=a|(c<<6)
            got=sum(raw[b*N+code] for b,s in enumerate(B) if s&m==s)%4
            need(got==value%4,'literal integer query')
            records.append([code,a,c,value])
    return dict(status='PASS',queries=len(records),literal_queries=records,full_bank_sha256=hashlib.sha256(raw).hexdigest(),basis_coordinates=len(B),tournaments=N)
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,separators=(',',':')))
