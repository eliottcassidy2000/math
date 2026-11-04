#!/usr/bin/env python3
"""Per-anchor endpoint criterion and explicit even quotient separator."""
from pathlib import Path
from tempfile import TemporaryDirectory
import hashlib,json,subprocess
R=(598,508,344,865,40,69,656,317,113,438)
ROOT=Path(__file__).resolve().parent

def need(ok,why):
    if not ok:raise RuntimeError(why)

def tournament(n,code):
    a=[0]*n;k=0
    for i in range(n):
        for j in range(i+1,n):
            if code>>k&1:a[i]|=1<<j
            else:a[j]|=1<<i
            k+=1
    return tuple(a)

def opposite(a):return tuple(((1<<len(a))-1)^r^(1<<i) for i,r in enumerate(a))
def ends(a):
    n=len(a);d=[0]*(1<<n)
    for i in range(n):d[1<<i]=1<<i
    for s in range(1,1<<n):
        m=d[s]
        while m:
            bit=m&-m;m-=bit;u=bit.bit_length()-1;unused=a[u]&~s
            while unused:
                v=unused&-unused;unused-=v;d[s|v]^=v
    return d

def integer_ends(a):
    """Independent integer predecessor DP, reducing only completed entries."""
    n=len(a);d=[[0]*n for _ in range(1<<n)];incoming=[0]*n
    for i,row in enumerate(a):
        for j in range(n):
            if row>>j&1:incoming[j]|=1<<i
        d[1<<i][i]=1
    for mask in range(1,1<<n):
        if mask&(mask-1)==0:continue
        for last in range(n):
            if mask>>last&1:
                rest=mask^(1<<last);preds=incoming[last]&rest;value=0
                while preds:
                    bit=preds&-preds;preds-=bit
                    value+=d[rest][bit.bit_length()-1]
                d[mask][last]=value
    return [sum((value&1)<<i for i,value in enumerate(row)) for row in d]

def detector(a):
    n=len(a);E=ends(a);S=ends(opposite(a));b=0;h=[]
    for mask in range(1,1<<n):
        k=mask.bit_count()
        if k>=2 and n-k>=2 and (n-k)%2==0:b^=E[mask]^S[mask]
    for j in range(n):
        col=0
        for mask in range(1,1<<n):
            if not mask>>j&1 and mask.bit_count()%2==0:col^=E[mask]^S[mask]
        h.append(col)
    return b,h,S,E

def direct_correction(S,E):
    full=len(S)-1;b=0;pairs=0
    for A in range(1,full+1):
        if A.bit_count()<2:continue
        rest=full^A;B=rest
        while B:
            if B.bit_count()>=2:
                if E[A].bit_count()%2:b^=S[B]^E[B]
                pairs+=1
            B=(B-1)&rest
    return b,pairs

def main():
    sources=['PROOF.md','verify.py','run_dp.cpp','literal_hostile.cpp'];pins={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in sources};small={};nchecked=0;support_pairs=0
    for n in range(3,6):
        for code in range(1<<(n*(n-1)//2)):
            T=tournament(n,code);b,h,S,E=detector(T);direct,pairs=direct_correction(S,E);need(direct==b,'support count correction');support_pairs+=pairs;nchecked+=1
            need(b.bit_count()%2==0,'total correction parity')
            if n%2:need(b!=(1<<n)-1,'odd order retains an anchor')
            if n==4:need(b==15,'four-name trace blind boundary')
            if n==5:need(b==sum((r.bit_count()%2)<<i for i,r in enumerate(T)),'five-name outdegree formula')
            small[n,code]=(b,h)
    b,h,S,E=detector(R);need(b==763 and (1023^b)==260,'ten-name anchor mask');need(sum(v.bit_count() for v in h)==44 and h[0]>>3&1,'ten-name shared cubic')
    need(S==integer_ends(opposite(R)) and E==integer_ends(R),'independent integer endpoint subsets')
    with TemporaryDirectory() as td:
        td=Path(td);outputs={}
        for name in ['run_dp','literal_hostile']:
            binary=td/name;subprocess.run(['g++','-std=c++17','-O2','-Wall','-Wextra','-Werror',str(ROOT/(name+'.cpp')),'-o',str(binary)],check=True)
            out=subprocess.check_output([str(binary)]);outputs[name]=json.loads(out);(ROOT/(name+'.json')).write_bytes(out)
    checks=0
    for row in outputs['run_dp']['small_trace_controls']:
        n,code=row[:2];b0,h0=small[n,code]
        need(row[2:2+n]==[4*(1^((b0>>j)&1)) for j in range(n)],'literal per-anchor trace')
        need(row[2+n:]==[0 if i==j else 4*((h0[j]>>i)&1) for j in range(n) for i in range(n)],'literal shared cubic');checks+=1
    need(outputs['run_dp']['coefficient_mod8']==4,'eight-query coefficient')
    literal=outputs['literal_hostile'];need((literal['source_mod8'],literal['target_mod8'])==(3,7),'literal physical separator')
    witness=next(v for v in outputs['run_dp']['actual_pair_queries'] if v['a']==3 and v['c']==0);need((witness['source_mod8'],witness['target_mod8'])==(3,7),'independent weighted separator')
    need(pins=={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in sources},'source frozen')
    return dict(status='PASS',support_formula_quotients=nchecked,disjoint_support_pair_checks=support_pairs,complete_small_literal_trace_controls=checks,hostile_rows=R,b_mask=b,zero_anchors=[j for j in range(10) if not b>>j&1],shared_cubic_columns=h,nonzero_shared_cubics=sum(v.bit_count() for v in h),literal_twenty_name_witness=literal,source_sha256=pins,data_sha256={n+'.json':hashlib.sha256((ROOT/(n+'.json')).read_bytes()).hexdigest() for n in ['run_dp','literal_hostile']},scope='Sufficient per-anchor and shared-cubic criteria, including an explicit even ten-name quotient; no universal even-order reconstruction claim.')
if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
