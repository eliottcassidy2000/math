"""Fast exact F1 branch, then optional finite certified minimal fallback."""
from math import isqrt
from pathlib import Path
import argparse
import json
from small_arm import g
from join_loads import analyze, analyze_bank, explicit_features
from selector import solve, profile

def affine_inverse(N,A,R,eps0,min_index=1):
    # Solve A*g(v)+R*(v+1)-eps0*[v even]=N.
    # v=Lq+r, L power of two; the appended zeros produce a coefficient of q.
    L=2
    while L*L*R<N:L*=2
    k=L.bit_length()-1
    for r in range(L):
        m=r;s=0
        for _ in range(k):s=(s<<1)+1-(m&1);m>>=1
        B=R*L+A*s
        C=N-R*(r+1)-A*(g(r)-1)+eps0*(r%2==0)
        lo=max(0,-(-(C-2*A)//(B+A)))
        hi=(C-A)//B
        for q in range(lo,hi+1):
            v=L*q+r
            if v+1>=min_index and R*(v+1)+A*g(v)-eps0*(v%2==0)==N:
                return v

def find_F1(N,b_limit=None):
    M=isqrt(N)
    if b_limit is not None:M=min(M,b_limit)
    for b in range(1,M+1):
        v=affine_inverse(N,g(b-1),b,b%2,b)
        if v is not None:return b,v+1

def archived(N):
    here=Path(__file__).resolve().parent/'targeted'
    for a in (32,35):
        lo,hi=2**a+1,2**a+1999
        if not lo<=N<=hi:continue
        if a==35 and N==34359739333:
            return analyze_bank((2,4,9,3,15),1012954554)
        p=here/f'band{a}.tsv'
        if not p.exists():return None
        for line in p.read_text().splitlines():
            fields=line.split()
            if len(fields)==3 and fields[0]==str(N):
                w=analyze(int(fields[1]),int(fields[2]))
                assert w['H']==N
                return w

if __name__=='__main__':
    cli=argparse.ArgumentParser()
    cli.add_argument('N',type=int)
    cli.add_argument('--no-archive',action='store_true',help='do not use verified frontier witness receipts')
    cli.add_argument('--bank',help='optional JSON array of small rooted outgoing-mask arrays (root last)')
    cli.add_argument('--b-limit',type=int,default=128,help='bound the fast family search; 0 means entire F1')
    cli.add_argument('--exhaustive',action='store_true',help='certify global presence or absence with finite fallback')
    cli.add_argument('--minimum',action='store_true',help='certify least order; may be impractical')
    cli.add_argument('--all-classes',action='store_true',help='enumerate every first-birth class (implies --minimum)')
    a=cli.parse_args();N=a.N
    if N<1 or N%2==0:cli.error('expected an odd positive integer')
    if N in {7,21}:print(json.dumps(dict(target=N,status='known hole')))
    else:
        witness=None if a.no_archive else archived(N)
        if witness is None:
            pair=find_F1(N,a.b_limit or None)
            witness=analyze(*pair) if pair is not None else None
        if witness is None:
            bank=[(2,4,9,3,15)] # certified near-frontier F2 repair
            if a.bank:
                with open(a.bank) as f:bank.extend(json.load(f))
            for rows in bank:
                feat=explicit_features(rows)
                v=affine_inverse(N,feat['P'][0],feat['C'][0],feat['E'][0])
                if v is not None:witness=analyze_bank(rows,v+1);break
        if witness is not None and not (a.minimum or a.all_classes):
            print(json.dumps(dict(status='structured witness, order not certified minimum', **witness)))
        elif a.exhaustive or a.minimum or a.all_classes:
            W=solve(N,a.all_classes,witness['order'] if witness else None)
            print(json.dumps(dict(target=N,status=('minimum proved' if W else 'globally absent'),
                                  witnesses=[profile(w) for w in W])))
        else:
            print(json.dumps(dict(target=N,status='not found in bounded construction search; no global inference')))
