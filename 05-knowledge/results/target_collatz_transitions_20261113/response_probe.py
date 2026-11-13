"""Exact target-equation probe using an inherited ordinary bank through eight."""
from pathlib import Path
from math import gcd
import gzip,json,subprocess

HERE=Path(__file__).resolve().parent
REF='06847e4c6e6709585bf1b4cfa6f84713081034ea'
BANK='05-knowledge/results/multicore_metagraph_20261113/rooted_growth.json.gz'

def need(ok,why):
    if not ok:raise RuntimeError(why)

def inventory(T):
    n=len(T);full=(1<<n)-1
    need(all(0<=row<=full and not(row>>i&1) for i,row in enumerate(T)), 'loop/range')
    need(all(((T[i]>>j)&1)+((T[j]>>i)&1)==1 for i in range(n) for j in range(i+1,n)), 'pair')
    dp=[[0]*n for _ in range(1<<n)]
    for i in range(n):dp[1<<i][i]=1
    for mask in range(1,1<<n):
        for v,value in enumerate(dp[mask]):
            rem=T[v]&~mask
            while value and rem:
                bit=rem&-rem;rem-=bit
                dp[mask|bit][bit.bit_length()-1]+=value
    H=[1]+[sum(row) for row in dp[1:]]
    return H,dp

def pair(T):
    H,_=inventory(T);full=len(H)-1
    return H[-1],sum(H[s]*H[full^s] for s in range(full+1))

def bank():
    raw=json.loads(gzip.decompress(subprocess.check_output(['git','show',REF+':'+BANK])))
    ordinary={}
    for r in raw['records']:
        if r['n']<=8:ordinary.setdefault((r['n'],r['ordinary_code']),tuple(r['rows']))
    need([sum(n==i for n,c in ordinary) for i in range(1,9)]==[1,1,2,4,12,56,456,6880], 'inherited bank universe')
    return ordinary

def main():
    profiles={}
    for (n,code),T in sorted(bank().items()):
        h,p=pair(T);profiles.setdefault((h,p),dict(n=n,rows=T,code=code,h=h,p=p))
    hits=[];small=[]
    data=list(profiles.values())
    for M in data:
        h,p=M['h'],M['p'];x=6*p+h
        for N in data:
            j,q=N['h'],N['p']
            if gcd(h,j)>1:continue
            e=x*j-6*h*q
            if abs(e)==1:
                H=12*p*q+6*h*q+2*p*j-h*j
                hit=dict(M=M,N=N,e=e,H=H,target=3*H+2*e)
                hits.append(hit)
                if h*j>1:small.append(hit)
    small.sort(key=lambda z:(z['M']['n']+z['N']['n'],z['H']))
    out=dict(status='FINITE-EXACT',source_ref=REF,scope='All response pairs represented in inherited ordinary types orders1..8; ordered module pairs tested for j(6p+h)-6hq=±1',profile_count=len(data),hit_count=len(hits),nontrivial_count=len(small),nontrivial=small[:40])
    (HERE/'response_probe.json').write_text(json.dumps(out,indent=2)+'\n')
    (HERE/'response_bank.json.gz').write_bytes(gzip.compress(json.dumps(data,separators=(',',':')).encode(),mtime=0))
    print(json.dumps({k:v for k,v in out.items() if k!='nontrivial'},indent=2))
    for x in small[:8]:print('H',x['H'],'target',x['target'],'e',x['e'],'M',(x['M']['n'],x['M']['h'],x['M']['p']),'N',(x['N']['n'],x['N']['h'],x['N']['p']))

if __name__=='__main__':main()
