"""Find actual one-edge interpolation witnesses, retaining Boolean endpoints."""
from collections import Counter
from itertools import combinations
from pathlib import Path
import json
HERE=Path(__file__).resolve().parent

def rows(n,mask):
    out=[0]*n
    for k,(i,j) in enumerate(combinations(range(n),2)):
        if mask>>k&1:out[i]|=1<<j
        else:out[j]|=1<<i
    return out

def h_count(out):
    n=len(out);size=1<<n
    dp=[[0]*n for _ in range(size)]
    for i in range(n):dp[1<<i][i]=1
    for mask in range(1,size):
        free=(size-1)^mask
        for last,value in enumerate(dp[mask]):
            if not value:continue
            nxt=out[last]&free
            while nxt:
                bit=nxt&-nxt;nxt-=bit
                dp[mask|bit][bit.bit_length()-1]+=value
    return sum(dp[-1])

def main():
    result=[]
    for n,target in ((5,7),(6,21)):
        N=n*(n-1)//2
        counts=[h_count(rows(n,x)) for x in range(1<<N)]
        pair=next((x,x^(1<<e),e) for x,a in enumerate(counts) if a==target-2
                  for e in range(N) if counts[x^(1<<e)]==target+2)
        x,y,e=pair
        row={'n':n,'target':target,'masks':[x,y],'edge_index':e,
             'edge':list(combinations(range(n),2))[e],'counts':[counts[x],counts[y]],
             'out_rows':[rows(n,x),rows(n,y)],
             'histogram':dict(sorted(Counter(counts).items()))}
        result.append(row)
        print(json.dumps({k:v for k,v in row.items() if k!='histogram'}),flush=True)
    (HERE/'scout.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n',encoding='utf-8')

if __name__=='__main__':main()
