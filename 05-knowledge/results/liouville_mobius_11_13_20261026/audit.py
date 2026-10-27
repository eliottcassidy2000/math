"""Independent arithmetic sieve and adjacency/permutation audit; no producer import."""
from collections import Counter
from functools import lru_cache
from itertools import permutations, product
from math import prod
from pathlib import Path
import json
import sys

sys.stdout.reconfigure(encoding='utf-8',newline='\n')
ROOT=Path(__file__).resolve().parent

def need(test,message):
    if not test: raise RuntimeError(message)

def automorphisms(rows):
    n=len(rows)
    arc=[[bool(rows[i]&(1<<j)) for j in range(n)] for i in range(n)]
    colour=[sum(row) for row in arc]
    while True:
        palette=sorted(set(colour))
        signatures=[(colour[i],tuple(sum(arc[i][j] and colour[j]==c for j in range(n))
                                      for c in palette)) for i in range(n)]
        keys={s:k for k,s in enumerate(sorted(set(signatures)))}
        new=[keys[s] for s in signatures]
        if len(set(new))==len(set(colour)): colour=new; break
        colour=new
    cells=[[i for i in range(n) if colour[i]==c] for c in sorted(set(colour))]
    tested=0; count=0
    for targets in product(*(permutations(cell) for cell in cells)):
        p=dict(zip(sum(cells,[]),sum((list(t) for t in targets),[])))
        tested+=1
        count+=all(arc[i][j]==arc[p[i]][p[j]] for i in range(n) for j in range(i))
    return {'order':count,'colour_cells':list(map(len,cells)),'permutations_tested':tested}

def path_count(rows):
    n=len(rows)
    predecessors=[tuple(j for j in range(n) if rows[j]&(1<<i)) for i in range(n)]
    @lru_cache(None)
    def ending(mask,last):
        remaining=mask^(1<<last)
        if not remaining: return 1
        return sum(ending(remaining,j) for j in predecessors[last] if remaining&(1<<j))
    return sum(ending((1<<n)-1,i) for i in range(n))

def main():
    bound=1000; mu=[1]*(bound+1); lam=[1]*(bound+1); prime=[True]*(bound+1)
    prime[0]=prime[1]=False
    for p in range(2,bound+1):
        if not prime[p]: continue
        for j in range(p,bound+1,p): mu[j]*=-1
        for j in range(p*p,bound+1,p*p): mu[j]=0
        power=p
        while power<=bound:
            for j in range(power,bound+1,power): lam[j]*=-1
            power*=p
        for j in range(p*p,bound+1,p): prime[j]=False
    div_mu=[0]*(bound+1); div_lam=[0]*(bound+1)
    for d in range(1,bound+1):
        for n in range(d,bound+1,d): div_mu[n]+=mu[d]; div_lam[n]+=lam[d]
    need(div_mu[1:]==[1]+[0]*(bound-1),'independent Mobius convolution')
    squares={d*d for d in range(1,32)}
    need(div_lam[1:]==[int(n in squares) for n in range(1,bound+1)],'independent Liouville convolution')
    source=json.loads((ROOT/'normal.json').read_text(encoding='utf-8'))
    audited=[]
    for item,expected in zip(source['mobius_Aut_order15'],[3,1,1,3]):
        rows=item['rows']; n=len(rows)
        need(n==15 and all(not (rows[i]>>i&1) for i in range(n)),'diagonal and size')
        need(all(((rows[i]>>j)&1)+((rows[j]>>i)&1)==1 for i in range(n) for j in range(i)),'tournament type')
        a=automorphisms(rows); H=path_count(rows)
        need(a['order']==expected==item['Aut'],'independent automorphism count')
        need(H==152377413==item['H'],'independent Hamilton count')
        audited.append({'Aut':a,'H':H})
    print(json.dumps({'status':'FINITE-EXACT independent implementation, same-agent proof review',
                      'sieve_bound':bound,'order15':audited},indent=2,sort_keys=True))

if __name__=='__main__': main()
