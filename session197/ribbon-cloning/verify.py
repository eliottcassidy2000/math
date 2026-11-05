"""Standalone graph-level checks for the literal cloning dichotomy."""
import hashlib, itertools, json
from pathlib import Path

def require(x, why):
    if not x: raise RuntimeError(why)

def tournament(n, mask):
    rows=[0]*n
    for k,(i,j) in enumerate(itertools.combinations(range(n),2)):
        if mask>>k&1: rows[j]|=1<<i
        else: rows[i]|=1<<j
    return rows

def clone(rows,v,m):
    old=[x for x in range(len(rows)) if x!=v]
    origins=old+[v]*m
    n=len(origins); out=[0]*n
    for i in range(n):
        for j in range(i+1,n):
            forward=(i>=len(old)) or bool(rows[origins[i]]>>origins[j]&1)
            out[i if forward else j]|=1<<(j if forward else i)
    return out

def paths(rows):
    n=len(rows)
    if not n: return 0
    dp=[[0]*n for _ in range(1<<n)]
    for i in range(n): dp[1<<i][i]=1
    full=(1<<n)-1
    for bits in range(1,full+1):
        for i in range(n):
            z=dp[bits][i]
            if not z: continue
            rest=rows[i]&(full^bits)
            while rest:
                b=rest&-rest;rest-=b
                dp[bits|b][b.bit_length()-1]+=z
    return sum(dp[full])

def triangle(rows,v):
    for x in range(len(rows)):
        if rows[v]>>x&1:
            for y in range(len(rows)):
                if rows[x]>>y&1 and rows[y]>>v&1: return x,y
    return None

def closure(rows):
    reach=[r|1<<i for i,r in enumerate(rows)]
    for k in range(len(rows)):
        for i in range(len(rows)):
            if reach[i]>>k&1: reach[i]|=reach[k]
    return reach

def span(rows,order):
    return max([0]+[j-i for i in range(len(order)) for j in range(i+1,len(order)) if rows[order[j]]>>order[i]&1])

def minspan(rows):
    return min(span(rows,p) for p in itertools.permutations(range(len(rows))))

def clock(m):
    # r=0, a=1, transitive B=2,...,m+1
    n=m+2;rows=[0]*n;rows[0]|=1<<1
    for b in range(2,n):
        rows[1]|=1<<b;rows[b]|=1
        for c in range(b+1,n):rows[b]|=1<<c
    return rows

def wordgraph(word):
    n=len(word);rows=[0]*n
    for j,c in enumerate(word):
        for i in range(j):
            back=i>=j-2 and (int(c)>>(i-j+2))&1
            rows[j if back else i]|=1<<(i if back else j)
    return rows

def main():
    totals=dict(tournaments=0,rooted_cases=0,clone_count_checks=0,induced_clocks=0,exact_clone_spans=0)
    digest=hashlib.sha256()
    for n in range(1,6):
        for mask in range(1<<(n*(n-1)//2)):
            rows=tournament(n,mask);reach=closure(rows);h=paths(rows)
            totals['tournaments']+=1
            for v in range(n):
                totals['rooted_cases']+=1
                tri=triangle(rows,v)
                nontrivial=any(x!=v and reach[v]>>x&1 and reach[x]>>v&1 for x in range(n))
                require(bool(tri)==nontrivial,('component',rows,v))
                for m in (1,2,3):
                    rr=clone(rows,v,m);hh=paths(rr)
                    totals['clone_count_checks']+=1
                    if tri is None:require(hh==h,('constant',rows,v,m))
                    else:
                        require(hh>=2**m+1,('growth',rows,v,m))
                        old=[x for x in range(n) if x!=v];x,y=tri
                        chosen=[old.index(x),old.index(y)]+list(range(n-1,n-1+m))
                        induced=[sum(((rr[a]>>b)&1)<<j for j,b in enumerate(chosen)) for a in chosen]
                        require(induced==clock(m),('induced clock',rows,v,m))
                        totals['induced_clocks']+=1
                    if n<=4:
                        bw=minspan(rr);totals['exact_clone_spans']+=1
                        if tri is not None:require(bw>=(m+3)//3,('span',rows,v,m,bw))
                    digest.update(f'{n},{mask},{v},{m},{hh};'.encode())
    clocks=[]
    for m in range(1,7):
        rr=clock(m);bw=minspan(rr);h=paths(rr)
        require(bw==(m+3)//3,('sharpness',m,bw))
        require(h==2**m+1,('clock H',m,h))
        clocks.append(dict(m=m,n=m+2,minimum_span=bw,H=h))
    witnesses=[]
    for word,target,cn,cspan in [('022333332',513,11,4),('003231331133133',65537,18,6)]:
        rows=wordgraph(word);h=paths(rows);reach=closure(rows)
        require(h==target,('witness H',word,h))
        require(span(rows,range(len(rows)))<=2,('witness span',word))
        require(all(r==(1<<len(rows))-1 for r in reach),('witness strong',word))
        require(2**(cn-2)+1==h and (cn+1)//3==cspan,('contrasting clock',word))
        witnesses.append(dict(word=word,rows=rows,n=len(rows),H=h,strong=True,displayed_span=span(rows,range(len(rows))),clock_order=cn,clock_minimum_span=cspan))
    print(json.dumps(dict(checks=totals,clone_data_sha256=digest.hexdigest(),sharp_clock_cases=clocks,witnesses=witnesses),sort_keys=True,indent=2))

if __name__=='__main__':main()
