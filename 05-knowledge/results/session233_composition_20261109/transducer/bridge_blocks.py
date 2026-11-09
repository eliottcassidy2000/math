"""Exact odd-cycle packing transfer for transitive blocks joined by spine bridges."""
import itertools,json,sys
from functools import lru_cache
@lru_cache(None)
def matrix(m,x=2):
    # states closed, even unfinished path, odd unfinished path; column convention
    A=[[0]*3 for _ in range(3)]; A[0][0]=1
    subsets=[0,1] if m==1 else [2**(m-2)]*2
    for q in (0,1):
        A[1+q][0]=subsets[q]
        A[0][1+q]=x*subsets[1-q]
    through=([0,1] if m==1 else [1,0] if m==2 else [(3**(m-2)+(-1)**(m-2))//2,(3**(m-2)-(-1)**(m-2))//2])
    for p in (0,1):
        for q in (0,1): A[1+(p+q)%2][1+p]+=through[q]
        if m>=2:
            for r in (0,1):
                q=1-p
                count=(3**(m-2)+(1 if (q-1)%2==0 else -1)+(1 if (r-1)%2==0 else -1)+(1 if (q+r+m-4)%2==0 else -1))//4
                A[1+r][1+p]+=x*int(count)
    return tuple(map(tuple,A))
def apply(A,v): return tuple(sum(a*b for a,b in zip(row,v)) for row in A)
def count(blocks,x=2):
    v=(1,0,0)
    for m in blocks: v=apply(matrix(m,x),v)
    return v[0]
def tournament(blocks):
    n=sum(blocks); owner=[]
    for i,m in enumerate(blocks): owner += [i]*m
    out=[0]*n
    for i in range(n):
        for j in range(i+1,n):
            a,b=(i,j) if owner[i]==owner[j] or j==i+1 else (j,i)
            out[a]|=1<<b
    return out
def paths(out):
    n=len(out); dp=[[0]*n for _ in range(1<<n)]
    for i in range(n): dp[1<<i][i]=1
    for S in range(1,1<<n):
        for i,v in enumerate(dp[S]):
            if not v: continue
            rem=out[i]&~S
            while rem:
                b=rem&-rem;rem-=b;dp[S|b][b.bit_length()-1]+=v
    return sum(dp[-1])
def compositions(n):
    if n==0: yield ();return
    for m in range(1,n+1):
        for tail in compositions(n-m):yield (m,)+tail

def main():
    if "--fail" in sys.argv: raise RuntimeError("deliberate sentinel")
    checked=0
    for n in range(1,10):
        for c in compositions(n):
            if count(c)!=paths(tournament(c)): raise RuntimeError(("path mismatch",c))
            checked+=1
    levels=[{(1,0,0):()}]; values={1}; rows=[]
    for n in range(1,23):
        states={}
        for m in range(1,n+1):
            for v,w in levels[n-m].items():
                new=apply(matrix(m),v)
                if new[0]<=10000:states.setdefault(new,w+(m,))
        levels.append(states);values.update(v[0] for v in states)
        rows.append({'n':n,'states':len(states),'attained_up_to_10000':len(values),'first_gaps':[i for i in range(1,501,2) if i not in values][:30]})
        if len(states)>600000:break
    result={'independent_path_checks_all_compositions_through_n9':checked,'levels':rows,'first_100_gaps_below_10000':[i for i in range(1,10001,2) if i not in values][:100], 'matrices_1_to_5':[matrix(m) for m in range(1,6)]}
    print(json.dumps(result,indent=2))
if __name__=='__main__':main()
