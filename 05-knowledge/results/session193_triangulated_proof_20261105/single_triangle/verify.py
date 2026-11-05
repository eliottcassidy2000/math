#!/usr/bin/env python3
"""Exact controls for the rooted response envelope. No third-party dependencies."""
from collections import Counter
from functools import lru_cache
from hashlib import sha256
from itertools import combinations, permutations
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent

def require(ok,msg):
    if not ok: raise RuntimeError(msg)

def B(n):
    if n==1:return 1
    if n<3:raise ValueError('No strong tournament order two.')
    q,r=divmod(n,3)
    return (3*5**(q-1),5**q,9*5**(q-1))[r]

def C(n):return 1 if n in (1,3) else 2*B(n-1)-B(n)
def L(n):return B(n)+C(n)

def valid(rows,p):return all(rows[i]>>j&1 for i,j in zip(p,p[1:]))

def strong(rows,S):
    if not S:return False
    b=S&-S
    for reverse in (False,True):
        Q=b
        while True:
            R=Q
            for i in range(len(rows)):
                if reverse:
                    if (S>>i&1) and rows[i]&Q:R|=1<<i
                elif Q>>i&1:R|=rows[i]&S
            if R==Q:break
            Q=R
        if Q!=S:return False
    return True

def near_transitive(rows):
    n=len(rows)
    for i in range(n):
        for j in range(i+1,n):
            r=list(rows);r[i]^=1<<j;r[j]^=1<<i
            if len(set(x.bit_count() for x in r))==n:return True
    return False

def profile(rows):
    """Independent exact subset DP for all start/end responses."""
    n=len(rows);N=1<<n;F=N-1
    incoming=[sum(1<<j for j in range(n) if rows[j]>>i&1) for i in range(n)]
    end=[[0]*n for _ in range(N)];start=[[0]*n for _ in range(N)];h=[0]*N;h[0]=1
    for S in range(1,N):
        bits=S
        while bits:
            b=bits&-bits;bits-=b;i=b.bit_length()-1;A=S^b
            if not A:end[S][i]=start[S][i]=1
            else:
                q=A&incoming[i]
                while q:
                    d=q&-q;q-=d;end[S][i]+=end[A][d.bit_length()-1]
                q=A&rows[i]
                while q:
                    d=q&-q;q-=d;start[S][i]+=start[A][d.bit_length()-1]
        h[S]=sum(end[S])
    R=[sum(start[S][a]*h[F^S] for S in range(N) if S>>a&1) for a in range(n)]
    return h,end,start,R

def labelled(n):
    pairs=list(combinations(range(n),2))
    for mask in range(1<<len(pairs)):
        rows=[0]*n
        for k,(i,j) in enumerate(pairs):
            if mask>>k&1:rows[i]|=1<<j
            else:rows[j]|=1<<i
        yield tuple(rows)

def literal_paths(rows):
    n=len(rows);answer={}
    for S in range(1<<n):
        answer[S]=tuple(p for p in permutations(i for i in range(n) if S>>i&1) if valid(rows,p))
    return answer

def literal_E(paths,S,a):
    result=[]
    for Qset in range(S+1):
        if Qset&~S or not(Qset>>a&1):continue
        for Q in paths[Qset]:
            if Q[0]!=a:continue
            for P in paths[S^Qset]:result.append(('R',P,Q))
    for P in paths[S]:
        if P[-1]==a:result.append(('e',P))
    return result

def check_injection(rows,A,Bset,a,v,paths):
    @lru_cache(None)
    def merge(P,Q):
        if not P:return Q
        if not Q:return P
        # Exhaustively choose a legal order-preserving interleaving.
        def rec(i,j,last):
            if i==len(P) and j==len(Q):return ()
            for useP in (True,False):
                if useP:
                    if i==len(P):continue
                    nxt=P[i]
                else:
                    if j==len(Q):continue
                    nxt=Q[j]
                if last!=-1 and not(rows[last]>>nxt&1):continue
                tail=rec(i+useP,j+(not useP),nxt)
                if tail is not None:return (nxt,)+tail
            return None
        out=rec(0,0,-1)
        require(out is not None,'Path merging failure')
        return out
    seen=set();input_count=0
    for obj in literal_E(paths,A,a):
        for word in paths[Bset]:
            i=word.index(v);Cm,Cp=word[:i],word[i+1:]
            if obj[0]=='R':
                _,P,Q=obj
                if v in P:
                    i=P.index(v);out=('R',merge(P[:i],Cm)+(v,)+merge(P[i+1:],Cp),Q)
                else:
                    i=Q.index(v);out=('R',merge(P,Cm),Q[:i]+(v,)+merge(Q[i+1:],Cp))
            else:
                W=obj[1];i=W.index(v)
                W=merge(W[:i],Cm)+(v,)+merge(W[i+1:],Cp)
                if W[-1]==a:out=('e',W)
                else:
                    i=W.index(a);out=('R',W[:i],W[i:])
            require(out not in seen,'Injection collision')
            seen.add(out);input_count+=1
            words=out[1:]
            require(all(valid(rows,p) for p in words),'Invalid output path')
            require(sorted(x for p in words for x in p)==list(range(len(rows))),'Bad output cover')
            require((out[0]=='e' and out[1][-1]==a) or (out[0]=='R' and out[2][0]==a),'Wrong root endpoint')
    return input_count

def main():
    literal_counts={};injection_pairs=0;injection_objects=0;small_min={}
    for n in range(1,6):
        count=0
        for rows in labelled(n):
            count+=1;F=(1<<n)-1;paths=literal_paths(rows);h,end,start,R=profile(rows)
            for S in range(1<<n):require(h[S]==len(paths[S]),'DP/literal H mismatch')
            for a in range(n):
                require(len(literal_E(paths,F,a))==R[a]+end[F][a],'DP/literal E mismatch')
            if strong(rows,F) and n!=2:
                for a in range(n):
                    require(R[a]>=B(n),'R floor failed')
                    require(R[a]+end[F][a]>=L(n),'E floor failed')
                    triangles=sum(1 for i,j in combinations([i for i in range(n) if i!=a],2) if ((rows[a]>>i&1) and (rows[i]>>j&1) and (rows[j]>>a&1)) or ((rows[a]>>j&1) and (rows[j]>>i&1) and (rows[i]>>a&1)))
                    if triangles==1 and n>=4:small_min[n]=min(small_min.get(n,10**10),h[F]-h[F^(1<<a)])
            # All arbitrary gluings with A,B each >=3; includes nonstrong factors.
            if n==5:
                for v in range(n):
                    for Arest in combinations([i for i in range(n) if i!=v],2):
                        A=(1<<v)|sum(1<<i for i in Arest);BB=(F^A)|(1<<v)
                        for a in Arest:
                            injection_pairs+=1;injection_objects+=check_injection(rows,A,BB,a,v,paths)
        literal_counts[n]=count
    require(small_min=={4:2,5:4},'Wrong small single-triangle minima')

    # Inherited complete unlabelled bank, independently evaluated here.
    file=HERE/'inherited_tournaments_through8.txt';raw=file.read_bytes();tokens=list(map(int,raw.split()));K=tokens[0];pos=1
    bank_counts=Counter();strong_counts=Counter();decomp_cases=Counter();root_counts=Counter();mins={};strong_masks_cache={}
    for identity in range(K):
        n=tokens[pos];pos+=1;rows=tuple(tokens[pos:pos+n]);pos+=n;F=(1<<n)-1;bank_counts[n]+=1
        require(all(not(rows[i]>>i&1) for i in range(n)),'Loop in bank')
        require(all(((rows[i]>>j&1)+(rows[j]>>i&1))==1 for i,j in combinations(range(n),2)),'Malformed tournament')
        if not strong(rows,F):continue
        strong_counts[n]+=1;h,end,start,R=profile(rows)
        st={S for S in range(1<<n) if S.bit_count()>=3 and strong(rows,S)}
        nt=near_transitive(rows) if n>=4 else True
        for a in range(n):
            root_counts[n]+=1;E=R[a]+end[F][a]
            require(R[a]>=B(n),'Bank R floor failed')
            require(E>=L(n),'Bank E floor failed')
            require(0<=end[F][a]<=h[F]<=R[a],'Weight interval failed')
            mins[n]=min(mins.get(n,10**30),E)
            if n<5:continue
            if nt:decomp_cases['nearly_transitive']+=1;continue
            if not any(not(S>>a&1) for S in st):decomp_cases['transitive_root_deletion']+=1;continue
            found=None
            for A in st:
                if A==F:continue
                for v in range(n):
                    if v==a or not(A>>v&1):continue
                    BB=(F^A)|(1<<v)
                    if BB in st:found=(A,BB,v);break
                if found:break
            require(found is not None,f'Root-avoiding decomposition failure: {rows}, root={a}')
            decomp_cases['root_avoiding_split']+=1
    require(pos==len(tokens),'Trailing bank data')
    require(dict(bank_counts)=={1:1,2:1,3:2,4:4,5:12,6:56,7:456,8:6880},'Unexpected inherited universe')
    require(all(mins[n]==L(n) for n in mins),'Bound not attained in bank')
    # All side orders through 120, rational weights, and residue formulas.
    for n in [1]+list(range(3,121)):
        require(0<C(n)<=B(n),'Invalid sharp endpoint coordinate')
    for p in range(3,81):
        for q in range(3,81):
            require(L(p)*B(q)>=L(p+q-1),'Induction arithmetic failed')
    result={'status':'PASS: finite controls, not the all-order proof','labelled_literal_universe':literal_counts,'tagged_merge_input_families':injection_pairs,'tagged_merge_objects':injection_objects,'small_single_triangle_minima':small_min,'inherited_bank_sha256':sha256(raw).hexdigest(),'inherited_unlabelled_counts':dict(bank_counts),'strong_counts':dict(strong_counts),'strong_root_counts':dict(root_counts),'decomposition_controls':dict(decomp_cases),'minimum_E':mins}
    print(json.dumps(result,indent=2,sort_keys=True))

if __name__=='__main__':main()
