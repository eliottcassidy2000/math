"""Independent F extraction and ordinal-block factorization controls."""
from pathlib import Path
from itertools import permutations,product
import json,gzip
import numpy as np
from frame_probe import ordinary_bank,need,hp,valid
from feedback_probe import all_profiles

HERE=Path(__file__).resolve().parent

def literal_F(T,r):
    n=len(T);count=0
    for w in permutations(range(n)):
        k=w.index(r)+1
        if all(T[a]>>b&1 for a,b in zip(w[:k],w[1:k])) and all(T[a]>>b&1 for a,b in zip(w[k:],w[k+1:])):count+=1
    return count

def join(*blocks):
    out=[];offset=0;N=sum(map(len,blocks))
    for T in blocks:
        future=((1<<N)-1)^((1<<(offset+len(T)))-1)
        out.extend((row<<offset)|future for row in T);offset+=len(T)
    return tuple(out)

def covers2(T):
    n=len(T);total=0
    for S in range(1<<n):
        vs=[v for v in range(n) if S>>v&1];ws=[v for v in range(n) if not(S>>v&1)]
        def sub(V):return tuple(sum(((T[v]>>w)&1)<<j for j,w in enumerate(V)) for v in V)
        total+=(hp(sub(vs)) if vs else 1)*(hp(sub(ws)) if ws else 1)
    return total

def triangle_context(T,r):
    n=len(T);allold=(1<<n)-1
    return tuple(row|((1<<n) if v==r else 0) for v,row in enumerate(T))+(allold^(1<<r)|(1<<(n+1)),allold)

def run():
    bank=ordinary_bank();literal=0
    for (n,code),T in sorted(bank.items()):
        if n>5:continue
        Z=all_profiles(np.array(T,dtype=np.int64))
        for r in range(n):
            need(literal_F(T,r)==int(Z[r,4]),'literal broken joins')
            child=triangle_context(T,r)
            need(valid(child) and hp(child)==int(Z[r,0]+2*Z[r,4]),'literal triangle context law')
            literal+=1
    blocks=[(),(0,),(2,0),(2,4,1),(14,12,8,0)]
    centers=[(T,r) for (n,code),T in sorted(bank.items()) if 2<=n<=4 for r in range(n)]
    factors=0
    for L,R in product(blocks,repeat=2):
        for C,r in centers:
            T=join(L,C,R);need(valid(T),'ordinal product tournament')
            actual=tuple(map(int,all_profiles(np.array(T,dtype=np.int64))[len(L)+r]))
            z=tuple(map(int,all_profiles(np.array(C,dtype=np.int64))[r]))
            hl=hp(L) if L else 1;hr=hp(R) if R else 1
            need(actual[0]==hl*z[0]*hr and actual[4]==covers2(L)*z[4]*hr,'exact preceding/following block factorization')
            factors+=1
    data=json.loads(gzip.decompress((HERE/'surplus.json.gz').read_bytes()))
    missing={str(d):[r['h'] for r in data['head'] if r['h']>=3 and str(d) not in r['witnesses']] for d in (1,2)}
    degree_pairs={}
    for r in data['small_surplus_occurrences']:
        key=str((r['delta'],r['indegree']));degree_pairs[key]=degree_pairs.get(key,0)+1
    result=dict(status='FINITE-EXACT independent controls; all-order small-head conclusion also uses the proved padding lemma and CITED Busch floor',literal_F_controls=literal,literal_triangle_context_controls=literal,ordinal_factorization_controls=factors,small_surplus_degree_pairs=degree_pairs,missing_small_surplus_by_h=missing)
    suffix='' if __debug__ else '.optimized'
    (HERE/('surplus_audit'+suffix+'.json')).write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print('LITERAL_F',literal,'BLOCK_CONTROLS',factors,'DEGREES',degree_pairs,'MISSING',missing,flush=True)

if __name__=='__main__':run()
