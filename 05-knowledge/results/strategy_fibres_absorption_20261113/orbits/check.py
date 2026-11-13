"""Exact capped orbit and renewal audit for the inherited thirteen matrices."""
from pathlib import Path
from collections import deque
from itertools import product
import hashlib
import json
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
SOURCE = ROOT/'05-knowledge/results/kps_s18_ternary_chains_digraph_weights_20261109/out/generators_base2_13.txt'

def need(ok,message):
    if not ok:
        raise RuntimeError(message)

def apply(m,s):
    h,e,o=s
    return tuple(m[i]*h+m[i+1]*e+m[i+2]*o for i in (0,3,6))

def orbit(matrices,limit):
    seen={(1,0,1)}
    todo=deque(seen)
    depths={(1,0,1):0}
    while todo:
        s=todo.popleft()
        for m in matrices:
            t=apply(m,s)
            need(t[0]>=s[0] and t[0]%2==1,'monotone odd observer')
            if t[0]<=limit and t not in seen:
                seen.add(t);todo.append(t);depths[t]=depths[s]+1
    return seen,max(depths.values())

def prime(n):
    return n>=2 and all(n%d for d in range(2,int(n**0.5)+1))

def compile_word(word,specs):
    rows=[0]
    root=0
    permutation=[0]
    for letter in word:
        k,bits,w,newport=specs[letter]
        n=len(rows)
        rows.extend([0]*k)
        for i in range(k):
            for j in range(k):
                if bits[k*i+j]=='1': rows[n+i]|=1<<(n+j)
        for i in range(n):
            for j in range(k):
                if i==root and j==w: rows[n+j]|=1<<i
                else: rows[i]|=1<<(n+j)
        root=n+newport
        permutation.extend([n+1,n] if letter==0 else range(n,n+k))
    return rows,root,permutation

def split_join(word,specs):
    pieces=[]
    excursion=[]
    for letter in word:
        if letter==0:
            pieces.append(excursion+[1])
            excursion=[]
        else: excursion.append(letter)
    pieces.append(excursion)
    rows=[]
    root=0
    for piece in pieces:
        block,port,_=compile_word(piece,specs)
        offset=len(rows)
        later=((1<<len(block))-1)<<offset
        rows=[row|later for row in rows]
        rows.extend(row<<offset for row in block)
        root=offset+port
    return rows,root

def direct_state(rows,root):
    """Independent subset path count, including parity-weighted root paths."""
    n=len(rows)
    size=1<<n
    dp=[[0]*n for _ in range(size)]
    rooted=[[0]*n for _ in range(size)]
    totals=[0]*size
    totals[0]=1
    for mask in range(1,size):
        for v in range(n):
            if not (mask>>v)&1: continue
            old=mask^(1<<v)
            if old==0:
                dp[mask][v]=1
                rooted[mask][v]=int(v==root)
            else:
                dp[mask][v]=sum(dp[old][u] for u in range(n) if (rows[u]>>v)&1)
                rooted[mask][v]=sum(rooted[old][u] for u in range(n) if (rows[u]>>v)&1)
        totals[mask]=sum(dp[mask])
    even=odd=0
    for mask in range(1,size):
        term=sum(rooted[mask])*totals[(size-1)^mask]
        if mask.bit_count()%2: odd+=term
        else: even+=term
    return (totals[-1],even,odd)

def main():
    matrices=[]
    specs=[]
    for line in SOURCE.read_text(encoding='utf-8').splitlines():
        if line and not line.startswith('#'):
            matrices.append(tuple(map(int,line.split('|')[0].split())))
            k,bits,w,port=line.split('|')[1].split()
            specs.append((int(k),bits,int(w),int(port)))
    need(len(matrices)==13,'alphabet size')
    reset=(1,2,0,0,0,0,1,2,0)
    need(matrices[0]==reset and matrices[1][:3]==reset[:3],'reset and closing observer')
    for m in matrices:
        need(min(m)>=0 and m[0]>=1 and m[1]>=2,'positive monotone matrix')
    limit=10001
    full,full_depth=orbit(matrices,limit)
    free,free_depth=orbit(matrices[1:],limit)
    values={s[0] for s in full}
    tails={s[0] for s in free}
    returns={s[0]+2*s[1] for s in free if s[0]+2*s[1]<=limit}
    need(returns<=tails,'closing observer substitutes a nonreset letter')
    renewal=set(tails)
    for h in range(1,limit+1,2):
        if h in renewal:
            for f in returns:
                if f*h<=limit: renewal.add(f*h)
    need(renewal==values,'full orbit equals reset renewal closure')
    prime_lane=[h for h in range(3,limit+1,2) if prime(h) and h!=7]
    seven_lane=[7*p for p in range(3,limit//7+1,2) if prime(p) and p!=3]
    need(all((h in values)==(h in tails) for h in prime_lane+seven_lane),'atomic lane preservation')
    extra=sorted(values-tails)
    need(all(not prime(h) and h not in seven_lane for h in extra),'reset only adds noncritical composites')
    # Symbolic renewal identity on all13^5 labelled words, retaining no
    # scalar state deduplication. Exact integer multiplication is independent
    # of the capped orbit search.
    word_checks=0
    for length in range(6):
        for word in product(range(13),repeat=length):
            state=(1,0,1)
            scale=1
            excursion=(1,0,1)
            for letter in word:
                state=apply(matrices[letter],state)
                if letter==0:
                    scale*=excursion[0]+2*excursion[1]
                    excursion=(1,0,1)
                else:
                    excursion=apply(matrices[letter],excursion)
            need(state==tuple(scale*x for x in excursion),'full state renewal identity')
            word_checks+=1
    graph_checks=direct_checks=0
    for length in range(4):
        for word in product(range(13),repeat=length):
            rows,root,permutation=compile_word(word,specs)
            joined,joined_root=split_join(word,specs)
            old_to_new={v:i for i,v in enumerate(permutation)}
            permuted=[sum(1<<old_to_new[v] for v in range(len(rows)) if (rows[u]>>v)&1)
                      for u in permutation]
            need(permuted==joined and old_to_new[root]==joined_root,'actual rooted graph renewal identity')
            graph_checks+=1
            if length<=2:
                matrix_state=(1,0,1)
                for letter in word: matrix_state=apply(matrices[letter],matrix_state)
                need(direct_state(rows,root)==matrix_state,'independent actual tournament path DP')
                direct_checks+=1
    report=dict(status='PASS',scope='Exact all-word reduction plus complete finite scalar caps; no coverage theorem',
                cap=limit,full_states=len(full),reset_free_states=len(free),
                full_shortest_depth_max=full_depth,reset_free_shortest_depth_max=free_depth,
                full_values=len(values),reset_free_values=len(tails),reset_return_values=len(returns),
                reset_added_values=extra,
                full_missing_allowed_odds=[h for h in range(1,limit+1,2) if h not in values and h not in (7,21)],
                prime_lane_misses=[h for h in prime_lane if h not in tails],
                seven_prime_lane_misses=[h for h in seven_lane if h not in tails],
                labelled_word_identity_checks=word_checks,
                actual_rooted_graph_identity_checks=graph_checks,
                independent_direct_state_checks=direct_checks,
                source_pins=[dict(path=str(SOURCE.relative_to(ROOT)).replace('\\','/'),
                                  sha256_lf=hashlib.sha256(SOURCE.read_bytes().replace(b'\r\n',b'\n')).hexdigest())])
    suffix='.optimized' if sys.flags.optimize else ''
    (HERE/('receipt'+suffix+'.json')).write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in report.items() if k not in ('source_pins','reset_added_values')}))
    print('reset-added count',len(extra),'first values',extra[:30])

if __name__=='__main__': main()
