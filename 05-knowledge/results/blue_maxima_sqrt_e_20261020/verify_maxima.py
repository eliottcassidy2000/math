#!/usr/bin/env python3
"""Independent orbit, moment, witness, completion, and search-record checks."""
from pathlib import Path
from itertools import combinations, permutations, product
from fractions import Fraction
from math import factorial, comb
import argparse, hashlib, json
HERE=Path(__file__).resolve().parent

def require(ok,msg):
    if not ok: raise AssertionError(msg)

def adjacency(k,A,b,matching=None):
    n=2*k+(matching is not None);a=[[0]*n for _ in range(n)]
    def edge(i,j,v): a[i][j]=int(v);a[j][i]=1-int(v)
    for z,(i,j) in enumerate(combinations(range(k),2)):
        v=A>>z&1;c=b>>z&1
        edge(2*i,2*j,v);edge(2*j+1,2*i+1,v)
        edge(2*i,2*j+1,c);edge(2*j,2*i+1,c)
    if matching is not None:
        for i in range(k):
            edge(2*i,2*i+1,matching>>i&1)
            edge(2*i,2*k,1);edge(2*k,2*i+1,1)
    return a

def h_dp(k,A,b):
    a=adjacency(k,A,b);dp=[[0]*(2*k) for _ in range(1<<k)]
    for i in range(2*k):dp[1<<(i//2)][i]=1
    for S in range(1,1<<k):
        for i,v in enumerate(dp[S]):
            if v:
                for j in range(2*k):
                    if not(S>>(j//2)&1) and a[i][j]: dp[S|1<<(j//2)][j]+=v
    return sum(dp[-1][::2])

def second_moment(k):
    ans=Fraction(0)
    for j in range(k//2+1):
        r=k-2*j+1
        c=1 if j==0 else sum(2**s*comb(r,s)*comb(j-1,s-1) for s in range(1,min(j,r)+1))
        ans+=Fraction(factorial(k-2*j)*c,4**j)
    return factorial(k)*ans

def transform(k,A,b,p):
    pairs=list(combinations(range(k),2));lookup={xy:z for z,xy in enumerate(pairs)}
    aa=bb=0
    for z,(i,j) in enumerate(pairs):
        u,v=p[i],p[j];q=lookup[tuple(sorted((u,v)))];val=(A>>z&1)^(u>v)
        aa|=val<<q;bb|=(b>>z&1)<<q
    return aa,bb

def witness(k,A,b,M):
    a=adjacency(k,A,b);auts=[]
    for p in permutations(range(k)):
        if all(a[2*i][2*j]==a[2*p[i]][2*p[j]] and a[2*i][2*j+1]==a[2*p[i]][2*p[j]+1] for i in range(k) for j in range(k) if i!=j):auts.append(p)
    require(len(auts)==5,'frame automorphism order')
    identity=tuple(range(k));generator=next(p for p in auts if p!=identity)
    power=identity;generated=set()
    for _ in range(5):generated.add(power);power=tuple(generator[power[i]] for i in range(k))
    require(generated==set(auts),'cyclic order-five group')
    cyc=[i for i in range(k) if generator[i]!=i];fixed=[i for i in range(k) if generator[i]==i]
    require(len(cyc)==5,'five moved indices')
    matching=1<<cyc[0]
    classes={}
    for vec in range(1<<k):
        orbit={sum(((vec>>i)&1)<<p[i] for i in range(k)) for p in auts}
        rep=min(orbit)
        stabilizer=sum(all((vec>>i&1)==(vec>>p[i]&1) for i in range(k)) for p in auts)
        require(stabilizer*len(orbit)==5,'matching orbit-stabilizer')
        classes[rep]=M//stabilizer
    fixed_count=sum(v==M//5 for v in classes.values());max_count=sum(v==M for v in classes.values())
    require(max_count==(2**k-2**(k-4))//5,'maximum class formula')
    require(fixed_count==2**(k-4),'fixed class formula')
    complete=adjacency(k,A,b,matching)
    rho=[i^1 for i in range(2*k)]+[2*k]
    require(all(complete[i][j]+complete[j][i]==1 for i in range(2*k+1) for j in range(i+1,2*k+1)),'full tournament')
    require(all(complete[i][j]==complete[rho[j]][rho[i]] for i in range(2*k+1) for j in range(2*k+1)),'anti-involution')
    regular_matching=sum((k-1-sum(a[2*i]))<<i for i in range(k))
    require(all(k-1-sum(a[2*i]) in (0,1) for i in range(k)),'regular completability')
    require(all(sum(row)==k for row in adjacency(k,A,b,regular_matching)),'regular completion degrees')
    require(h_dp(k,A,b)==M,'independent witness DP')
    count=0;digest=hashlib.sha256()
    for p in permutations(range(k)):
        for bits in range(1<<(k-1)):
            prefix=[2*p[t]+(bits>>t&1) for t in range(k)]
            if all(a[u][v] for u,v in zip(prefix,prefix[1:])):
                word=prefix+[2*k]+[v^1 for v in reversed(prefix)]
                require(len(set(word))==2*k+1,'symmetric Hamiltonian word')
                require(all(complete[u][v] for u,v in zip(word,word[1:])),'full word edges')
                count+=1;digest.update(bytes(word))
    require(count==M,'direct prefix enumeration')
    return {'k':k,'order':2*k+1,'A':A,'b':b,'h':M,'matching':matching,
            'adjacency_row_masks':[sum(bit<<j for j,bit in enumerate(row)) for row in complete],
            'frame_automorphisms':[list(p) for p in auts],'group_generator':list(generator),
            'moved_indices':cyc,'fixed_indices':fixed,'matching_stabilizer_order':1,'tournament_regular_completion_matching':regular_matching,
            'maximum_blue_classes':max_count,'other_classes_in_same_hollow_type':fixed_count,
            'other_blue_mass':M//5,'maximum_classes_per_triangle_parity_shore':max_count//2,
            'enumerated_symmetric_words':count,'word_stream_sha256':digest.hexdigest()}

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--skip-orbit-coverage',action='store_true');args=parser.parse_args()
    result={'finite_exact':[],'heuristic_records':[]}
    for k in (6,7):
        x=json.loads((HERE/f'exact_k{k}.json').read_text());N=1<<(k*(k-1));M=x.get('maximum')
        require(x['full_labelled_frames']==N,'frame count')
        require(x['sum_h']==N*factorial(k),'mean')
        require(x['sum_h_squared']==N*second_moment(k),'second moment')
        hist={int(v):n for v,n in x['histogram'].items()}
        require(sum(hist.values())==N and sum(v*n for v,n in hist.items())==x['sum_h'],'histogram count/mean')
        require(sum(v*v*n for v,n in hist.items())==x['sum_h_squared'],'histogram second moment')
        require(max(hist)==M and all(v>0 and v%2 for v in hist),'positive odd values and maximum')
        reps=x['tournament_representatives'];ps=list(permutations(range(k)))
        text_reps=[int(v) for v in (HERE/f'A_reps_k{k}.txt').read_text().split()]
        require(text_reps==[r['A'] for r in reps],'replay representative list')
        covered=set()
        if not args.skip_orbit_coverage:
            for r in reps:
                orb={transform(k,r['A'],0,p)[0] for p in ps}
                require(len(orb)==r['orbit_size'],'A orbit size')
                require(len(orb)*r['aut_order']==factorial(k),'A stabilizer')
                require(not covered.intersection(orb),'disjoint A orbits');covered.update(orb)
            require(len(covered)==1<<(k*(k-1)//2),'all A orientations covered')
        require(len(x['maximizing_reduced_frames'])==1,'unique reduced maximizing frame')
        A,b=x['maximizing_reduced_frames'][0];w=witness(k,A,b,M)
        require(hist[M]==factorial(k)//5,'labelled maximizing frame orbit')
        replay=json.loads((HERE/f'replay_k{k}.json').read_text())
        require(replay['all_equal'],'independent replay status')
        require(replay.get('pointwise_checks',replay.get('pointwise_comparisons'))==x['reduced_frames_checked'],'independent pointwise count')
        (HERE/f'witness_n{2*k+1}.json').write_text(json.dumps(w,indent=2)+'\n')
        result['finite_exact'].append({'k':k,'M':M,'labelled_maximizers':hist[M],
            'A_types':len(reps),'A_orientations_independently_covered':len(covered),
            'reduced_frames_pointwise_replayed':x['reduced_frames_checked'],
            'mean':factorial(k),'variance':str(second_moment(k)-factorial(k)**2),
            'max_blue_classes':w['maximum_blue_classes']})
    for k,name in [(8,'local_k8_deep.json'),(9,'local_k9.json'),(10,'local_k10.json')]:
        x=json.loads((HERE/name).read_text());h=h_dp(k,x['A'],x['b']);require(h==x['max_h_found'],'heuristic witness DP')
        result['heuristic_records'].append({'k':k,'certified_lower_bound':h,'A':x['A'],'b':x['b'],'global_optimality_claimed':False})
    (HERE/'verification.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
