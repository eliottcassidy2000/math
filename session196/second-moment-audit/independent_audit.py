#!/usr/bin/env python3
"""Independent exact controls: no import or execution of the audited verifier.
Uses Hamilton-word completion accumulation rather than endpoint DP for the census,
bitstring runs and odd-block recurrences for W, polynomial division for asymptotics,
and direct averages of compatible completions to audit small witness traces.
"""
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from itertools import permutations
from math import comb, factorial
from pathlib import Path
import hashlib, json

ROOT = Path(__file__).resolve().parent

def check(ok, message):
    if not ok:
        raise RuntimeError(message)


def formula(n):
    return factorial(n) + sum(2**r * comb(j-1,r-1) * comb(n-2*j,r) * factorial(n-2*j)
        for j in range(1,(n-1)//2+1) for r in range(1,min(j,n-2*j)+1)) if n else 1


def all_completions(n):
    edges = [(i,j) for i in range(n) for j in range(i+1,n)]
    index = {e:b for b,e in enumerate(edges)}
    full = (1<<len(edges))-1
    counts = [0]*(full+1)
    for word in permutations(range(n)):
        present = absent = 0
        for a,b in zip(word,word[1:]):
            bit = 1 << index[tuple(sorted((a,b)))]
            if a < b: present |= bit
            else: absent |= bit
        free = full ^ (present | absent)
        sub = free
        while True:
            counts[present | sub] += 1
            if not sub: break
            sub = (sub-1)&free
    return edges, counts


def orbit_audit(n, edges, hs):
    unseen = set(range(len(hs)))
    masses = [F(0),F(0)]
    reps = []
    pms = list(permutations(range(n)))
    idx = {e:b for b,e in enumerate(edges)}
    while unseen:
        mask = min(unseen)
        orbit = []
        for p in pms:
            image = 0
            for b,(i,j) in enumerate(edges):
                x,y = p[i],p[j]
                sign = (mask >> idx[tuple(sorted((x,y)))]) & 1
                if sign == (x<y): image |= 1<<b
            orbit.append(image)
        aut = orbit.count(mask)
        distinct = set(orbit)
        check(len(distinct)*aut == factorial(n), ('orbit stabilizer',n,mask))
        check(all(hs[x] == hs[mask] for x in distinct), ('H invariant',n,mask))
        check(hs[mask]%aut == 0, ('free path action',n,mask))
        unseen.difference_update(distinct)
        h = hs[mask]
        masses[0] += F(h,aut)
        masses[1] += F(h*h,aut)
        reps.append({'mask':mask,'aut':aut,'H':h})
    check(len(reps)==[1,1,1,2,4,12,56][n], ('A000568',n))
    check(masses[0] == 2**((n-1)*(n-2)//2), ('orbit first',n))
    check(masses[1]*4**(n-1) == 2**len(edges)*formula(n), ('orbit second',n))
    return {'classes':len(reps),'first_mass':str(masses[0]),'second_mass':str(masses[1])}


def bitstring_run_mass(n):
    total = 0
    for bits in range(1<<(n-1)):
        runs = [len(s) for s in format(bits, f'0{n-1}b').split('0') if s]
        if any(k%2 for k in runs): continue
        total += 2**len(runs) * factorial(n-sum(runs))
    return total


def block_recurrence(N):
    # B[n][m] is coefficient of z^n in (z+2z^3+2z^5+...)^m.
    B = [[0]*(N+1) for _ in range(N+1)]
    B[0][0] = 1
    for n in range(1,N+1):
        for m in range(1,n+1):
            B[n][m] = B[n-1][m-1] + 2*sum(B[n-k][m-1] for k in range(3,n+1,2))
    return [sum(B[n][m]*factorial(m) for m in range(n+1)) for n in range(N+1)]


def mul(a,b):
    ans=[F(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b): ans[i+j] += x*y
    return ans


def expansion(K):
    ans=[F(0)]*(K+1)
    for j in range(1,K+1):
        for r in range(1,j+1):
            shift=2*j-r
            if shift>K: continue
            num=[F(2**r*comb(j-1,r-1),factorial(r))]
            den=[F(1)]
            for i in range(r): num=mul(num,[1,-2*j-i])
            for i in range(2*j): den=mul(den,[1,-i])
            q=[]
            for k in range(K-shift+1):
                q.append((num[k] if k<len(num) else 0)-sum(den[i]*q[k-i] for i in range(1,min(k,len(den)-1)+1)))
            for k,x in enumerate(q): ans[k+shift]+=x
    return [str(x) for x in ans]


def direct_weighted(weights):
    n=len(weights)
    @lru_cache(None)
    def extend(remaining,last):
        if not remaining: return 1
        return sum(weights[last][v]*extend(remaining^(1<<v),v)
                   for v in range(n) if remaining>>v&1)
    return sum(extend(((1<<n)-1)^(1<<v),v) for v in range(n))


def check_trace(wit, census):
    n=wit['n']
    weights=[[int(i!=j) for j in range(n)] for i in range(n)]
    for i in range(n-1): weights[i][i+1],weights[i+1][i]=2,0
    check(direct_weighted(weights)==formula(n), ('initial mean',n))
    if n<=6:
        edges,hs=census[n]
        candidates=[m for m in range(len(hs)) if all(m>>edges.index((i,i+1))&1 for i in range(n-1))]
    for step in wit['decision_trace']:
        i,j=step['pair']
        check(direct_weighted(weights)==step['before_scaled_mean'], ('trace before',n,i,j))
        for forwards,key in [(True,'forward_scaled_mean'),(False,'backward_scaled_mean')]:
            weights[i][j],weights[j][i]=(2,0) if forwards else (0,2)
            check(direct_weighted(weights)==step[key], ('trace branch',n,i,j,forwards))
            if n<=6:
                b=edges.index((i,j))
                sub=[m for m in candidates if bool(m>>b&1)==forwards]
                check(F(sum(hs[m] for m in sub),len(sub))*2**(n-1)==step[key], ('literal completion average',n,i,j,forwards))
        chosen=step['chosen_forward']
        check(chosen==(step['forward_scaled_mean']>=step['backward_scaled_mean']), ('tie rule',n,i,j))
        weights[i][j],weights[j][i]=(2,0) if chosen else (0,2)
        if n<=6: candidates=[m for m in candidates if bool(m>>edges.index((i,j))&1)==chosen]
    rows=wit['outgoing_masks']
    check(all(weights[i][j]==2*int(bool(rows[i]>>j&1)) for i in range(n) for j in range(n)), ('output adjacency',n))
    check(direct_weighted(weights)==wit['H']*2**(n-1), ('output count',n))
    check(wit['H']*2**(n-1)>=formula(n), ('witness guarantee',n))
    return {'n':n,'H':wit['H'],'decisions':len(wit['decision_trace']),'completion_averages_checked':n<=6}


def main():
    report={}
    census={}
    census_report=[]
    for n in range(1,7):
        edges,hs=all_completions(n)
        census[n]=edges,hs
        check(sum(hs)*2**(n-1)==len(hs)*factorial(n), ('first moment',n))
        check(sum(x*x for x in hs)*4**(n-1)==len(hs)*factorial(n)*formula(n), ('second moment',n))
        census_report.append({'n':n,'sum_H':sum(hs),'sum_H2':sum(x*x for x in hs),'orbit':orbit_audit(n,edges,hs)})
    report['independent_completion_census_n1_to_n6']=census_report
    for n in range(1,19): check(bitstring_run_mass(n)==formula(n), ('bitstring runs',n))
    report['bitstring_run_enumeration_through_n']=18
    B=block_recurrence(80)
    check(B==[formula(n) for n in range(81)], 'odd-block recurrence')
    report['odd_block_recurrence_through_n']=80
    coeff=expansion(10)
    check(coeff[:8]==['0','2','0','-14/3','-58/3','-916/15','-6986/45','-13010/63'],'series coefficients')
    report['variance_coefficients_x0_to_x10']=coeff
    count=0
    for n in range(12,121):
        levels=[F(0)]+[F(sum(2**r*comb(j-1,r-1)*comb(n-2*j,r)*factorial(n-2*j)
                     for r in range(1,min(j,n-2*j)+1)),factorial(n)) for j in range(1,(n-1)//2+1)]
        tail=F(0)
        for J in range(len(levels)-1,-1,-1):
            bound=F(2,3)*F(6,n)**(J+1)/(1-F(6,n))
            check(tail<=bound, ('tail',n,J))
            tail+=levels[J]
            count+=1
    report['all_J_tail_checks_n12_to_n120']=count
    original=json.loads((ROOT/'source_snapshot/receipt.json').read_text())
    report['independent_trace_checks']=[check_trace(w,census) for w in original['deterministic_witnesses']]
    report['status']='PASS'
    text=json.dumps(report,indent=2,sort_keys=True)+'\n'
    (ROOT/'independent_receipt.json').write_text(text)
    print(text,end='')

if __name__=='__main__': main()
