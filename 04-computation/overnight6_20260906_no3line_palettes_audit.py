#!/usr/bin/env python3
"""Independent palette audit by row-neighbor board enumeration and graph BFS.

No producer imports, palette-board generator, or full shore-label census.
Small exact probability kernels use partial injections into abstract cycles.
"""
from collections import Counter, defaultdict
from fractions import Fraction as Q
from functools import lru_cache
from hashlib import sha256
from itertools import combinations, combinations_with_replacement, permutations, product
from math import comb, factorial
from pathlib import Path
import sys

sys.stdout.reconfigure(newline="\n")
GATES=0
BASE=Path(__file__).resolve().parent


def need(ok,label):
    global GATES
    GATES+=1
    if not ok:
        raise RuntimeError(label)


def falling(n,r):
    return factorial(n)//factorial(n-r)


def events(n):
    return [frozenset(((a,x),(b,y),(c,z)))
            for a,b,c in combinations(range(n),3)
            for x,y,z in permutations(range(n),3)
            if (b-a)*(z-y)==(c-b)*(y-x)]


def components(board,n):
    adjacency=[[] for _ in range(2*n)]
    for r,c in board:
        adjacency[r].append(n+c)
        adjacency[n+c].append(r)
    need(all(len(v)==2 for v in adjacency),'literal row/column degree-two object')
    unseen=set(range(2*n))
    out=[]
    while unseen:
        seed=min(unseen)
        seen={seed}
        stack=[seed]
        while stack:
            x=stack.pop()
            for y in adjacency[x]:
                if y not in seen:
                    seen.add(y)
                    stack.append(y)
        unseen-=seen
        rows=tuple(sorted(x for x in seen if x<n))
        cols=tuple(sorted(x-n for x in seen if x>=n))
        need(len(rows)==len(cols)>=2,'simple bipartite cycle component')
        out.append((rows,cols))
    return tuple(out)


def board_census(n):
    parts=defaultdict(list)
    rowchoices=tuple(combinations(range(n),2))
    for choice in product(rowchoices,repeat=n):
        counts=Counter(c for pair in choice for c in pair)
        if any(counts[c]!=2 for c in range(n)):
            continue
        board=frozenset((r,c) for r,pair in enumerate(choice) for c in pair)
        cc=components(board,n)
        key=tuple(sorted(len(rows) for rows,cols in cc))
        parts[key].append((board,cc))
    return parts


def central_stats(values):
    mu=sum(map(Q,values),Q(0))/len(values)
    return mu,sum((Q(x)-mu)**2 for x in values)/len(values),sum((Q(x)-mu)**3 for x in values)/len(values)


def normalized_edges(F):
    F=tuple(F)
    rr={x:i for i,x in enumerate(sorted({r for r,c in F}))}
    cc={x:i for i,x in enumerate(sorted({c for r,c in F}))}
    return tuple(sorted((rr[r],cc[c]) for r,c in F))


@lru_cache(None)
def injection_probability(k,F):
    if not F:
        return Q(1)
    nr=1+max(r for r,c in F)
    nc=1+max(c for r,c in F)
    if max(nr,nc)>k:
        return Q(0)
    cycle={(i,i) for i in range(k)}|{(i,(i+1)%k) for i in range(k)}
    count=0
    for rr in permutations(range(k),nr):
        for cc in permutations(range(k),nc):
            count+=all((rr[r],cc[c]) in cycle for r,c in F)
    return Q(count,falling(k,nr)*falling(k,nc))


def local_kernel(F,rowsets,colsets):
    ro={x:i for i,rs in enumerate(rowsets) for x in rs}
    co={x:i for i,cs in enumerate(colsets) for x in cs}
    if any(ro[r]!=co[c] for r,c in F):
        return Q(0)
    out=Q(1)
    for i,rs in enumerate(rowsets):
        block=normalized_edges((r,c) for r,c in F if ro[r]==i)
        out*=injection_probability(len(rs),block)
    return out


def palette_head(n):
    bank=board_census(n)
    expected={(2,2):18,(4,):72} if n==4 else {(2,3):600,(5,):1440}
    need({key:len(rows) for key,rows in bank.items()}==expected,('full binary-board cycle-type universe',n))
    ev=events(n)
    need(len(ev)==(12 if n==4 else 52),'complete nonaxis line-event universe')
    target=(2,n-2)
    palettes=defaultdict(list)
    for board,cc in bank[target]:
        # For two equal C4s, either actual component may be the first named one.
        for rs,cs in cc:
            if len(rs)==2:
                palettes[rs,cs].append(board)
    need(len(palettes)==comb(n,2)**2 and {len(bs) for bs in palettes.values()}==({1} if n==4 else {6}),
         ('all named palettes and conditional board multiplicities',n))
    counts={board:sum(T<=board for T in ev) for board,cc in bank[target]}
    vals=[counts[b] for boards in palettes.values() for b in boards]
    mus=[]
    variances=[]
    thirds=[]
    for (rs,cs),boards in palettes.items():
        rowsets=(rs,tuple(x for x in range(n) if x not in rs))
        colsets=(cs,tuple(x for x in range(n) if x not in cs))
        X=[counts[b] for b in boards]
        a,b,c=central_stats(X)
        mus.append(a)
        variances.append(b)
        thirds.append(c)
        predicted=Q(0)
        for T in ev:
            prob=local_kernel(T,rowsets,colsets)
            need(prob==Q(sum(T<=board for board in boards),len(boards)),('every single-event probability',n,rs,cs,T))
            predicted+=prob
        need(predicted==a,('all palette means',n,rs,cs))
    avg=lambda xs:sum(xs,Q(0))/len(xs)
    overall=central_stats(vals)
    terms=(avg(thirds),central_stats(mus)[2],3*(avg([x*y for x,y in zip(mus,variances)])-avg(mus)*avg(variances)))
    need(sum(terms)==overall[2] and avg(variances)+central_stats(mus)[1]==overall[1],('total cumulance and total variance',n))
    if n==4:
        need(Counter(vals)=={0:18,2:8,4:4,6:4,8:2},'named n4 histogram')
        need(overall==(Q(2),Q(56,9),Q(16)) and terms==(0,16,0),'n4 outer-only ordinary cumulant')
        factorial_third=overall[2]-3*overall[1]+2*overall[0]
        need(factorial_third==Q(4,3) and 2*overall[0]==4,'ordinary/factorial conditional mismatch retained')
    else:
        need(overall==(Q(13,3),Q(1769,225),Q(26159,675)),'n5 unconditional exact moments')
        need(terms==(Q(-871,675),Q(6467,360),Q(7949,360)),'n5 all three total-cumulance terms')
        # Complete lists of up to three events, including repeated events,
        # at two named palettes; their complete union controls presence.
        for key in (((0,1),(0,1)),((0,2),(1,4))):
            boards=palettes[key]
            rowsets=(key[0],tuple(x for x in range(n) if x not in key[0]))
            colsets=(key[1],tuple(x for x in range(n) if x not in key[1]))
            for r in (1,2,3):
                for family in combinations_with_replacement(ev,r):
                    union=frozenset().union(*family)
                    need(local_kernel(union,rowsets,colsets)==Q(sum(union<=b for b in boards),len(boards)),
                         ('complete repeated-event union kernel',key,r))
    aut=32 if n==4 else 24
    need(len(bank[target])*aut==factorial(n)**2,'fixed-skeleton label multiplicity matches board law')
    print('HEAD',n,'all_types',{key:len(rows) for key,rows in bank.items()},'palettes',len(palettes))
    print('  HIST',sorted(Counter(vals).items()))
    print('  CENTRAL',overall,'TOTAL_CUMULANCE',terms)


def conditional_six():
    palettes=((0,1,2),(3,4,5))
    boards=[]
    choices=[tuple(combinations(palettes[0 if r<3 else 1],2)) for r in range(6)]
    for rows in product(*choices):
        counts=Counter(c for pair in rows for c in pair)
        if all(counts[c]==2 for c in range(6)):
            boards.append(frozenset((r,c) for r,pair in enumerate(rows) for c in pair))
    need(len(boards)==36,'n6 independent conditional row-neighbor universe')
    A=frozenset((i,i) for i in (0,1,2))
    B=frozenset((i,i) for i in (3,4,5))
    C=frozenset((i,i) for i in (1,2,3))
    def moment(*family):
        union=frozenset().union(*family)
        actual=Q(sum(union<=board for board in boards),36)
        need(actual==local_kernel(union,palettes,palettes),'n6 connected-union injection kernel')
        return actual
    def third(a,b,c):
        return moment(a,b,c)-moment(a,b)*moment(c)-moment(a,c)*moment(b)-moment(b,c)*moment(a)+2*moment(a)*moment(b)*moment(c)
    probs=tuple(moment(*f) for f in ((A,),(B,),(C,),(A,B),(A,C),(B,C),(A,B,C)))
    need(probs==(Q(1,3),Q(1,3),Q(1,3),Q(1,9),Q(2,9),Q(1,6),Q(1,9)),'n6 exact seven probabilities')
    need(third(A,A,B)==0 and third(A,B,C)==Q(1,54),'conditional disconnected and connected actual cumulants')
    print('CONDITIONAL_N6',probs,'disconnected0 connected1/54')


def paired_partitions(labels):
    if not labels:
        yield ()
        return
    first=labels[0]
    for j in range(1,len(labels)):
        other=labels[j]
        for rest in paired_partitions(labels[1:j]+labels[j+1:]):
            yield ((first,other),)+rest


def det(a,b,c):
    return (b[0]-a[0])*(c[1]-a[1])==(c[0]-a[0])*(b[1]-a[1])


def all_C4_six():
    n=6
    pairings=tuple(paired_partitions(tuple(range(n))))
    need(len(pairings)==15,'unordered row-pair palette count')
    boards=set()
    ev=events(n)
    for rr in pairings:
        for rawcc in pairings:
            for cc in permutations(rawcc):
                blocks=[frozenset(product(r,c)) for r,c in zip(rr,cc)]
                board=frozenset().union(*blocks)
                need(board not in boards,'unique paired-label representation')
                boards.add(board)
                pairterm=0
                for i,(r,c) in enumerate(zip(rr,cc)):
                    diagonals=(((r[0],c[0]),(r[1],c[1])),((r[0],c[1]),(r[1],c[0])))
                    for j in range(3):
                        if i!=j:
                            pairterm+=sum(det(a,b,p) for a,b in diagonals for p in blocks[j])
                tripleterm=sum(det(*points) for points in product(*blocks))
                need(sum(T<=board for T in ev)==pairterm+tripleterm,'all1350 three-C4 geometric decompositions')
    need(len(boards)==1350==factorial(n)**2//(2**n*factorial(n//2)),'all-C4 board cardinality')
    center=(1,1)
    need(det((0,0),(2,2),center) and det((0,2),(2,0),center),'two diagonals retain two distinct triples through one outside point')
    print('ALL_C4_N6 boards1350 unique representations; all two-block/three-block decompositions PASS')


def cylinder_gram():
    n=3
    full=tuple(permutations(range(n)))
    partial=[frozenset(zip(dom,im)) for r in range(n+1) for dom in combinations(range(n),r) for im in permutations(range(n),r)]
    need(len(partial)==34,'complete one-shore partial injection basis')
    rows={a:tuple(int(all(p[i]==x for i,x in a)) for p in full) for a in partial}
    for a in partial:
        need(Q(sum(rows[a]),6)==Q(1,falling(n,len(a))),'uniform cylinder expectation')
        for b in partial:
            union=a|b
            valid=len({i for i,x in union})==len(union)==len({x for i,x in union})
            expected=Q(1,falling(n,len(union))) if valid else Q(0)
            need(Q(sum(x*y for x,y in zip(rows[a],rows[b])),6)==expected,'exact cylinder Gram/product equality')
    for p in full:
        need(rows[frozenset(enumerate(p))].count(1)==1,'full assignment orthogonal idempotent')
    print('CYLINDER_GRAM 34x34 exactly A*Atranspose/6; two-shore positivity follows by tensor product')


def matching_and_type_controls():
    for k in range(2,10):
        L=2*k
        for s in range(1,4):
            actual=sum(all((j+1)%L not in chosen for j in chosen)
                       for chosen in combinations(range(L),s))
            expected=(2*k if s==1 else k*(2*k-3) if s==2 else
                      2*k*(k-2)*(2*k-5)//3 if k>=3 else 0)
            need(actual==expected,('complete small matching formula control',k,s))
            if s<=k:
                probability=Q(actual*factorial(s),falling(k,s)**2)
                closed=(Q(2,k) if s==1 else Q(2*(2*k-3),k*(k-1)**2) if s==2 else
                        Q(4*(2*k-5),k*(k-1)**2*(k-2)))
                need(probability==closed,('specified matching edges have s-factorial and no shore swaps',k,s))
    need(injection_probability(2,((0,0),(1,1)))==1,'non-induced matching inside K22')
    need(injection_probability(3,((0,0),(0,1),(0,2)))==0,'degree-three union impossible')
    need(injection_probability(3,((0,0),(0,1),(1,0),(1,1)))==0,'closed C4 cannot embed in C6')
    need(injection_probability(3,((0,0),(1,1),(2,2)))==Q(1,3),'three specified matching edges normalization')
    print('MATCHING_KERNEL k2..9 s1..3; non-induced, degree-three, closed-cycle and s-factorial controls')


def main():
    producer=BASE/'overnight6_20260906_no3line_palettes.py'
    digest=sha256(producer.read_bytes()).hexdigest()
    need(digest=='cc0fc2774d8475766da8f90ab4a5bb12d070170bfc5ba2b4d3bbf0c74e05f80f','frozen producer bytes')
    print('PRODUCER_SHA256',digest)
    palette_head(4)
    palette_head(5)
    conditional_six()
    all_C4_six()
    cylinder_gram()
    matching_and_type_controls()
    print('PASS_OPTIMIZATION_LIVE_GATES',GATES)


if __name__=='__main__':
    main()
