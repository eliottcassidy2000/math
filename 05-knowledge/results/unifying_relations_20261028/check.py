"""Exact controls for the arithmetic, incidence and descent synthesis."""
from collections import Counter, deque
from fractions import Fraction
from itertools import combinations, combinations_with_replacement, product
from math import gcd, isqrt, prod
import json
import sys
from sympy import Matrix, ZZ
from sympy.matrices.normalforms import smith_normal_form

sys.stdout.reconfigure(encoding='utf-8', newline='\n')

def need(ok, message):
    if not ok:
        raise RuntimeError(message)

def factor(n):
    result = {}; p = 2
    while p*p <= n:
        while n % p == 0:
            result[p] = result.get(p, 0)+1; n //= p
        p += 1
    if n > 1: result[n] = 1
    return result

def divisors(n):
    out = [1]
    for p, a in factor(n).items():
        out = [d*p**j for d in out for j in range(a+1)]
    return sorted(out)

def sigma(n):
    return prod(sum(p**i for i in range(a+1)) for p,a in factor(n).items())

def squarefree_part(n):
    return prod(p for p,a in factor(n).items() if a % 2)

def square(n):
    return n >= 0 and isqrt(n)**2 == n

def fus(n):
    ds = divisors(n)[1:-1]
    return [len(ds), sum(len(factor(d))==1 and next(iter(factor(d).values()))==1 for d in ds),
            sum(all(a==1 for a in factor(d).values()) for d in ds)]

def mask_of_squareclass(n, coordinates):
    f = factor(n)
    need(all(p in coordinates for p,a in f.items() if a % 2), 'coordinate coverage')
    return sum((f.get(p,0)%2)<<i for i,p in enumerate(coordinates))

def bank(primes, coordinates):
    need(all(factor(p)=={p:1} for p in primes),'prime bank guard')
    cols = [mask_of_squareclass(p+1, coordinates) for p in primes]
    solutions = []
    for bits in product((0,1), repeat=len(primes)):
        indices = [i for i,b in enumerate(bits) if b]
        total = prod(primes[i]+1 for i in indices)
        source=prod(primes[i] for i in indices)
        need(total==sum(divisors(source)), 'divisor-addition control')
        parity = 0
        for i in indices: parity ^= cols[i]
        need((parity==0)==square(total), 'parity versus integer square')
        if parity==0:
            solutions.append({'support':[primes[i] for i in indices],
                              'source':prod(primes[i] for i in indices), 'sigma':total})
    supports = [frozenset(row['support']) for row in solutions if row['support']]
    circuits = [sorted(a) for a in supports if not any(b<a for b in supports)]
    return {'primes':primes, 'coordinates':coordinates, 'columns':cols,
            'weight_enumerator':dict(sorted(Counter(len(row['support']) for row in solutions).items())),
            'circuits':sorted(circuits), 'square_sigma_sources':solutions}

def incidence(n, edges):
    return Matrix([[int(v==a or v==b) for v in range(n)] for a,b in edges])

def connected_bipartite(n, edges):
    adj=[[] for _ in range(n)]
    for a,b in edges: adj[a].append(b); adj[b].append(a)
    colours={0:0}; todo=[0]; bip=True
    for u in todo:
        for v in adj[u]:
            if v not in colours: colours[v]=1-colours[u]; todo.append(v)
            elif colours[v]==colours[u]: bip=False
    return len(colours)==n, bip

def smith_receipt(n, edges):
    B=incidence(n,edges)
    D=smith_normal_form(B,domain=ZZ)
    diagonal=[abs(int(D[i,i])) for i in range(min(D.rows,D.cols)) if D[i,i]]
    connected,bip=connected_bipartite(n,edges)
    need(connected,'connected universe')
    expected=[1]*(n-1)+([] if bip else [2])
    need(diagonal==expected,'integer incidence Smith law')
    return {'vertices':n,'edges':len(edges),'bipartite':bip,'smith_nonzero':diagonal,
            'cokernel_free_rank':len(edges)-len(diagonal)}

def parent(pair):
    m,n=pair
    if pair==(2,1): return None
    if m<2*n: return (n,2*n-m),'U'
    if m<3*n: return (n,m-2*n),'A'
    return (m-2*n,n),'D'

def word_and_character(pair):
    letters=[]; x=pair
    while x!=(2,1):
        result=parent(x)
        need(result is not None, 'nonroot parent')
        y,letter=result
        need(y[0]>y[1]>0 and gcd(*y)==1 and (sum(y)%2)==1 and sum(y)<sum(x), 'Berggren descent')
        x=y; letters.append(letter)
    return ''.join(reversed(letters)), (-1)**letters.count('A')

def apply_word(word):
    m,n=2,1
    for c in word:
        if c=='U': m,n=2*m-n,m
        elif c=='A': m,n=2*m+n,m
        else: m,n=m+2*n,n
    return m,n

def main():
    small=bank([2,3,5,7,11,13],[2,3,7])
    need(small['circuits']==[[2,5,7],[2,11],[3],[5,7,11]],'small-bank minimal dependencies')
    need([r['source'] for r in small['square_sigma_sources'] if 13 in r['support']]==[], '13 is a bridge in this bank')
    fano=bank([7,2,19,5,89,59,29],[2,3,5])
    need(set(fano['columns'])==set(range(1,8)),'all Fano points occur')
    need(fano['weight_enumerator']=={0:1,3:7,4:7,7:1},'Hamming weight polynomial')
    need(len(fano['circuits'])==14,'seven lines and seven four-element circuits')

    tower=[]
    for a in range(41):
        values=[(2**a)*m for m in (33,35,47)]
        totals=[sigma(n) for n in values]
        need(totals==[48*(2**(a+1)-1)]*3,'common sigma doubling towers')
        need(square(totals[0])==(a==1),'unique square layer')
        for n in values[:2]:
            F,U,S=fus(n)
            need((F==U+S)==(a==2),'unique FUS layer')
        if a<=4: tower.append({'a':a,'sources':values,'sigma':totals[0], 'F_U_S':fus(values[0])})
    fibres={str(t):[n for n in range(1,t) if sigma(n)==t] for t in (48,144,336)}
    for t in (48,144,336):
        need(fibres[str(t)]==[n for n in range(1,t) if sum(d for d in range(1,n+1) if n%d==0)==t],
             'complete fibres by independent direct divisor addition')
    need(fibres['48']==[33,35,47], 'complete sigma48 fibre')
    need([n for n in fibres['336'] if (lambda v:v[0]==v[1]+v[2])(fus(n))]==[132,140], 'FUS sources in complete sigma336 fibre')
    need(fus(650)==[10,3,7] and sigma(650)-650==652 and len(divisors(650)[:-1])==11,
         'inherited balanced excess-two positive control')
    profile_count=0
    for r in range(1,7):
        for aa in combinations_with_replacement(range(1,7),r):
            F=prod(a+1 for a in aa)-2
            U=r-int(aa==(1,))
            S=2**r-1-int(all(a==1 for a in aa))
            need((F==U+S)==(aa in ((1,),(3,),(1,1,2))), 'inherited FUS shape formula')
            profile_count+=1

    exponents={2:5,3:7,5:3,7:4,11:1,13:2}
    N=prod(p**a for p,a in exponents.items())
    heights={}
    for b in divisors(N)[1:]:
        h=min(exponents[p]//a for p,a in factor(b).items())
        k=0; power=1
        while N%(power*b)==0: power*=b; k+=1
        need(h==k,'valuation minimum versus integer division')
        heights[b]=h
    tail=[prod(a//k+1 for a in exponents.values())-1 for k in range(1,9)]
    need(tail==[5759,143,23,7,3,1,1,0],'inherited all-divisor tails')
    cycle=[4,2,5,3,7]
    need([heights[b] for b in cycle]==cycle[1:]+cycle[:1] and heights[13]==2,'literal repaired cycle')
    need([min(factor(N//2).get(p,0)//a for p,a in factor(b).items()) for b in (2,4)]==[4,2], 'normalization changes cycle')
    for a,b,k in product(range(1,41),range(1,41),range(1,5)):
        root=prod(p**(e//k) for p,e in factor(a).items())
        need((a%(b**k)==0)==(root%b==0),'power/root divisibility adjunction')

    graph_count=0
    for n in range(2,6):
        possible=list(combinations(range(n),2))
        for mask in range(1<<len(possible)):
            edges=[e for i,e in enumerate(possible) if mask>>i&1]
            if connected_bipartite(n,edges)[0]: smith_receipt(n,edges); graph_count+=1
    need(graph_count==771,'complete connected graph universe through five vertices')
    k5=list(combinations(range(5),2))
    k33=[(a,b) for a in range(3) for b in range(3,6)]
    petersen=sorted({tuple(sorted(e)) for i in range(5) for e in ((i,(i+1)%5),(i,i+5),(i+5,(i+2)%5+5))})
    graph_results={'K5':smith_receipt(5,k5),'K33':smith_receipt(6,k33),'Petersen':smith_receipt(10,petersen)}
    contracted={tuple(sorted((a%5,b%5))) for a,b in petersen if a%5!=b%5}
    need(contracted==set(k5),'Petersen spoke contraction gives K5')
    need(all(sum(1 for w in range(10) if tuple(sorted((u,w))) in petersen and tuple(sorted((v,w))) in petersen)<=1 for u,v in combinations(range(10),2)), 'Petersen common-neighbour control')
    k33_labels=[0,1,2,3,6,9]
    need(sorted(abs(k33_labels[a]-k33_labels[b]) for a,b in k33)==list(range(1,10)), 'explicit graceful K33')
    A=[1,193,673]; B=[3,168,483]
    need(all(square(a+b) for a in A for b in B),'inherited positive square-sum K33')
    need(not any(square(a+b) for shore in (A,B) for a,b in combinations(shore,2)), 'induced K33 guard')
    triangle=incidence(3,[(0,1),(1,2),(0,2)])
    need(triangle.inv()*Matrix([1,1,1])==Matrix([Fraction(1,2)]*3), 'odd-cycle integrality hostile')
    need(sorted(triangle*Matrix([6,19,30]))==[25,36,49], 'positive square-sum triangle')
    need([abs(a-b) for a,b in ((0,1),(1,2))]==[1,1], 'tree absolute-difference collision')
    L=Matrix([[2,-1,-1],[-1,2,-1],[-1,-1,2]])
    LD=smith_normal_form(L,domain=ZZ)
    need([abs(int(LD[i,i])) for i in range(3)]==[1,3,0], 'Laplacian versus incidence torsion')

    berggren_count=0
    for m in range(2,101):
        for n in range(1,m):
            if gcd(m,n)==1 and (m+n)%2:
                word,chi=word_and_character((m,n))
                need(apply_word(word)==(m,n),'exact inverse branch reconstruction')
                berggren_count+=1
    for modulus in range(1,201):
        m,n=16*modulus+2,4*modulus+1
        word,chi=word_and_character((m,n))
        need(chi==-1 and all((a-b)%modulus==0 for a,b in zip((m*m-n*n,2*m*n,m*m+n*n),(3,4,5))), 'every-modulus character hostile')
    ray_prime_indices=[]
    for k in range(1,101):
        q=12*k*k-1; m,n=apply_word('D'*(k-1))
        need((m,n)==(2*k,1), 'Berggren D-ray coordinates')
        need((q-2,12*k,q+4)==tuple(3*x for x in (m*m-n*n,2*m*n,m*m+n*n)),
             'scale-three centered triangle')
        if factor(q)=={q:1}:
            need(sigma(2*q)==36*k*k and sigma(6*q)==144*k*k, 'prime-guarded arithmetic ray')
            ray_prime_indices.append(k)
    need(factor(299)=={13:1,23:1} and sigma(598)==1008, 'ray primality hostile')
    bracket_sources=[]
    for n in range(2,1001):
        m=1
        while n>(2*m+1)**2: m+=1
        if 2*n<=(2*m+1)**2: bracket_sources.append(n)
    need(bracket_sources==[2,3,4,10,11,12],'inherited bracket sources')

    print(json.dumps({'small_prime_bank':small,'fano_prime_bank':fano,'doubling_towers':tower,
          'complete_sigma_fibres':fibres,'FUS_profile_universe':profile_count,
          'divisor_height':{'N':N,'tail_counts':tail,'exact_height_counts':dict(sorted(Counter(heights.values()).items())), 'cycle':cycle},
          'graphs':graph_results,'connected_graphs_through5':graph_count,
          'Petersen_K5_minor_branch_sets':[[i,i+5] for i in range(5)],
          'Berggren_pair_count_m_at_most100':berggren_count,'Berggren_modulus_controls':200,
          'Berggren_semiprime_ray_prime_indices_k_at_most100':ray_prime_indices,
          'odd_square_bracket_sources_through1000':bracket_sources},indent=2,sort_keys=True))

if __name__=='__main__': main()
