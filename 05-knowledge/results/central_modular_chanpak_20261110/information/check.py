#!/usr/bin/env python3
"""Exact row fibres and globally primitive golden content contractions."""
from collections import Counter, deque
from fractions import Fraction
from functools import lru_cache, reduce
from hashlib import sha256
from itertools import combinations, permutations, product
from math import gcd
from pathlib import Path
import ast
import json

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
OLD = ROOT/'05-knowledge/results/chanpak_55_105_golden_20261110/golden'
I = ((1,0,0),(0,1,0),(0,0,1))
M = {
 's':((0,1,1),(0,0,1),(1,0,1)),
 'c':((1,0,4),(1,0,3),(3,1,6)),
 'q0':((1,0,8),(1,0,7),(7,1,18)),
 'q1':((3,0,6),(3,0,5),(9,3,14)),
 'q2':((4,1,7),(3,1,6),(6,1,14)),
 'q4':((3,0,12),(3,0,9),(5,1,14)),
 'q5':((4,1,9),(4,0,7),(7,2,11)),
}
WORDS = {
 5:('q2','c','c','q2','s','s','s','s','q2'),
 9:('q1','s','s','q1','s','s','q1'),
 21:('q5','s','q0','c','q0','q5','q0','c','s','s','s','q5'),
}

def require(q, message):
    if not q: raise RuntimeError(message)

def mul(A, B):
    return tuple(tuple(sum(A[i][k]*B[k][j] for k in range(3)) for j in range(3))
                 for i in range(len(A)))

def row(v, A):
    return tuple(sum(x*y for x,y in zip(v,col)) for col in zip(*A))

def word_matrix(word):
    A = I
    for letter in word: A = mul(A, M[letter])
    return A

def determinant(A):
    return (A[0][0]*(A[1][1]*A[2][2]-A[1][2]*A[2][1])
            -A[0][1]*(A[1][0]*A[2][2]-A[1][2]*A[2][0])
            +A[0][2]*(A[1][0]*A[2][1]-A[1][1]*A[2][0]))

def permutation_det(A):
    return sum((-1)**sum(p[i]>p[j] for i in range(3) for j in range(i+1,3))
               *A[0][p[0]]*A[1][p[1]]*A[2][p[2]] for p in permutations(range(3)))

def content(v): return reduce(gcd, v, 0)

def smith(A):
    d1 = content(x for a in A for x in a)
    d2 = content(A[i][k]*A[j][l]-A[i][l]*A[j][k]
                 for i,j in combinations(range(3),2) for k,l in combinations(range(3),2))
    return (d1,d2//d1,abs(determinant(A))//d2)

def inverse_unimodular(A):
    require(determinant(A)==1,'positive unimodularity')
    out = []
    for i in range(3):
        line = []
        for j in range(3):
            rs = [r for r in range(3) if r!=j]
            cs = [c for c in range(3) if c!=i]
            line.append((-1)**(i+j)*(A[rs[0]][cs[0]]*A[rs[1]][cs[1]]
                                     -A[rs[0]][cs[1]]*A[rs[1]][cs[0]]))
        out.append(tuple(line))
    B = tuple(out)
    require(mul(A,B)==mul(B,A)==I,'integer inverse')
    return B

def rref(A, p):
    A = [[x%p for x in a] for a in A]; r=0
    for col in range(3):
        pivot = next((i for i in range(r,len(A)) if A[i][col]),None)
        if pivot is None: continue
        A[r],A[pivot]=A[pivot],A[r]
        unit=pow(A[r][col],-1,p); A[r]=[unit*x%p for x in A[r]]
        for i in range(len(A)):
            if i!=r:
                q=A[i][col]; A[i]=[(x-q*y)%p for x,y in zip(A[i],A[r])]
        r+=1
        if r==len(A): break
    return tuple(tuple(a) for a in A[:r])

def alignment(A, primes, singular, desired_ranks):
    start=tuple(rref(A,p) for p in primes)
    queue=deque([start]); paths={start:()}
    while queue:
        state=queue.popleft()
        candidate=tuple(rref(mul(S,M[singular]),p) for S,p in zip(state,primes))
        if tuple(map(len,candidate))==desired_ranks:
            return paths[state],len(paths)
        for letter in ('s','c','q0'):
            nxt=tuple(rref(mul(S,M[letter]),p) for S,p in zip(state,primes))
            if nxt not in paths:
                paths[nxt]=paths[state]+(letter,); queue.append(nxt)
    raise RuntimeError('alignment unavailable')

def factors(n):
    out={};p=2
    while p*p<=n:
        while n%p==0: out[p]=out.get(p,0)+1;n//=p
        p+=1
    if n>1:out[n]=out.get(n,0)+1
    return out

def valuation(n,p):
    k=0
    while n%p==0:n//=p;k+=1
    return k

def atom(name):
    n=1 if name=='s' else 3 if name=='c' else 4
    bits=0 if name in ('s','c') else int(name[1:])
    pairs=[(i,j) for i in range(n) for j in range(i+2,n)]
    flags={pair:(bits>>k)&1 for k,pair in enumerate(pairs)}
    # Complement every nonpath edge of the chart, then reverse vertex labels.
    graph=[[False]*n for _ in range(n)]
    for i,j in combinations(range(n),2):
        forward=True if j==i+1 else not flags[i,j]
        a,b=(i,j) if forward else (j,i)
        graph[n-1-a][n-1-b]=True
    return tuple(sum(int(bit)<<j for j,bit in enumerate(a)) for a in graph)

def literal(word):
    rows=()
    for letter in word:
        block=atom(letter);a,b=len(rows),len(block)
        new=[r|(((1<<b)-1)<<a) for r in rows]+[r<<a for r in block]
        if a:
            new[a-1]^=1<<a; new[a]|=1<<(a-1)
        rows=tuple(new)
    return rows

def literal_state(rows):
    n=len(rows);full=(1<<n)-1;port=n-1
    @lru_cache(None)
    def end(mask,v):
        if mask==1<<v:return 1
        rest=mask^(1<<v)
        return sum(end(rest,u) for u in range(n) if rest>>u&1 and rows[u]>>v&1)
    @lru_cache(None)
    def start(mask,v):
        if mask==1<<v:return 1
        rest=mask^(1<<v)
        return sum(start(rest,u) for u in range(n) if rest>>u&1 and rows[v]>>u&1)
    @lru_cache(None)
    def H(mask):return 1 if not mask else sum(end(mask,v) for v in range(n) if mask>>v&1)
    h,e=H(full),end(full,port)
    c=sum(start(mask,port)*H(full^mask) for mask in range(1<<n) if mask>>port&1)
    return h-e,e,c

def main():
    assignments={}
    for node in ast.parse((OLD/'check.py').read_text()).body:
        if isinstance(node,ast.Assign) and len(node.targets)==1 and isinstance(node.targets[0],ast.Name):
            if node.targets[0].id in ('M','WORDS'):
                assignments[node.targets[0].id]=ast.literal_eval(node.value)
    require(M==assignments['M'],'pinned inherited matrices')
    require(WORDS[5]==assignments['WORDS'][5],'inherited unimodular5 word')
    dets={name:determinant(A) for name,A in M.items()}
    require(dets=={'s':1,'c':1,'q0':1,'q1':9,'q2':5,'q4':9,'q5':21},'determinant alphabet')
    rows_checked=0; central_checked=0; row_examples={}
    for modulus in range(2,62):
        for a in range(modulus):
            histogram=Counter(a*x%modulus for x in range(modulus));d=gcd(a,modulus)
            require(len(histogram)==modulus//d and set(histogram.values())=={d},'row fibres')
            rows_checked+=1
        ff=factors(modulus)
        for p,e in ff.items():
            actual=Fraction(sum(valuation(gcd(a,modulus),p) for a in range(modulus)),modulus)
            predicted=sum((Fraction(1,p**j) for j in range(1,e+1)),Fraction())
            require(actual==predicted,'exact mean information-loss exponent')
        nonzero_lossless=sum(gcd(a,modulus)==1 for a in range(1,modulus))
        require((nonzero_lossless==modulus-1)==(len(ff)==1 and next(iter(ff.values()))==1),'field row criterion')
        if modulus in (5,7,9,15):
            row_examples[modulus]={'lossless_nonzero_rows':nonzero_lossless,
                                  'fibre_sizes':[gcd(a,modulus) for a in range(1,modulus)]}
    for modulus in range(3,102,2):
        k=(modulus-1)//2;u=k*k%modulus;v=k*(k+1)%modulus
        require(u+v==modulus and gcd(u,modulus)==gcd(v,modulus)==1,'central units and integer sums')
        outputs=Counter(((u*x+v*y)%modulus,(v*x+u*y)%modulus)
                        for x in range(modulus) for y in range(modulus))
        require(len(outputs)==modulus and set(outputs.values())=={modulus},'central information loss')
        require(set(outputs)=={(t,-t%modulus) for t in range(modulus)},'central image')
        central_checked+=1
    primitive={}; inverses={}
    for scalar,word in WORDS.items():
        A=word_matrix(word)
        replay=[]
        for v in I:
            for letter in word:v=row(v,M[letter])
            replay.append(v)
        require(A==tuple(replay),'independent row-word replay')
        require(all(x%scalar==0 for a in A for x in a),'scalar divisibility')
        U=tuple(tuple(x//scalar for x in a) for a in A)
        require(determinant(A)==permutation_det(A)==scalar**3,'two determinant calculations')
        require(smith(A)==(scalar,)*3,'scalar Smith partition')
        inverse=inverse_unimodular(U);inverses[scalar]=inverse
        primitive[scalar]={'word':word,'vertex_cost':sum(1 if l=='s' else 3 if l=='c' else 4 for l in word),
                           'matrix':A,'normalized':U,'normalized_inverse':inverse,'normalized_det':1}
    alignment_trace=[];A=I;word=()
    for rank in (2,1,0):
        units,seen=alignment(A,(3,7),'q5',(rank,rank));word+=units+('q5',)
        for letter in units+('q5',):A=mul(A,M[letter])
        alignment_trace.append({'target_ranks':[rank,rank],'unit_word':units,'visited_states':seen})
    require(word==WORDS[21],'simultaneous3/7 construction')
    A=M['q1'];word=('q1',);nine_trace=[]
    for _ in range(2):
        units,seen=alignment(A,(3,),'q1',(0,));word+=units+('q1',)
        for letter in units+('q1',):A=mul(A,M[letter])
        require(all(x%3==0 for a in A for x in a),'remove exact3 valuation')
        A=tuple(tuple(x//3 for x in a) for a in A)
        nine_trace.append({'unit_word':units,'visited_states':seen})
    require(word==WORDS[9] and determinant(A)==1,'primitive9 construction')
    content_checks=0
    for b,k,c in product(range(3),repeat=3):
        scalar=5**b*9**k*21**c
        A=word_matrix(WORDS[5]*b+WORDS[9]*k+WORDS[21]*c)
        U=tuple(tuple(x//scalar for x in a) for a in A)
        require(determinant(U)==1,'unbounded-family formula finite controls')
        for v in product(range(-2,3),repeat=3):
            require(content(row(v,A))==scalar*content(v),'exact content for signed rows')
            content_checks+=1
    matrix_fibres=[]
    for scalar in (5,9,21,105):
        A=word_matrix(WORDS[scalar] if scalar in WORDS else WORDS[5]+WORDS[21])
        for modulus in (3,5,7,9,11,15):
            counts=Counter(tuple(x%modulus for x in row(v,A)) for v in product(range(modulus),repeat=3))
            fibres=gcd(scalar,modulus)**3
            require(len(counts)==modulus**3//fibres and set(counts.values())=={fibres},'all modular fibres')
            matrix_fibres.append({'scalar':scalar,'modulus':modulus,'kernel':fibres,'image':len(counts)})
    physical=[]
    for w in (('q1','s','s'),('q5','s','q0')):
        rows=literal(w);st=literal_state(rows)
        require(st==word_matrix(w)[0],'actual small bridge prefix')
        physical.append({'word':w,'rows':rows,'state':st})
    oldword=assignments['WORDS'][3]+assignments['WORDS'][5]+assignments['WORDS'][7]
    oldmatrix=word_matrix(oldword)
    oldnormal=tuple(tuple(x//105 for x in a) for a in oldmatrix)
    require(determinant(oldnormal)==105**3,'locally primitive product is not globally primitive')
    sources=[OLD/'check.py',OLD/'NOTE.md',ROOT/'05-knowledge/results/golden_hamilton_transfer_20261103/RECURRENCE.md',
             ROOT/'05-knowledge/results/session91_packing_golden_carriers_20261029/zeckendorf/NORM_CONTROLLED_OBSERVERS.md']
    print(json.dumps({'status':'PASS','scope':'exact row/central fibre controls and literal golden word certificates',
                      'row_maps_checked':rows_checked,'row_examples':row_examples,'central_moduli_checked':central_checked,
                      'primitive_words':primitive,'alignment21':alignment_trace,'alignment9':nine_trace,
                      'signed_content_checks':content_checks,'matrix_fibre_controls':matrix_fibres,
                      'physical_prefix_controls':physical,'nonprimitive_old_product_normalized_det':determinant(oldnormal),
                      'source_sha256':{str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in sources},
                      'script_sha256':sha256(Path(__file__).read_bytes()).hexdigest()},indent=2,sort_keys=True))

if __name__=='__main__':main()
