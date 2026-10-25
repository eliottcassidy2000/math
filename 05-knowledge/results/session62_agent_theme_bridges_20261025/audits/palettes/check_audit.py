#!/usr/bin/env python3
"""Independent audit: marked-block recurrence, literal witness joins and direct response counts.
Does not import or execute producer code. All tests remain active with -O.
"""
from itertools import product, permutations
from math import comb
from fractions import Fraction
import json

checks = 0
counts = {}
def need(value, label):
    global checks
    checks += 1
    if not value:
        raise RuntimeError(label)

def signature(values, k):
    # No first-occurrence canonicalization: named pins and the complete equality graph.
    return (tuple(v if v < k else None for v in values),
            tuple(values[i] == values[j] for i in range(len(values)) for j in range(i)))

def classes(q, k, n):
    result = {}
    for x in product(range(q), repeat=n):
        result.setdefault(signature(x,k), set()).add(x)
    return result

def block_recurrence(q,k,w):
    # Grow one variable: any of k names, any existing unnamed block, or one new block.
    t = [1] + [0] * (q-k)
    for _ in range(w):
        t = [(k+r)*t[r] + (t[r-1] if r else 0) for r in range(q-k+1)]
    return sum(t)

def stirling(n,b):
    # Onto functions / b!, independent of the recurrence under test.
    from math import factorial
    return sum((-1)**(b-j)*comb(b,j)*j**n for j in range(b+1))//factorial(b)

def stated(q,k,w):
    return sum(comb(w,j)*k**(w-j)*sum(stirling(j,b) for b in range(min(j,q-k)+1))
               for j in range(w+1))

for q in range(2,10):
    for k in range(q+1):
        for w in range(10):
            need(block_recurrence(q,k,w)==stated(q,k,w), 'marked recurrence versus claimed formula')
            if k==0:
                need(stated(q,k,w)==sum(stirling(w,b) for b in range(min(w,q)+1)), 'k=0')
            if k==q:
                need(stated(q,k,w)==q**w,'k=q')
counts['recurrence_parameter_cases'] = sum(q+1 for q in range(2,10))*10

explicit = 0
for q in range(2,5):
    for k in range(q+1):
        group = [tuple(range(k))+tail for tail in permutations(range(k,q))]
        for n in range(5):
            cs = classes(q,k,n)
            need(len(cs)==stated(q,k,n), 'equality graphs count the claimed orbits')
            for assignments in cs.values():
                representative = next(iter(assignments))
                orbit = {tuple(g[c] for c in representative) for g in group}
                need(orbit==assignments, 'equality graph is exactly one group orbit')
            # Burnside directly on combined exterior plus dead/free/forced state carrier.
            fixed_sum = sum(sum(g[c]==c for c in range(q))**n *
                            (2+sum(g[c]==c for c in range(q))) for g in group)
            need(fixed_sum//len(group)==2*stated(q,k,n)+stated(q,k,n+1), 'actual action row budget')
            explicit += 1
counts['explicit_group_parameter_cases'] = explicit
need(2*stated(6,0,6)+stated(6,0,7)==1282,'published path budget')
need(2*stated(6,0,7)+stated(6,0,8)==5863,'published anchored budget')

def subsets(keys):
    for bits in range(1 << len(keys)):
        yield {keys[i] for i in range(len(keys)) if (bits>>i)&1}

joins = 0
for q,k in ((2,0),(2,1),(2,2),(3,0),(3,1)):
    binary=classes(q,k,2)
    ternary=classes(q,k,3)
    keys=list(binary)
    for R in subsets(keys):
        left=set().union(*(binary[x] for x in R)) if R else set()
        for S in subsets(keys):
            right=set().union(*(binary[x] for x in S)) if S else set()
            # Literal witness matching, not evaluation of the orbit join formula.
            glued={(a,b,c) for a,b in left for bb,c in right if b==bb}
            literal={signature(x,k) for x in glued}
            # A union orbit qualifies by restriction of a representative.
            candidate={key for key,values in ternary.items()
                       if any(signature(x[:2],k) in R and signature(x[1:],k) in S
                              for x in [next(iter(values))])}
            need(candidate==literal,'literal witness joins versus partition amalgamation')
            projected={signature((x[0],x[2]),k) for x in glued}
            orbit_projected={signature((next(iter(ternary[key]))[0],
                                       next(iter(ternary[key]))[2]),k) for key in candidate}
            need(projected==orbit_projected,'literal witness projection versus orbit projection')
            joins+=1
counts['exhaustive_relation_pairs'] = joins

# Minimal independent-representative hostile and missing symmetry hostile.
need({signature((a,b),0) for a,b in product(range(2),repeat=2)} !=
     {signature((0,0),0)}, 'independent unary representatives lose unequal alignment')
need(signature((0,),0)==signature((1,),0) and
     len({0}-{0})==0 and len({0}-{1})==1,
     'unnamed quotient cannot consume a fixed named singleton letter')

# Literal proper colourings; all matrix entries below are obtained from enumeration.
def full_count(q, word):
    return sum(all(c in L for c,L in zip(colours,word)) and
               all(colours[i]!=colours[i+1] for i in range(len(word)-1))
               for colours in product(range(q),repeat=len(word)))

def pinned(q,word):
    result=[[0]*q for _ in range(q)]
    if not word:
        return [[int(i==j) for j in range(q)] for i in range(q)]
    for a in range(q):
        for colours in product(range(q),repeat=len(word)):
            if (colours[0]!=a and all(c in L for c,L in zip(colours,word)) and
                all(colours[i]!=colours[i+1] for i in range(len(word)-1))):
                result[colours[-1]][a]+=1
    return result

def rank(rows, p=0):
    rows = [[x%p if p else Fraction(x) for x in row] for row in rows]
    r=0
    for j in range(len(rows[0])):
        i=next((i for i in range(r,len(rows)) if rows[i][j]),None)
        if i is None:
            continue
        rows[r],rows[i]=rows[i],rows[r]
        factor=pow(rows[r][j],-1,p) if p else 1/rows[r][j]
        rows[r]=[(x*factor)%p if p else x*factor for x in rows[r]]
        for i in range(r+1,len(rows)):
            multiplier=rows[i][j]
            rows[i]=[(x-multiplier*y)%p if p else x-multiplier*y
                     for x,y in zip(rows[i],rows[r])]
        r+=1
        if r==len(rows):
            break
    return r

rank_receipts=[]
for q in range(2,9):
    basis=[]
    for i,j in product(range(q),repeat=2):
        word=({i},) if i==j else ({j},{i})
        T=pinned(q,word)
        need(T==[[int(a==i and b!=j) for b in range(q)] for a in range(q)],
             'literal word realizes the asserted outer product')
        basis.append(sum(T,[]))
    need(rank(basis)==q*q,'direct-count named two-port basis over Q')
    small_words=[()]+[({i},) for i in range(q)]
    H=[[full_count(q,u+v) for v in small_words] for u in small_words]
    # Subtract the first row from every other row: [1 1^T; 0 -I].
    reduced=[H[0]]+[[x-y for x,y in zip(row,H[0])] for row in H[1:]]
    need(reduced==[[1]*(q+1)]+[[0]+[-int(i==j) for j in range(q)] for i in range(q)],
         'unimodular triangular Hankel determinant proof')
    need(rank(H)==q+1,'direct-count full series minor over Q')
    for p in (2,3,5,7):
        need(rank(H,p)==q+1,'Hankel minor remains full in every tested characteristic')
        A=[[int(i!=j) for j in range(q)] for i in range(q)]
        need(rank(A+[[1]*q],p)==q,'prefix rank survives with empty continuation')
        exceptional=(q-1)%p==0
        expected=q*(q-1) if exceptional else q*q
        need(rank(basis,p)==expected,'positive-word span characteristic boundary')
        I=sum([[int(i==j) for j in range(q)] for i in range(q)],[])
        need(rank(basis+[I],p)==expected+int(exceptional),'whole two-port response characteristic boundary')
    rank_receipts.append({'q':q,'prefix_Q':q,'whole_series_Q':q+1,'two_port_Q':q*q})
counts['direct_basis_words'] = sum(q*q for q in range(2,9))

# A larger nontrivial Hankel block, beyond only the decisive minor.
for q,depth in ((2,2),(3,1)):
    alphabet=[frozenset(i for i in range(q) if (mask>>i)&1) for mask in range(1<<q)]
    words=[word for n in range(depth+1) for word in product(alphabet,repeat=n)]
    block=[[full_count(q,u+v) for v in words] for u in words]
    need(rank(block)==q+1,'larger literal Hankel block')

# Nonnegative integer counts are not signed spans and are not Boolean support mod p.
need(full_count(3,(frozenset(range(3)),frozenset(range(3))))==6,
     'six positive counts become zero modulo two and three')
need(full_count(2,())==1 and full_count(2,({0,1},))==2,
     'empty prefix and free nonempty prefix have different weights')

# Higher-order relation via permutations of the inner triangle, independently of 3^6 enumeration.
weights={tips:sum(all(triangle[i]!=tips[i] for i in range(3))
                  for triangle in permutations(range(3))) for tips in product(range(3),repeat=3)}
need(all(w==(0 if len(set(tips))==1 else 2) for tips,w in weights.items()),'three-port permanent hostile')
for i,j in ((0,1),(0,2),(1,2)):
    need({(tips[i],tips[j]) for tips,w in weights.items() if w}==set(product(range(3),repeat=2)),
         'all pair marginals full')
need(sum(w>0 for w in weights.values())==24 and sum(weights.values())==48,'hostile support and multiplicity')

# A scalar attached to the current deterministic Boolean state cannot represent
# nonproportional continuation rows, even though a general q+2-space realization can.
u,v=({0,1},),({1,2},)
x=tuple(full_count(3,u+({i},)) for i in (0,2))
y=tuple(full_count(3,v+({i},)) for i in (0,2))
need(x==(1,2) and y==(2,1) and x[0]*y[1]!=x[1]*y[0],
     'nonproportional free-state responses refute one scalar sidecar')

# Integral row lattice: exact constructive preimages, then the identity quotient.
integral_vectors=0
residue_vectors=0
for q in range(2,7):
    m=q-1
    for v in product(range(-2,3),repeat=q):
        if sum(v)%m:
            continue
        t=sum(v)//m
        coefficients=tuple(t-vj for vj in v)
        recovered=tuple(sum(coefficients[j] for j in range(q) if j!=i) for i in range(q))
        need(recovered==v,'constructive integral row-lattice inverse')
        integral_vectors+=1
    fibers={}
    for residue in product(range(m),repeat=q):
        image=tuple((residue[i]-residue[-1])%m for i in range(q-1))
        fibers[image]=fibers.get(image,0)+1
        residue_vectors+=1
    need(len(fibers)==m**(q-1) and set(fibers.values())=={m},
         'identity kills exactly the diagonal cyclic residue subgroup')
counts['integral_row_vectors'] = integral_vectors
counts['integral_residue_vectors'] = residue_vectors

print(json.dumps({'status':'PASS','checks':checks,'universes':counts,'ranks':rank_receipts,
                 'exceptional_characteristic_example':{'q':3,'p':2,'positive_word_span':6,
                                                      'including_empty_word_span':7,
                                                      'whole_series_rank':4}},indent=2))
