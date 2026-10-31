"""Independent digital audits; never used as a substitute for the proofs."""
from itertools import combinations,product
import json
from pathlib import Path
from fractions import Fraction
import digital

BASE=Path(__file__).resolve().parent
ROOTS=json.loads((BASE/'actual_seed_roots.json').read_text())
states=[(r['state'][0],r['state'][1],r['state'][2],r['state'][0]-r['state'][3]) for r in ROOTS]
def require(v,m):
    if not v:raise RuntimeError(m)

# Direct independent update loop, every actual seed, every word up to length 10.
count=0
for initial in states:
    h,F0,q0,g0=initial
    for k in range(11):
        for word in product((0,1),repeat=k):
            F,q,g=F0,q0,g0
            for b in word:
                if not b:F+=q+g
                q=2*q+h*b;g=h*b
            require(F==digital.output(initial,word),('digital output',initial,word))
            if k:
                require(g==h*(q%2),('canonical g',initial,word))
                kappa=F0-digital.phi(h,q0)+(g0-h*(q0%2) if word[0]==0 else 0)
                require(F-digital.phi(h,q)==kappa,('potential invariant',initial,word))
                first=digital.singleton(initial,word[0])[2]
                back=digital.ancestor_word(h,first,q)
                require(back is not None,('backward ancestry',initial,word))
            count+=1
print('all_actual_273_seeds_all_words_length_0_to_10',count,'PASS')

for n in range(1<<15):
    require(digital.f(n)==digital.f_pair_formula(n),('pair formula',n))
    if n<1<<12:
        v=digital.radical_inverse(n)
        require(digital.f(n)-n*(1-v)==digital.error_runs(n),('exact radical-inverse error',n))
        if n:require(0<digital.error_runs(n)<=1+n.bit_length(),('asymptotic bound',n))
print('32768_recurrence_pair_formula_and_4096_exact_fraction_error_checks_PASS')

for h in sorted({r[0] for r in states}):
    for q in range(512):
        p=digital.phi(h,q)
        require(digital.phi(3*h,3*q)==3*p,('scaling',h,q))
        for j in range(4):
            expect=(p+3*q+h*(q%2),p+q+h*(q%2),p+2*q+2*h,p)[j]
            require(digital.phi(h,4*q+j*h)==expect,('macro potential',h,q,j))
print('all_seed_h_512_q_macro_potential_and_scaling_PASS')

# Operator identities tested on all actual seed states and nonempty-word states.
macro=0
for initial in states:
    for word in product((0,1),repeat=4):
        state=initial
        for b in word:state=digital.singleton(state,b)
        h,F,q,g=state
        for j in range(4):
            shadow=digital.singleton(digital.singleton(state,j//2),j%2)
            actual=digital.triangle(state,j)
            defect=[a-3*b for a,b in zip(actual,shadow)]
            expected=((0,0,0,0),(0,3*q-g+2*h,0,-2*h),(0,-3*q+g-4*h,0,2*h),(0,0,0,0))[j]
            require(tuple(defect)==expected,('two-bit defect',state,j))
            macro+=1
print('actual_state_two_bit_triangle_operator_identities',macro,'PASS')

# Literal adjacency and forward path DP: no profile or potential is used to count.
def actual_graph(seed,blocks):
    n,code,cut=seed;A=[[False]*n for _ in range(n)]
    for k,(i,j) in enumerate(combinations(range(n),2)):
        A[i][j]=bool(code>>k&1);A[j][i]=not A[i][j]
    bits=[bool(cut>>i&1) for i in range(n)]
    for block in blocks:
        old=len(A);size=3 if isinstance(block,str) else 1
        A=[row+[True]*size for row in A]+[[False]*(old+size) for _ in range(size)]
        if size==3:
            j=int(block[1])
            for i in range(3):A[old+i][old+(i+1)%3]=True
            bits.extend(i<j for i in range(3))
        else:bits.append(bool(block))
    return [row+[not bit] for row,bit in zip(A,bits)]+[bits+[False]]
def count_paths(A):
    n=len(A);dp=[[0]*n for _ in range(1<<n)]
    for i in range(n):dp[1<<i][i]=1
    for mask in range(1,1<<n):
        for v,c in enumerate(dp[mask]):
            if c:
                for u in range(n):
                    if not(mask>>u&1) and A[v][u]:dp[mask|1<<u][u]+=c
    return sum(dp[-1])
literal=0
for record,initial in zip(ROOTS,states):
    for blocks in ((0,1,0),(1,0,1),('C1',),('C2',),(0,'C1',1),(1,'C2',0)):
        state=initial
        for block in blocks:state=digital.triangle(state,int(block[1])) if isinstance(block,str) else digital.singleton(state,block)
        require(count_paths(actual_graph(record['witness'],blocks))==state[1],('literal adjacency',record,blocks))
        literal+=1
print('actual_seed_literal_adjacency_DP_checks',literal,'PASS')

# Global exclusion in the one-seed transitive grammar follows after stripping
# trailing binary ones from B: 13 = f(2n) forces 1<=n<=12.
require(all(digital.f(n)+n+n%2!=13 for n in range(1,13)),'f image contains 13')
print('global_one_transitive_seed_f_image_excludes_13_by_finite_inverse_reduction_PASS')
print('PASS_no_global_spectrum_completeness_claim')

inverse=0
for state in states:
    for k in range(5):
        for word in product((0,1),repeat=k):
            target=digital.output(state,word)
            solutions=digital.inverse_addresses(target,state)
            B=digital.address(word)
            if state[2]==0:
                while B%2:B//=2
            require(B in solutions,('complete canonical inverse',state,word,target,B,solutions))
            for C in solutions:
                recovered=tuple(1-((C>>i)&1) for i in range(C.bit_length()))
                require(digital.output(state,recovered)==target,('inverse soundness',state,C))
            inverse+=1
print('actual_seed_scalar_inverse_checks',inverse,'PASS')

# The fiber table is not a union of formal offsets: validate recurrence and
# frozen boundary inheritance for each actual initial h and all finite sources.
for h in sorted({s[0] for s in states}):
    sources={}
    for state in states:
        if state[0]==h:
            for q,kappa,b in digital.kappa_sources(state):sources.setdefault(q,set()).add(kappa)
    L=max(sources)
    for q in range(2*L+h+1):
        kappa={v-digital.phi(h,q) for v in digital.fiber(h,q,states)}
        previous=digital.predecessor(h,q)
        inherited=set() if previous is None else {v-digital.phi(h,previous[0]) for v in digital.fiber(h,previous[0],states)}
        require(kappa==set(sources.get(q,()))|inherited,('fiber ancestry recurrence',h,q))
        if q>L:require(kappa==inherited,('eventual frozen fiber',h,q))
print('exact_actual_seed_fiber_ancestry_and_eventual_freezing_PASS')
