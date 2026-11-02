#!/usr/bin/env python3
"""Standalone integer referee for finite-bank H-spectrum order compression.
No external packages or repository imports. The 80265 witness bank is inherited,
not regenerated here; all new inequalities/batching rules are independently checked.
"""
import argparse, itertools, json, math, random

B, D, L = 80265, 80267, 285

def demand(ok, detail):
    if not ok:
        raise RuntimeError(detail)

def floor_strong(n):
    if n < 3: raise ValueError(n)
    k, r = divmod(n, 3)
    return (3*5**(k-1), 5**k, 9*5**(k-1))[r]

def cutoff(h):
    """Least n >= 0 with D**n >= h**22; exact ceil(22 log_D h)."""
    n, power, target = 0, 1, h**22
    while power < target:
        n += 1; power *= D
    return n

def small_cost(h):
    if h == 1: return 0
    if h in (7, 21) or h <= 0 or h % 2 == 0 or h >= L:
        raise ValueError(h)
    if h == 3: return 3
    if h == 5: return 4
    if h <= 15: return 5
    if h <= 33: return 6
    if h <= 61: return 7
    return 8

def certify_repack(factors):
    """Factors are (actual H, order), where supplied floor is a necessary guard.
    This returns an arithmetic recipe using inherited small/11-vertex banks,
    not synthetic tournament adjacencies. For tests all factors come from
    actual enumerated small strong tournaments, Moon minimum chains, or both.
    """
    product, order, pending = 1, 0, 1
    pieces=[]
    for h, n in factors:
        demand(h >= 3 and h % 2 and h not in (7,21), ('factor',h,n))
        demand(h >= floor_strong(n), ('floor',h,n))
        product *= h
        if h < L:
            pending *= h
            if pending >= L:
                demand(pending <= (L-2)**2 <= B, ('batch',pending))
                pieces.append(('bank', pending, 11)); order += 11
                pending = 1
        elif h <= B:
            pieces.append(('bank', h, 11)); order += 11
        else:
            demand(h**22 >= D**n, ('large',h,n))
            pieces.append(('original', h, n)); order += n
    if pending > 1:
        n=small_cost(pending)
        pieces.append(('residual',pending,n)); order += n
    demand(math.prod(p[1] for p in pieces)==product, 'product mismatch')
    demand(order <= cutoff(product), ('cutoff', factors, order, cutoff(product)))
    return {'input_order':sum(n for h,n in factors),'count':product,
            'recipe_order':order,'universal_cutoff':cutoff(product),'pieces':pieces}

def adjacency(n, code):
    a=[0]*n
    for bit,(i,j) in enumerate(itertools.combinations(range(n),2)):
        u,v=(i,j) if (code>>bit)&1 else (j,i)
        a[u] |= 1<<v
    return a

def hcount(a):
    n=len(a)
    if not n: return 1
    dp=[[0]*n for _ in range(1<<n)]
    for v in range(n): dp[1<<v][v]=1
    for mask in range(1,1<<n):
        for v,c in enumerate(dp[mask]):
            if c:
                todo=a[v]&~mask
                while todo:
                    wbit=todo&-todo; todo-=wbit
                    dp[mask|wbit][wbit.bit_length()-1] += c
    return sum(dp[-1])

def literal(a):
    return sum(all((a[u]>>v)&1 for u,v in zip(p,p[1:]))
               for p in itertools.permutations(range(len(a))))

def sccs(a):
    n=len(a); reach=[x|(1<<v) for v,x in enumerate(a)]
    for v in range(n):
        for u in range(n):
            if (reach[u]>>v)&1: reach[u] |= reach[v]
    unseen=set(range(n)); result=[]
    while unseen:
        u=min(unseen)
        c=[v for v in unseen if (reach[u]>>v)&1 and (reach[v]>>u)&1]
        unseen.difference_update(c); result.append(c)
    return result

def induced(a, vs):
    return [sum(((a[u]>>v)&1)<<j for j,v in enumerate(vs)) for u in vs]

def chain(parts):
    n=1+sum(parts); a=[sum(1<<j for j in range(i+1,n)) for i in range(n)]
    start=0
    for p in parts:
        end=start+p; a[start]^=1<<end; a[end]|=1<<start; start=end
    return a

def main():
    parser=argparse.ArgumentParser(); parser.add_argument('--out'); args=parser.parse_args()
    demand((L-2)**2 <= B < D < L*L, 'batch sandwich')
    demand(5**22 >= D**3, 'strong-floor induction multiplier')
    base_checks={str(n): floor_strong(n)**22 >= D**n for n in (23,24,25)}
    demand(all(base_checks.values()), 'strong-floor bases')
    for n in range(23,301):
        demand(floor_strong(n)**22 >= D**n, ('floor extended',n))
    allowed=[h for h in range(3,L,2) if h not in (7,21)]
    for h in allowed:
        # Strictly n-1 < 22 log_D(h), with no floating arithmetic.
        demand(D**(small_cost(h)-1) < h**22, ('residual',h))
    # All possible products at the crossing step, even without realizability.
    crossing=0
    for a in [1]+allowed:
        for b in allowed:
            if a*b>=L:
                crossing += 1
                demand(L<=a*b<=B and D**11 <= (a*b)**22, ('crossing',a,b))
    # Reproduce all labelled tournaments <=6; independent permutation DP <=5.
    strong_values={}; actual_min={1:0}; count=0; literal_count=0
    for n in range(1,7):
        values=set()
        for code in range(1<<(n*(n-1)//2)):
            a=adjacency(n,code); h=hcount(a); cs=sccs(a); count+=1
            actual_min[h]=min(actual_min.get(h,math.inf),n)
            if n<=5:
                demand(h==literal(a), ('literal',n,code)); literal_count+=1
            demand(h%2==1 and h not in (7,21), ('parity/gap',n,code))
            demand(h==math.prod(hcount(induced(a,c)) for c in cs),'SCC multiplication')
            demand(cutoff(h)>=min(n,small_cost(h)) if 1<h<L else True,'small cutoff')
            if len(cs)==1 and n>=3:
                values.add(h); demand(h>=floor_strong(n),'strong floor')
        strong_values[str(n)]=sorted(values)
    # Verify SCC-cost closure against exact all-tournament order supports <=6.
    sigma={h:n for n in range(6,2,-1) for h in strong_values[str(n)]}
    mu={1:0}
    for h in range(3,46,2):
        candidates=[sigma.get(h,math.inf)]
        for d in range(3,h,2):
            if h%d==0: candidates.append(mu.get(d,math.inf)+mu.get(h//d,math.inf))
        mu[h]=min(candidates)
    for h in range(1,46,2):
        value=mu.get(h,math.inf)
        demand((value if value<=6 else math.inf)==actual_min.get(h,math.inf), ('tropical replay',h))
    demand(mu[9]==5 and mu[3]+mu[3]==6, 'cost-primitive hostile')
    demand(mu[15]==5 and mu[3]+mu[5]==7, '15 cost primitive')
    # Actual chain controls: no claim that scalar repacking preserves inclusion.
    controls=[]
    for k in range(1,7):
        a=chain([2]*k)
        demand(hcount(a)==3**k and len(sccs(a))==1, ('Moon-chain',k))
        controls.append({'count':3**k,'ordered_triangle_core_order':3*k,
                         'strong_chain_order':len(a),'cutoff':cutoff(3**k)})
    # Actual n<=6 factor bank and Moon-minimum witnesses at arbitrary n.
    bank=[(h,int(n)) for n,vals in strong_values.items() for h in vals]
    bank += [(floor_strong(n),n) for n in range(7,61)]
    rng=random.Random(136)
    tested=0
    for length in range(1,61):
        for repeat in range(20):
            fs=[rng.choice(bank) for _ in range(length)]
            certify_repack(fs); tested+=1
    # Two-value relative-mesh amplifier, with exact arithmetic.
    mesh=[]
    for a,b in [(3,5),(9,11),(65,67)]:
        m0=0
        while b**m0 < a**(m0+1): m0+=1
        checked=0
        for m in range(m0,m0+11):
            for fraction in range(11):
                x=a**m+(a**(m+1)-a**m-1)*fraction//10
                y=a**m; j=0
                while y<x:
                    y=y//a*b; j+=1
                demand(j<=m and x<=y and a*y<b*x, ('relative mesh',a,b,m,x,y))
                checked+=1
        mesh.append({'a':a,'b':b,'m0':m0,'integer_samples':checked})
    demand(mesh[-1]['m0']==138, '65/67 threshold')
    # Exact boundary hostile: increasing the batch threshold to287 makes
    # (285,285) a possible crossing >B, losing the inherited bank guarantee.
    demand(285**2>B, 'too-large threshold hostile')
    receipt={'status':'PASS','scope':'New integer inequalities, recipe packing, and labelled n<=6 controls; large witness bank inherited, not regenerated',
      'constants':{'B':B,'D':D,'L':L,'coefficient_display_only':22/math.log(D),
                   'old_core_coefficient_display_only':3/math.log(3),
                   'batch_max':(L-2)**2},
      'floor_base_checks':base_checks,'floor_orders_checked':[23,300],
      'residual_values_checked':len(allowed),'crossings_checked':crossing,
      'labelled_tournaments':count,'literal_controls':literal_count,
      'strong_values_n_le_6':strong_values,'minimum_order_convolution_values_checked':23,
      'relative_mesh_controls':mesh,'random_actual_factor_recipes':tested,
      'scope_hostile_triangle_chains':controls,
      'examples':[certify_repack([(3,3)]*20),certify_repack([(5,4)]*20),
                  certify_repack([(3,3),(5,4),(11,5),(283,8),(floor_strong(30),30)])],
      'claims_not_made':['Every allowed odd h is attained','Output witness is induced in input',
                          'Output is strong','Universal cutoff is exact minimal order',
                          'New proof of inherited bank','Literature priority']}
    text=json.dumps(receipt,sort_keys=True,indent=2)+'\n'
    if args.out: open(args.out,'w').write(text)
    else: print(text,end='')
if __name__=='__main__': main()
