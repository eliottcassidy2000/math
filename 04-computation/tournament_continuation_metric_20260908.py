"""Exact audit of homogeneous-pair contacts and pinned quotient envelopes.

All checks remain active under -O. Canonical references: THM-2221/2249/2256.
No full unpinned reversal kernel, H>=disc result, or LRC transfer is asserted.
"""
from itertools import permutations, product, combinations_with_replacement
from collections import Counter
from hashlib import sha256


def check(ok, label):
    if not ok:
        raise RuntimeError(label)


def tournament(n, mask):
    A = [[0]*n for _ in range(n)]
    bit = 0
    for i in range(n):
        for j in range(i+1,n):
            A[i][j] = mask>>bit&1
            A[j][i] = 1-A[i][j]
            bit += 1
    return tuple(map(tuple,A))


def interaction(A,s,t):
    n=len(A)
    return sum(A[i][k] and A[t[k]][s[i]] for i in range(n) for k in range(n))


def mismatch(A,p):
    return sum(A[i][j] != A[p[i]][p[j]] for i in range(len(A)) for j in range(i+1,len(A)))


def separation(A,u,v):
    if u==v:
        return 0
    return sum(A[u][x] != A[v][x] for x in range(len(A)) if x not in (u,v))


def pairs(A):
    return tuple((u,v) for u in range(len(A)) for v in range(u+1,len(A))
                 if separation(A,u,v)==0)


def matchings(edges):
    ans=[()]
    for edge in edges:
        ans += [m+(edge,) for m in ans if not any(set(edge)&set(e) for e in m)]
    return tuple(ans)


def matching_perm(n, matching):
    p=list(range(n))
    for u,v in matching:
        p[u],p[v]=p[v],p[u]
    return tuple(p)


def components(n, edges):
    graph=[set() for _ in range(n)]
    for u,v in edges:
        graph[u].add(v)
        graph[v].add(u)
    check(all(len(row)<=2 for row in graph), "homogeneous-pair degree <=2")
    remaining=set(range(n))
    ans=[]
    while remaining:
        seed=min(remaining)
        stack=[seed]
        comp=set()
        while stack:
            v=stack.pop()
            if v in comp:
                continue
            comp.add(v)
            stack.extend(graph[v]-comp)
        remaining-=comp
        check(sum(len(graph[v]) for v in comp)==2*(len(comp)-1), "homogeneous-pair component is path")
        ans.append(tuple(sorted(comp)))
    return tuple(ans)


def fibonacci(n):
    a,b=0,1
    for _ in range(n):
        a,b=b,a+b
    return a


def strong_components(A):
    n=len(A)
    reach=[set() for _ in range(n)]
    for start in range(n):
        stack=[start]
        while stack:
            v=stack.pop()
            if v in reach[start]:
                continue
            reach[start].add(v)
            stack.extend(u for u in range(n) if A[v][u] and u not in reach[start])
    unseen=set(range(n))
    comps=[]
    while unseen:
        v=min(unseen)
        comp={u for u in unseen if v in reach[u] and u in reach[v]}
        comps.append(tuple(sorted(comp)))
        unseen-=comp
    return tuple(sorted(comps,key=lambda comp:-sum(A[comp[0]][u] for u in range(n) if u not in comp)))


def classify_contacts():
    ledger=[]
    checked=0
    all_contact_pairs=0
    for n in range(2,6):
        sym=tuple(permutations(range(n)))
        identity=tuple(range(n))
        bounded=linear=contacts=0
        structural=Counter()
        for mask in range(1<<(n*(n-1)//2)):
            A=tournament(n,mask)
            edges=pairs(A)
            scc=strong_components(A)
            expected=set()
            for block in scc:
                induced=tuple(tuple(A[u][v] for v in block) for u in block)
                expected.update(tuple(sorted((block[u],block[v]))) for u,v in pairs(induced))
            expected.update(tuple(sorted((a[0],b[0]))) for a,b in zip(scc,scc[1:]) if len(a)==len(b)==1)
            check(set(edges)==expected,"strong-component homogeneous-pair decomposition")
            comps=components(n,edges)
            # Every component is a transitive module; graph edges are consecutive.
            for comp in comps:
                order=sorted(comp,key=lambda u:-sum(A[u][v] for v in comp))
                check(all(A[order[i]][order[j]] for i in range(len(order)) for j in range(i+1,len(order))), "component transitive")
                check(all(A[x][u]==A[x][comp[0]] for x in range(n) if x not in comp for u in comp), "component is module")
            matching_words=matchings(edges)
            predicted={matching_perm(n,m):len(m) for m in matching_words}
            fib=1
            for comp in comps:
                fib*=fibonacci(len(comp)+1)
            check(len(predicted)==fib, "Fibonacci product contact count")
            autos=[s for s in sym if mismatch(A,s)==0]
            auto_set=set(autos)
            actual=set()
            for tau in sym:
                w=interaction(A,identity,tau)+interaction(A,tau,identity)
                check(w==sum(separation(A,v,tau[v]) for v in range(n)), "contact as external row separation")
                check(all(sum(A[u][x]!=A[v][x] for x in range(n))==separation(A,u,v)+int(u!=v)
                          for u in range(n) for v in range(n)), "restored endpoint coordinate gives full row Hamming metric")
                if w==0:
                    actual.add(tau)
                    check(mismatch(A,tau)==predicted.get(tau,-1), "matching diagonal equals swap count")
                checked+=1
            check(actual==set(predicted), "all zero contacts classified")
            if edges:
                bounded+=1
            else:
                linear+=1
            structural['strong' if len(scc)==1 else 'nonstrong','bounded' if edges else 'linear']+=1
            # Independent all-auto contact enumeration at n<=4, and all n5
            # automorphism/nonautomorphism pairs; no quotient representatives.
            count=0
            for sigma in autos:
                inverse=[0]*n
                for i,j in enumerate(sigma):
                    inverse[j]=i
                for tau in sym:
                    w=interaction(A,sigma,tau)+interaction(A,tau,sigma)
                    rho=tuple(inverse[tau[v]] for v in range(n))
                    check((w==0)==(rho in predicted), "automorphism-relative contact classification")
                    if tau in auto_set and tau!=sigma:
                        moved=sum(v!=rho[v] for v in range(n))
                        check(w>=moved>=3, "sharp automorphism contact floor three")
                    if tau!=sigma and w==0:
                        check(mismatch(A,tau)>0, "nonempty matching is nonautomorphism")
                        count+=1
            check(count==len(autos)*(fib-1), "complete contact-pair count")
            contacts+=count
            all_contact_pairs+=count
            ledger.append((n,mask,len(edges),tuple(map(len,comps)),count))
        print(f"CONTACT n={n}: bounded={bounded}, linear={linear}, ordered_zero_contacts={contacts} PASS")
        print(f"CONTACT n={n} strong-component census={dict(sorted(structural.items()))}")
    semantic="\n".join(map(str,ledger)).encode()
    print(f"CONTACT permutation_probes={checked}; total_ordered_zero_contacts={all_contact_pairs}")
    print(f"CONTACT semantic_sha256={sha256(semantic).hexdigest()}")


def cut(word, mass=1):
    return tuple(tuple(mass*int(a!=b) for b in word) for a in word)


def hall_minimum(A,D,N):
    n=len(A)
    sym=tuple(permutations(range(n)))
    diag=tuple(mismatch(A,s) for s in sym)
    ext=tuple(sum(D[i][s[i]] for i in range(n)) for s in sym)
    W=tuple(tuple(interaction(A,s,t)+interaction(A,t,s) for t in sym) for s in sym)
    free={i for i in range(len(sym)) if diag[i]==ext[i]==0}
    best=None
    witness=None
    count=0
    for layers in combinations_with_replacement(range(len(sym)),N):
        if layers[0]==layers[-1] and layers[0] in free:
            continue
        frequencies=Counter(layers)
        active=sorted(frequencies)
        value=sum(diag[i]*frequencies[i]**2+ext[i]*frequencies[i] for i in active)
        value+=sum(W[i][j]*frequencies[i]*frequencies[j]
                   for pos,i in enumerate(active) for j in active[pos+1:])
        if best is None or value<best:
            best=value
            witness=tuple((sym[i],frequencies[i]) for i in active)
        count+=1
    check(best is not None,"nonfree transport exists")
    return best,witness,count


def pin_controls():
    check(len(strong_components(tournament(4,43)))==1 and pairs(tournament(4,43))==((0,1),),
          "strong bounded-contact hostile")
    check(len(strong_components(tournament(4,61)))==2 and pairs(tournament(4,61))==(),
          "nonstrong linear-contact hostile")
    print("SCC hostile: C3[T2,1,1] is strong and bounded; C3 join1 is nonstrong and linear PASS")
    # T3 with one alternating pinned observer: both homogeneous pairs are
    # distinguished, but one-vertex exchanges cost exactly three at every N.
    A=tournament(3,7)
    D=cut((0,1,0))
    check(pairs(A)==((0,1),(1,2)), "transitive-three homogeneous pairs")
    check(all(D[u][v]==1 for u,v in pairs(A)), "one pin separates every homogeneous pair")
    identity=(0,1,2)
    p=(1,0,2)
    q=(0,2,1)
    qp=tuple(q[p[i]] for i in range(3))
    check(interaction(A,identity,p)+interaction(A,p,identity)==0, "first zero-contact generator")
    check(interaction(A,identity,q)+interaction(A,q,identity)==0, "second zero-contact generator")
    check(interaction(A,identity,qp)+interaction(A,qp,identity)==1, "zero contacts not closed under composition")
    check(separation(A,0,1)==separation(A,1,2)==0 and separation(A,0,2)==1,
          "pair-dependent deletion is not a pseudometric")
    print("CONTACT hostile: T3 adjacent zero-contact swaps compose to contact cost1; pairwise zero relation is not transitive PASS")
    tested=0
    for N in range(1,8):
        minimum,_,count=hall_minimum(A,D,N)
        check(minimum==3, "single pin does not create macroscopic floor")
        tested+=count
    print("PIN hostile: transitive3 + alternating singleton pin separates both pairs, yet floor=3 for N=1..7 PASS")
    # Complete order-three quotient/context bank, including constant pins and
    # arbitrary multiplicity two, tests min(N,c) and exactness at N>=c.
    bank=Counter()
    for mask in range(8):
        A=tournament(3,mask)
        edges=pairs(A)
        for word in product((0,1),repeat=3):
            for mass in (0,1,2):
                D=cut(word,mass)
                c=min((1+2*D[u][v] for u,v in edges),default=None)
                for N in range(1,8):
                    value,_,count=hall_minimum(A,D,N)
                    tested+=count
                    if c is not None:
                        check(value<=c,"homogeneous swap pinned upper bound")
                        check(value>=min(N,c),"uniform pinned lower bound")
                        if N>=c:
                            check(value==c,"eventual pinned exact floor")
                        bank['bounded']+=1
                    else:
                        check(value>=N,"homogeneous-pair-free linear floor")
                        bank['linear']+=1
    print(f"PIN complete n=3 bank: masks=8, words=8, multiplicities=0,1,2, N=1..7; cases={dict(bank)}")
    print(f"PIN total_nonfree_Hall_multisets={tested} PASS")
    # The scale-bearing repair: replicate the alternating pin N times.
    values=[]
    for N in range(1,8):
        value,_,_=hall_minimum(tournament(3,7),cut((0,1,0),N),N)
        values.append(value)
    print(f"PIN growing-mass control: transitive3, pin multiplicity=N, minima N=1..7={values}")


if __name__=='__main__':
    classify_contacts()
    pin_controls()
    print("PASS: quotient forced envelope only; internal response and unpinned exchanges retained as separate obligations")
