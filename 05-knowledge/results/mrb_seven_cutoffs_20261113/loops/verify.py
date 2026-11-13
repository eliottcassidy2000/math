"""Independent exact replay: full uniform-context TT_n profiles modulo q."""
from math import factorial, comb, lcm
import hashlib, json

def check(ok, msg):
    if not ok: raise RuntimeError(msg)

def factors(q):
    out={};p=2
    while p*p<=q:
        while q%p==0:out[p]=out.get(p,0)+1;q//=p
        p+=1
    if q>1:out[q]=out.get(q,0)+1
    return out

def cutoff(q):
    k=1;w=1
    while w%q:
        k+=1;w*=k*(k-1)
    return k

def predicted(q):
    fs=factors(q); periods=[]
    for p,a in fs.items():
        if p==2:periods.append(1 if a<=2 else 2**max(1,a-3))
        else:periods.append((p-1)*p**(a-1))
    return max(fs.values()),lcm(*periods)

def next_state(state,q):
    return tuple((k*state[k-1]+(k*(k-1)*state[k-2] if k>1 else 0))%q for k in range(1,len(state)+1))

def direct(n,q,kmax):
    return tuple(factorial(k-1)*sum((-1)**(k-j)*comb(k,j)*pow(j,n,q) for j in range(k+1))%q for k in range(1,kmax+1))

def orbit(q):
    k=cutoff(q);s=(1,)+(0,)*(k-2);seen={};n=1;states=[]
    while s not in seen:
        seen[s]=n;states.append(s);s=next_state(s,q);n+=1
    got=(seen[s],n-seen[s]);check(got==predicted(q),f'orbit {q}: {got} vs {predicted(q)}')
    for n in sorted({1,2,3,got[0],got[0]+got[1]-1,got[0]+got[1]}):
        idx=n-1 if n<=len(states) else got[0]-1+(n-got[0])%got[1]
        check(direct(n,q,k-1)==states[idx],f'finite-difference {q} n={n}')
    return {'q':q,'K':k,'first_periodic_n':got[0],'least_period':got[1],'states':len(states)}

def main():
    universe=set(range(2,257))|{2**a for a in range(1,15)}|{3**a for a in range(1,8)}|{5**a for a in range(1,6)}|{7**a for a in range(1,5)}|{11**a for a in range(1,4)}|{14**2,14**3,77**2}
    rows=[orbit(q) for q in sorted(universe)]
    primes=[p for p in range(2,98) if len(factors(p))==1 and factors(p).get(p)==1]
    wheel=lcm(*(p-1 for p in primes))
    check(len(primes)==25,'25 primes')
    selected=[r for r in rows if r['q'] in {2,4,8,16,32,64,14,77,49,121,196,5929}]
    hostile={'mod14_gate_counts':[3,1+2**7], 'mod77_gate_counts':[3,1+2**31], 'scalar_H':1,'first_profile_coordinates':'D1=1; D2=2^n-2'}
    check((129-3)%14==0 and 129!=3,'14 hostile')
    check((2**31+1-3)%77==0,'77 hostile')
    # Explicit direct path enumeration for the gate of TT_n, n=1..6.
    from itertools import permutations
    gate=[]
    for n in range(1,7):
        N=n+2;a=n;b=n+1
        edges={(i,j) for i in range(n) for j in range(i+1,n)}|{(a,b)}|{(b,i) for i in range(n)}|{(i,a) for i in range(n)}
        h=sum(all((perm[i],perm[i+1]) in edges for i in range(N-1)) for perm in permutations(range(N)))
        check(h==1+2**n,f'gate H n={n}');gate.append(h)
    # Independent physical inflation check for the finite search-box theorem.
    from itertools import product
    from functools import lru_cache
    def h_inflation(sizes):
        labels=tuple(i for i,size in enumerate(sizes) for _ in range(size));N=len(labels)
        pred=[sum(1<<u for u in range(N) if u!=v and ((labels[u]==labels[v] and u<v) or (labels[v]-labels[u])%3==1)) for v in range(N)]
        @lru_cache(None)
        def dp(mask,v):
            rest=mask^(1<<v)
            if not rest:return 1
            prev=rest&pred[v];total=0
            while prev:
                bit=prev&-prev;prev-=bit;total+=dp(rest,bit.bit_length()-1)
            return total
        return sum(dp((1<<N)-1,v) for v in range(N))
    inflation={sizes:h_inflation(sizes) for sizes in product(range(1,5),repeat=3)}
    residues={}
    for q in [3,4]:
        tau,L=predicted(q)
        def red(n):return n if n<tau else tau+(n-tau)%L
        for sizes,h in inflation.items():check(h%q==inflation[tuple(map(red,sizes))]%q,f'inflation q{q} sizes{sizes}')
        residues[q]=sorted({h%q for h in inflation.values()})
    small=[]
    for sizes in product((1,2),repeat=3):
        labels=tuple(i for i,size in enumerate(sizes) for _ in range(size));N=len(labels)
        edges={(u,v) for u in range(N) for v in range(N) if u!=v and ((labels[u]==labels[v] and u<v) or (labels[v]-labels[u])%3==1)}
        h=sum(all((u,v) in edges for u,v in zip(perm,perm[1:])) for perm in permutations(range(N)))
        check(h==inflation[sizes],f'independent C3 count {sizes}')
        small.append({'sizes':sizes,'H':h})
    for sizes,h in inflation.items():
        check(h%3==(0 if len({n%2 for n in sizes})==1 else 2),f'parity cubic {sizes}')
        check(h%4==(3-2*sum(n>1 for n in sizes))%4,f'singleton cutoff {sizes}')
    thirteen_rows=(18,8,3,5,14)
    thirteen=sum(all(thirteen_rows[u]>>v&1 for u,v in zip(perm,perm[1:])) for perm in permutations(range(5)))
    check(thirteen==13,'five-vertex escape from triangle inflation residue obstruction')
    # False stronger period phi(2^a): it overstates actual full-profile period.
    check(predicted(4)[1]==1 and predicted(8)[1]==2 and predicted(16)[1]==2,'two-adic boundary')
    print(json.dumps({'status':'PASS','moduli_tested':len(rows),'universe':'2..256 plus listed prime powers and 14^2,14^3,77^2','selected':selected,'first25_primes':primes,'first25_full_context_period':wheel,'first25_run_cutoff':97,'finite_orbit_sha256':hashlib.sha256(json.dumps(rows,sort_keys=True).encode()).hexdigest(),'gate_counts_direct':gate,'C3_inflations_checked':len(inflation),'C3_inflation_residues':residues,'C3_small_representatives':small,'H13_escape_rows':thirteen_rows,'hostile':hostile},indent=2,sort_keys=True))

if __name__=='__main__':main()
