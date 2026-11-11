"""Exact all-order series evaluations and a separate small-class birth audit."""
from collections import Counter
from pathlib import Path
from math import factorial, gcd, prod
import csv, gzip, hashlib, json, time

HERE = Path(__file__).resolve().parent
PREV = HERE.parent / 'metagraph_aspect_sequences_20261111'
LIMIT = 64

def require(ok, msg):
    if not ok:
        raise RuntimeError(msg)

def mul(a, b):
    return [sum(a[j]*b[i-j] for j in range(i+1)) for i in range(len(a))]

def inv(a):
    require(a[0] == 1, 'series inverse constant')
    b = [1] + [0]*(len(a)-1)
    for i in range(1, len(a)):
        b[i] = -sum(a[j]*b[i-j] for j in range(1, i+1))
    return b

def parts(n, lo=1, prefix=(), exponent=0):
    if n == 0:
        yield prefix, exponent
    for p in range(lo, n+1, 2):
        yield from parts(n-p, p, prefix+(p,), exponent+(p-1)//2+sum(gcd(p,q) for q in prefix))

def burnside(n):
    f = factorial(n)
    a = q = 0
    for pp, e in parts(n):
        z = prod(p**m*factorial(m) for p,m in Counter(pp).items())
        mass = (f//z) << e
        a += mass
        q += mass << len(pp)
    require(a % f == q % f == 0, 'Burnside integrality')
    return a//f, q//f

def scc(rows):
    n = len(rows)
    reach = [rows[i] | (1<<i) for i in range(n)]
    for k in range(n):
        for i in range(n):
            if reach[i]>>k & 1:
                reach[i] |= reach[k]
    todo = (1<<n)-1
    components = []
    while todo:
        bit = todo & -todo
        i = bit.bit_length()-1
        c = sum(1<<j for j in range(n) if reach[i]>>j & 1 and reach[j]>>i & 1)
        components.append(c)
        todo &= ~c
    return sorted(components, key=lambda c: -reach[(c & -c).bit_length()-1].bit_count())

def hcount(rows, mask=None):
    if mask is None:
        mask = (1<<len(rows))-1
    vs = [i for i in range(len(rows)) if mask>>i & 1]
    n = len(vs)
    if n == 0:
        return 1
    dp = [[0]*n for _ in range(1<<n)]
    for j in range(n):
        dp[1<<j][j] = 1
    for s in range(1, 1<<n):
        for j in range(n):
            if s>>j & 1:
                old = s ^ (1<<j)
                if old:
                    dp[s][j] = sum(dp[old][i] for i in range(n) if rows[vs[i]]>>vs[j] & 1)
    return sum(dp[-1])

def child(rows, out):
    n = len(rows)
    return [r | (0 if out>>i & 1 else 1<<n) for i,r in enumerate(rows)] + [out]

def interval(comps, out):
    left = next((i for i,c in enumerate(comps) if out & c), len(comps))
    right = next((i for i in range(len(comps)-1,-1,-1) if comps[i] & ~out), -1)
    return (left, right) if left <= right else None

def submasks(mask):
    s = mask
    while True:
        yield s
        if s == 0:
            break
        s = (s-1) & mask

def run():
    start = time.perf_counter()
    aq = [burnside(n) for n in range(LIMIT+1)]
    A, Q = [list(t) for t in zip(*aq)]
    Ai, Qi = inv(A), inv(Q)
    S, R = [-v for v in Ai], [-v for v in Qi]
    S[0] = R[0] = 0
    A2 = mul(A,A)
    mixed = [r-2*s for r,s in zip(R,S)]
    nonempty = [r-s for r,s in zip(R,S)]
    widths = [A2, mul(A2,mixed), mul(A2,mul(nonempty,nonempty))]
    for k in range(3,LIMIT+1):
        widths.append(mul(widths[-1],R))
    require(all(sum(w[n] for w in widths) == Q[n] for n in range(LIMIT+1)), 'width partition')
    rooted_strong = mul(Q,mul(Ai,Ai))
    rooted_strong[0] = 0
    with (HERE/'birth_sequences.tsv').open('w', newline='', encoding='utf-8') as f:
        w = csv.writer(f,delimiter='\t')
        w.writerow(['source_n','A000568','subset_orbits_Q','strong_classes_S','strong_subset_orbits_R',
                    'strong_pointed_children','singleton_births','one_component_births','multi_component_births'])
        for n in range(LIMIT+1):
            w.writerow([n,A[n],Q[n],S[n],R[n],rooted_strong[n],widths[0][n],widths[1][n],
                        sum(w[n] for w in widths[2:])])
    (HERE/'birth_width_triangle.json').write_text(json.dumps(
        {str(n):[widths[k][n] for k in range(n+1)] for n in range(LIMIT+1)},indent=2)+'\n',encoding='utf-8')
    observed = {n:Counter() for n in range(1,9)}
    obs_strong = Counter()
    obs_s = Counter()
    obs_r = Counter()
    nclasses = Counter()
    direct_scc = direct_h = factor_checks = ternary_checks = orbits = 0
    witnesses = {}
    kernels = {}
    core_witnesses = []
    transitions = {}
    by_h = {1:[1]+[0]*8}
    with gzip.open(PREV/'records.jsonl.gz','rt',encoding='utf-8') as f:
        for line in f:
            r = json.loads(line)
            n, rows = r['n'], r['rows']
            comps = scc(rows)
            ch = [hcount(rows,c) for c in comps]
            require(prod(ch)==r['H'], 'component H product')
            nclasses[n] += 1
            by_h.setdefault(r['H'],[0]*9)[n] += 1
            orbit = r['attachment_orbit']
            alphas = [len({orbit[m] for m in submasks(c)}) for c in comps]
            require(prod(alphas)==len(set(orbit)), 'component orbit product')
            if len(comps)==1:
                obs_s[n] += 1
                obs_r[n] += len(set(orbit))
            local = Counter({0:len(comps)+1, 1:sum(a-2 for a in alphas)})
            for i in range(len(comps)):
                for j in range(i+1,len(comps)):
                    local[j-i+1] += (alphas[i]-1)*(alphas[j]-1)*prod(alphas[i+1:j])
            actual = Counter()
            for out in range(1<<n):
                if orbit[out] != out:
                    continue
                orbits += 1
                span = interval(comps,out)
                width = 0 if span is None else span[1]-span[0]+1
                actual[width] += 1
                active = 0 if span is None else sum(comps[span[0]:span[1]+1])
                outside = prod(ch) if span is None else prod(ch[:span[0]])*prod(ch[span[1]+1:])
                require(r['extension'][out] % outside == 0, 'outside divisor')
                if r['extension'][out] in (3*r['H']+2, 3*r['H']+4):
                    ternary_checks += 1
                    require(outside == 1, 'ternary outside singleton restriction')
                if width == len(comps) and span is not None:
                    obs_strong[n] += 1
                    g = r['extension'][out]
                    if g < 125:
                        kernels.setdefault((r['H'],g),[0]*9)[n] += 1
                        core_witnesses.append([n,r['rep'],out,r['H'],g,','.join(map(str,rows))])
                g = r['extension'][out]
                if g < 125:
                    transitions.setdefault((r['H'],g),[0]*9)[n] += 1
                if n <= 7:
                    cc = scc(child(rows,out))
                    containing = next(c for c in cc if c>>n & 1)
                    require(containing == active | (1<<n), 'literal SCC fusion')
                    direct_scc += 1
                if n <= 6:
                    hc = hcount(child(rows,out))
                    require(hc == r['extension'][out], 'independent child H')
                    require(hc == outside*hcount(child(rows,out),active | (1<<n)), 'literal factorization')
                    direct_h += 1
                    factor_checks += 1
                if n <= 7 and width >= 2 and len(witnesses) < 12:
                    witnesses.setdefault(str(width), dict(n=n,rep=r['rep'],rows=rows,out=out,
                        source_H=r['H'],target_H=r['extension'][out],components=comps,active=active))
            require(actual == +local, 'source width orbit partition')
            observed[n].update(actual)
    for n in range(1,9):
        require(nclasses[n]==A[n], 'complete class count')
        require([observed[n][k] for k in range(n+1)]==[w[n] for w in widths[:n+1]], 'independent width census')
        require(obs_strong[n]==rooted_strong[n] and obs_s[n]==S[n] and obs_r[n]==R[n], 'strong series census')
    # Independent assembly of the fixed-(source H,target H) transfer identity.
    divisors = lambda m: [d for d in range(1,m+1) if m%d==0]
    zero = [0]*9
    D = {d:[sum(mul(by_h.get(a,zero),by_h.get(d//a,zero))[n]
                   for a in divisors(d)) for n in range(9)] for d in range(1,125,2)}
    transition_checks = 0
    for h in range(1,125,2):
        for g in range(h,125,2):
            predicted = list(D[h]) if h==g else list(zero)
            for d in divisors(gcd(h,g)):
                term = mul(D[d],kernels.get((h//d,g//d),zero))
                predicted = [a+b for a,b in zip(predicted,term)]
            actual = transitions.get((h,g),zero)
            require(predicted[1:]==actual[1:], 'fixed H transfer identity '+str((h,g)))
            transition_checks += 8
    exported = {str(h)+','+str(g):c for (h,g),c in sorted(kernels.items())}
    with (HERE/'primitive_births.tsv').open('w',newline='',encoding='utf-8') as f:
        w=csv.writer(f,delimiter='\t')
        w.writerow(['source_n','source_rep','new_vertex_out_mask','source_H','target_H','source_rows'])
        w.writerows(core_witnesses)
    (HERE/'strong_birth_kernels_below125.json').write_text(json.dumps(
        {'status':'FINITE-EXACT complete at every order using CITED Busch f(10)=125',
         'coefficient_n':'source order; child order n+1; ordinary pointed classes',
         'kernels':exported},indent=2)+'\n',encoding='utf-8')
    ternary = {str(h):{str(g):kernels.get((h,g),zero) for g in (3*h+2,3*h+4)}
               for h in range(1,40,2)}
    (HERE/'ternary_kernels_through39.json').write_text(json.dumps(ternary,indent=2)+'\n',encoding='utf-8')
    result = dict(status='PASS',series_max_order=LIMIT,ordinary_classes=sum(nclasses.values()),
        attachment_orbits=orbits,direct_child_SCC_checks=direct_scc,direct_child_H_checks=direct_h,
        H_factorization_checks=factor_checks,ternary_target_checks=ternary_checks,
        fixed_H_coefficient_checks=transition_checks,nonzero_strong_birth_kernels=len(kernels),
        input_sha256=hashlib.sha256((PREV/'records.jsonl.gz').read_bytes()).hexdigest(),
        seconds=round(time.perf_counter()-start,3),witnesses=witnesses)
    (HERE/'birth_checks.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in result.items() if k!='witnesses'}))

if __name__ == '__main__':
    run()
