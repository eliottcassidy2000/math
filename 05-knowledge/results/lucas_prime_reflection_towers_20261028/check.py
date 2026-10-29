#!/usr/bin/env python3
"""Exact controls, standard library only; no assertion-dependent checks."""
import itertools
import json
import math
import hashlib
from collections import Counter, deque

checks = 0


def ck(ok, message):
    global checks
    checks += 1
    if not ok:
        raise RuntimeError(message)


def fibluc(n):
    f, l = [0, 1], [2, 1]
    for _ in range(2, n + 1):
        f.append(f[-1] + f[-2])
        l.append(l[-1] + l[-2])
    return f, l


def prime(n):
    return n >= 2 and all(n % d for d in range(2, math.isqrt(n) + 1))


def comp(a, b):
    return tuple(a[x] for x in b)


def inv(a):
    b = [0] * len(a)
    for i, j in enumerate(a):
        b[j] = i
    return tuple(b)


def power(a, n):
    b = tuple(range(len(a)))
    while n:
        if n & 1:
            b = comp(b, a)
        a = comp(a, a)
        n //= 2
    return b


def root(a):
    seen, b = set(), [0] * len(a)
    for x in range(len(a)):
        if x in seen:
            continue
        cyc, y = [], x
        while y not in seen:
            seen.add(y)
            cyc.append(y)
            y = a[y]
        ck(len(cyc) % 2 == 1, 'root must have odd cycles')
        step = (len(cyc) + 1) // 2
        for i, v in enumerate(cyc):
            b[v] = cyc[(i + step) % len(cyc)]
    return tuple(b)


def closure(gens, n):
    ident = tuple(range(n))
    seen, todo = {ident}, deque([ident])
    while todo:
        x = todo.popleft()
        for g in gens:
            y = comp(g, x)
            if y not in seen:
                seen.add(y)
                todo.append(y)
    return seen


def words(n):
    return [w for w in range(1 << n)
            if not (w & (w << 1)) and not ((w & 1) and (w >> (n - 1)))]


def reflect(w, n, edge=False):
    return sum(((w >> i) & 1) << (((1 if edge else 0) - i) % n)
               for i in range(n))


def reflection_controls():
    f, l = fibluc(1002)
    for m in range(2, 500):
        ck(l[2*m+1] - f[m+2]**2 == 2*f[m]*f[m-1], 'odd fixed gap')
        n = 2*m+1
        if n % 3:
            ck(((l[n]-f[m+2])//2) % 2 == (n % 12 == 11), 'odd sign')
    rows = []
    for n in range(3, 20):
        ws = words(n)
        ck(len(ws) == l[n], 'Lucas word count')
        fixed = sum(reflect(w, n) == w for w in ws)
        ck(fixed == (f[(n+3)//2] if n % 2 else f[n//2+3]), 'vertex mirror')
        if n % 2 == 0:
            fe = sum(reflect(w, n, True) == w for w in ws)
            ck(fe == f[n//2], 'edge mirror')
        if n % 2 and prime(l[n]):
            p = l[n]
            a = next(x for x in range(p) if (x*x-x-1) % p == 0 and pow(x,n,p)==1)
            labels = {w: sum(((w >> i) & 1)*pow(a,i,p) for i in range(n)) % p for w in ws}
            ck(len(set(labels.values())) == p, 'golden reader bijective')
            r = [0]*p
            for w in ws:
                r[labels[w]] = labels[reflect(w,n)]
            r = tuple(r)
            t = tuple((x+1) % p for x in range(p))
            ck(comp(r,r) == tuple(range(p)), 'reader involution')
            group = 'S' if ((p-fixed)//2) % 2 else 'A'
            if p == 11:
                ck(len(closure([t,r],p)) == 660, 'exceptional group')
                group = 'PSL(2,11)'
            if p == 29:
                g = tuple(range(p))
                for k in [19,17,22,23,18,1,18,21,6,1,15,17]:
                    g = comp(r,comp(power(t,k),g))
                witness = power(g,25)
                ck([i for i in range(p) if witness[i] != i] == [0,20,23], 'A29 3-cycle')
            rows.append(dict(n=n,p=p,fixed=fixed,golden_root=a,group=group))
    for k in range(3,10):
        m=2**(k-1)
        ck(l[2*m] == 5*f[m]**2+2, 'even gap')
        ck(((l[2*m]-f[m])//2) % 2 == (k % 2 == 0), 'even sign')
    ck(prime(l[8]) and prime(l[16]),'even Lucas-prime examples')
    # Hostile: one-fixed-point involution and a p-cycle can generate only D_p.
    for p in [7,13,29]:
        ck(len(closure([tuple((i+1)%p for i in range(p)), tuple((-i)%p for i in range(p))],p))==2*p,
           'dihedral boundary')
    return rows


def matmul(a,b,p):
    return tuple(sum(a[i*3+k]*b[k*3+j] for k in range(3)) % p
                 for i in range(3) for j in range(3))


def projective_controls():
    report=[]
    for q in [2,3]:
        pts=[]
        for v in itertools.product(range(q),repeat=3):
            if any(v) and next(x for x in v if x) == 1:
                pts.append(v)
        counts=Counter()
        for a in itertools.product(range(q),repeat=9):
            if not any(a) or next(x for x in a if x) != 1:
                continue
            sq=matmul(a,a,q)
            c=sq[0]
            if not c or sq != (c,0,0,0,c,0,0,0,c):
                continue
            if a == (1,0,0,0,1,0,0,0,1):
                continue
            fixed=0
            for v in pts:
                av=tuple(sum(a[3*i+j]*v[j] for j in range(3))%q for i in range(3))
                first=next(x for x in av if x)
                norm=tuple(x*pow(first,-1,q)%q for x in av)
                fixed += norm == v
            ck(fixed*fixed > len(pts), 'PGL3 involution bound')
            counts[fixed]+=1
        report.append(dict(q=q,degree=len(pts),involutions_by_fixed=dict(counts)))
    # All semilinear projective maps on P1(F4), including genuine field twists.
    def mul(a,b):
        c=0
        while b:
            if b&1: c ^= a
            b >>= 1
            a <<= 1
            if a&4: a ^= 7
        return c
    def reciprocal(a):
        return mul(a,a)
    perms=set()
    for a,b,c,d in itertools.product(range(4),repeat=4):
        if not (mul(a,d)^mul(b,c)):
            continue
        for sigma in [0,1]:
            vals=[]
            for x in range(5):
                if x==4:
                    num,den=a,c
                else:
                    z=mul(x,x) if sigma else x
                    num,den=mul(a,z)^b,mul(c,z)^d
                vals.append(mul(num,reciprocal(den)) if den else 4)
            perms.add(tuple(vals))
    ck(len(perms)==120,'PGLamma2(4) order')
    fixeds=Counter()
    for r in perms:
        if r != tuple(range(5)) and comp(r,r)==tuple(range(5)):
            f=sum(i==v for i,v in enumerate(r))
            ck(f==1 or f*f>5,'semilinear involution gap')
            fixeds[f]+=1
    for q in range(2,33):
        for k in range(1,9):
            count=sum(q**i for i in range(2*k+1))
            lower=sum(q**i for i in range(k+1))
            ck(lower**2>count,'odd dimension polynomial inequality')
            for a in range(2*k+2):
                f=sum(q**i for i in range(a))+sum(q**i for i in range(2*k+1-a))
                ck(f>=lower,'eigenspace count')
        for d in range(2,10):
            ck(sum(q**i for i in range(d))**2>sum(q**(2*i) for i in range(d)), 'semilinear square gap')
    return dict(linear=report,semilinear_degree5=dict(fixeds))


def qmul(x,y,p):
    a,b=x; c,d=y
    return ((a*c+b*d)%p,(a*d+b*c+b*d)%p)


def qpower(x,n,p):
    y=(1,0)
    while n:
        if n&1: y=qmul(x,y,p)
        x=qmul(x,x,p); n//=2
    return y


def factors(n):
    result=[]
    for d in range(2,math.isqrt(n)+1):
        if n%d==0:
            result.append(d)
            while n%d==0: n//=d
    if n>1: result.append(n)
    return result


def val(x,p):
    e=0
    while x and x%p==0:
        x//=p; e+=1
    return e


def arithmetic_controls():
    f,l=fibluc(5002)
    table={1:(1,1),5:(1,1),7:(-1,1),11:(-1,1),13:(-1,-1),17:(-1,-1),19:(1,-1),23:(1,-1)}
    bins=Counter(); rows=[]; eligible=0; lucas_divisors=0
    for p in range(3,5001,2):
        if not prime(p) or p==5:
            continue
        chi=1 if pow(5,(p-1)//2,p)==1 else -1
        d=p-chi
        theta=(-1%p,-1%p)
        ck(qpower(theta,d,p)==(1,0),'golden projective order bound')
        for r in factors(d):
            while d%r==0 and qpower(theta,d//r,p)==(1,0):
                d//=r
        # Independent Fibonacci recurrence for the exact projective order.
        a,b=0,1
        for rank in range(1,p-chi+1):
            a,b=b,(a+b)%p
            if a==0:
                break
        ck(rank==d,'rank versus finite algebra order')
        bins[val(d,2)]+=1; eligible+=1; lucas_divisors += d%2==0
        for n in range(1,121):
            predicted=d%2==0 and n%(d)==d//2
            ck((l[n]%p==0)==predicted,'Lucas midpoint divisibility')
            if predicted:
                h=d//2
                ck(val(l[n],p)==val(l[h],p)+val(n//h,p),'valuation lift')
                if n%2:
                    ck(chi==1,'odd Lucas divisors split in golden field')
        if p>3:
            s6=1 if pow(6,(p-1)//2,p)==1 else -1
            sm6=1 if pow(-6,(p-1)//2,p)==1 else -1
            ck((s6,sm6)==table[p%24],'mirror squareclass table')
        if p<=100:
            rows.append(dict(p=p,chi5=chi,projective_order=d,first_lucas_index=d//2 if d%2==0 else None))
    for ell in [5,7,11,13,17,19,23,29,31]:
        for r in factors(l[ell]):
            ck(prime(r) and r%(2*ell)==1,'prime-index Lucas factor congruence')
    # Golden projective half-turn versus the twelve-frame mirror at eleven.
    p,n=11,5
    a,b,c,d=f[n+1]%p,f[n]%p,f[n]%p,f[n-1]%p
    def mobius(A,x):
        a,b,c,d=A
        num,den=(a,c) if x==p else ((a*x+b)%p,(c*x+d)%p)
        return num*pow(den,-1,p)%p if den else p
    fixQ=[x for x in range(p+1) if mobius((a,b,c,d),x)==x]
    fix6=[x for x in range(p+1) if mobius((0,6,1,0),x)==x]
    ck(len(fixQ)==2 and len(fix6)==0,'two involution classes at eleven')
    group_rows=[]
    for p in [5,7,11,13]:
        t=tuple((x+1)%p for x in range(p))+(p,)
        j=tuple(p if x==0 else 6*pow(x,-1,p)%p for x in range(p))+(0,)
        actual=len(closure([t,j],p))
        square=pow(-6,(p-1)//2,p)==1
        expected=p*(p*p-1)//(2 if square else 1)
        ck(actual==expected,'uniform projective six-mirror group')
        group_rows.append(dict(p=p,order=actual,group='PSL' if square else 'PGL'))
    return dict(prime_bound=5000,eligible_primes=eligible,lucas_divisor_primes=lucas_divisors,
                rank_v2_counts=dict(sorted(bins.items())),small_rows=rows,
                eleven=dict(golden_half_turn_fixed=fixQ,frame_mirror_fixed=fix6),projective_groups=group_rows)


def wreath_controls(k):
    n=3**k
    digits=list(itertools.product(range(3),repeat=k))
    index={a:i for i,a in enumerate(digits)}
    t0=tuple(index[tuple((-a)%3 for a in d)] for d in digits)
    gens=[]
    for depth in range(k):
        g=[]
        for d in digits:
            e=list(d)
            if all(v==0 for v in d[:depth]): e[depth]=(e[depth]+1)%3
            g.append(index[tuple(e)])
        g=tuple(g); gens.append(g)
        ck(comp(t0,comp(g,t0))==inv(g),'zero-prefix generator in K')
    def forward(i,j):
        for a,b in zip(digits[i],digits[j]):
            if a!=b: return (b-a)%3==1
        return False
    for g in gens:
        for i in range(n):
            for j in range(i):
                ck(forward(i,j)==forward(g[i],g[j]),'native automorphism')
    marks={t0}; todo=deque([t0])
    while todo:
        s=todo.popleft()
        for a in gens:
            v=comp(a,comp(s,inv(a)))
            if v not in marks: marks.add(v);todo.append(v)
    cexp=(3**k-2*k-1)//4
    ck(len(marks)==n*3**cexp,'complete mark orbit count')
    bases={s:root(comp(s,t0)) for s in marks}
    baseinv={s:inv(b) for s,b in bases.items()}
    corrections=set()
    for a in gens:
        t=comp(a,comp(t0,inv(a)))
        ck(bases[t]==a,'generator is exact square-root transversal')
        for s in marks:
            v=comp(a,comp(s,inv(a)))
            correction=comp(baseinv[v],comp(a,bases[s]))
            direct=comp(baseinv[v],comp(root(comp(v,t)),bases[t]))
            ck(correction==direct,'Schreier generator is an actual triangle')
            ck(comp(correction,t0)==comp(t0,correction),'correction centralizes mark')
            corrections.add(correction)
    h=closure(list(corrections),n)
    ck(len(h)==3**cexp,'full centralizer holonomy')
    # Decompose every retained correction by native wreath coordinates, then
    # replace its zero-prefix generator word by actual based transport triangles.
    def local_coordinates(g,path=()):
        if len(g)==1: return []
        size=len(g)//3; shift=g[0]//size
        result=[(path,shift)] if shift else []
        for block in range(3):
            child=tuple(g[block*size+j]%size for j in range(size))
            result+=local_coordinates(child,path+(block,))
        return result
    signed={(i,1):g for i,g in enumerate(gens)}
    signed.update({(i,-1):inv(g) for i,g in enumerate(gens)})
    longest=0
    for target in h:
        word=[]
        for path,r in local_coordinates(target):
            prefix=[(i,1 if digit==1 else -1) for i,digit in enumerate(path) if digit]
            word+=prefix+[(len(path),1 if r==1 else -1)]+[(i,-v) for i,v in reversed(prefix)]
        longest=max(longest,len(word))
        ck(len(word)<=1+(k-1)*n,'uniform generator length bound')
        value=tuple(range(n))
        for letter in word: value=comp(value,signed[letter])
        ck(value==target,'native wreath word reconstructs correction')
        s=t0; transported=tuple(range(n))
        for letter in reversed(word):
            a=signed[letter]; v=comp(a,comp(s,inv(a))); t=comp(a,comp(t0,inv(a)))
            triangle=comp(baseinv[v],comp(root(comp(v,t)),bases[t]))
            transported=comp(triangle,transported);s=v
        ck(s==t0 and transported==target,'paid marking loop realizes exact correction')
    remaining=set(range(n)); sizes=[]
    while remaining:
        x=min(remaining); orb={g[x] for g in h}
        sizes.append(len(orb)); remaining-=orb
        signatures={next(((j,d) for j,d in enumerate(digits[y]) if d),(k,0)) for y in orb}
        ck(len(signatures)==1,'first nonzero digit invariant')
    expected=sorted([1]+[3**j for j in range(k) for _ in range(2)])
    ck(sorted(sizes)==expected,'all single-port orbits')
    if k==2:
        pairs={(j,l) for j in range(3) for l in range(3)}
        while pairs:
            j,l=min(pairs)
            orb={(g[3+j]-3,g[6+l]-6) for g in h}
            ck(orb=={(a,b) for a in range(3) for b in range(3) if (a+b-j-l)%3==0},'joint port charge')
            pairs-=orb
    nonabelian=any(comp(a,b)!=comp(b,a) for a in corrections for b in corrections)
    if k==3: ck(nonabelian,'depth-three holonomy is nonabelian')
    moments=port_moments(k,8)[1][k]
    for r in range(9):
        total=sum(sum(i==v for i,v in enumerate(g))**r for g in h)
        ck(total==len(h)*moments[r],'ordered-port recurrence versus explicit centralizer Burnside')
    return dict(depth=k,vertices=n,mark_count=len(marks),triangle_generator_count=len(corrections),
                holonomy_order=len(h),port_orbit_sizes=sorted(sizes),nonabelian=nonabelian,
                maximum_constructed_generator_length=longest,uniform_transport_edge_bound=3*(1+(k-1)*n),
                holonomy_sha256=hashlib.sha256('\n'.join(','.join(map(str,g)) for g in sorted(h)).encode()).hexdigest())


def port_moments(depth,max_r):
    a=[[1]*(max_r+1)];c=[[1]*(max_r+1)]
    for k in range(1,depth+1):
        aa=[1];cc=[1]
        for r in range(1,max_r+1):
            raw=sum(math.comb(r,i)*math.comb(r-i,j)*a[-1][i]*a[-1][j]*a[-1][r-i-j]
                    for i in range(r+1) for j in range(r-i+1))
            ck(raw%3==0,'tuple orbit division integral')
            aa.append(raw//3)
            cc.append(sum(math.comb(r,s)*2**(r-s)*c[-1][s]*a[-1][r-s] for s in range(r+1)))
        a.append(aa);c.append(cc)
    return a,c


def port_polynomial_controls():
    a,c=port_moments(12,8)
    for k in range(13):
        ck(a[k][2]==2*k+1 and a[k][3]==6*k*k+2*k+1,'automorphism low moments')
        ck(c[k][1]==2*k+1 and c[k][2]==8*k*k+1,'centralizer low moments')
    for r in range(1,9):
        for table,degree,lead in [(a,r-1,math.factorial(r)),(c,r,2**r*math.factorial(r))]:
            row=[v[r] for v in table]
            for _ in range(degree):row=[y-x for x,y in zip(row,row[1:])]
            ck(set(row)=={lead*math.factorial(degree)},'exact polynomial leading finite difference')
    return dict(universe='ordered tuples with repetitions',max_r=8,depth=12,
                rows=[dict(k=k,automorphism=a[k][:5],centralizer=c[k][:5]) for k in range(5)])


def composite_control():
    ws=words(9)
    def labels(w):
        a=0
        for i in range(9):
            if (w>>i)&1: a ^= [1,2,3][i%3]
        return a,sum(((w>>i)&1)*pow(5,i,19) for i in range(9))%19
    mirror={labels(w):labels(reflect(w,9)) for w in ws}
    ck(len(mirror)==76,'period-nine CRT bijection')
    square=[0,1,3,2]
    for (a,b),(c,d) in mirror.items():
        ck(c==square[a],'F4 Frobenius reflection')
    ambiguity={b:len({mirror[a,b][1] for a in range(4)}) for b in range(19)}
    ck(ambiguity=={b:4 if b else 1 for b in range(19)},'F19 reflection needs full F4 tag')
    def omega(a,e):
        for _ in range(e%3): a=[0,2,3,1][a]
        return a
    for (a,b),image in mirror.items():
        if not b: continue
        sign,i=next((s,i) for s in [1,-1] for i in range(9) if s*pow(5,i,19)%19==b)
        c=omega(a,-i); aa,bb=mirror[c,sign%19]
        outsign,j=next((s,j) for s in [1,-1] for j in range(9) if s*pow(5,j,19)%19==bb)
        predicted=(omega(square[c],-i),outsign*pow(5,(j-i)%9,19)%19)
        ck(predicted==image,'eight-chart phase reflection formula')
    # A full 76-cycle and an eight-fixed-point involution in an imprimitive group.
    cycle=tuple(4*((b+1)%19)+(j+(b==18))%4 for b in range(19) for j in range(4))
    r=tuple(4*b+(j if b<2 else j^1) for b in range(19) for j in range(4))
    ck(len({power(cycle,i)[0] for i in range(76)})==76,'composite full cycle')
    ck(sum(i==v for i,v in enumerate(r))==8,'composite low-fixed involution')
    for g in [cycle,r]:
        ck(all(len({g[4*b+j]//4 for j in range(4)})==1 for b in range(19)), 'imprimitive partition survives')
    return dict(ring='F4 x F19',mirror_F4=square,F19_output_counts=ambiguity,
                eight_seed_images=[dict(sign=s,tag=a,image=list(mirror[a,s%19])) for s in [1,-1] for a in range(4)])


def main():
    report={
        'scope':'Exact finite controls for written all-size proofs; same-agent audit',
        'lucas_reflections':reflection_controls(),
        'projective_involutions':projective_controls(),
        'prime_arithmetic':arithmetic_controls(),
        'marking_holonomy':[wreath_controls(k) for k in [1,2,3]],
        'composite_period9':composite_control(),
        'port_moments':port_polynomial_controls(),
    }
    report['active_checks']=checks
    print(json.dumps(report,indent=2,sort_keys=True))


if __name__=='__main__':
    main()
