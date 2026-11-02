#!/usr/bin/env python3
"""Actual width-three local universality: independent profile and controller audit.
Standard library only; all assertions are explicit and survive python -O.
No imports from a repository checkout. Graphs are adjacency bitmask tuples.
"""
from functools import lru_cache
from itertools import product, permutations
from math import factorial, gcd, comb
import json


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def divisor(c):
    return factorial(c) * factorial(c - 1)


def cutoff(m):
    c = 1
    while divisor(c) % m:
        c += 1
    return c


def graph(n, mask):
    out = [0] * n
    k = 0
    for i in range(n):
        for j in range(i + 1, n):
            u, v = (i, j) if mask >> k & 1 else (j, i)
            out[u] |= 1 << v
            k += 1
    return tuple(out)


def linear(n):
    return tuple(sum(1 << j for j in range(i + 1, n)) for i in range(n))


def substitute(q, blocks):
    offsets = [0]
    for b in blocks:
        offsets.append(offsets[-1] + len(b))
    out = [0] * offsets[-1]
    for i, block in enumerate(blocks):
        for u, row in enumerate(block):
            out[offsets[i] + u] = row << offsets[i]
            for j, b in enumerate(blocks):
                if q[i] >> j & 1:
                    out[offsets[i] + u] |= ((1 << len(b)) - 1) << offsets[j]
    return tuple(out)


ONE = (0,)
JOIN = (2, 0)
CYCLE = (2, 4, 1)


def all_h(a):
    n = len(a)
    dp = [[0] * n for _ in range(1 << n)]
    out = [0] * (1 << n)
    out[0] = 1
    for v in range(n):
        dp[1 << v][v] = 1
    for s in range(1, 1 << n):
        for v, value in enumerate(dp[s]):
            if value:
                out[s] += value
                rest = a[v] & ~s
                while rest:
                    bit = rest & -rest
                    rest -= bit
                    dp[s | bit][bit.bit_length() - 1] += value
    return out


def hp(a):
    return all_h(a)[-1]


@lru_cache(None)
def profile(a):
    """Works also for digraphs with missing arcs, needed by the universal kernel."""
    n = len(a)
    h = all_h(a)
    @lru_cache(None)
    def go(s):
        if not s:
            return (1,) + (0,) * n
        bit = s & -s
        tail = s ^ bit
        ans = [0] * (n + 1)
        u = tail
        while True:
            block = u | bit
            old = go(s ^ block)
            for c in range(1, n + 1):
                ans[c] += h[block] * old[c - 1]
            if not u:
                break
            u = (u - 1) & tail
        return tuple(ans)
    return go((1 << n) - 1)


def literal_profile(a):
    n = len(a)
    out = [0] * (n + 1)
    for word in permutations(range(n)):
        bad = sum(not (a[word[i]] >> word[i + 1] & 1) for i in range(n - 1))
        for c in range(bad + 1, n + 1):
            out[c] += comb(n - 1 - bad, c - 1 - bad)
    for c in range(1, n + 1):
        require(out[c] % factorial(c) == 0, 'ordered cover division')
        out[c] //= factorial(c)
    return tuple(out)


from collections import deque
from math import lcm


def factor(n):
    out = {}
    p = 2
    while p * p <= n:
        while n % p == 0:
            out[p] = out.get(p, 0) + 1
            n //= p
        p += 1
    if n > 1:
        out[n] = out.get(n, 0) + 1
    return out


def phi(n):
    ans = n
    for p in factor(n):
        ans = ans // p * (p - 1)
    return ans


def inv(a, m):
    return 0 if m == 1 else pow(a, -1, m)


def crt(rows):
    ans, modulus = 0, 1
    for residue, m in rows:
        if m == 1:
            continue
        require(gcd(modulus, m) == 1, 'coprime CRT factors')
        ans += modulus * ((residue - ans) * inv(modulus, m) % m)
        modulus *= m
    return ans % modulus


@lru_cache(None)
def unit_recipes(m):
    start = 1 % m
    found = {start: (0, 0, 0)}
    pending = [start]
    for x in pending:
        for i, gen in enumerate((5, 17, 23)):
            y = x * gen % m
            if y not in found:
                exponents = list(found[x])
                exponents[i] += 1
                found[y] = tuple(exponents)
                pending.append(y)
    require(len(found) == phi(m), '2/3-primary unit bank generation')
    return found


def mod5_word(start, target):
    pending = [start]
    words = {start: ''}
    for h, g in pending:
        if h == target and g:
            return words[(h, g)]
        for name, nxt in (('R', (g, h)), ('E', ((2*h+4*g)%5, (3*h+3*g)%5))):
            if nxt not in words:
                words[nxt] = words[(h,g)] + name
                pending.append(nxt)
    raise RuntimeError('mod-five primitive orbit target absent')


def construct(m, target):
    """Compact positive actual-graph recipe, represented by proved gate words.
    Inverse primitive gates in R/E/S compile to positive finite permutation powers.
    The returned scalar state is verified; no integer spectrum claim is made.
    """
    target %= m
    require(m % 2 or target % 2, 'parity-compatible target')
    f = factor(m)
    a, b, c = f.get(2,0), f.get(3,0), f.get(5,0)
    two, three = 2**a, 3**b
    n = m // two // three
    five, n0 = 5**c, n // 5**c
    v = 0
    while v < b and target % (3**(v+1)) == 0:
        v += 1
    seed_modulus = two * 3**(b-v)
    u = crt([(target * inv(3**v, two) % two, two),
             ((target // 3**v) % (3**(b-v)), 3**(b-v))])
    k, ell, e = unit_recipes(seed_modulus)[u]
    if not k:
        k += phi(seed_modulus)
    seed_h = pow(3,v,m)*pow(5,k,m)*pow(17,ell,m)*pow(23,e,m) % m
    seed_g = pow(12,v,m)*pow(32,k,m)*pow(224,ell,m)*pow(296,e,m) % m
    require(seed_h % two == target % two and seed_h % three == target % three,
            'anchor seed matches both primary targets')
    require(gcd(gcd(seed_h,seed_g),n) == 1, 'primitive non-three seed pair')
    if five > 1:
        require(seed_h % 5 == 0 and seed_g % 5 != 0, 'forced zero/nonzero mod-five seed')
    h, g = seed_h % n, seed_g % n
    r = two * three
    tau = 15 * r * inv(4,n) % n
    steps = []
    def shear(power):
        nonlocal h
        h = (h + tau * power * g) % n
        if power:
            steps.append(('S',power))
    def R():
        nonlocal h,g
        h,g = ((5*h-g)*inv(4,n)%n,h)
        steps.append(('R',1))
    def E():
        nonlocal h,g
        h,g = ((11*h-3*g)*inv(8,n)%n,(3*h-2*g)%n)
        steps.append(('E',1))
    if n > 1:
        # Ensure h is a unit at every p|n0 before the initial swap-like R.
        rows=[]
        for p, exponent in factor(n0).items():
            modulus=p**exponent
            power=0 if g%p==0 else (1-h)*inv(tau*g,modulus)%modulus
            rows.append((power,modulus))
        shear(crt(rows))
        require(gcd(h,n0)==1, 'initial preconditioning makes h a unit')
        R()
        require(gcd(g,n0)==1, 'initial R makes g a unit')
        word=mod5_word((h%5,g%5),target%5) if five>1 else ''
        for op in word:
            desired = g if op == 'R' else 0
            shear((desired-h)*inv(tau*g,n0)%n0)
            R() if op == 'R' else E()
            require(gcd(g,n0)==1, 'safe intermediate g at non-five primes')
        require(gcd(g,n)==1, 'terminal g is a unit at every odd non-three prime')
        divisor_gcd=gcd(tau*g,n)
        require((target-h)%divisor_gcd==0, 'final shear solvability')
        reduced=n//divisor_gcd
        power=((target-h)//divisor_gcd)*inv((tau*g)//divisor_gcd,reduced)%reduced
        shear(power+n)  # Positive identity period additionally guarantees a strong last operation.
    else:
        # S still has an actual positive C3 context; n=1 is the identity on its empty n-side.
        shear(1)
    require(h == target % n, 'final exact n-coordinate')
    result = crt([(seed_h%two,two),(seed_h%three,three),(h,n)])
    require(result == target, 'simultaneous CRT landing, not separate witnesses')
    return {'M':m,'target':target,'seed_exponents_3_5_17_23':[v,k,ell,e],
            'seed_H_mod_M':seed_h,'seed_Q2_mod_M':seed_g,'shear_base_power':r,
            'n':n,'tau_mod_n':tau,'steps':steps,'terminal_Q2_mod_n':g,'result':result}


def evaluate_polynomial(p, z, m):
    value, fall = 0, 1
    for c in range(1,len(p)):
        fall *= z-c+1
        value += p[c]*fall
    return value % m


def pf(p,m):
    return tuple(evaluate_polynomial(p,z,m) for z in range(m))


def pad(F, exponent, m):
    return tuple(v*pow(z,exponent,m)%m for z,v in enumerate(F))


def join_functions(F,G,m):
    return tuple(x*y%m for x,y in zip(F,G))


def branch(F,G,m):
    return tuple(z*(F[z]*G[z]+(F[(z+1)%m]-F[z])*(G[(z+1)%m]-G[z])
                    +(F[(z-1)%m]-F[z])*(G[(z-1)%m]-G[z]))%m for z in range(m))


def matrix_mul(A,B,m):
    return tuple(sum(A[2*i+k]*B[2*k+j] for k in range(2))%m for i in range(2) for j in range(2))


def matrix_order(A,m):
    if m==1:return 1
    ident=(1,0,0,1);cur=ident
    for k in range(1,2*m**3+1):
        cur=matrix_mul(A,cur,m)
        if cur==ident:return k
    raise RuntimeError('small audit matrix-order guard exceeded')


def unit_order(x,m):
    if m==1:return 1
    cur=1
    for k in range(1,phi(m)+1):
        cur=cur*x%m
        if cur==1:return k
    raise RuntimeError('unit order')


def actual_gate_audit(m):
    """Construct gates through actual full polynomial-function operations, not pair matrices."""
    f=factor(m);a,b=f.get(2,0),f.get(3,0)
    two,three=2**a,3**b;n=m//two//three;odd=m//two
    N=(max(a,b)+4)*phi(odd)*phi(n)
    U0=tuple(pow(z,N,m) for z in range(m))
    def normalized(F,G):return pad(branch(F,G,m),N-1,m)
    def D1(F):return normalized(F,U0)
    U2=D1(U0);U4=D1(U2)
    def A(F):return normalized(F,U2)
    def A4(F):return normalized(F,U4)
    L2=tuple(z*z%m for z in range(m))
    V=pad(branch(U2,L2,m),N-3,m)
    def Vgate(F):return join_functions(F,V,m)
    orders={
      'D1':lcm(unit_order(3,two),unit_order(2,three),unit_order(2,n)),
      'A':lcm(unit_order(9,two),unit_order(4,three),matrix_order((5%n,(-1)%n,1%n,0),n)),
      'V':lcm(unit_order(-3,two),unit_order(4,three),unit_order(-inv(2,three),three),unit_order(inv(4,n),n))}
    def power(F,op,count):
        for _ in range(count):F=op(F)
        return F
    def inverse(F,which):return power(F,{'D1':D1,'A':A,'V':Vgate}[which],orders[which]-1)
    def U(F):return A(D1(D1(inverse(Vgate(F),'A'))))
    def R(F):return inverse(inverse(A(F),'D1'),'D1')
    def E(F):return inverse(inverse(inverse(A4(F),'D1'),'D1'),'D1')
    def state(F):return (F[1%m]%two,F[1%m]%three,F[2%m]%three,F[1%m]%n,F[2%m]%n)
    checked=0
    seeds=[ONE,CYCLE,linear(3),substitute(CYCLE,[linear(2),ONE,ONE]),graph(5,41)]
    for graph_seed in seeds:
        F=pad(pf(profile(graph_seed),m),N,m)
        h2,h3,g3,h,g=state(F)
        require(F[2%m]%two==0 and F[3%m]%three==0,'reset sidecars')
        expected={
          'D1':(3*h2%two,2*h3%three,2*g3%three,2*h%n,g),
          'A':(9*h2%two,4*h3%three,4*g3%three,(5*h-g)%n,h),
          'A4':(27*h2%two,8*h3%three,8*g3%three,(11*h-3*g)%n,(3*h-2*g)%n),
          'V':(-3*h2%two,4*h3%three,-inv(2,three)*g3%three,h,inv(4,n)*g%n),
          'U':(-27*h2%two,16*h3%three,-2*g3%three,(h+15*inv(4,n)*g)%n,g),
          'R':(h2,h3,g3,(5*h-g)*inv(4,n)%n,h),
          'E':(h2,h3,g3,(11*h-3*g)*inv(8,n)%n,(3*h-2*g)%n)}
        for name,op in [('D1',D1),('A',A),('A4',A4),('V',Vgate),('U',U),('R',R),('E',E)]:
            actual=op(F)
            require(state(actual)==expected[name],'actual polynomial gate '+name)
            require(actual[2%m]%two==0 and actual[3%m]%three==0,'preserved gate sidecars')
            checked+=1
        for name,op in [('D1',D1),('A',A),('V',Vgate)]:
            require(state(inverse(op(F),name))==state(F),'positive-power controlled inverse')
            checked+=1
    return checked


def actual_recipe_audit(m,target):
    """Construct gates through actual full polynomial-function operations, not pair matrices."""
    f=factor(m);a,b=f.get(2,0),f.get(3,0)
    two,three=2**a,3**b;n=m//two//three;odd=m//two
    N=(max(a,b)+4)*phi(odd)*phi(n)
    U0=tuple(pow(z,N,m) for z in range(m))
    def normalized(F,G):return pad(branch(F,G,m),N-1,m)
    def D1(F):return normalized(F,U0)
    U2=D1(U0);U4=D1(U2)
    def A(F):return normalized(F,U2)
    def A4(F):return normalized(F,U4)
    L2=tuple(z*z%m for z in range(m))
    V=pad(branch(U2,L2,m),N-3,m)
    def Vgate(F):return join_functions(F,V,m)
    orders={
      'D1':lcm(unit_order(3,two),unit_order(2,three),unit_order(2,n)),
      'A':lcm(unit_order(9,two),unit_order(4,three),matrix_order((5%n,(-1)%n,1%n,0),n)),
      'V':lcm(unit_order(-3,two),unit_order(4,three),unit_order(-inv(2,three),three),unit_order(inv(4,n),n))}
    def power(F,op,count):
        for _ in range(count):F=op(F)
        return F
    def inverse(F,which):return power(F,{'D1':D1,'A':A,'V':Vgate}[which],orders[which]-1)
    def U(F):return A(D1(D1(inverse(Vgate(F),'A'))))
    def R(F):return inverse(inverse(A(F),'D1'),'D1')
    def E(F):return inverse(inverse(inverse(A4(F),'D1'),'D1'),'D1')
    def state(F):return (F[1%m]%two,F[1%m]%three,F[2%m]%three,F[1%m]%n,F[2%m]%n)
    recipe=construct(m,target)
    v,k,ell,e=recipe['seed_exponents_3_5_17_23']
    F3=tuple((z**3+2*z)%m for z in range(m))
    F5=tuple((z**4+4*z*z)%m for z in range(m))
    F17=tuple((z**6+8*z**4+8*z*z)%m for z in range(m))
    F23=branch(tuple(pow(z,3,m) for z in range(m)),tuple(pow(z,2,m) for z in range(m)),m)
    F=tuple(pow(F3[z],v,m)*pow(F5[z],k,m)*pow(F17[z],ell,m)*pow(F23[z],e,m)%m for z in range(m))
    F=pad(F,N,m)
    for op,count in recipe['steps']:
        if op=='S':F=power(F,U,recipe['shear_base_power']*count)
        elif op=='R':F=R(F)
        else:F=E(F)
    require(F[1%m]==target%m,'end-to-end actual positive polynomial recipe')
    return 1


def main():
    # Foundation: actual vertex-level covers verify the bilinear discrete-energy identity.
    small=[ONE,linear(2),linear(3),CYCLE]
    foundation=0
    for T,U in product(small,repeat=2):
        a=substitute(CYCLE,[T,U,ONE]);p=profile(a)
        require(p==literal_profile(a),'direct vs literal branch covers')
        F,G=profile(T),profile(U)
        for z in range(-4,7):
            eval_exact=lambda q,x:sum(q[c]*__import__('math').prod(x-j for j in range(c)) for c in range(1,len(q)))
            ft,gu=eval_exact(F,z),eval_exact(G,z)
            value=z*(ft*gu+(eval_exact(F,z+1)-ft)*(eval_exact(G,z+1)-gu)
                        +(eval_exact(F,z-1)-ft)*(eval_exact(G,z-1)-gu))
            require(value==eval_exact(p,z),'bilinear polynomial kernel')
            foundation+=1
    seed_values=[]
    for T,H,Q2 in [(CYCLE,3,12),(substitute(CYCLE,[linear(2),ONE,ONE]),5,32),
                    (substitute(CYCLE,[linear(4),ONE,ONE]),17,224),
                    (substitute(CYCLE,[linear(3),linear(2),ONE]),23,296)]:
        p=profile(T);require(p[1]==H and 2*(p[1]+p[2])==Q2,'actual anchor seed')
        require(p==literal_profile(T),'literal anchor seed')
        seed_values.append({'H':H,'Q2':Q2,'adjacency':T})
    require(len(unit_recipes(72))==24,'full unit base72')
    gate_checks=sum(actual_gate_audit(m) for m in range(2,41))
    end_to_end=sum(actual_recipe_audit(m,z) for m in range(2,25) for z in range(m) if m%2 or z%2)
    end_to_end+=sum(actual_recipe_audit(m,z) for m,z in [(35,7),(65,13),(91,13),(120,7),(144,13)])
    require(((11*2-3*3)*inv(8,13)%13,(3*2-2*3)%13)==(0,0),'E prime13 kernel hostile')
    require((2+4)%5==1 and (3+3)%5==1,'fixed diagonal mod5 hostile')
    checks=0;examples=[];max_steps=0
    for m in range(1,257):
        for target in range(m):
            if m%2==0 and target%2==0:continue
            result=construct(m,target);checks+=1;max_steps=max(max_steps,len(result['steps']))
            if (m,target) in [(91,13),(120,7),(144,13),(210,21),(255,7),(256,21)]:examples.append(result)
    for m in (729,1024,2187,3125,4096,5**4*7**2,2**8*3**4*5**3*7*13,2**4*3**3*5**2*7*37):
        for target in (1,7,13,21,m-1,m//2|1):
            result=construct(m,target);checks+=1;max_steps=max(max_steps,len(result['steps']))
            examples.append(result)
    orbit_words={str(t):{str(z):mod5_word((t,0),z) for z in range(5)} for t in range(1,5)}
    print(json.dumps({'status':'PASS','foundation_direct_checks':foundation,'actual_polynomial_gate_checks':gate_checks,
                     'end_to_end_positive_recipe_checks':end_to_end,'E13_kernel_hostile':[2,3],
                     'mod5_fixed_diagonal_hostile':[1,1],'simultaneous_constructor_checks':checks,'max_recipe_steps':max_steps,
                     'actual_anchor_seeds':seed_values,'units_mod72':sorted(unit_recipes(72)),
                     'mod5_words_from_t_zero':orbit_words,'examples':examples,
                     'scope':'All-modulus theorem is PROOF.md; receipts verify finite implementation and actual gate semantics.'},indent=2,sort_keys=True))

if __name__=='__main__':main()
