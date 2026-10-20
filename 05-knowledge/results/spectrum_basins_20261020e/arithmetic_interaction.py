"""Local interaction/transfer audit; no imports from earlier research checkers."""
from collections import Counter
from functools import lru_cache
from hashlib import sha256
from itertools import permutations
import json
from pathlib import Path


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


@lru_cache(None)
def words(n):
    if not n:
        return ((),)
    return tuple((a,)+w for a in (1,3,4) if a <= n for w in words(n-a))


def occupancy(w):
    n = sum(w)
    chi = [0]*max(n-2,0)
    i = 0
    for a in w:
        if a >= 3:
            chi[i] = 1
        if a == 4:
            chi[i+1] = 1
        i += a
    return tuple(chi)


def join(a,b):
    m = len(a)
    return tuple(x | (((1 << len(b))-1) << m) for x in a) + tuple(x << m for x in b)


def graph(w):
    atoms = {1:(0,),3:(2,4,1),4:(6,4,8,3)}
    t = ()
    for a in w:
        t = join(t,atoms[a])
    return t


def literal_h(t):
    n = len(t)
    if n == 0:
        return 1
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1,1 << n):
        for v in range(n):
            if dp[mask][v]:
                available = t[v] & ~mask
                while available:
                    bit = available & -available
                    available -= bit
                    dp[mask | bit][bit.bit_length()-1] += dp[mask][v]
    return sum(dp[-1])


def polynomial_add(a,b,scale=1):
    r = dict(a)
    for monomial,value in b.items():
        r[monomial] = r.get(monomial,0) + scale*value
    return {k:v for k,v in r.items() if v}


def polynomial_mul(a,b):
    r = {}
    for x,u in a.items():
        for y,v in b.items():
            key = tuple(i+j for i,j in zip(x,y))
            r[key] = r.get(key,0) + u*v
    return {k:v for k,v in r.items() if v}


def determinant_control():
    # Monomial indices are (z degree, edge-weight u degree, interaction v degree).
    one = {(0,0,0):1}
    t = [[one,{(0,1,0):1},{},{}],
         [{},{},one,{(0,1,1):1}],
         [one,{},{},{}],
         [{},{},one,{}]]
    z = {(1,0,0):1}
    a = [[polynomial_add(one if i==j else {},polynomial_mul(z,t[i][j]),-1)
          for j in range(4)] for i in range(4)]
    det = {}
    for p in permutations(range(4)):
        term = one
        for i,j in enumerate(p):
            term = polynomial_mul(term,a[i][j])
        sign = (-1)**sum(p[i]>p[j] for i in range(4) for j in range(i+1,4))
        det = polynomial_add(det,term,sign)
    expected = {(0,0,0):1,(1,0,0):-1,(3,1,0):-1,(4,2,1):-1}
    need(det == expected, ('determinant',det))
    return [{'z_degree':k[0],'u_degree':k[1],'v_degree':k[2],'coefficient':v}
            for k,v in sorted(det.items())]


def main():
    # Initial state has two virtual zero occupation bits.
    states = {(0,0):Counter({(0,0):1})}
    transfer = {0:Counter({(0,0):1}),1:Counter({(0,0):1}),2:Counter({(0,0):1})}
    for length in range(1,17):
        target = {}
        for (older,last),poly in states.items():
            for new in (0,1):
                if older and new:
                    continue
                dest = target.setdefault((last,new),Counter())
                for (e,c),count in poly.items():
                    dest[(e+new,c+last*new)] += count
        states = target
        poly = Counter()
        for value in states.values():
            poly.update(value)
        transfer[length+2] = poly
    hcoeff = [1]
    for n in range(1,19):
        hcoeff.append(sum(weight*hcoeff[n-d] for d,weight in ((1,1),(3,3),(4,5)) if n>=d))
    integer_transfer = ((1,1,0,0),(0,0,3,5),(1,0,0,0),(0,0,1,0))
    finish = (1,3,1,1)
    state = (1,0,0,0)
    for length in range(17):
        need(sum(a*b for a,b in zip(state,finish)) == hcoeff[length+2], ('integer H transfer',length))
        state = tuple(sum(state[i]*integer_transfer[i][j] for i in range(4)) for j in range(4))
    rows = []
    literal_checks = 0
    for n in range(19):
        distribution = Counter()
        for w in words(n):
            chi = occupancy(w)
            need(all(not (chi[i] and chi[i+2]) for i in range(len(chi)-2)), ('exclusion',w))
            e = sum(chi)
            c = sum(a*b for a,b in zip(chi,chi[1:]))
            need(c == w.count(4), ('interaction',w,chi,c))
            need(e == w.count(3)+2*w.count(4), ('triangles',w))
            need(e-c == w.count(3)+w.count(4), ('feedback',w))
            need(e-2*c >= 0, ('disjoint interactions',w))
            h = 3**(e-2*c)*5**c
            need(h == 3**w.count(3)*5**w.count(4), ('H weight',w))
            if n <= 8:
                need(h == literal_h(graph(w)), ('literal H',w))
                literal_checks += 1
            distribution[(e,c)] += 1
        need(distribution == transfer[n], ('weighted transfer',n))
        weighted = sum(count*3**(e-2*c)*5**c for (e,c),count in distribution.items())
        need(weighted == hcoeff[n], ('H GF',n))
        rows.append({'n':n,'words':sum(distribution.values()),'sum_H':weighted,
                     'edge_interaction_distribution':[{'E':e,'C':c,'count':count}
                                                       for (e,c),count in sorted(distribution.items())]})
    # Exact shear composition on (q,p), independent of any knot or metagraph claim.
    L=((1,1),(0,1)); R=((1,0),(1,1))
    M=tuple(tuple(sum(L[i][k]*R[k][j] for k in range(2)) for j in range(2)) for i in range(2))
    need(M == ((2,1),(1,1)), 'LR composition')
    fib=[0,1]
    for _ in range(105):
        fib.append(fib[-1]+fib[-2])
    for k in range(101):
        q,p = fib[k+2],fib[k+1]
        need((2*q+p,q+p) == (fib[k+4],fib[k+3]), ('two-step address',k))
    examples=[]
    for w in ((3,3),(4,1,1),(1,1,1,3),(3,3)):
        chi=occupancy(w);e=sum(chi);c=sum(a*b for a,b in zip(chi,chi[1:]))
        examples.append({'word':w,'chi':chi,'E':e,'C':c,'feedback':e-c,'H':3**(e-2*c)*5**c})
    payload={'status':'FINITE-EXACT supporting the uniform local proof',
             'universe':'All7920golden SCC words through vertex order18; independently generated four-state occupancy strings',
             'source_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),
             'rows':rows,'literal_H_checks_through_n8':literal_checks,
             'integer_H_transfer':integer_transfer,'integer_H_terminal_vector':finish,
             'determinant_terms':determinant_control(),'examples':examples,'LR_matrix':M,
             'limitations':['Golden SCC language only','Endpoint weights do not constrain intermediate black moves',
                            'Shear action changes order; it is not a fixed-order black move']}
    dest=Path(__file__).with_suffix('.json')
    dest.write_text(json.dumps(payload,indent=2,sort_keys=True)+'\n',encoding='utf-8')
    print(json.dumps({'all_checks':'PASS','words':sum(r['words'] for r in rows),
                      'literal_H_checks':literal_checks,'H_GF_initial_coefficients':hcoeff[:14],
                      'receipt_sha256':sha256(dest.read_bytes()).hexdigest()},sort_keys=True))


if __name__ == '__main__':
    main()
