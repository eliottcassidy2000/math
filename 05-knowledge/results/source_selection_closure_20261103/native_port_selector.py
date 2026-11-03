#!/usr/bin/env python3
"""Literal port selectors and guarded balanced renewal. No repository imports."""
from functools import lru_cache
from itertools import combinations, permutations
from math import gcd, isqrt
import json


def need(ok, label):
    if not ok:
        raise RuntimeError(label)


def endpoint_dp(T):
    n = len(T)
    D = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        D[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v, count in enumerate(D[mask]):
            if count:
                rest = T[v] & ~mask
                while rest:
                    bit = rest & -rest
                    rest -= bit
                    D[mask | bit][bit.bit_length()-1] += count
    return D


@lru_cache(None)
def H(T):
    return sum(endpoint_dp(T)[-1]) if T else 1


def induced(T, names):
    return tuple(sum(((T[u] >> v) & 1) << j for j, v in enumerate(names)) for u in names)


def profile(T, root):
    D = endpoint_dp(T)
    n = len(T)
    all_mask = (1 << n)-1
    hs = [sum(row) for row in D]
    hs[0] = 1
    F = sum(D[mask][root]*hs[all_mask ^ mask] for mask in range(1, 1 << n) if mask >> root & 1)
    # Independent initial-root DP on the graph with that first vertex removed.
    names = tuple(i for i in range(n) if i != root)
    R = induced(T, names)
    if not R:
        s = 1
    else:
        rev = tuple(((1 << len(R))-1) ^ row ^ (1 << i) for i, row in enumerate(R))
        starts = endpoint_dp(rev)[-1]
        s = sum(starts[i] for i, v in enumerate(names) if T[root] >> v & 1)
    return hs[-1], s, H(R), F


def literal_profile(T, root):
    h = s = F = 0
    for p in permutations(range(len(T))):
        good = [bool(T[a] >> b & 1) for a, b in zip(p, p[1:])]
        if all(good):
            h += 1
            s += p[0] == root
        cut = p.index(root)
        F += all(value for i, value in enumerate(good) if i != cut)
    return h, s, H(induced(T, tuple(i for i in range(len(T)) if i != root))), F


def ear(T, incoming):
    n = len(T)
    return tuple(row | ((1 << n) if incoming >> i & 1 else 0) for i, row in enumerate(T)) + (((1 << n)-1) ^ incoming,)


def birth(T, root):
    return ear(T, 1 << root)


def history(word):
    T = (2, 4, 1)
    for root in word:
        T = birth(T, root)
    return T


def R(z):
    h, s, q, F = z
    return h+F-s, q, q, 2*F-q


def N(z):
    h, s, q, F = z
    return h+F-s, h-s, h, h+F


def block_profile(z, A, B):
    h, s, q, F = z
    a, d = h-s, F-q
    k = A+B+1
    return (a+q+(2**k-1)*d,
            2**B*d if A else a+(2**B-1)*d,
            h if k == 1 else a+q+(2**(k-1)-1)*d,
            2**A*a+(2**(B+1)*3**A-2**A)*d+2**(A+1)*q+(s-q if B == 0 else 0))


QUERIES = ((), (1,), (2,), (3,), (0,1), (1,1), (2,1), (3,1),
           (0,2), (1,2), (2,2), (3,2), (2,0,1), (3,0,1), (0,2,1), (2,2,1), (3,2,1))


def X(T, ports):
    out = []
    for query in QUERIES:
        U = T
        for mask in query:
            U = ear(U, sum(1 << ports[i] for i in range(2) if mask >> i & 1))
        out.append(H(U))
    return tuple(out)


def set_module_order(T, order):
    V = list(T)
    for i, u in enumerate(order):
        for v in order[i+1:]:
            V[u] |= 1 << v
            V[v] &= ~(1 << u)
    return tuple(V)


def native(T, p):
    need(sorted(p) == list(range(len(T))) and all(T[u] >> v & 1 for u, v in zip(p,p[1:])), 'native path guard')
    edges = {frozenset((u,v)) for u,v in zip(p,p[1:])}
    V = list(T)
    for u in range(len(T)):
        for v in range(u):
            if frozenset((u,v)) not in edges:
                V[u] ^= 1 << v
                V[v] ^= 1 << u
    return tuple(V)


def module_path(T, order):
    # Redei insertion on the quotient, then expansion of the transitive module.
    representative = order[0]
    outside = [v for v in range(len(T)) if v not in order]
    p = []
    for v in outside+[representative]:
        i = next((i for i,u in enumerate(p) if T[v] >> u & 1), len(p))
        p.insert(i, v)
    p = tuple(v for u in p for v in (order if u == representative else (u,)))
    need(all(T[u] >> v & 1 for u,v in zip(p,p[1:])), 'module expansion Hamilton path')
    return p


def adjacent_native_swap(T, order, index, tag):
    u, v = order[index:index+2]
    P = module_path(T, order)
    n = len(T)
    if v != P[-1]:
        # Temporary vertex beats only u.
        E = ear(T, ((1 << n)-1) ^ (1 << u))
        P1 = P+(n,)
        branch = 'append'
    else:
        # Temporary vertex loses only to v.
        E = ear(T, 1 << v)
        P1 = (n,)+P
        branch = 'prepend'
    i = P.index(u)
    need(P[i+1] == v, 'traced adjacent clone arc')
    Q = P[:i+1]+(n,)+P[i+1:]
    U = native(E, P1)
    V = native(U, Q)
    result = tuple(row & ((1 << n)-1) for row in V[:n])
    expected = list(T)
    expected[u] ^= 1 << v
    expected[v] ^= 1 << u
    need(result == tuple(expected), 'exact single-old-arc output')
    new_order = list(order)
    new_order[index:index+2] = [v,u]
    need(result == set_module_order(T, new_order), 'module permutation output')
    return result, tuple(new_order), {'temporary_name': tag, 'temporary_index': n,
        'branch': branch, 'swapped_old_arc': (u,v), 'temporary_incoming_mask': sum(1 << i for i in range(n) if E[i] >> n & 1),
        'first_path': P1, 'second_path': Q, 'result_rows': result}


def compile_history(births, final_ports):
    states = {(0,1): ([], (0,1)), (1,2): ([], (1,2)), (2,0): ([], (2,0))}
    T = [2,4,1]
    available = {0,1,2}
    for name, root in births:
        need(root in available and name not in available, 'named chronological birth')
        available.add(name)
        nxt = {}
        for ports,(word,initial) in states.items():
            if root not in ports:
                continue
            slot = ports.index(root)
            for drop in range(3):
                pp = list(ports)
                if drop < 2:
                    pp[drop] = name
                nxt.setdefault(tuple(pp), (word+[3*slot+drop], initial))
        states = nxt
    need(final_ports in states, 'terminal marked reentry')
    word, initial = states[final_ports]
    ports = list(initial)
    for (name,root), code in zip(births,word):
        slot, drop = divmod(code,3)
        need(ports[slot] == root, 'literal program root')
        if drop < 2:
            ports[drop] = name
    need(tuple(ports) == final_ports, 'literal terminal slots')
    return {'initial_ports': initial, 'commands': word, 'final_ports': final_ports}


def replay_history(births):
    n = len(births)+3
    T = [0]*n
    T[0],T[1],T[2] = 2,4,1
    available = {0,1,2}
    for x, root in births:
        for u in available:
            if u == root:
                T[u] |= 1 << x
            else:
                T[x] |= 1 << u
        available.add(x)
    need(available == set(range(n)), 'history covers literal names')
    return tuple(T)


def twin_controls():
    count = 0
    examples = []
    for n in range(1,5):
        pairs = list(combinations(range(n),2))
        for mask in range(1 << len(pairs)):
            Q = [0]*n
            for j,(u,v) in enumerate(pairs):
                a,b = (u,v) if mask >> j & 1 else (v,u)
                Q[a] |= 1 << b
            Q = tuple(Q)
            for u in range(n):
                z = profile(Q,u)
                need(z == literal_profile(Q,u), 'independent literal rooted profile')
                T = tuple(row | ((1 << n) if i == u or row >> u & 1 else 0)
                          for i,row in enumerate(Q)) + (Q[u],)
                zu,zv = profile(T,u),profile(T,n)
                need(zu[1]-zv[1] == z[1] and zv[3]-zu[3] == z[3], 'both twin injections')
                gap = H(birth(T,n))-H(birth(T,u))
                need(gap == z[1]+z[3] > 0, 'strict exact twin birth gap')
                count += 1
                if Q == (2,4,1) and u == 0:
                    examples.append({'quotient_rows':Q,'root':u,'expanded_rows':T,'profile':z,'birth_gap':gap})
    return {'rooted_quotients_all_labelled_orders1_to4':count,'examples':examples}


def selector_controls():
    sources = []
    for n in range(1,4):
        pairs = list(combinations(range(n),2))
        for mask in range(1 << len(pairs)):
            T = [0]*n
            for j,(u,v) in enumerate(pairs):
                a,b = (u,v) if mask >> j & 1 else (v,u)
                T[a] |= 1 << b
            sources += [(tuple(T),r,3) for r in range(n)]
    sources += [(history((0,)),3,4), (history((0,3)),3,3), ((6,12,24,17,3),0,3)]
    checked = 0
    selected = []
    for T,v,limit in sources:
        z = profile(T,v)
        old_g = gcd(*z)
        later_g = gcd(z[2],z[3],z[0]-z[1])
        for k in range(1,limit+1):
            U = T
            for _ in range(k):
                U = birth(U,v)
            leaf = len(U)-1
            original_order = tuple(range(len(U)-1,len(T)-1,-1))
            previous_child = None
            for A in range(k):
                B = k-A-1
                order = list(original_order[1:])
                order.insert(A,leaf)
                V = set_module_order(U,order)
                actual = profile(V,leaf)
                need(actual == block_profile(z,A,B), 'full marked block formula')
                xx = X(V,(v,leaf))
                need(gcd(*xx) == (old_g if B == 0 else later_g), 'exact oldest versus later content')
                child = H(birth(V,leaf))
                need(child == actual[0]+actual[3]-actual[1], 'actual root birth')
                if previous_child is not None:
                    need(child > previous_child, 'strict ordered selector')
                previous_child = child
                checked += 1
                if (T,v) in ((history((0,)),3),(history((0,3)),3)):
                    selected.append({'source_rows':T,'root':v,'source_profile':z,'k':k,'A':A,'B':B,
                                     'rows':V,'ports':(v,leaf),'profile':actual,'queries':xx,'D':gcd(*xx),'child_H':child})
    return {'actual_position_query_controls':checked,'selected':selected}


def family_menu():
    rows = []
    for k in range(2,13):
        h = 3*(2**k+1)
        menu = []
        for A in range(k):
            B = k-A-1
            p = 3*2**k+3+6*2**B*3**A+6*2**A-3*2**B-(B == 0)
            menu.append(p)
        need(all(a<b for a,b in zip(menu,menu[1:])), 'family ordered menu')
        need(not set(menu+[3*(2**(k+1)+1)]) & {3*h,3*h+2,3*h+4}, 'all-position ternary hostile')
        rows.append({'k':k,'H':h,'ordered_leaf_child_menu':menu,'root_child_H':3*(2**(k+1)+1)})
    return rows


def renewal_cycle(T, root, births, b, t, k, cycle_tag, literal=True):
    need(k in (3,5), 'controller alphabet')
    if literal:
        need(profile(T,root) == (b,t,t,b), 'actual incoming balanced profile')
    stages = []
    new_root = len(T)
    T = birth(T,root)
    births = births+[(new_root,root)]
    stages.append({'kind':'N','rows':T,'H':2*b-t})
    leaves = []
    for i in range(k):
        x = len(T)
        T = birth(T,new_root)
        births.append((x,new_root))
        leaves.append(x)
        stages.append({'kind':'root-block','rows':T,'H':b*(2**(i+1)+1)})
    z = leaves[-1]
    order = tuple(reversed(leaves))
    pre_rows = T
    trace = []
    for i in range(k-1):
        T,order,row = adjacent_native_swap(T,order,i,f'{cycle_tag}_temporary_{i}')
        births[-k:] = [(u,new_root) for u in reversed(order)]
        need(replay_history(births) == T, 'reset endpoint full adjacency reentry')
        row['marked_reentry'] = compile_history(births,(new_root,z))
        trace.append(row)
    need(order[-1] == z, 'selected physical leaf oldest')
    oldest = block_profile(N((b,t,t,b)),k-1,0)
    need(oldest[0] == b*(2**k+1), 'oldest count')
    if literal:
        need(profile(T,z) == oldest, 'literal oldest rooted profile')
        if len(T) <= 9:
            need(gcd(*X(pre_rows,(new_root,z))) == b, 'literal newly introduced factor')
            need(gcd(*X(T,(new_root,z))) == gcd(b,t), 'literal erased factor')
    y = len(T)
    T = birth(T,z)
    births.append((y,z))
    exit_h = 2*b*(2**k+3**(k-1))-t
    stages.append({'kind':'oldest-exit','rows':T,'H':exit_h})
    last = len(T)
    T = birth(T,y)
    births.append((last,y))
    bp = b*(3*2**k+4*3**(k-1)+1)-2*t
    tp = b*(2**k+1)
    stages.append({'kind':'balanced-return','rows':T,'H':bp})
    need(R(N(oldest)) == (bp,tp,tp,bp), 'symbolic balanced renewal')
    need(replay_history(births) == T, 'whole renewed source graph')
    final_program = compile_history(births,(y,last))
    if literal:
        need(profile(T,y) == (bp,tp,tp,bp), 'literal final balanced profile')
        for stage in stages:
            need(H(stage['rows']) == stage['H'], 'literal cycle Hamilton stage')
    return T,y,births,bp,tp, {'k':k,'incoming_balanced_pair':(b,t),'outgoing_balanced_pair':(bp,tp),
        'pre_reset_rows':pre_rows,'retained_reset_ports':(new_root,z),'native_trace':trace,'stages':stages,
        'final_ports':(y,last),'final_root':y,'birth_history':births,'marked_program':final_program,
        'verification':'all stage counts and rooted profiles by subset DP' if literal else 'literal graph/path/reentry, count by proved transfer',
        'clocks':{'adaptive_births':k+3,'temporary_births':k-1,'temporary_deletions':k-1,'native_moves':2*(k-1)}}


def prime(n):
    return n >= 2 and all(n%d for d in range(2,isqrt(n)+1))


def renewals():
    seed = ((2,4,1),0,[],3,1)
    c3 = renewal_cycle(*seed,3,'first_3')
    c5 = renewal_cycle(*seed,5,'first_5')
    c33 = renewal_cycle(*c3[:5],3,'second_3')
    c35 = renewal_cycle(*c3[:5],5,'second_5',literal=False)
    # One larger selected source is counted literally; no support census.
    T,r,bs,b,t = c35[:5]
    seven_graph = birth(T,r)
    need((b,t) == (76147,5973), 'selected branching source')
    need(H(seven_graph) == 146321 == 7*20903 and prime(20903), 'literal 18-vertex seven-prime')
    c355 = renewal_cycle(*c35[:5],5,'third_5',literal=False)
    exit_stage = next(s for s in c355[-1]['stages'] if s['kind'] == 'oldest-exit')
    need(exit_stage['H'] == 17203249 == 7*2457607 and prime(2457607), 'larger exact transfer seven-prime')
    arithmetic = []
    b,t = 3,1
    for j in range(12):
        need(gcd(b,t) == 1, 'fixed3 primitive recursion')
        arithmetic.append((b,t))
        b,t = 61*b-2*t,9*b
    return {'literal_cycles':[c3[-1],c5[-1],c33[-1]],'branching_cycles':[c35[-1],c355[-1]],
            'first_seven_prime':{'rows':seven_graph,'H':146321,'quotient':20903,'order':len(seven_graph),
                                 'count_method':'independent literal graph subset DP'},
            'second_seven_prime':{'rows':exit_stage['rows'],'H':17203249,'quotient':2457607,
                                  'order':len(exit_stage['rows']),'count_method':'proved actual-source transfer; not subset DP here'},
            'fixed3_balanced_pairs':arithmetic}


def projective(b,t,p):
    return 'infinity' if b%p == 0 else (t*pow(b,-1,p))%p


def controller():
    states = (0,1,3,4,5,6,7,'infinity')
    transitions = []
    for r in states:
        b,t = (0,1) if r == 'infinity' else (1,r)
        r3 = projective(61*b-2*t,9*b,11)
        need(r3 in states, 'three-command state closure')
        row = {'state':r,'3':r3,'5':None if r == 7 else 0}
        if r != 7:
            need(projective(421*b-2*t,33*b,11) == 0, 'five-command reset to zero')
        transitions.append(row)
    mod7 = []
    b,t = 3,1
    for _ in range(3):
        counts = ((2*b-t)%7,3*b%7,5*b%7,9*b%7,(34*b-t)%7,(61*b-2*t)%7)
        need(all(counts), 'fixed-three seven-lane hostile')
        mod7.append({'balanced_pair_mod7':(b,t),'six_stage_residues':counts})
        b,t = (61*b-2*t)%7,9*b%7
    need((b,t) == (3,1), 'exact period-three cycle')
    return {'projective_modulus':11,'states':states,'transitions':transitions,
            'fixed3_mod7_cycle':mod7}


def main():
    print(json.dumps({'scope':'actual module port transport and guarded 3/5 balanced renewal; no all-target coverage',
                     'static_query_words':QUERIES,'twin_controls':twin_controls(),
                     'selector_controls':selector_controls(),'special_family_menu':family_menu(),
                     'controller':controller(),'renewal':renewals()},sort_keys=True,indent=2))


if __name__ == '__main__':
    main()
