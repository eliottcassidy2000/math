# Two-loop (excess-1) kernels of the linear-forest cluster expansion.
# For a quotient graph G with fixed reference edge orientations, K(G) = sum over labelled run decompositions
# (each degree-2 vertex: pass-through, mu = 1, or junction of two run ends, mu = -1; each core vertex: a set of
# disjoint half-edge pairs = pass-throughs, the rest run ends, mu = (-1)^(s-1)(s-1)! with s = block size)
# of mu * prod over runs of [(-1)^alpha + (-1)^(l - alpha)], alpha = edges run against reference.  Closed runs
# are forbidden.  Then (vertex-transitive T):
#   figure-eight contribution to log Phi = (1/8) sum_(a,b) K8(a,b) t_a t_b / n,         t_a = tr Shat^a,
#   theta contribution                   = (1/12) sum_(a,b,c) Kth(a,b,c) Theta(a,b,c),
#   Theta(a,b,c) = sum_(u,v) (Shat^a)_uv (Shat^b)_uv (Shat^c)_uv  (all three paths oriented u -> v).
import itertools, sys
from math import factorial

def mob(s): return (-1) ** (s - 1) * factorial(s - 1)

def center_configs(halves):
    """all sets of disjoint pairs of the given half-edges"""
    halves = list(halves)
    def rec(rem):
        if not rem: yield []; return
        h, rest = rem[0], rem[1:]
        for cfg in rec(rest): yield cfg                       # h unpaired
        for i, g in enumerate(rest):
            for cfg in rec(rest[:i] + rest[i + 1:]): yield [(h, g)] + cfg
    yield from rec(halves)

def kernel(paths, core_of):
    """paths: list of (length, start_core, end_core). core_of: list of core vertex names.
    Edge (p, i), i = 1..len, goes (reference) from node (p, i-1) to node (p, i); node (p,0) = start core,
    node (p, len) = end core.  Returns K."""
    # half-edges at core vertices: ('s', p) = first edge of path p at its start core, ('f', p) = last edge at end core
    halves = {c: [] for c in core_of}
    for p, (l, a, b) in enumerate(paths):
        halves[a].append(('s', p)); halves[b].append(('f', p))
    internal = [(p, i) for p, (l, a, b) in enumerate(paths) for i in range(1, l)]
    total = 0
    cfgs_per_core = [list(center_configs(halves[c])) for c in core_of]
    for choice in itertools.product([0, 1], repeat=len(internal)):          # 1 = pass-through
        passint = {internal[k] for k in range(len(internal)) if choice[k]}
        mu_int = (-1) ** (len(internal) - len(passint))
        for cc in itertools.product(*cfgs_per_core):
            mu = mu_int
            link = {}
            for ci, cfg in enumerate(cc):
                h = halves[core_of[ci]]
                npairs = len(cfg); s = len(h) - npairs              # block size: endpoints + pass-throughs
                mu *= mob(s)
                for x, y in cfg: link[x] = y; link[y] = x
            # build runs: edges (p, i); neighbour of edge across node (p, i) (between edge i and i+1)
            def across(p, i, side):
                """edge reached from edge (p,i) through its 'side' end ('hi' = node (p,i), 'lo' = node (p,i-1))"""
                l = paths[p][0]
                if side == 'hi':
                    if i < l: return ((p, i + 1), 'lo') if (p, i) in passint else None
                    h = ('f', p)
                else:
                    if i > 1: return ((p, i - 1), 'hi') if (p, i - 1) in passint else None
                    h = ('s', p)
                if h not in link: return None
                g = link[h]
                if g[0] == 's': return ((g[1], 1), 'lo')
                return ((g[1], paths[g[1]][0]), 'hi')
            seen = set(); prod = 1; ok = True
            alledges = [(p, i) for p, (l, a, b) in enumerate(paths) for i in range(1, l + 1)]
            for e0 in alledges:
                if e0 in seen: continue
                # find one end of the run containing e0: walk 'lo'-wards
                cur, entered = e0, 'hi'      # pretend we entered e0 from its hi side, walk out the lo side
                steps = 0; closed = False
                while True:
                    out = 'lo' if entered == 'hi' else 'hi'
                    nxt = across(cur[0], cur[1], out)
                    if nxt is None: break
                    cur, entered = nxt[0], nxt[1]
                    steps += 1
                    if cur == e0: closed = True; break
                    if steps > 10000: raise RuntimeError
                if closed: ok = False; break
                # now cur is an end edge; traverse forward from the free side
                # the free side of cur is 'out' computed last: we leave through the opposite side
                start, free_side = cur, ('lo' if entered == 'hi' else 'hi')
                # traverse: enter start through free_side
                l_run = 0; alpha = 0; cur, ent = start, free_side
                while True:
                    seen.add(cur); l_run += 1
                    # traversing edge cur from side ent to the other: reference goes lo -> hi
                    if ent == 'hi': alpha += 1
                    out = 'hi' if ent == 'lo' else 'lo'
                    nxt = across(cur[0], cur[1], out)
                    if nxt is None: break
                    cur, ent = nxt
                prod *= (-1) ** alpha + (-1) ** (l_run - alpha)
                if prod == 0: break
            if not ok or prod == 0: continue
            total += mu * prod
    return total

def K8(a, b): return kernel([(a, 'v', 'v'), (b, 'v', 'v')], ['v'])
def Kth(a, b, c): return kernel([(a, 'u', 'v'), (b, 'u', 'v'), (c, 'u', 'v')], ['u', 'v'])
def Kdb(a, b, c):   # dumbbell: loop a at u, loop b at w, path c from u to w
    return kernel([(a, 'u', 'u'), (b, 'w', 'w'), (c, 'u', 'w')], ['u', 'w'])

if __name__ == '__main__':
    print('K8(2,2) =', K8(2, 2), '(cluster_shapes.py needs -28)')
    print('Kth(2,1,1) =', Kth(2, 1, 1), ' Kth(1,2,1) =', Kth(1, 2, 1), '(cluster_shapes.py needs 4)')
    M = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    print('\nK8(a,b), a,b = 1..M:')
    for a in range(1, M + 1): print(' ', [K8(a, b) for b in range(1, M + 1)])
    print('\nKth(a,b,c) for a <= b <= c, a+b+c <= M+2:')
    for a in range(1, M + 1):
        for b in range(a, M + 1):
            for c in range(b, M + 1):
                if a + b + c <= M + 2: print(f'  ({a},{b},{c}): {Kth(a, b, c)}')
