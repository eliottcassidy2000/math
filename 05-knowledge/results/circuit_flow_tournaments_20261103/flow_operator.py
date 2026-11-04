"""Exact one-module circuit compiler; no numerical interpolation."""
from collections import defaultdict
from itertools import permutations
from math import comb


def need(ok, why):
    if not ok:
        raise RuntimeError(why)


def patterns(seed, root=0):
    """Unordered families of nonempty exterior words, canonically ordered."""
    outside = [i for i in range(len(seed)) if i != root]
    if not outside:
        yield (0, 0, 0)
        return
    for order in permutations(outside):
        for cut in range(1 << (len(outside)-1)):
            words = []
            start = 0
            for i in range(len(outside)-1):
                if cut >> i & 1:
                    words.append(order[start:i+1]); start = i+1
            words.append(order[start:])
            if [min(x) for x in words] != sorted(min(x) for x in words):
                continue
            allowed = forced = 0
            valid = True
            for word in words:
                allowed += (seed[root] >> word[0]) & 1
                allowed += (seed[word[-1]] >> root) & 1
                for a,b in zip(word,word[1:]):
                    bridge = bool(seed[a] >> root & 1 and seed[root] >> b & 1)
                    direct = bool(seed[a] >> b & 1)
                    allowed += bridge
                    if not direct:
                        if not bridge:
                            valid = False; break
                        forced += 1
                if not valid:
                    break
            if valid:
                yield (len(words), allowed, forced)


def compile_operator(seed, root=0):
    """Map (shift, falling-factorial degree) to integer coefficient."""
    out = defaultdict(int)
    for ell, allowed, forced in patterns(seed, root):
        for i in range(forced+1):
            out[allowed-ell-i, ell] += (-1)**i * comb(forced,i)
    return {k:v for k,v in sorted(out.items()) if v}


def falling(r, k):
    t = 1
    for i in range(k): t *= r-i
    return t


def apply(op, r, f):
    return sum(c*falling(r,k)*f(r+s) for (s,k),c in op.items()
               if falling(r,k))


def symbolic(op):
    import sympy as sp
    r=sp.Symbol('r'); out={}
    for s in sorted({s for s,k in op}):
        out[s]=str(sp.factor(sum(c*sp.prod(r-i for i in range(k))
                               for (ss,k),c in op.items() if ss==s)))
    return out


if __name__ == '__main__':
    import json
    lo=(16,17,3,7,12); hi=(16,25,3,5,12)
    L=compile_operator(lo); U=compile_operator(hi)
    D={k:U.get(k,0)-L.get(k,0) for k in U.keys()|L.keys()}
    D={k:v for k,v in D.items() if v}
    print(json.dumps({'L':symbolic(L),'U':symbolic(U),'difference':symbolic(D)},indent=2))
