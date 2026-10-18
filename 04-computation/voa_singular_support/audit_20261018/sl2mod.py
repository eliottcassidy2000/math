"""Independent implementation (audit) of the universal affine vertex algebra V^k(sl2).

Basis e, h, f with [h,e]=2e, [h,f]=-2f, [e,f]=h;  (e|f)=(f|e)=1, (h|h)=2.
[a_m, b_n] = [a,b]_{m+n} + m k (a|b) delta_{m+n,0}.
PBW monomial: tuple of modes (n, t) with n >= 1 (meaning t_{-n}), t in {0:e, 1:h, 2:f},
sorted by key (-n, t) ascending  (largest n first, ties e<h<f).
Polynomial side (gr_F): the same tuple read as the monomial prod x_{t, n-1}; d x_{t,j} = (j+1) x_{t,j+1}.
"""
import sys
import flint

sys.setrecursionlimit(100000)

E, H, F = 0, 1, 2
NAMES = 'ehf'
CHARGE = {E: 2, H: 0, F: -2}
# bracket table: BR[(a,b)] = dict t -> coeff
BR = {
    (E, E): {}, (H, H): {}, (F, F): {},
    (H, E): {E: 2}, (E, H): {E: -2},
    (H, F): {F: -2}, (F, H): {F: 2},
    (E, F): {H: 1}, (F, E): {H: -1},
}
FORM = {(E, F): 1, (F, E): 1, (H, H): 2}


def modekey(md):
    return (-md[0], md[1])


def colored_partitions(D, maxkey=None):
    """All PBW monomials of weight D (tuples of (n,t) sorted by modekey)."""
    out = []
    modes = [(n, t) for n in range(D, 0, -1) for t in (E, H, F)]  # already in modekey order

    def rec(rem, start, cur):
        if rem == 0:
            out.append(tuple(cur))
            return
        for i in range(start, len(modes)):
            n, t = modes[i]
            if n <= rem:
                cur.append((n, t))
                rec(rem - n, i, cur)
                cur.pop()
    rec(D, 0, [])
    return out


def weight(mono):
    return sum(n for n, _ in mono)


def charge(mono):
    return sum(CHARGE[t] for _, t in mono)


class AffSl2:
    def __init__(self, k, one):
        self.k = k * one
        self.one = one
        self.cache = {}

    def act(self, t, m, mono):
        """t_m applied to PBW monomial mono."""
        key = (t, m, mono)
        r = self.cache.get(key)
        if r is not None:
            return r
        one = self.one
        if not mono:
            res = {((-m, t),): one} if m < 0 else {}
        else:
            n1, t1 = mono[0]
            rest = mono[1:]
            if m < 0 and (m, t) <= (-n1, t1):
                res = {((-m, t),) + mono: one}
            else:
                res = {}
                # t_m t1_{-n1} rest = t1_{-n1} (t_m rest) + [t_m, t1_{-n1}] rest
                for m2, c2 in self.act(t, m, rest).items():
                    for m3, c3 in self.act(t1, -n1, m2).items():
                        v = res.get(m3, 0) + c2 * c3
                        if v == 0:
                            res.pop(m3, None)
                        else:
                            res[m3] = v
                for t3, cb in BR[(t, t1)].items():
                    for m3, c3 in self.act(t3, m - n1, rest).items():
                        v = res.get(m3, 0) + cb * c3
                        if v == 0:
                            res.pop(m3, None)
                        else:
                            res[m3] = v
                if m == n1 and (t, t1) in FORM:
                    cc = m * self.k * FORM[(t, t1)]
                    if cc != 0:
                        v = res.get(rest, 0) + cc
                        if v == 0:
                            res.pop(rest, None)
                        else:
                            res[rest] = v
        self.cache[key] = res
        return res

    def actv(self, t, m, vec):
        out = {}
        for mono, c in vec.items():
            for m2, c2 in self.act(t, m, mono).items():
                v = out.get(m2, 0) + c * c2
                if v == 0:
                    out.pop(m2, None)
                else:
                    out[m2] = v
        return out


# ------------------------- polynomial side -------------------------------------------
def pmono_mul(a, b):
    return tuple(sorted(a + b, key=modekey))


def pmul(p1, p2):
    out = {}
    for a, ca in p1.items():
        for b, cb in p2.items():
            m = pmono_mul(a, b)
            v = out.get(m, 0) + ca * cb
            if v == 0:
                out.pop(m, None)
            else:
                out[m] = v
    return out


def pd(p):
    """d x_{t,j} = (j+1) x_{t,j+1}; variable (n,t) <-> x_{t,n-1}, so d(n,t) = n * (n+1,t)."""
    out = {}
    for m, c in p.items():
        seen = set()
        for i, md in enumerate(m):
            if md in seen:
                continue
            seen.add(md)
            mult = m.count(md)
            lst = list(m)
            lst.remove(md)
            lst.append((md[0] + 1, md[1]))
            nm = tuple(sorted(lst, key=modekey))
            v = out.get(nm, 0) + c * mult * md[0]
            if v == 0:
                out.pop(nm, None)
            else:
                out[nm] = v
    return out


def lidegree(mono):
    """Li degree = weight - number of modes."""
    return weight(mono) - len(mono)


def pretty(poly):
    terms = []
    for m, c in sorted(poly.items(), key=lambda kv: kv[0]):
        mon = '*'.join(f"{NAMES[t]}{n - 1}" for n, t in m)
        terms.append(f"({c})*{mon}")
    return ' + '.join(terms)
