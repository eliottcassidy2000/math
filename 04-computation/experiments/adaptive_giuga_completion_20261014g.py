"""Complete fixed-count signed prime completion; no chosen prime cutoff.

The remaining count is a termination rank, not a claim about unbounded count.
Every accepted leaf retains its index and all old/new native prime guards.
"""
from dataclasses import dataclass, asdict, replace
from fractions import Fraction
from itertools import combinations
from math import gcd, lcm, prod
import json

import giuga_korselt_defects_20261014f as one
import giuga_block_lattice_20261014f as two


def integer(n, minimum=None):
    return one.integer(n, minimum)


def sign(mu):
    return one.sheet(mu)


def boolean(x):
    if type(x) is not bool:
        raise ValueError("exact Boolean required")
    return x


def allowed_primes(m, last, count):
    """The least count primes above last not dividing m; exact finite search."""
    integer(m, 2)
    integer(last, 1)
    integer(count, 1)
    out = []
    p = last + 1
    while len(out) < count:
        if m % p and one.prime(p):
            out.append(p)
        p += 1
    return tuple(out)


def capacity(m, last, count, mu):
    sign(mu)
    ps = allowed_primes(m, last, count)
    return sum((Fraction(1, p) for p in ps), Fraction()) + Fraction(mu, m*prod(ps))


def divisors(n):
    integer(n, 1)
    out = [1]
    for p, e in one.factor(n):
        out = [a*p**j for a in out for j in range(e+1)]
    return tuple(sorted(out))


@dataclass(frozen=True)
class Pending:
    core: one.Core
    mu: int
    joint: bool
    residue: int
    modulus: int


def pending(c, mu, joint=False):
    """Condition on the WHOLE remaining product, never on one partial choice."""
    one.audit_core(c)
    sign(mu)
    boolean(joint)
    m, lam = c.source, c.period
    r = (-mu*pow(c.cofactor_sum, -1, m)) % m
    if not joint:
        return Pending(c, mu, False, r, m)
    if gcd(m, lam) != 1:
        return None
    s = pow(m, -1, lam) if lam > 1 else 0
    t = ((s-r)*pow(m, -1, lam)) % lam if lam > 1 else 0
    return Pending(c, mu, True, r+m*t, m*lam)


def audit_pending(p):
    if type(p) is not Pending:
        raise ValueError("Pending required")
    one.audit_core(p.core)
    sign(p.mu)
    boolean(p.joint)
    integer(p.residue, 0)
    integer(p.modulus, 1)
    if p != pending(p.core, p.mu, p.joint):
        raise ValueError("pending guard/source mismatch")
    return p


def accepts_pending(p, remaining_product):
    audit_pending(p)
    integer(remaining_product, 1)
    return (remaining_product-p.residue) % p.modulus == 0


@dataclass(frozen=True)
class Completion:
    core: one.Core
    mu: int
    index: int
    added: tuple
    target: int
    joint: bool


def audit_completion(row):
    if type(row) is not Completion or type(row.added) is not tuple:
        raise ValueError("exact Completion required")
    one.audit_core(row.core)
    sign(row.mu)
    boolean(row.joint)
    integer(row.index, 0)
    integer(row.target, 2)
    for p in row.added:
        integer(p, 2)
    if row.added and (tuple(sorted(set(row.added))) != row.added
                      or not all(one.prime(p) for p in row.added)):
        raise ValueError("added primes must be distinct, sorted and prime")
    if any(row.core.source % p == 0 for p in row.added):
        raise ValueError("new-prime guard failed")
    ps = tuple(sorted(row.core.primes+row.added))
    n = row.core.source*prod(row.added)
    if n != row.target:
        raise ValueError("target/source mismatch")
    if any((n//p+row.mu) % p for p in ps):
        raise ValueError("old or new signed local guard failed")
    if sum(n//p for p in ps)+row.mu != row.index*n:
        raise ValueError("wrong exact integer index")
    if row.joint and (len(ps) < 2 or any((n-1) % (p-1) for p in ps)):
        raise ValueError("full composite Korselt target failed")
    return row


def native_guards(row):
    audit_completion(row)
    return tuple((p, (row.target//p+row.mu) % p,
                  (row.target-1) % (p-1))
                 for p in sorted(row.core.primes+row.added))


def compile_completion(c, count, mu, joint=False, target_index=None):
    """Complete for this exact added-prime count; ordered prime tuple once only."""
    one.audit_core(c)
    integer(count, 0)
    sign(mu)
    boolean(joint)
    if target_index is not None:
        integer(target_index, 0)
    n, co = c.source, c.cofactor_sum
    S = Fraction(co, n)
    lower = S+Fraction(mu,n) if count == 0 else S
    upper = lower if count == 0 else S+capacity(n, 1, count, mu)
    stats = {"nodes": 0, "nonpositive_gap": 0, "capacity_prunes": 0,
             "old_korselt_prunes": 0, "prime_candidates": 0,
             "one_term_nodes": 0, "two_term_nodes": 0,
             "two_term_divisors": 0, "max_next_bound": 0,
             "max_factor_target": 0, "leaf_korselt_rejections": 0}
    rows = []

    def emit(k, added):
        ps = tuple(sorted(c.primes+added))
        target = n*prod(added)
        plain = Completion(c, mu, k, added, target, False)
        audit_completion(plain)
        if joint and (len(ps) < 2 or any((target-1) % (p-1) for p in ps)):
            stats["leaf_korselt_rejections"] += 1
            return
        rows.append(audit_completion(replace(plain, joint=joint)))

    def visit(m, A, s, last, added, k):
        stats["nodes"] += 1
        if A <= 0:
            stats["nonpositive_gap"] += 1
            return
        ps = tuple(sorted(c.primes+added))
        if joint and gcd(m, lcm(*(p-1 for p in ps))) != 1:
            stats["old_korselt_prunes"] += 1
            return
        if Fraction(A, m) > capacity(m, last, s, mu):
            stats["capacity_prunes"] += 1
            return
        if s == 1:
            stats["one_term_nodes"] += 1
            q, rem = divmod(m+mu, A)
            if not rem and q > last and m % q and one.prime(q):
                emit(k, added+(q,))
            return
        if s == 2:
            stats["two_term_nodes"] += 1
            D = m*m+mu*A
            stats["max_factor_target"] = max(stats["max_factor_target"], D)
            if D <= 0:
                return
            for u in divisors(D):
                v = D//u
                if u >= v:
                    break
                stats["two_term_divisors"] += 1
                if (u+m) % A or (v+m) % A:
                    continue
                q, r = (u+m)//A, (v+m)//A
                if last < q < r and m % q and m % r and one.prime(q) and one.prime(r):
                    emit(k, added+(q, r))
            return
        bound = (s*m-1)//A if mu == -1 else (s*m+1)//A
        stats["max_next_bound"] = max(stats["max_next_bound"], bound)
        for q in range(last+1, bound+1):
            if m % q and one.prime(q):
                stats["prime_candidates"] += 1
                visit(m*q, A*q-m, s-1, q, added+(q,), k)

    if count == 0:
        k, rem = divmod(co+mu, n)
        if not rem and (target_index is None or target_index == k):
            emit(k, ())
        indices = [] if rem or (target_index is not None and target_index != k) else [k]
    else:
        indices = list(range(S.numerator//S.denominator+1,
                             upper.numerator//upper.denominator+1))
        if target_index is not None:
            indices = [k for k in indices if k == target_index]
        for k in indices:
            visit(n, k*n-co, count, 1, (), k)
    rows.sort(key=lambda r: (r.index, r.added))
    return {"core": n, "count": count, "mu": mu, "joint": joint,
            "requested_index": target_index,
            "core_support": list(c.primes), "index_interval": [str(lower), str(upper)],
            "index_lower_strict": count > 0,
            "indices": indices, "statistics": stats, "completions": tuple(rows)}


def json_catalog(cat):
    return {**cat, "completions": [{"index": r.index, "added": list(r.added),
                                    "target": r.target,
                                    "all_prime_guards": native_guards(r)}
                                   for r in cat["completions"]]}


def marked_trace(row):
    """Each unfinished node carries the same target index and remaining guards."""
    audit_completion(row)
    out = []
    for i in range(len(row.added)+1):
        m = row.core.source*prod(row.added[:i])
        c = one.core(m)
        A = row.index*m-c.cofactor_sum
        g = pending(c, row.mu, row.joint)
        remaining = prod(row.added[i:])
        if g is None or not accepts_pending(g, remaining):
            raise ValueError("completed target violates a retained prefix guard")
        if i < len(row.added) and A <= 0:
            raise ValueError("unfinished positive deficit missing")
        if i == len(row.added) and A != row.mu:
            raise ValueError("terminal signed unit defect missing")
        out.append({"selected": list(row.added[:i]), "partial_product": m,
                    "remaining_product": remaining, "remaining_count": len(row.added)-i,
                    "target_index": row.index, "A": A,
                    "pending_residue": g.residue, "pending_modulus": g.modulus})
    return out


def main():
    checks = 0
    def check(ok):
        nonlocal checks
        checks += 1
        if not ok:
            raise RuntimeError("exact control failed")
    def reject(fn):
        try:
            fn()
        except (ValueError, TypeError):
            check(True)
        else:
            raise RuntimeError("malformed packet accepted")

    # Independent local-prime census: bounded control universe, unbounded compiler.
    bank = tuple(p for p in range(2, 24) if one.prime(p))
    core_count = comparisons = 0
    for n in range(2, 81):
        if any(e > 1 for _, e in one.factor(n)):
            continue
        core_count += 1
        c = one.core(n)
        for mu in (-1, 1):
            for count in range(4):
                cat = compile_completion(c, count, mu)
                seen = {r.added for r in cat["completions"]}
                joint = compile_completion(c, count, mu, True)
                joint_seen = {r.added for r in joint["completions"]}
                check(joint_seen == {r.added for r in cat["completions"]
                                    if len(c.primes+r.added) > 1
                                    and all((r.target-1) % (p-1) == 0 for p in c.primes+r.added)})
                for added in combinations(tuple(p for p in bank if n % p), count):
                    target = n*prod(added)
                    ps = c.primes+added
                    signed = all((target//p+mu) % p == 0 for p in ps)
                    check((added in seen) == signed)
                    korselt = len(ps) > 1 and all((target-1) % (p-1) == 0 for p in ps)
                    check((added in joint_seen) == (signed and korselt))
                    comparisons += 1
                for row in cat["completions"]:
                    check(audit_completion(row) == row)
                    check(sum((Fraction(1,p) for p in c.primes+row.added), Fraction())
                          + Fraction(mu,row.target) == row.index)
                    marked_trace(row)
                if count == 1:
                    old = one.complete_one(c, -mu)
                    check(seen == (set() if old is None else {(old.prime,)}))
                if count == 2:
                    old = two.compile_two(n, mu)
                    check(seen == {tuple(r["primes"]) for r in old["completions"]})

    cats = [compile_completion(one.core(n), count, mu, joint)
            for n in (30,561,1105,1729) for count in (3,4)
            for mu in (-1,1) for joint in (False,True)]
    plain = [c for c in cats if not c["joint"]]
    # Independent Fraction tree: no factorization, capacity pruning, or integer-A update.
    # The coarse index interval and looser next-prime bound still cover every solution.
    independent_nodes = 0
    def fraction_tree(n, count, mu):
        S = sum((Fraction(1,p) for p in one.core(n).primes),Fraction())
        answers = []
        def step(m, gap, left, last, added, k):
            nonlocal independent_nodes
            independent_nodes += 1
            if gap <= 0:
                return
            if left == 1:
                q = (1+Fraction(mu,m))/gap
                if q.denominator == 1 and q > last and m % int(q) and one.prime(int(q)):
                    answers.append((k,added+(int(q),)))
                return
            bound = int((left+Fraction(1,m))/gap)
            for q in range(last+1,bound+1):
                if m % q and one.prime(q):
                    step(m*q,gap-Fraction(1,q),left-1,q,added+(q,),k)
        for k in range(int(S)+1,int(S+Fraction(count,2)+Fraction(1,2))+1):
            step(n,Fraction(k)-S,count,1,(),k)
        return set(answers)
    for cat in plain:
        check(fraction_tree(cat["core"],cat["count"],cat["mu"])
              == {(r.index,r.added) for r in cat["completions"]})
    hits = [r for cat in plain for r in cat["completions"]]
    check([(r.core.source,r.mu,r.added,r.target,r.index) for r in hits]
          == [(561,1,(2,101,149,3109),52495396602,1)])
    check(all(not cat["completions"] for cat in cats if cat["joint"]))
    chosen = hits[0]
    check(compile_completion(chosen.core,4,1,target_index=1)["completions"] == (chosen,))
    check(not compile_completion(chosen.core,4,1,target_index=2)["completions"])
    trace = marked_trace(chosen)
    check(all(one.signed(one.core(step["partial_product"]),-1) is None
              for step in trace[:-1]))
    # Pending means remaining PRODUCT: testing one chosen prime would prune a success.
    c6 = one.core(6)
    g = pending(c6,-1)
    check(not accepts_pending(g,7) and accepts_pending(g,7*41))
    row = next(r for r in compile_completion(c6,2,-1)["completions"] if r.added==(7,41))
    check(marked_trace(row)[1]["remaining_product"] == 41)
    # A successful signed terminal need not be a full joint terminal.
    r561 = compile_completion(one.core(561),2,-1)["completions"][0]
    check(r561.added == (2,59) and r561.target == 66198)
    reject(lambda: audit_completion(replace(r561,joint=True)))
    # Any unfinished real suffix has positive gap; terminal A is mu, not always positive.
    check(marked_trace(row)[-1]["A"] == -1)
    check(trace[-1]["A"] == 1)
    # Closed target on count0 remains legal; no invented strict S<k boundary there.
    check(len(compile_completion(one.core(30),0,-1)["completions"]) == 1)
    check(len(compile_completion(one.core(6),0,1)["completions"]) == 1)
    check(not compile_completion(one.core(2),0,-1,True)["completions"])
    check(compile_completion(one.core(2),0,-1)["index_interval"] == ["0","0"])
    check(not compile_completion(one.core(2),0,-1)["index_lower_strict"])
    for bad in (True,1.0,-1):
        reject(lambda bad=bad: compile_completion(c6,bad,-1))
    for bad in (True,-1.0,0,2):
        reject(lambda bad=bad: compile_completion(c6,2,bad))
    reject(lambda: compile_completion(c6,2,-1,1))
    reject(lambda: compile_completion(c6,2,-1,target_index=True))
    reject(lambda: compile_completion(c6,2,-1,target_index=1.0))
    reject(lambda: compile_completion(replace(c6,source=6.0),2,-1))
    for bad in (replace(chosen,index=True),replace(chosen,mu=True),
                replace(chosen,target=float(chosen.target)),replace(chosen,joint=1),
                replace(chosen,added=(2.0,101,149,3109)),replace(chosen,added=(2,2,149,3109)),
                replace(chosen,added=(2,101,149,3108)),replace(chosen,added=(3,101,149,3109)),
                replace(chosen,added=chosen.added[::-1]),replace(chosen,core=one.core(1105)),
                replace(chosen,index=2),replace(chosen,target=chosen.target+1)):
        reject(lambda bad=bad:audit_completion(bad))
    for bad in (replace(g,residue=float(g.residue)),replace(g,modulus=True),
                replace(g,joint=1),replace(g,core=one.core(30))):
        reject(lambda bad=bad:audit_pending(bad))
    reject(lambda:accepts_pending(g,True))
    reject(lambda:accepts_pending(g,0))
    print(json.dumps({"status":"PROVED fixed-count algorithm / FINITE-EXACT catalog",
                      "checks":checks,"control_universe":{"squarefree_cores_2_to_80":core_count,
                      "added_count_range":[0,3],"independent_local_prime_bank":list(bank),
                      "independent_tuple_comparisons":comparisons,
                      "independent_fraction_tree_nodes":independent_nodes,
                      "catalog_prime_cutoff":None},
                      "requested_catalogs":[json_catalog(c) for c in cats],
                      "positive_completion_trace":trace,
                      "scope":"Fixed added count terminates. No uniform count or prime bound, global Giuga/Carmichael resolution, or unsigned-prefix rejection is asserted."},
                     indent=2,sort_keys=True))


if __name__ == "__main__":
    main()
