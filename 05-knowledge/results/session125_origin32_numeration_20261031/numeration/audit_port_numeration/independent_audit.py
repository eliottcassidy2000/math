"""Independent Session 125 arithmetic audit; explicit failures survive python -O.

Only the unit-under-test is imported. Oracle Fibonacci weights, regex legality,
coin DP, all local rewrites, and DFA refinement are independently constructed.
No source or remote files are changed by this program.
"""
import hashlib
import importlib.util
import itertools
import json
import random
import re
import sys
from collections import deque
from functools import lru_cache
from pathlib import Path

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent / 'port_numeration.py'
spec = importlib.util.spec_from_file_location('under_test', SOURCE)
p = importlib.util.module_from_spec(spec)
spec.loader.exec_module(p)
checks = 0

def require(ok, why):
    global checks
    checks += 1
    if not ok:
        raise RuntimeError(why)

# Generate all Fibonacci terms, rather than the producer's even-subsequence recurrence.
fib = [0, 1]
for _ in range(2050):
    fib.append(sum(fib[-2:]))
A = fib[2::2]

def direct(ds):
    return sum(a*d for a,d in zip(A,ds))

def legal(w, locked=False):
    return set(w) <= set('012') and re.search('21*2', ('2' if locked else '') + w) is None

def trim(ds):
    ds = list(ds)
    while ds and ds[-1] == 0:
        ds.pop()
    return tuple(ds)

def rewrite_choices(ds):
    """All currently enabled positive rewrites, with padded tuples canonicalized."""
    for i,d in enumerate(ds):
        if d >= 3:
            out = list(ds) + [0]
            out[i] -= 3
            out[i+1] += 1
            if i:
                out[i-1] += 1
            yield ('three',i), trim(out)
    for lo in range(len(ds)):
        if ds[lo] != 2:
            continue
        for hi in range(lo+1,len(ds)):
            if ds[hi] == 2:
                out = list(ds) + [0]
                out[lo] -= 2
                out[hi] -= 2
                for j in range(lo+1,hi):
                    out[j] -= 1
                out[hi+1] += 1
                if lo:
                    out[lo-1] += 1
                yield ('block',lo,hi), trim(out)
            if ds[hi] != 1:
                break

@lru_cache(None)
def terminal_forms(ds):
    targets = list(rewrite_choices(ds))
    if not targets:
        word = ''.join(map(str,reversed(ds))) or '0'
        require(legal(word), 'irreducible noncanonical vector')
        return frozenset([ds])
    result = set()
    for _,out in targets:
        require(min(out, default=0) >= 0, 'rewrite nonnegativity')
        require(direct(out) == direct(ds), 'rewrite value preservation')
        require(sum(out) < sum(ds), 'rewrite strict potential')
        result.update(terminal_forms(out))
    return frozenset(result)

# Independent complete DFA, including one absorbing rejection state.
SINK = ('sink',)
START = (0,0,0)  # grammar lock, current numerical parity, shifted parity

def step(q,d):
    if q == SINK or (q[0] and d == 2):
        return SINK
    lock,n,pair = q
    return (int(d == 2 or (lock and d == 1)), (pair+d)%2, (n+pair+d)%2)

def run(q,w):
    for c in w:
        q = step(q,int(c))
    return q

def accept(q):
    return q != SINK and q[1] == 0

def label(q):
    return 'sink' if q == SINK else ('L' if q[0] else 'F') + str(q[1]) + str(q[2])

def main():
    require(A[:7] == [1,3,8,21,55,144,377], 'Fibonacci oracle seed')
    enumerated = 0
    for k in range(12):
        by_state = [[],[]]
        for digits in itertools.product('012', repeat=k):
            w = ''.join(digits)
            enumerated += 1
            v = direct(map(int,reversed(w)))
            require(p.value(w) == v, 'direct vs producer evaluator')
            require(p.valid(w) == legal(w), 'regex vs producer grammar')
            require(p.valid(w,True) == legal(w,True), 'locked regex vs grammar')
            require(accept(run(START,w)) == (legal(w) and v%2 == 0), 'DFA language')
            for lock in (0,1):
                if legal(w,lock):
                    by_state[lock].append(v)
        require(by_state[0] == list(range(A[k])), 'free interval/rank')
        require(by_state[1] == list(range(A[k]-(A[k-1] if k else 0))), 'locked interval/rank')

    limit = 250000
    # Full 2-dimensional coin-type DP. At level j, alternatives have coefficient
    # d_j=0 or d_j>0 and are disjoint. Thus counts refer to unordered coefficient
    # arrays, never permutations/ordered compositions of a summand list.
    best = [[0] + [limit+1]*limit]
    ways = [[1] + [0]*limit]
    coins = [a for a in A if a <= limit]
    for coin in coins:
        prev,pc = best[-1],ways[-1]
        row,cnt = prev.copy(),pc.copy()
        for n in range(coin,limit+1):
            candidate = row[n-coin]+1
            if candidate < row[n]:
                row[n],cnt[n] = candidate,cnt[n-coin]
            elif candidate == row[n]:
                cnt[n] += cnt[n-coin]
        best.append(row)
        ways.append(cnt)
    for n in range(limit+1):
        w = p.encode(n)
        require(legal(w), 'greedy legality')
        require(direct(map(int,reversed(w))) == n, 'greedy value')
        require(sum(map(int,w)) == best[-1][n], 'actual coin minimum')
        require(ways[-1][n] == 1, 'unique minimizing coefficient vector')
        require(p.successor(w) == p.encode(n+1), 'numerical successor')
        require(p.parity(w) == n%2, 'parity value')
        require((n+32)%2 == n%2 and (n+31)%2 != n%2, 'origin conventions')

    # Exhaustive coefficient-array enumeration separately validates small DP
    # counts. Each vector is emitted once, regardless of summand permutations.
    small = 500
    sums = [small+1]*(small+1)
    counts = [0]*(small+1)
    small_coins = [a for a in A if a <= small]
    def enumerate_vectors(i,value,cost):
        if i == len(small_coins):
            if cost < sums[value]:
                sums[value],counts[value] = cost,1
            elif cost == sums[value]:
                counts[value] += 1
            return
        c = small_coins[i]
        for d in range((small-value)//c+1):
            enumerate_vectors(i+1,value+d*c,cost+d)
    enumerate_vectors(0,0,0)
    require(sums == best[-1][:small+1], 'enumerated coefficient-array minima')
    require(counts == ways[-1][:small+1], 'enumerated coefficient-array counts')

    confluence_inputs = 0
    for k in range(6):
        for ds in itertools.product(range(4),repeat=k):
            forms = terminal_forms(trim(ds))
            require(len(forms) == 1, 'all positive rewrite orders agree')
            expected = ''.join(map(str,reversed(next(iter(forms))))) or '0'
            require(p.normalize(ds) == expected, 'producer direct normalization')
            require(expected == p.encode(direct(ds)), 'normalization vs interval normal form')
            confluence_inputs += 1
    rng = random.Random(12513809)
    random_rewrites = 0
    for _ in range(1000):
        ds = trim(rng.randrange(12) for _ in range(rng.randrange(1,50)))
        n = direct(ds)
        while True:
            choices = list(rewrite_choices(ds))
            if not choices:
                break
            _,target = rng.choice(choices)
            require(direct(ds) == direct(target), 'random rewrite value')
            require(sum(target) < sum(ds), 'random rewrite potential')
            ds = target
            random_rewrites += 1
        out = ''.join(map(str,reversed(ds))) or '0'
        require(out == p.encode(n), 'random rewrite strategy normal form')

    for n in range(100001):
        w,z = p.encode(n),p.zeckendorf(n)
        require(set(z) <= set('01') and '11' not in z, 'Zeckendorf grammar')
        require(sum(int(d)*fib[i+2] for i,d in enumerate(reversed(z))) == n, 'Zeckendorf value')
        zlsd = [0,0] + [int(d) for d in reversed(z)] + [0,0,0]
        signed = [zlsd[2*j+2]+zlsd[2*j+1]-zlsd[2*j+3] for j in range((len(z)+2)//2)]
        require(all(c in (-1,0,1) for c in signed), 'signed coefficient range')
        require(direct(signed) == n, 'signed coefficient bridge')
        x,y = p.split_into_zeckendorf(w)
        require('11' not in x and '11' not in y, 'each split is individually canonical')
        require(sum(int(d)*fib[i+2] for i,d in enumerate(reversed(x))) +
                sum(int(d)*fib[i+2] for i,d in enumerate(reversed(y))) == n, 'two-summand bridge')
    for i in range(1000):
        a,b,c = A[i:i+3]
        require(a*a-3*a*b+b*b == 1, 'Pell-type invariant')
        require(a*b+1 == (b-a)**2 and a*c+1 == b*b and b*c+1 == (c-b)**2, 'D(1) triple')
    require(3+8 == 11 and 3*8+1 == 25, 'square-sum vs D(1) separation')

    # Reachability, partition refinement, and all Myhill--Nerode witnesses.
    witness = {START:''}
    todo = deque([START])
    while todo:
        q = todo.popleft()
        for d in range(3):
            target = step(q,d)
            if target not in witness:
                witness[target] = witness[q]+str(d)
                todo.append(target)
    states = list(witness)
    require(len(states) == 9, 'nine reachable complete states')
    blocks = {q:int(accept(q)) for q in states}
    refinement_sizes = [len(set(blocks.values()))]
    while True:
        sigs = {q:(blocks[q],*(blocks[step(q,d)] for d in range(3))) for q in states}
        catalog = {s:i for i,s in enumerate(dict.fromkeys(sigs.values()))}
        refined = {q:catalog[sigs[q]] for q in states}
        refinement_sizes.append(len(catalog))
        if all((blocks[a] == blocks[b]) == (refined[a] == refined[b]) for a in states for b in states):
            break
        blocks = refined
    require(len(catalog) == 9, 'minimal complete DFA')
    separators = []
    for a,b in itertools.combinations(states,2):
        found = None
        for k in range(4):
            for letters in itertools.product('012',repeat=k):
                s = ''.join(letters)
                if accept(run(a,s)) != accept(run(b,s)):
                    found = s
                    break
            if found is not None:
                break
        require(found is not None, 'short distinguishing suffix')
        separators.append({'left':label(a),'right':label(b),'suffix':found})
    certificate = {
        'states':[{'state':label(q),'access_word':witness[q], 'accepting':accept(q),
                   'transitions':{str(d):label(step(q,d)) for d in range(3)}} for q in states],
        'partition_refinement_class_counts':refinement_sizes,
        'pairwise_distinguishing_suffixes':separators,
        'minimal_complete_state_count':9,
        'minimal_live_state_count':8,
        'alphabet':'012','leading_zeros_allowed':True,'empty_word_accepted':True,
    }
    (HERE/'minimal_dfa_certificate.json').write_text(json.dumps(certificate,indent=2)+'\n')
    report = {
        'status':'PASS', 'optimized_python':not __debug__, 'explicit_checks':checks,
        'producer_source_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
        'ternary_words_through_length':11,'ternary_words_enumerated':enumerated,
        'independent_unordered_coin_minimum_and_unique_minimizer_cases':limit+1,
        'direct_coefficient_array_enumeration_value_limit':small,
        'all_rewrite_orders_inputs':confluence_inputs,
        'rewrite_DAG_states':terminal_forms.cache_info().currsize,
        'random_rewrite_strategy_cases':1000,'random_rewrite_steps':random_rewrites,
        'Zeckendorf_bridge_cases':100001,'D1_triples':1000,
        'complete_DFA_states':9,'Myhill_Nerode_distinguished_pairs':len(separators),
        'scope':'Arithmetic/automata audit only; no graph Hamiltonicity inference; no source/remote edits.',
    }
    out = json.dumps(report,indent=2)+'\n'
    print(out,end='')
    if len(sys.argv)>1:
        Path(sys.argv[1]).write_text(out)

if __name__ == '__main__':
    main()
