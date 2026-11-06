"""Independent symbolic word expansion, binary-pair correction and inversion controls."""
from fractions import Fraction
import json, random

def need(q, why):
    if not q:
        raise RuntimeError(why)

def g(u):
    total = 1
    while u > 1:
        if u & 1:
            u = (u + 1) // 2
            total += u - (u & 1)
        else:
            u //= 2
    return total

def pair_sum(k):
    ones = [i for i in range(k.bit_length()) if (k >> i) & 1]
    return sum((1 << (j-i-1)) + (j == i+1)
               for i in ones for j in ones if i < j)

def value(seed, u):
    X, e, V, d = seed
    return X-V+(V-d)*u+d*g(u)+(d-e)*(1-u%2)

def literal(seed, word):
    X, e, V, d = seed
    for bit in word:
        if bit:
            X, e, V = X+V-e, d, 2*V-d
        else:
            e, V = 0, 2*V
    return X

def add(a,b): return tuple(x+y for x,y in zip(a,b))
def sub(a,b): return tuple(x-y for x,y in zip(a,b))
def twice(a): return tuple(2*x for x in a)
unit_d = (0,0,0,1)
qs = [pair_sum(k) for k in range(65536)]
for k, q in enumerate(qs):
    need(q == k+1-k%2-g(k+1), 'binary-pair correction equals g expression')
    need(g(k+1) >= 1 and g(k+1)%2 == 1, 'positive odd g')
stack = [(0,0,(1,0,0,0),(0,1,0,0),(0,0,1,0))]
symbolic = 0
while stack:
    depth,k,X,V,e = stack.pop()
    need(X == (1,k,-(k&1),-qs[k]), 'full symbolic word expansion')
    symbolic += 1
    if depth < 16:
        stack.append((depth+1,k,X,twice(V),(0,0,0,0)))
        stack.append((depth+1,k+(1<<depth),sub(add(X,V),e),sub(twice(V),unit_d),unit_d))
seeds = [(677,27,1061,49),(759,29,927,29),(1547,101,1691,109),(821,33,1079,33)]
controls = seeds + [(0,0,2,1),(0,1,2,1),(0,2,6,5),(0,5,6,5)]
rng = random.Random(203)
long_cases = 0
for length in (0,1,2,3,7,16,31,64,127,256):
    words = [[0]*length,[1]*length,[i%2 for i in range(length)]]
    words += [[rng.randrange(2) for _ in range(length)] for _ in range(8)]
    for w in words:
        k = sum(b<<i for i,b in enumerate(w))
        need(pair_sum(k) == k+1-k%2-g(k+1), 'long correction identity')
        for row in controls:
            X = literal(row,w)
            need(X == value(row,k+1), 'long direct word')
            need(X == literal(row,w+[0,0,0]), 'trailing O neutrality')
            need(X >= row[0]+(row[2]-row[3])*k, 'lower growth bound')
            long_cases += 1
inverse_checks = 0
for row in seeds:
    X0,e,V,d = row
    need(value(row,1) == X0, 'empty word')
    need(literal(row,[0,1]) != literal(row,[1,0]), 'chronology hostile')
    # V0 is not a safe replacement for V0-d, already at the first Z.
    need(1+(value(row,2)-X0)//V < 2, 'unsafe V0 denominator hostile')
    values = {value(row,u) for u in range(1,2049)}
    for target in range(X0-3,X0+20001):
        upper = max(0,1+(target-X0)//(V-d))
        found = any(value(row,u)==target for u in range(1,upper+1))
        need(found == (target in values), 'complete bounded inverse')
        inverse_checks += 1
    for N in (0,X0-1,X0,X0+V,10000,100000,1000000):
        actual = sum(x<=N for x in values)
        upper = max(0,1+(N-X0)//(V-d))
        need(actual<=upper, 'finite counting bound')
    for u in range(1,65537):
        need(value(row,u)%2 == 1, 'odd physical bases')
# The lower bound requires e0<=d; dropping this hypothesis really can break it.
need(value((0,2,2,1),2) < 1, 'e0>d hostile')
bound = 2*sum((Fraction(1,V-d) for X,e,V,d in seeds),Fraction())
need(bound == Fraction(693609731,93988567442) < Fraction(3,400), 'exact odd-relative density')
print(json.dumps(dict(status='PASS',symbolic_words=symbolic,maximum_symbolic_length=16,
    binary_pair_identities=65536,long_direct_cases=long_cases,maximum_long_length=256,
    complete_inverse_targets=inverse_checks,odd_base_formula_checks=4*65536,
    odd_relative_upper_density=str(bound),scope='Four O/Z carrier orbits only; no global occupancy.'),sort_keys=True,indent=2))
