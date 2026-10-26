"""Small exact controls for the recovered cross-thread connections."""
from fractions import Fraction
from math import gcd
from itertools import combinations, permutations
import json
import sys

sys.stdout.reconfigure(encoding='utf-8', newline='\n')

def need(test, message):
    if not test:
        raise RuntimeError(message)

a = [2]
product = 2
for _ in range(7):
    a.append(product + 1)
    product *= a[-1]
need(all(gcd(x,y)==1 for i,x in enumerate(a) for y in a[i+1:]), 'coprimality')
need(sum((Fraction(1,x) for x in a),Fraction()) + Fraction(1,product)==1, 'Egyptian remainder')
orbit = [4]
for _ in range(3):
    orbit.append(orbit[-1]**2-2)
need(orbit == [4,14,194,37634], 'Lucas-Lehmer seed')
cyclic = (2**25 + 4*2**5 + 20*2)//25
elementary = (2**25 + 24*2**5)//25
need((cyclic,elementary)==(1342184,1342208), 'cycle-index observer')
mass = Fraction(2,81)+Fraction(2,27)
need(mass==Fraction(128,1296)==Fraction(8,81), 'nonuniform stabilizer mass')
# Parent's independent five-vertex check: direct matrices and all 120 words,
# with no imports from the producer metagraph code.
def presentation(n, bits):
    pairs=list(combinations(range(n),2))
    free=[(i,j) for i,j in pairs if j>i+1]
    arcs=set()
    for i,j in pairs:
        forward=j==i+1 or bool(bits & (1<<free.index((i,j))))
        arcs.add((i,j) if forward else (j,i))
    return arcs

def complemented_prefix(arcs):
    joined={(0,j) for j in range(1,5)} | {(i+1,j+1) for i,j in arcs}
    return {(i,j) if abs(i-j)==1 else (j,i) for i,j in joined}

def hp(arcs,n):
    return sum(all((word[i],word[i+1]) in arcs for i in range(n-1))
               for word in permutations(range(n)))

port_h=[hp(complemented_prefix(presentation(4,bits)),5) for bits in (2,5)]
need(port_h==[11,13], 'independent rooted composition hostile')
print(json.dumps({'sylvester_terms':a, 'lucas_lehmer_from_4':orbit,
                  'subset_orbits_C25':cyclic, 'subset_orbits_C5_squared':elementary,
                  'observer_fibre_mass':str(mass),
                  'independent_order5_hamilton_counts':port_h}, indent=2))
