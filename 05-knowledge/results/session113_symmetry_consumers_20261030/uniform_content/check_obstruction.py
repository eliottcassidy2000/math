"""Exact all-integer/Fraction check of an epoch-envelope obstruction; no search."""
from fractions import Fraction
from math import gcd
import json
from pathlib import Path

def v5(n):
    if not n != 0:
        raise AssertionError('verification failed')
    h = 0
    while n % 5 == 0:
        h += 1
        n //= 5
    return h

def audit(word, E=None):
    lam = Fraction(5, 8) * (1 - Fraction(1, 5 ** E)) if E else Fraction(5, 8)
    P = Q = 1
    B = 0
    clock = 0
    score = Fraction(0)
    actual = 0
    events = []
    for i, a in enumerate(word, 1):
        P, B, Q = (3 * P, 3 * B + Q, Q * 2 ** a)
        clock += a
        Z = B + Q - P
        h = v5(Z)
        ht = min(h, E) if E else h
        active = (clock + i) % 2 == 0
        if a > 1:
            k = v5(2 ** a - 2)
            kt = min(k, E) if E else k
            score += (max(ht, kt) if active else 0) - lam * a
            events.append({'prefix': i, 'clock': clock, 'a': a, 'height': h, 'active': active, 'return_height': v5(gcd(abs(P - Q), B))})
        else:
            score -= lam - (Fraction(5, 8) * (1 - Fraction(1, 5 ** ht)) if active else 0)
        actual += min(v5(gcd(abs(P - Q), B)), E) if E else v5(gcd(abs(P - Q), B))
    return {'E': E, 'length': len(word), 'clock': clock, 'envelope_excess': str(score), 'actual_content': actual, 'sharp_bound': str(lam * clock), 'events': events}
W = [5] + [1] * 8 + [5] + [1] * 233 + [5] + [1] * 8 + [5]
U = W + [1] * 8 + W
if not (len(W) == 253 and sum(W) == 269):
    raise AssertionError('verification failed')
if not v5(3 ** 9 + 2 ** 13) == 3:
    raise AssertionError('verification failed')
if not v5(2 ** 251 + 3 ** 243) == 4:
    raise AssertionError('verification failed')
if not v5(2 ** (269 + 8) + 3 ** (253 + 8)) == 3:
    raise AssertionError('verification failed')
results = [audit(U, E) for E in (8, 9, 10, 11, None)]
if not results[2]['envelope_excess'] == '1043012/1953125':
    raise AssertionError('verification failed')
if not results[-1]['envelope_excess'] == '479371/312500':
    raise AssertionError('verification failed')
if not all((r['actual_content'] <= Fraction(r['sharp_bound']) for r in results)):
    raise AssertionError('verification failed')
if not Fraction(results[2]['envelope_excess']) > 0:
    raise AssertionError('verification failed')
if __name__ == '__main__':
    out = Path(__file__).with_name('obstruction.json')
    out.write_text(json.dumps(results, indent=2) + '\n')
    for r in results:
        print({k: v for k, v in r.items() if k != 'events'})
