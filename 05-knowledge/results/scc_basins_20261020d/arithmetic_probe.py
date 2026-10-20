"""Exact SCC/Fibonacci codec, deficit hostiles, and marked golden junctions."""
from collections import Counter
from functools import lru_cache
import hashlib
from itertools import combinations, product
import json
from math import comb, factorial, gcd
from pathlib import Path
import sympy as sp

HERE = Path(__file__).resolve().parent


def need(condition, message):
    if not condition:
        raise RuntimeError(message)


@lru_cache(None)
def words(n, alphabet):
    if n == 0:
        return ((),)
    return tuple((a,) + tail for a in alphabet if a <= n
                 for tail in words(n-a, alphabet))


def parse_fibonacci(word):
    result, at = [], 0
    while at < len(word):
        if word[at] == 1:
            result.append(1)
            at += 1
        elif at + 1 == len(word):
            return tuple(result), 1
        else:
            result.append(3 if word[at+1] == 1 else 4)
            at += 2
    return tuple(result), 0


def decode_fibonacci(word, flag):
    token = {1:(1,), 3:(2,1), 4:(2,2)}
    return tuple(v for a in word for v in token[a]) + ((2,) if flag else ())


def tournament(word):
    atoms = {1:(0,), 3:(2,4,1), 4:(6,4,8,3)}
    out, offset = [], 0
    n = sum(word)
    for a in word:
        later = ((1 << n)-1) ^ ((1 << (offset+a))-1)
        out.extend((row << offset) | later for row in atoms[a])
        offset += a
    return tuple(out)


def path_inventory(out):
    n, total, blue = len(out), 0, 0
    full = (1 << n)-1

    def visit(path, used):
        nonlocal total, blue
        if used == full:
            total += 1
            r = dict(zip(path, path[::-1]))
            blue += all(bool(out[i] >> j & 1) != bool(out[r[i]] >> r[j] & 1)
                        for i in range(n) for j in range(i+1,n))
            return
        allowed = (out[path[-1]] if path else full) & (full ^ used)
        while allowed:
            bit = allowed & -allowed
            allowed -= bit
            visit(path+(bit.bit_length()-1,), used | bit)
    visit((), 0)
    return total, blue


def canonical_word(word):
    rename = {}
    return tuple(rename.setdefault(c, len(rename)) for c in word)


def main():
    z,y = sp.symbols('z y')
    q = 1-z-z**3
    fden = 1-z-z**3-z**4
    need(sp.expand((1+z**2)*(1-z-z**2)) == fden, 'golden factorization')
    need(sp.expand((1+y*z**2)*(1-z-y*z**2)) == 1-z-y*z**3-y*y*z**4,
         'triangle-potential-preserving factorization')
    need(sp.cancel((z+3)/(5*(1-z-z**2))+(z+2)/(5*(1+z**2))-1/fden) == 0,
         'partial fractions')
    counts = [0]*201
    counts[0] = 1
    for n in range(1,len(counts)):
        counts[n] = sum(counts[n-a] for a in (1,3,4) if n >= a)
    fib, lucas = [0,1], [2,1]
    for n in range(2,204):
        fib.append(fib[-1]+fib[-2])
        lucas.append(lucas[-1]+lucas[-2])
    c = lambda n: 0 if n < 0 or n % 2 else (-1)**(n//2)
    for n, count in enumerate(counts):
        need(count+(counts[n-2] if n >= 2 else 0) == fib[n+1], 'Fibonacci count')
        need(5*count == lucas[n+2]+2*c(n)+c(n-1), 'Lucas exact formula')
    rigid = [0]*201
    rigid[0] = 1
    for n in range(1,len(rigid)):
        rigid[n] = sum(rigid[n-a] for a in (1,3) if n >= a)
    terminal_counts = []
    for n in range(17):
        palindromes = sum(counts[(n-a)//2] for a in (0,1,3,4)
                          if n >= a and (n-a)%2 == 0)
        pure = sum(rigid[(n-a)//2] for a in (0,1,3)
                   if n >= a and (n-a)%2 == 0)
        expected = (counts[n]+palindromes)//2-pure-int(n == 4)
        reps = {min(word,word[::-1]) for word in words(n,(1,3,4))}
        actual = sum(len(word) >= 2 and
                     (4 in word or word != word[::-1]) for word in reps)
        need(actual == expected, 'restricted terminal Burnside formula')
        terminal_counts.append(actual)
    codec_tests = 0
    for n in range(17):
        images = set()
        for word in words(n,(1,2)):
            image = parse_fibonacci(word)
            need(decode_fibonacci(*image) == word, 'literal codec inverse')
            need(sum(image[0])+2*image[1] == n, 'codec weight')
            need(word.count(2) == image[0].count(3)+2*image[0].count(4)+image[1],
                 'codec triangle potential and terminal flag')
            images.add(image)
            codec_tests += 1
        expected = {(word,0) for word in words(n,(1,3,4))}
        if n >= 2:
            expected |= {(word,1) for word in words(n-2,(1,3,4))}
        need(images == expected, 'codec surjectivity')

    coefficient_rows = []
    for r in range(5):
        series = [0]*81
        for n in range(4*r,81):
            for t in range((n-4*r)//3+1):
                singles = n-4*r-3*t
                series[n] += factorial(singles+t+r)//(
                    factorial(singles)*factorial(t)*factorial(r))
        # Independent recurrence from multiplication by the exact denominator.
        denominator = sp.Poly(q**(r+1),z)
        for n in range(81):
            lhs = sum(int(coefficient)*series[n-int(exponent[0])]
                      for exponent,coefficient in denominator.terms()
                      if exponent[0] <= n)
            need(lhs == int(n == 4*r), 'fixed insertion generating function')
        coefficient_rows.append(dict(insertions=r, coefficients=series))

    inventory_tests, deficit_examples = 0, []
    for n in range(1,13):
        for word in words(n,(1,3,4)):
            out = tournament(word)
            total, blue = path_inventory(out)
            r, triangles = word.count(4), word.count(3)
            auto = 3**triangles
            h = 5**r
            palindrome = word == word[::-1]
            b = 0 if not palindrome else 5**(r//2)*(
                3 if len(word)%2 and word[len(word)//2] == 4 else 1)
            need(total == auto*h and blue == auto*b, 'literal SCC h,b formula')
            need(comb(n,3)-sum(comb(row.bit_count(),2) for row in out) == triangles+2*r,
                 'literal directed-triangle potential')
            inventory_tests += 1
            if all(a == 4 for a in word):
                deficit_examples.append(dict(word=word,h=h,b=b,beta=h-b))
    # Four R4 blocks are still cheap but outside the complete order<=12 universe.
    total, blue = path_inventory(tournament((4,4,4,4)))
    need((total,blue) == (625,25), 'quartic deficit hostile')
    deficit_examples.append(dict(word=(4,4,4,4),h=total,b=blue,beta=total-blue))

    boundary_words = sorted({canonical_word(w) for w in product(range(3),repeat=7)
                             if all(w.count(c)%2 for c in range(3))})
    rows, flags = [], []
    for direction in range(1,8):
        matching = sorted({tuple(sorted((v,v^direction)))
                           for v in range(1,8) if v != direction})
        for removed in matching:
            triple = tuple(sorted((direction,)+removed))
            wires = [e for e in matching if e != removed]
            need(triple[0]^triple[1]^triple[2] == 0, 'Fano line')
            flags.append(dict(direction=direction,triple=triple,wires=wires))
            rows.append([int(len({w[i-1] for i in triple}) == 3 and
                             all(w[i-1] == w[j-1] for i,j in wires))
                         for w in boundary_words])
    j = sp.Matrix(rows)
    p = sp.Matrix([[sum(rows[3*d+k][c] for k in range(3))
                    for c in range(91)] for d in range(7)])
    w = sp.Matrix([[(1 if sorted(Counter(word).values()) == [1,1,5] else -1)
                    if word.count(word[i]) == 1 else 0 for word in boundary_words]
                   for i in range(7)])
    need(j.rank() == 21 and p.rank() == 7, 'marked Fano bank ranks')
    need(p*p.T == 20*sp.eye(7)+7*sp.ones(7), 'parallel bank Gram')
    need(j*w.T == sp.zeros(21,7), 'golden bank misses every Kempe direction')
    radical = [1,3,5,15,7,21,35,105]
    golden_flags = [dict(triple=[radical[v] for v in f['triple']],
                         wires=[[radical[v] for v in e] for e in f['wires']])
                    for f in flags if f['direction'] == 2]

    # Replay the inherited distinction: scalar modulus versus cyclotomic ideal.
    pair, scalar_period = (1,0), 0
    while True:
        pair = (pair[1] % 105, (pair[0]+pair[1]) % 105)
        scalar_period += 1
        if pair == (1,0):
            break
    need(scalar_period == 80, 'golden scalar modulus105 period')
    remainder = sp.rem(sp.cyclotomic_poly(105,z),z*z-z-1,z)
    aa, bb = int(remainder.coeff(z,0)), int(remainder.coeff(z,1))
    norm = aa*aa+aa*bb-bb*bb
    need((aa,bb,norm) == (5224431949,8453308464,16271615641), 'inherited ideal')
    need(gcd(bb,norm) == 1, 'ideal quotient generator')
    root = (-aa*pow(bb,-1,norm)) % norm
    need((root*root-root-1) % norm == 0 and pow(root,105,norm) == 1,
         'ideal golden relation and order divides105')
    need(all(pow(root,105//prime,norm) != 1 for prime in (3,5,7)), 'exact order105')
    conjugate = (1-root) % norm
    need(pow(conjugate,210,norm) == 1 and
         all(pow(conjugate,210//prime,norm) != 1 for prime in (2,3,5,7)),
         'exact conjugate order210')

    receipt = dict(
        universe=dict(codec_orders=list(range(17)),codec_words=codec_tests,
                      literal_scc_orders=list(range(1,13)),scc_words=inventory_tests,
                      extra_hostile_word=[4,4,4,4],boundary_colour_orbits=91,
                      fano_flags=21),
        word_counts=counts, fixed_insertions=coefficient_rows,
        restricted_merged_terminal_counts=terminal_counts,
        R4_deficits=deficit_examples, golden_character_flags=golden_flags,
        bank_ranks=dict(fano_junctions=21,parallel_sums=7),
        junction_gram_polynomial=str(sp.factor((j*j.T).charpoly().as_expr())),
        parallel_gram='20 I + 7 J; eigenvalues69 once,20 six times',
        killed_Kempe_dimension=7,
        bivariate_identity='1/(1-z-y*z^3-y^2*z^4)=1/((1+y*z^2)*(1-z-y*z^2)); y counts directed triangles',
        inherited_golden_clocks=dict(scalar_modulus=105,scalar_period=80,
                                     ideal_coefficients=[aa,bb],norm=norm,
                                     ideal_root=root,ideal_period=105,
                                     conjugate_period=210),
        source_sha256_LF=hashlib.sha256(Path(__file__).read_bytes().replace(b'\r\n',b'\n')).hexdigest())
    Path(__file__).with_suffix('.json').write_text(json.dumps(receipt,indent=2)+'\n',
                                                  encoding='utf-8',newline='\n')
    print('codec tests:',codec_tests,'literal SCC tests:',inventory_tests)
    print('word coefficients n0..16:',counts[:17])
    print('R4 deficits:',deficit_examples)
    print('golden character flags:',golden_flags)
    print('ranks:',receipt['bank_ranks'],'Kempe blind:',7)
    print('golden periods: scalar105 ->80; cyclotomic ideal ->105; conjugate ->210')
    print('PASS exact codecs, recurrences, boundary ranks, and distinct arithmetic clocks.')


if __name__ == '__main__':
    main()
