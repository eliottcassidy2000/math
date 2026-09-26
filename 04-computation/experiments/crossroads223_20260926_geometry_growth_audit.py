"""Independent exact audit of local carry banks and the no-dip cutoff."""
from fractions import Fraction
from itertools import combinations
from math import comb, factorial


def check(ok,message):
    if not ok:
        raise RuntimeError(message)


def words(length,budget):
    if length == 0:
        yield ()
    else:
        for k in range(1,budget-length+2):
            for tail in words(length-1,budget-k):
                yield (k,)+tail


def factor(n):
    answer = {}
    p = 2
    while p*p <= n:
        while n%p == 0:
            answer[p] = answer.get(p,0)+1
            n //= p
        p = 3 if p == 2 else p+2
    if n>1:
        answer[n] = answer.get(n,0)+1
    return answer


def odd_steps(n,length):
    nodes,exponents = [n],[]
    for _ in range(length):
        numerator = 3*n+1
        k = (numerator & -numerator).bit_length()-1
        n = numerator//2**k
        exponents.append(k)
        nodes.append(n)
    return nodes,tuple(exponents)


def coefficients(word):
    carry,scale = 0,1
    for k in word:
        carry,scale = 3*carry+scale,scale*2**k
    return carry,scale


def main():
    universes,realizations,pairs,collisions = 0,0,0,0
    for length in (2,3,4):
        seen = 0
        bound = 96**comb(length,2)
        a,fact = 0,1
        while fact*(a+1) < bound:
            a,fact = a+1,fact*(a+1)
        check(factorial(a)<bound<=factorial(a+1),"strict factorial inverse")
        for word in words(length,3*length):
            seen += 1
            total = sum(word)
            carry,scale = coefficients(word)
            modulus = 2*scale
            source = ((scale-carry)*pow(3**length,-1,modulus))%modulus
            check(source>0 and source%2 == 1,"exact valuation cylinder source")
            base,obtained = odd_steps(source,length)
            check(obtained == word,"exact cylinder valuations")
            check(odd_steps(source+scale,length)[1] != word,"missing final parity bit hostile")
            partial = [0]
            for k in word:
                partial.append(partial[-1]+k)
            slopes = [3**i*2**(total+1-partial[i]) for i in range(length)]
            banks = [{2,3} for _ in range(length)]
            determinants = {}
            for i,j in combinations(range(length),2):
                subcarry,_ = coefficients(word[i:j])
                determinant = slopes[j]*base[i]-slopes[i]*base[j]
                check(determinant == -3**i*2**(total+1-partial[j])*subcarry,
                      "incident determinant identity")
                check(determinant != 0 and subcarry%2 and subcarry%3,"nonzero coprime carry")
                primes = set(factor(subcarry))
                banks[i] |= primes
                banks[j] |= primes
                determinants[i,j] = determinant
            for t in (0,1,17):
                nodes,actual = odd_steps(source+modulus*t,length)
                check(actual == word,"high cylinder realization")
                node_primes = [set(factor(n)) for n in nodes[:-1]]
                for i in range(length):
                    check(nodes[i] == slopes[i]*t+base[i],"affine node")
                    for p in node_primes[i]-banks[i]:
                        check(all(p not in node_primes[j] for j in range(length) if j != i),
                              "outside incident bank prime is not private")
                for i,j in combinations(range(length),2):
                    determinant = determinants[i,j]
                    x = Fraction(slopes[j]*nodes[i],determinant)
                    y = Fraction(-slopes[i]*nodes[j],determinant)
                    check(x+y == 1,"unit equation identity")
                    if node_primes[i] <= banks[i] and node_primes[j] <= banks[j]:
                        check(set(factor(abs(x.numerator)))|set(factor(x.denominator)) <= banks[i],
                              "first coordinate bank")
                        check(set(factor(abs(y.numerator)))|set(factor(y.denominator)) <= banks[j],
                              "second coordinate bank")
                    if nodes[i] == nodes[j]:
                        check(determinant%nodes[i] == 0,"collision determinant divisor")
                        collisions += 1
                    pairs += 1
                realizations += 1
        check(seen == comb(3*length,length),"bounded composition count")
        universes += seen
    print("FINITE-EXACT cylinders",universes,"realizations",realizations,
          "pair identities",pairs,"collision controls",collisions)

    intervals = 0
    for n in (1,2,3,4):
        length = 3*n
        modulus = 2**(length+1)
        tail = sum(comb(length,j) for j in range(n))
        observed = sum(sum(odd_steps(x,n)[1])>length for x in range(1,modulus,2))
        check(observed == tail,"odd source L+1-bit tail")
        for start in (1,17,2**20+223):
            for width in (17,113,1000):
                count = sum(sum(odd_steps(x,n)[1])>length
                            for x in range(start,start+width) if x%2)
                check(count*modulus <= (width+modulus)*tail,"uniform interval residue bound")
                intervals += 1
    print("FINITE-EXACT Terras normalizations 4; translated interval gates",intervals)

    hard,cutoffs = 0,0
    for source in range(3,65536,2):
        horizon = source.bit_length()-1
        value = source
        for _ in range(horizon):
            value = (3*value+1)//2 if value%2 else value//2
            if value<source:
                break
        else:
            hard += 1
            for length in range(2,9):
                ell = (3**length).bit_length()+1
                if source >= 2**ell:
                    check(sum(odd_steps(source,length)[1]) <= ell,"conditioned no-dip cutoff")
                    cutoffs += 1
    check(sum(odd_steps(1,5)[1]) > (3**5).bit_length()+1,"small-source hostile control")
    print("FINITE-EXACT no-dip sources below65536",hard,"eligible cutoff gates",cutoffs)
    print("Hostile n=1,N=5 violates cutoff without the source-size condition")
    print("PASS: local banks, factorial inverse, odd normalization, uniformity, no-dip direction")


if __name__ == '__main__':
    main()
