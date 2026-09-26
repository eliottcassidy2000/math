"""Independent TP2 trace and rational-bridge audit; standard library only."""
from collections import Counter
from itertools import combinations, product
from math import ceil, exp, floor, log, pi, sin


def check(ok, message):
    if not ok:
        raise RuntimeError(message)


def multiply(a, b):
    return [[sum(a[i][k]*b[k][j] for k in range(len(b)))
             for j in range(len(b[0]))] for i in range(len(a))]


def identity(n):
    return [[int(i == j) for j in range(n)] for i in range(n)]


def trace(a):
    return sum(a[i][i] for i in range(len(a)))


def tp2(a):
    return all(a[i][j]*a[k][l] >= a[i][l]*a[k][j]
               for i, k in combinations(range(len(a)), 2)
               for j, l in combinations(range(len(a[0])), 2))


def ceil_count(q, length):
    e, power = 0, 1
    while power < 2**length:
        e, power = e+1, q*power
    return e


def bridge_matrix(length, ones, width):
    grids = [[n for n in range(1, length*width) if (n+ones*j)%length == 0]
             for j in range(length+1)]
    result = identity(width-1)
    for j in range(length):
        nxt = {n: i for i, n in enumerate(grids[j+1])}
        step = [[0]*len(nxt) for _ in grids[j]]
        for i, n in enumerate(grids[j]):
            for delta in (-ones, length-ones):
                if n+delta in nxt:
                    step[i][nxt[n+delta]] += 1
        check(tp2(step), "rectangular one-step kernel not TP2")
        result = multiply(result, step)
        check(tp2(result), "kernel product not TP2")
    check(grids[0] == grids[-1], "bridge grids fail to close")
    return result


def direct_bridge(q, length, ones, width):
    total, images = 0, Counter()
    for positions in combinations(range(length), ones):
        word = tuple(int(j in positions) for j in range(length))
        centered, prefix = 0, [0]
        for bit in word:
            centered += length*bit-ones
            prefix.append(centered)
        check(centered == 0, "not a rational bridge")
        for start in range(1, width):
            if not all(0 < length*start+x < length*width for x in prefix):
                continue
            total += 1
            cut = min(range(length), key=lambda j: prefix[j])
            rotated = word[cut:]+word[:cut]
            images[rotated] += 1
            e = 0
            for j, bit in enumerate(rotated, 1):
                e += bit
                check(q**e > 2**j, "rotation does not give actual strict survival")
                check(q**e < q**(width+1)*2**j, "rotation peak exceeds width+1")
    check(not images or max(images.values()) <= length*(width-1), "rotation multiplicity")
    check(total <= len(images)*length*(width-1), "image cardinality lower bound")
    return total, len(images)


def centered_count(length, ones, width, start, nonnegative=False):
    states = {start*length: 1}
    for _ in range(length):
        newer = Counter()
        for x, count in states.items():
            for delta in (-ones, length-ones):
                y = x+delta
                if (0 <= y if nonnegative else 0 < y) and y < width*length:
                    newer[y] += count
        states = newer
    return states.get(start*length,0)


def mechanical_word(length, ones, positive):
    prefix = [(ones*j+length-1)//length if positive else ones*j//length
              for j in range(length+1)]
    return [prefix[j+1]-prefix[j] for j in range(length)]


def connector_word(s, t, width, half_time, c):
    a = round(c*half_time+width/2-s)
    u = s+a-c*half_time
    b = round(c*half_time+t-u)
    check(0<a<half_time and 0<b<half_time, "connector degenerate slope")
    check(abs(u+b-c*half_time-t) < 1e-8, "connector phase mismatch")
    first = mechanical_word(half_time,a,u>=s)
    reflected_start, reflected_end = width-t,width-u
    reflected = mechanical_word(half_time,b,reflected_end>=reflected_start)
    word = first+list(reversed(reflected))
    ones, minimum, maximum = 0, width, 0
    for j, bit in enumerate(word,1):
        ones += bit
        x = s+ones-c*j
        check(0 < x < width, "connector numerical geometry")
        minimum, maximum = min(minimum,x),max(maximum,x)
    check(abs(s+ones-2*c*half_time-t)<1e-8, "connector numerical endpoint")
    h = width//8
    log_probability = 0.0
    for odd_count in (a,b):
        count = centered_count(half_time,odd_count,h,0,True)
        check(count>0, "empty one-sided connector family")
        log_probability += log(count)+odd_count*log(c)+(half_time-odd_count)*log(1-c)
    half_bound = (-pi*pi*half_time/(8*h*h)-pi**4*half_time/(135*h**4)
                  -width*width/(c*(1-c)*half_time)-log(half_time*(h-1)))
    check(log_probability >= 2*half_bound-1e-8, "connector probability bound")
    return log_probability-2*half_bound


def connector_checks():
    exact, numeric = 0, 0
    min_probability_margin = float('inf')
    for width in (16,32,64):
        length = floor(width*width/log(width))
        h = width//8
        for numerator,denominator in ((1,4),(1,3),(1,2),(2,3),(3,4)):
            ones = length*numerator//denominator
            positive = centered_count(length,ones,h,0,True)
            tr = sum(centered_count(length,ones,h,start) for start in range(1,h))
            check(positive*length*(h-1) >= tr, "one-sided trace rotation count")
            check(positive>0 and tr>0, "nonempty connector bridge")
            exact += 1
    for q in (3,5,7,31):
        c = log(2)/log(q)
        for width in (32,64):
            length = width**3
            k = floor(width*width/log(width))
            middle_time, middle_width = length-4*k,width-4
            e = ceil(c*middle_time)
            check(q**(e-1)<2**middle_time<q**e, "exact central ceil threshold")
            delta,alpha = e-c*middle_time,(-2*c*k)%1
            final = floor(c*length+width-1)-c*length
            for i in (1,(middle_width-1)//2,middle_width-1):
                x = alpha+1+i
                y = x+delta
                min_probability_margin = min(min_probability_margin,
                    connector_word(0,x,width,k,c),connector_word(y,final,width,k,c))
                check(1 < alpha+1 and alpha+1+middle_width+delta < width-1,
                      "central shifted strip guard")
                numeric += 2
    print("FINITE-EXACT one-sided connector trace-count gates",exact)
    print("VERIFIED phase-compatible connector geometries",numeric,
          "(q=3,5,7,31; M=32,64; three central heights)")
    print("VERIFIED connector probability lower bounds",numeric,
          "minimum log margin",min_probability_margin)


def main():
    accepted, powers = 0, 0
    for flat in product((0,1), repeat=9):
        matrix = [list(flat[3*i:3*i+3]) for i in range(3)]
        if not tp2(matrix):
            continue
        accepted += 1
        current = identity(3)
        for n in range(1,9):
            current = multiply(current, matrix)
            check(trace(current) <= trace(matrix)**n, "TP2 trace-power inequality")
            powers += 1
    hostile = [[0,1],[1,0]]
    check(not tp2(hostile) and trace(multiply(hostile,hostile)) > trace(hostile)**2,
          "periodic non-TP2 hostile control")
    print("FINITE-EXACT binary 3x3 TP2 matrices", accepted, "trace-power gates", powers)
    print("Non-TP2 hostile [[0,1],[1,0]] has trace 0 and spectral radius 1")

    cases, counted, distinct = 0, 0, 0
    minimum_margin = float('inf')
    for q in (3,5,7):
        for length in range(3,12):
            ones = ceil_count(q,length)
            check(0 < ones < length, "degenerate Bernoulli parameter")
            for width in range(2,6):
                matrix = bridge_matrix(length,ones,width)
                direct, images = direct_bridge(q,length,ones,width)
                check(trace(matrix) == direct, "matrix/direct bridge counts differ")
                p, d, theta = ones/length, 1-ones/length, pi/width
                beta = log(d*sin(p*theta)/(p*sin(d*theta)))
                mu = d*exp(-p*beta)*sin(theta)/sin(d*theta)
                weighted_trace = direct*p**ones*d**(length-ones)
                gap = weighted_trace-mu**length
                minimum_margin = min(minimum_margin,gap)
                check(gap >= -1e-12, "numerical sine/trace lower comparison")
                cases, counted, distinct = cases+1, counted+direct, distinct+images
    print("FINITE-EXACT rational bridges", cases, "closed paths", counted,
          "distinct rotated images summed per universe", distinct)
    print("VERIFIED sine/trace comparisons", cases, "smallest absolute margin", minimum_margin)
    connector_checks()
    print("PASS: zero patterns, strict survival, grid closure, rotation multiplicity audited")


if __name__ == '__main__':
    main()
