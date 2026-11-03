#!/usr/bin/env python3
"""Exact-rational MRB precision audit; only the Python standard library.

Analytic input: the inherited Euler/Cauchy remainder bound stated in README.md.
Every evaluated log/exp is enclosed by rational series bounds. No numerical
library's directed-rounding behavior is assumed. All checks survive python -O.
"""
from fractions import Fraction as F
from math import comb, factorial
from pathlib import Path
import argparse
import json
import time


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def validate_weights(weights, N):
    require(len(weights) == N, 'wrong number of Euler weights')
    require(all(0 < w < 2**N for w in weights), 'weight range failed')
    require(sum(weights) == N*2**(N-1), 'weight L1 mass failed')
    for k in range(N):
        direct = sum(comb(j,k)*2**(N-j-1) for j in range(k,N))
        require(weights[k] == direct, 'independent weight identity failed')


def validate_term(flo, fhi, d, q):
    require(flo <= fhi, 'reversed elementary-function enclosure')
    fd, eps = F(d, 2**q), F(1, 2**q)
    require(abs(fd-flo) < eps and abs(fd-fhi) < eps,
            'dyadic term approximation failed its exact certificate')


def hostile_controls():
    controls = [
        ('nonpositive_precision', lambda: certificate(0)),
        ('negative_atanh_argument', lambda: atanh_log_lower(F(-1,10), 3)),
        ('large_atanh_argument', lambda: atanh_log_lower(F(1,2), 3)),
        ('negative_exp_argument', lambda: exp_lower(F(-1,10), 3)),
        ('large_exp_argument', lambda: exp_lower(F(2), 3)),
        ('corrupted_weight', lambda: validate_weights([3,2], 2)),
        ('corrupted_rounding', lambda: validate_term(F(1,3), F(1,3), 3, 2)),
        ('undersized_precision', lambda: validate_term(F(0), F(1), 2, 2)),
    ]
    passed = []
    for name, operation in controls:
        try:
            operation()
        except RuntimeError:
            passed.append(name)
        else:
            raise RuntimeError('hostile input was accepted: '+name)
    return passed


def ceil_log2(n):
    require(isinstance(n, int) and n >= 1, 'ceil_log2 requires positive integer')
    return (n - 1).bit_length()


def atanh_log_lower(z, J):
    """2 sum_{j<J} z^(2j+1)/(2j+1), for rational 0<=z<=1/3."""
    require(F(0) <= z <= F(1, 3), 'atanh argument out of certified range')
    power, square, total = z, z*z, F(0)
    for j in range(J):
        total += 2*power/(2*j+1)
        power *= square
    return total


def exp_lower(x, T):
    require(F(0) <= x <= F(1), 'exp argument out of certified range')
    term, total = F(1), F(1)
    for j in range(1, T):
        term *= x/j
        total += term
    return total


def nearest_integer(x):
    """Round a rational to nearest integer; a tie goes upward."""
    return (2*x.numerator + x.denominator)//(2*x.denominator)


def certified_terms(max_n, q):
    """Yield integer d_n with |d_n/2^q - (n^(1/n)-1)| < 2^-q."""
    max_k = max_n.bit_length()-1
    J = q + ceil_log2(max_k+1) + 4
    T = q + 6
    log2_lo = atanh_log_lower(F(1, 3), J)
    log_remainder = F(1, 9**J)
    exp_remainder = F(2, factorial(T))
    eps = F(1, 2**q)
    largest_width = F(0)
    max_intermediate_bits = 0
    values = []
    for n in range(1, max_n+1):
        if n == 1:
            values.append(0)
            continue
        k = n.bit_length()-1
        scale = 2**k
        z = F(n-scale, n+scale)
        log_lo = k*log2_lo + atanh_log_lower(z, J)
        log_hi = log_lo + (k+1)*log_remainder
        xlo, xhi = log_lo/n, log_hi/n
        flo = exp_lower(xlo, T)-1
        fhi = exp_lower(xhi, T)-1+exp_remainder
        d = nearest_integer((flo+fhi)*2**q/2)
        validate_term(flo, fhi, d, q)
        values.append(d)
        largest_width = max(largest_width, fhi-flo)
        max_intermediate_bits = max(max_intermediate_bits,
            flo.numerator.bit_length(), flo.denominator.bit_length(),
            fhi.numerator.bit_length(), fhi.denominator.bit_length())
    return values, J, T, largest_width, max_intermediate_bits


def certificate(p):
    require(isinstance(p, int) and p >= 1, 'p must be a positive integer')
    start = time.perf_counter()
    N, M = p+4, 4*(p+4)
    q = p+ceil_log2(M+N)+3
    # W_k = sum_{j=k+1}^N binom(N,j). O(N) exact binomial recurrence.
    weights = []
    binomial, running = 1, 2**N
    for k in range(N):
        running -= binomial
        weights.append(running)
        numerator = binomial*(N-k)
        require(numerator % (k+1) == 0, 'nonintegral binomial recurrence')
        binomial = numerator//(k+1)
    validate_weights(weights, N)
    values, J, T, widest_term, max_bits = certified_terms(M+N-1, q)
    require(widest_term < F(1, 2**(q+3)), 'elementary interval width bound failed')
    prefix = sum(((-1)**n)*values[n-1] for n in range(1,M))
    tail = sum(((-1)**k)*weights[k]*values[M+k-1] for k in range(N))
    central = F(2**N*prefix + (-1)**M*tail, 2**(q+N))
    # Conservative exact bound deduced from the inherited analytic bound.
    analytic_radius = F(64, 4**N)
    amplification = F(M-1)+F(N,2)
    arithmetic_radius = amplification/F(2**q)
    output_integer = nearest_integer(central*2**(p+2))
    output = F(output_integer, 2**(p+2))
    radius = analytic_radius+arithmetic_radius+abs(output-central)
    require(analytic_radius <= F(1, 2**(p+3)), 'universal analytic bound failed')
    require(arithmetic_radius < F(1, 2**(p+3)), 'arithmetic amplification failed')
    require(abs(output-central) <= F(1, 2**(p+3)), 'final rounding bound failed')
    require(radius < F(1, 2**p), 'final precision target failed')
    require(radius < F(3, 2**(p+3)), 'stronger total radius failed')
    # Report exact fractions: the [lo, hi] interval is an analytic certificate.
    return {
        'p': p, 'M': M, 'N': N, 'q': q,
        'log_series_terms': J, 'exp_series_terms': T,
        'output_dyadic_numerator': str(output_integer),
        'output_dyadic_denominator_power': p+2,
        'C_lower': str(output-radius), 'C_upper': str(output+radius),
        'certified_radius': str(radius),
        'target_error': str(F(1, 2**p)),
        'radius_times_2_to_p': float(radius*2**p),
        'analytic_radius': str(analytic_radius),
        'term_error_amplification': str(amplification),
        'max_elementary_interval_width_less_than': str(F(1, 2**(q+3))),
        'max_final_elementary_fraction_bits': max_bits,
        'seconds': round(time.perf_counter()-start, 6),
    }


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--precisions', type=int, nargs='+', default=[1,8,16,32])
    parser.add_argument('--output', type=Path)
    parser.add_argument('--hostile-only', action='store_true')
    args=parser.parse_args()
    hostile=hostile_controls()
    print('PASS: '+str(len(hostile))+' hostile inputs rejected', flush=True)
    if args.hostile_only:
        if args.output:
            args.output.write_text(json.dumps({'status':'PASS', 'hostile_controls':hostile}, indent=2)+'\n')
        return
    rows=[]
    previous=None
    for p in args.precisions:
        row=certificate(p)
        lo,hi=F(row['C_lower']),F(row['C_upper'])
        if previous is not None:
            require(max(previous[0],lo) <= min(previous[1],hi),
                    'independent precision enclosures do not overlap')
        previous=lo,hi
        rows.append(row)
        print('PASS p={p} M={M} N={N} q={q}; radius/2^-p={radius_times_2_to_p:.8f}; '
              'seconds={seconds}'.format(**row), flush=True)
    result={
        'status':'PASS: exact rational elementary enclosures + inherited analytic Euler bound',
        'optimized_python':not __debug__,
        'no_floating_point_used_for_certification':True,
        'hostile_controls':hostile,
        'certificates':rows,
    }
    if args.output:
        args.output.write_text(json.dumps(result,indent=2)+'\n')
    print('ALL EXPLICIT CHECKS PASSED', flush=True)


if __name__ == '__main__':
    main()
