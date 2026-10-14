#!/usr/bin/env python3
"""Independent finite audit of the break clock, factor four, and gap sidecar."""
from fractions import Fraction
import argparse
import itertools
import json
from pathlib import Path


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def W(n):
    return (0, -1, 1)[n % 3]


def direct_energy(t, x, y):
    h, energy = 0, 0
    for n in range(1, y):
        h += (2 if n % 3 == 2 else -1)*t[n]
        if n >= x:
            energy += h*h
    return energy


def clock_energy(t, x, y):
    # Half-integral v permits zero-valued, as well as binary, labels.
    breaks = [n for n in range(1, len(t)-1) if t[n] != t[n+1]]
    j, energy, start = Fraction(0), Fraction(0), 1
    for stop in breaks + [len(t)-1]:
        left, right = max(x, start), min(y, stop+1)
        if left < right:
            length = right-left
            sw = sum(W(n) for n in range(left, right))
            sw2 = sum(W(n)**2 for n in range(left, right))
            sign = t[start]
            energy += 4*length*j*j-4*j*sign*sw+sign*sign*sw2
        if stop < len(t)-1:
            j += Fraction(W(stop)*(t[stop+1]-t[stop]), 2)
        start = stop+1
    require(energy.denominator == 1, 'integer source energy')
    return int(energy)


def from_breaks(bs, size):
    t = [0]
    sign = 1
    for n in range(1, size+1):
        t.append(sign)
        if n in bs:
            sign = -sign
    return t


def summary(steps):
    j = Fraction(0)
    u = v = Fraction(0)
    for step in steps:
        u += j
        v += j*j
        j += step
    return (len(steps), j, u, v)


def concat(a, b):
    la, da, ua, va = a
    lb, db, ub, vb = b
    return (la+lb, da+db, ua+ub+lb*da, va+vb+2*da*ub+lb*da*da)


def centered(a):
    length, _, u, v = a
    return v-u*u/length


def grammar_controls():
    cuts = triples = 0
    for word in itertools.product((Fraction(-1, 2), Fraction(0), Fraction(1, 2)), repeat=6):
        whole = summary(word)
        require(centered(whole) >= 0, 'positive clock matrix')
        for cut in range(1, 6):
            a, b = summary(word[:cut]), summary(word[cut:])
            require(concat(a, b) == whole, 'exact source-port concatenation')
            la, da, ua, _ = a
            lb, _, ub, _ = b
            innovation = Fraction(la*lb, la+lb)*(ua/la-da-ub/lb)**2
            require(centered(whole) == centered(a)+centered(b)+innovation, 'positive variance refinement')
            cuts += 1
        for left, right in itertools.combinations(range(1, 6), 2):
            a, b, c = summary(word[:left]), summary(word[left:right]), summary(word[right:])
            require(concat(concat(a, b), c) == concat(a, concat(b, c)) == whole, 'ordered associativity')
            triples += 1
    return {'half_integral_words_length_6': 3**6, 'cut_and_variance_checks': cuts,
            'three_block_associativity_checks': triples}


def producer_controls():
    # Exercise the original Python function bodies independently of the optional
    # Numba runtime. This is a finite functional audit, not a full Numba rerun.
    import importlib.util
    import sys
    import types
    import numpy as np
    producer = Path(__file__).with_name('collatz_dyadic_energy_breaks_20261013.py')
    old_numba, old_argv = sys.modules.get('numba'), sys.argv
    shim = types.ModuleType('numba')
    shim.njit = lambda **kwargs: lambda function: function
    sys.modules['numba'] = shim
    sys.argv = [str(producer), '14']
    try:
        spec = importlib.util.spec_from_file_location('breaks_audit_target', producer)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
    finally:
        sys.argv = old_argv
        if old_numba is None:
            del sys.modules['numba']
        else:
            sys.modules['numba'] = old_numba
    nmax = 1 << 14
    p = module.parity_table(nmax)
    breaks = [n for n in range(1024, nmax) if p[n+1] != p[n]]
    gaps = [b-a for a, b in zip(breaks, breaks[1:])]
    short, mods, large, classes, trans, twist, direct = module.gap_census(p, nmax, 1024, nmax)
    require(int(classes.sum()) == len(breaks), 'producer full break count')
    require(int(mods.sum()) == len(gaps), 'producer full gap count')
    require(list(mods) == [sum(g % 3 == r for g in gaps) for r in range(3)], 'producer all gap residues')
    require(list(short) == [gaps.count(g) for g in range(64)], 'producer short histogram')
    require(large == sum(g >= 64 for g in gaps), 'producer full long-gap count')
    synthetic = np.zeros(202, dtype=np.int8)
    synthetic[2:101] = 1  # Artificial binary labels, not a Collatz claim.
    sh, md, lg, cl, *_ = module.gap_census(synthetic, 200, 1, 200)
    require(int(sh.sum()) == 0 and list(md) == [1, 0, 0] and lg == 1 and int(cl.sum()) == 2,
            'artificial 99-long gap survives a short histogram')
    require(np.array_equal(twist, direct), 'producer twist versus direct moment')
    require(int(trans.sum()) == len(gaps), 'producer transition count')
    require(module.identity6_check(p, nmax) == 0, 'producer summation by parts')
    bad, checked = module.check_identities(p, nmax)
    require(bad == 0, 'producer first-step identities')
    for k in range(1, 12):
        x, y = module.transfer_check(p, nmax, k)
        require(x == y, 'producer complete transfer with boundaries')
    native_path = Path(__file__).parents[2]/'05-knowledge/results/collatz_break_clock_20261013d.json'
    native = json.loads(native_path.read_text())
    j = 0
    sums = [[0, 0, 0, 0] for _ in range(14)]
    for n in range(1, nmax):
        k = n.bit_length()-1
        if n == 1 << k:
            sums[k][0] = j
        h = W(n)*(1-2*int(p[n]))-2*j
        sums[k][1] += j
        sums[k][2] += j*j
        sums[k][3] += h*h
        j += W(n)*(int(p[n])-int(p[n+1]))
    for row in native['octaves'][:14]:
        require(sums[row['k']] == [row[key] for key in
                ('incoming_current', 'current_sum', 'current_square_sum', 'H_energy')],
                'independent Python/native octave clock moments')
    # Independent circulant diagonalization tests the all-break denominator.
    markov_checks = 0
    omega = complex(-.5, 3**.5/2)
    for a, b, c in ((1/3, 1/3, 1/3), (.2, .5, .3), (.6, .1, .3)):
        matrix = np.array([[a, b, c], [c, a, b], [b, c, a]])
        spectral = a+b*omega+c*omega**2
        for lag, rho in enumerate(module.renewal_rho(matrix, 12), 1):
            require(abs(rho-(-1)**lag*(spectral**lag).real) < 1e-12, 'Markov normalization versus Fourier eigenvalue')
            markov_checks += 1
    return {'scope': 'Python function bodies, optional Numba decorator bypassed; no full compiled rerun claimed',
            'source_limit': nmax, 'breaks': len(breaks), 'gaps_at_least_64': int(large),
            'artificial_long_gap_control': 99,
            'first_step_indices': checked, 'octave_transfers': 11,
            'native_octave_moment_rows_checked': 14,
            'circulant_markov_checks': markov_checks}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    controls = 0
    for word in itertools.product((-1, 0, 1), repeat=8):
        t = [0] + list(word)
        for x, y in ((1, 9), (4, 8)):
            require(direct_energy(t, x, y) == clock_energy(t, x, y), 'complete plateau identity')
            controls += 1
    a, b = [34, 35, 59, 61], [55, 56, 59, 61]
    ta, tb = from_breaks(a, 64), from_breaks(b, 64)
    currents_a = [W(n)*ta[n+1] for n in a]
    currents_b = [W(n)*tb[n+1] for n in b]
    require(currents_a == currents_b == [1, 1, -1, -1], 'same current sequence')
    require([n % 3 for n in a] == [n % 3 for n in b], 'same residue sequence')
    ea, eb = direct_energy(ta, 32, 64), direct_energy(tb, 32, 64)
    require(ea == clock_energy(ta, 32, 64) and eb == clock_energy(tb, 32, 64), 'hostile clock audit')
    require(ea-eb == 336, 'same event word different source energy')
    iid = []
    for m in range(1, 9):
        total = sum(4*sum(word)**2 for word in itertools.product((-1, 0, 1), repeat=m))
        expectation = Fraction(total, 3**m)
        require(expectation == Fraction(8*m, 3), 'factor-four model law')
        iid.append({'m': m, 'expectation_H_tilde_square': str(expectation)})
    # The source jump identity holds with zero values; the binary shortcut does not.
    require(0-1 != 2*0*(0 != 1), 'ROOT-support type hostile')
    result = {'status': 'FINITE-EXACT controls on arbitrary typed labels; artificial hostiles are not Collatz trajectories',
              'ternary_words_length_8': 3**8, 'plateau_energy_checks': controls,
              'clock_hostile': {'X': 32, 'breaks_A': a, 'breaks_B': b, 'current': currents_a,
                                'energy_A': ea, 'energy_B': eb, 'difference': ea-eb},
              'iid_model': iid, 'source_clock_grammar': grammar_controls(),
              'incoming_producer_audit': producer_controls()}
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result['clock_hostile']))


if __name__ == '__main__':
    main()
