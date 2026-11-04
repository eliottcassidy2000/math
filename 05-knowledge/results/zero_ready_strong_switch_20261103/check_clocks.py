"""Complete two-clock solver, with an independent bounded arithmetic control."""
import json
from pathlib import Path


def v2(x):
    x = abs(x)
    assert x
    return (x & -x).bit_length() - 1


def log_power_two(x):
    return x.bit_length() - 1 if x > 0 and x & (x - 1) == 0 else None


def solve(a, b, d):
    assert a > 0 and b > 0
    alpha, beta = v2(a), v2(b)
    a0, b0 = a >> alpha, b >> beta
    if d == 0:
        return {'infinite': a0 == b0,
                'equation_if_infinite': {'t_minus_s': beta-alpha},
                'solutions': []}
    e = v2(d)
    candidates = set()
    # u=e, so t=e-alpha. No exponent search.
    t = e-alpha
    if t >= 1:
        rhs = a*(1 << t)-d
        s = log_power_two(rhs//b) if rhs % b == 0 else None
        if s is not None and s >= 1:
            candidates.add((s, t))
    # v=e, so s=e-beta.
    s = e-beta
    if s >= 1:
        rhs = d+b*(1 << s)
        t = log_power_two(rhs//a) if rhs % a == 0 else None
        if t is not None and t >= 1:
            candidates.add((s, t))
    # u=v.
    diff = a0-b0
    if diff and d % diff == 0:
        u = log_power_two(d//diff)
        if u is not None:
            s, t = u-beta, u-alpha
            if min(s, t) >= 1:
                candidates.add((s, t))
    for s, t in candidates:
        assert a*(1 << t)-b*(1 << s) == d
    return {'infinite': False, 'solutions': sorted(candidates)}


def main():
    count = 0
    for a in range(1, 10):
        for b in range(1, 10):
            for d in range(-24, 25):
                actual = solve(a, b, d)
                brute = {(s, t) for s in range(1, 9) for t in range(1, 9)
                         if a*(1 << t)-b*(1 << s) == d}
                if actual['infinite']:
                    delta = actual['equation_if_infinite']['t_minus_s']
                    predicted = {(s, t) for s in range(1, 9) for t in range(1, 9)
                                 if t-s == delta}
                else:
                    predicted = {(s, t) for s, t in actual['solutions']
                                 if s <= 8 and t <= 8}
                assert brute == predicted, (a, b, d, brute, actual)
                count += 1
    folder = Path(__file__).resolve().parent
    input_file = 'literal/realized_suffix_depth3_moving1_all_prefixes.json'
    witnesses = json.loads((folder / input_file).read_text())
    results = []
    for index, w in enumerate(witnesses):
        h, f, c = w['suffix_h'], w['suffix_F'], w['offset']
        solved = solve(f, 3*42153, 3*2195+c-h)
        assert solved['solutions'] == [(2, len(w['prefix_names']))]
        k = len(w['prefix_names'])-2
        results.append({'index': index, 'offset': c, 'h': h, 'F': f,
                        'all_positive_clock_solutions': solved,
                        'fixed_shift_k': k,
                        'fixed_shift_offset_constant': h-3*2195,
                        'fixed_shift_offset_coefficient_of_2_to_s': (1 << k)*f-3*42153})
    out = {'status': 'PROVED_SOLVER_WITH_FINITE_EXACT_INPUTS',
           'input_file': input_file,
           'bounded_control_equations': count, 'results': results}
    (folder / 'clock_check.json').write_text(json.dumps(out, indent=2)+'\n')
    print(json.dumps(out, indent=2))


if __name__ == '__main__':
    main()
