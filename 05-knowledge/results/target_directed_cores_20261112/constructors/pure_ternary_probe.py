"""Search for a certified ready fixed-module K*3^m-1 family, not a dyadic fit."""
from itertools import product
from pathlib import Path
import json
import subprocess
from paired_scale_probe import evaluate, inflate

HERE = Path(__file__).resolve().parent
COMPILER = '05-knowledge/results/session296_multicore_metagraph_20261113/cores/verify_ready.py'


def main():
    ns = {'__name__': 'inherited_funnel_compiler'}
    exec(compile(subprocess.check_output(['git', 'show', 'HEAD:' + COMPILER], text=True),
                 COMPILER, 'exec'), ns)
    bank = {1: [(ns['SINGLE'], ())]}
    cases, metadata = [], []
    for n in range(3, 10):
        bank[n] = []
        for m in (1, *range(3, n)):
            t = n - m - 1
            words = [''] if t == 0 else [''.join(x) + '1' for x in product('01', repeat=t - 1)]
            for S, address in bank[m]:
                for word in words:
                    N = ns['advance'](S, word)
                    new_address = address + (word,)
                    bank[n].append((N, new_address))
                    h = N[0] + N[4] - 1
                    if h % 3 != 2:
                        continue
                    rows = ns['make'](new_address)
                    for vertex in range(n - 1):
                        sizes = [1] * len(rows); sizes[vertex] = 2
                        child, groups = inflate(rows, sizes)
                        cases.append((child, -1))
                        metadata.append({'parent_order': n, 'address': new_address, 'H1': h,
                                         'rows': rows, 'module': vertex})
    candidates = []
    for data, (h2, _, _) in zip(metadata, evaluate(cases)):
        if h2 == 3 * data['H1'] + 2:
            candidates.append({**data, 'H2': h2})
    profiles = []
    for data in candidates:
        for size in (2, 3, 4):
            sizes = [1] * len(data['rows']); sizes[data['module']] = size
            child, groups = inflate(data['rows'], sizes)
            profiles.append((child, groups[-1][0]))
    answers = iter(evaluate(profiles))
    for data in candidates:
        data['profiles_sizes2_3_4'] = [next(answers) for _ in range(3)]
        data['pure_ternary_ready_certificate'] = all(
            h == (data['H1'] + 1) * 3**(size - 1) - 1 and q - s == 1
            for size, (h, q, s) in zip((2, 3, 4), data['profiles_sizes2_3_4']))
        h1=data['H1'];h2,h3,h4=[p[0] for p in data['profiles_sizes2_3_4']]
        a1=h1;a2=(h2-h1)//2;a3=(h3-a1-6*a2)//6;a4=(h4-a1-14*a2-36*a3)//24
        assert h2==a1+2*a2 and h3==a1+6*a2+6*a3 and h4==a1+14*a2+36*a3+24*a4
        data['nonempty_run_coefficients1_4']=[a1,a2,a3,a4]
        if a4==0:
            data['all_size_exponential_coefficients_base3_base2_base1']=[a3,a2-3*a3,a1-2*a2+3*a3]
            data['all_size_reason']='M1 positive-run support descends: A4=0 forces every later coefficient0'
        data['ready_at_sizes1_2_3_4']=[1]+[q-s for _,q,s in data['profiles_sizes2_3_4']]
    out = {'scope': 'Every strong unique-start parent throughorder9; exterior root and predecessor; H1=2mod3; eligible vertex modules',
           'compiler': COMPILER, 'parents_by_order': {n: len(v) for n, v in bank.items()},
           'module_tests_at_size2': len(cases), 'candidates': candidates}
    (HERE / 'pure_ternary_receipt.json').write_text(json.dumps(out, indent=2) + '\n')
    print(json.dumps({'tests': len(cases), 'first_step_hits': len(candidates),
                      'certified': [c for c in candidates if c['pure_ternary_ready_certificate']]}))


if __name__ == '__main__':
    main()
