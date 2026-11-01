"""Build fixed, call-free arithmetic sources from the shared odd-cycle kernel.

This is an artifact-generation step, never a runtime oracle. Every emitted source
is fully contained in its own *_source.py file and has just one natural input n.
"""
from pathlib import Path
from textwrap import indent
ROOT = Path(__file__).resolve().parent


def kernel(base=2, color_exponent=None):
    original = (ROOT / 'a000568_source.py').read_text()
    body = original[original.index('    factorial = 1'):original.index('    return total // factorial')]
    body = '\n'.join(line[4:] for line in body.splitlines()) + '\n'
    if color_exponent is not None:
        body = body.replace('exponent = 0\n', 'exponent = 0\n        length = 0\n')
        body = body.replace('exponent += multiplicity * ((d - 1) // 2)',
                            'length += multiplicity\n            exponent += multiplicity * ((d - 1) // 2)')
        power = f'(2 ** ({2 if base == 4 else 1} * exponent + {color_exponent} * length))'
    else:
        power = f'({base} ** exponent)'
    body = body.replace('(2 ** exponent)', power)
    return body


def build():
    sources = {}
    sources['self_converse'] = ('def count(n):\n    parity = n % 2\n    n = n // 2\n' +
                               indent(kernel(4, 'parity'), '    ') + '    return total // factorial\n')
    # Use one kernel twice, with mode-dependent edge-orbit exponent multiplier.
    merged_kernel = kernel(2, 'parity').replace('1 * exponent + parity * length',
                                               '(1 + phase) * exponent + parity * length')
    sources['converse_merged'] = ('def count(n):\n    order = n\n    phase = 0\n    accumulated = 0\n'
        '    while phase < 2:\n        n = order\n        parity = 0\n'
        '        if phase == 1:\n            parity = n % 2\n            n = n // 2\n' +
        indent(merged_kernel, '        ') +
        '        accumulated += total // factorial\n        phase += 1\n    return accumulated // 2\n')
    sources['pair_tournaments'] = 'def count(n):\n' + indent(kernel(4), '    ') + '    return total // factorial\n'
    sources['colored_pair_tournaments'] = 'def count(n):\n' + indent(kernel(4, '1'), '    ') + '    return total // factorial\n'
    sources['colored_tournaments'] = 'def count(n):\n' + indent(kernel(2, '1'), '    ') + '    return total // factorial\n'
    sources['strong'] = ('def count(n):\n    if n == 0:\n        return 0\n'
        '    order = n\n    positive = 0\n    negative = 0\n    cut_code = 0\n'
        '    cut_limit = 2 ** (n - 1)\n    while cut_code < cut_limit:\n'
        '        cut_scan = cut_code\n        position = 1\n        component_size = 0\n'
        '        components = 0\n        component_product = 1\n'
        '        while position <= order:\n            component_size += 1\n            cut = 1\n'
        '            if position < order:\n                cut = cut_scan % 2\n                cut_scan = cut_scan // 2\n'
        '            if cut == 1:\n                n = component_size\n' + indent(kernel(), '                ') +
        '                component_product = component_product * (total // factorial)\n'
        '                component_size = 0\n                components += 1\n            position += 1\n'
        '        if components % 2 == 1:\n            positive += component_product\n'
        '        else:\n            negative += component_product\n        cut_code += 1\n'
        '    return positive - negative\n')
    sources['extremizer_tail'] = ('def count(n):\n    residue = n % 3\n    k = n // 3 + 2\n'
        '    if residue == 0:\n        return 1\n    if residue == 1:\n        return k\n'
        '    return (k * (k - 1)) // 2\n')
    for name, source in sources.items():
        compile(source, name, 'exec')
        (ROOT / f'{name}_source.py').write_text('# Fixed source: no calls, arrays, recursion, or sequence tables.\n' + source)
    return list(sources)


if __name__ == '__main__':
    print('\n'.join(build()))
