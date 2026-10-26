"""Fast, deterministic comparison of completed exhaustive replay receipts."""
import argparse
import hashlib
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def semantics(receipt):
    require(receipt['status'] == 'PASS', 'receipt failed')
    keys = ('n', 'universe', 'unique_literal_certificates',
            'certificate_stream_lf_sha256', 'producer_arguments', 'status')
    orders = [{key: item[key] for key in keys} for item in receipt['orders']]
    require([item['n'] for item in orders] == [8, 9, 10], 'not the complete finite head')
    for item, total in zip(orders, (640, 10240, 143360)):
        require(item['status'] == 'PASS', 'order failed')
        require(item['unique_literal_certificates'] == item['universe']['sources'] == total,
                'wrong complete-universe count')
        require(item['producer_arguments'] == [str(item['n']), '4'], 'filtered invocation')
    for item in receipt['orders']:
        require(any(line.startswith('# decoder controls:')
                    for line in item['producer_stdout_without_certificates']),
                'decoder controls absent')
        require(not any(line.startswith('NO ')
                        for line in item['producer_stdout_without_certificates']),
                'failed source in output')
    return {'status': receipt['status'], 'compiler': receipt['compiler'],
            'build_options': receipt['build_options'],
            'corrected_source_lf_sha256': receipt['corrected_source_lf_sha256'],
            'independent_hostiles': receipt['independent_hostiles'], 'orders': orders}


def paired_summary(directory=HERE):
    normal = semantics(json.loads((directory / 'normal.json').read_text()))
    optimized = semantics(json.loads((directory / 'optimized.json').read_text()))
    require(normal == optimized, 'mathematical results or certificate streams differ')
    require(hashlib.sha256((directory / 'blue_line_black_detours_corrected.cpp').read_bytes()).hexdigest()
            == normal['corrected_source_lf_sha256'], 'corrected source changed since replay')
    return {'status': 'PASS', 'comparison': 'normal versus optimized complete census',
            'not_reexecuted': 'This fast comparison reads completed receipts; use check.py --census to regenerate them.',
            'excluded_observational_fields': [
                'compiler_stderr (includes absolute source paths)',
                'orders[].elapsed_seconds',
                'orders[].producer_stdout_without_certificates (timings; decoder control and no-failure lines still checked)',
                'orders[].producer_stderr (progress times)'],
            'canonical_mathematical_summary': normal}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    payload = json.dumps(paired_summary(), indent=2) + '\n'
    if args.output:
        args.output.write_bytes(payload.encode())
    else:
        print(payload, end='')


if __name__ == '__main__':
    main()
