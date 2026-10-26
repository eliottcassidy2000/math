"""Read-only source/receipt validation; --refresh explicitly writes pins."""
import ast
import hashlib
import json
from pathlib import Path
import re
import sys
from replay import JOBS

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
INHERITED = [
    '01-canon/theorems/THM-305-two-adic-valuation-tournament-count.md',
    '01-canon/theorems/THM-2228-mahler-three-halves-carry-tail-and-integral-stabilization.md',
    '01-canon/theorems/THM-3227-selected-root-residue-contact-trie-primitive-carry-and-delayed-resplitting.md',
    '05-knowledge/results/endpoint_four_repairs_20261021/README.md',
    '05-knowledge/results/endpoint_four_repairs_20261021/audit/AUDIT.md',
    '05-knowledge/results/endpoint_four_repairs_20261021/audit/independent_audit.py',
    '05-knowledge/results/endpoint_four_repairs_20261021/audit/audit_core.py',
    '05-knowledge/results/ports_owners_f1_20261025/endpoint.md',
    '05-knowledge/results/unprotected_ports_cancellation_20261025/ports/README.md',
    '05-knowledge/results/unprotected_ports_cancellation_20261025/ports/ENDCAP_NEXT.md',
    '05-knowledge/results/unprotected_ports_cancellation_20261025/ports/INDEPENDENT_REVIEW.md',
    '05-knowledge/results/port_residual_integral_20261025/PORT_FRONTIER.md',
    '05-knowledge/results/port_residual_integral_20261025/LOW_ORDERS.md',
    '05-knowledge/results/port_residual_integral_20261025/LOW_ORDERS_REVIEW.md',
    '05-knowledge/results/port_residual_integral_20261025/h3/README.md',
    '05-knowledge/results/port_residual_integral_20261025/h3/INDEPENDENT_REVIEW.md',
    '05-knowledge/results/port_residual_integral_20261025/h4/README.md',
    '05-knowledge/results/port_residual_integral_20261025/h4/INDEPENDENT_REVIEW.md',
    '05-knowledge/results/port_residual_integral_20261025/counting/README.md',
    '05-knowledge/results/port_residual_integral_20261025/counting/check.py',
    '05-knowledge/results/port_residual_integral_20261025/counting/REVERSE_POTENTIAL.md',
    '05-knowledge/results/port_residual_integral_20261025/counting/REVERSE_REVIEW.md',
    '05-knowledge/results/port_residual_integral_20261025/precision/README.md',
    '05-knowledge/results/port_residual_integral_20261025/precision/check.py',
    '05-knowledge/results/session56_rmq_lca_connections_20261025/README.md',
    '05-knowledge/results/session56_rmq_lca_connections_20261025/valuation-lca/REPORT.md',
    '05-knowledge/results/session56_rmq_lca_connections_20261025/debt-range/RESULT.md',
    '05-knowledge/results/normseven_boundary_20261025/PAPERS_AND_TRANSPORT.md',
    '01-canon/theorems/THM-4732-black-detours-of-closed-blue-lines-order-9-census-order-10-detours-pure-blue-layers-and-the-klein-dictionary.md',
    '04-computation/experiments/blue_line_black_detours_20261020.cpp',
    '04-computation/experiments/blue_line_black_detours_20261020_n8.out',
    '04-computation/experiments/blue_line_black_detours_20261020_n9.out',
    '04-computation/experiments/blue_line_black_detours_20261020_n10.out',
    '05-knowledge/results/session64_lca_reverse_valuation_20261026/session64-repairs/next11/README.md',
    '05-knowledge/results/session64_lca_reverse_valuation_20261026/session64-audit-eleven/AUDIT.md',
    '05-knowledge/results/session64_lca_reverse_valuation_20261026/session64-counting/REVERSE_IMPLEMENTATION.md',
    '05-knowledge/results/session64_lca_reverse_valuation_20261026/session64-counting/SKIP_CONTRACTION.md',
    '05-knowledge/results/session64_lca_reverse_valuation_20261026/session64-counting/reverse_ordinary.py',
    '05-knowledge/results/session64_lca_reverse_valuation_20261026/session64-counting/contracted_ordinary.py',
    '05-knowledge/results/session64_lca_reverse_valuation_20261026/session64-audit-counting/REVERSE_AUDIT.md',
    '05-knowledge/results/session64_lca_reverse_valuation_20261026/session64-audit-counting/CONTRACTION_AUDIT.md',
    '05-knowledge/results/session64_lca_reverse_valuation_20261026/session64-reverse-valuation/REPORT.md',
    '05-knowledge/results/session64_lca_reverse_valuation_20261026/session64-audit-dynamics/REPORT.md',
]


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def digest(path):
    return hashlib.sha256(path.read_bytes().replace(b'\r\n', b'\n')).hexdigest()


def main():
    excluded = {'manifest.json', 'verification.json', 'verification-optimized.json'}
    files = sorted(p for p in HERE.rglob('*') if p.is_file()
                   and '__pycache__' not in p.parts and p.name not in excluded)
    manifest = HERE / 'manifest.json'
    if '--refresh' in sys.argv:
        data = {'hash_semantics':'SHA256 with CRLF normalized to LF',
                'packet':{p.relative_to(HERE).as_posix():digest(p) for p in files},
                'inherited':{name:digest(ROOT/name) for name in INHERITED}}
        manifest.write_bytes((json.dumps(data,indent=2,sort_keys=True)+'\n').encode('utf-8'))
    data = json.loads(manifest.read_text(encoding='utf-8'))
    need(set(data['packet']) == {p.relative_to(HERE).as_posix() for p in files}, 'packet file set changed')
    for name, expected in data['packet'].items():
        need(digest(HERE/name) == expected, ('packet hash', name))
    need(set(data['inherited']) == set(INHERITED), 'inherited file set changed')
    for name, expected in data['inherited'].items():
        need(digest(ROOT/name) == expected, ('inherited hash', name))
    for source, normal, optimized in JOBS:
        a = (HERE/normal).read_text(encoding='utf-8')
        b = (HERE/optimized).read_text(encoding='utf-8')
        need(a == b, ('receipt mismatch', source))
        need(json.loads(a).get('status') == 'PASS', ('receipt failed', source))
    links = parsed = 0
    for path in files:
        if path.suffix == '.py':
            ast.parse(path.read_text(encoding='utf-8'),filename=str(path))
            parsed += 1
        if path.suffix != '.md':
            continue
        prose = re.sub(r'```.*?```', '', path.read_text(encoding='utf-8'), flags=re.S)
        for url in re.findall(r'\]\(([^)]+)\)', prose):
            if '://' in url or url.startswith('#'):
                continue
            target = url.split('#',1)[0].strip('<>')
            need((path.parent/target).exists(), ('broken link', str(path), url))
            links += 1
    print(json.dumps({'status':'PASS','packet_files':len(files),
                      'inherited_files':len(INHERITED),'receipt_pairs':len(JOBS),
                      'python_files':parsed,'local_links':links},indent=2,sort_keys=True))


if __name__ == '__main__':
    main()
