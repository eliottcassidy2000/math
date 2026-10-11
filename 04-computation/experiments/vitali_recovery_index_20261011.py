#!/usr/bin/env python3
"""Reproducible discovery index, not a proof-status audit of every match."""
import argparse
import json
import re
import subprocess
from collections import Counter
from pathlib import Path

BASE = '078ac6b72c7efe8d933c3ee9cbf9aeab7ab6f343'


def git(*args):
    return subprocess.check_output(['git', *args], text=True)


def lane(path, title):
    name = (path+' '+title).lower()
    if any(w in name for w in ['collatz', 'selector_certificate']):
        return 'Source measures, proof receipts, and tail equivalence'
    if any(w in name for w in ['lrc', 'lonel', 'miss-zeta', 'covering-depth',
                               'equioscillation', 'cyclotomic', 'saddle-index',
                               'antipodal', 'two-indices', 'hard-core', 'brouwer']):
        return 'LRC, moment atoms, compatibility, and boundary witnesses'
    if any(w in name for w in ['chromatic', 'hadwiger', 'lattice', 'seven-col',
                               'graph-color', 'graph-colour']):
        return 'Colouring, lattice geometry, and integrality'
    if any(w in name for w in ['lambda', 'vitali-two-channel', 'vitali-atom-character',
                               'c5-', 'tournament', 'walsh', 'hadamard', 'hamilton']):
        return 'Finite tournament observers and completion data'
    if any(w in name for w in ['chemistry', 'protein', 'attention', 'eigenvectors',
                               'towers', 'seventeen', 'three-prime', 'vitali-atom.md']):
        return 'Speculative towers, physics, and information analogies'
    return 'Cross-topic recovery and other companion ideas'


def make(ref):
    names = git('grep', '-il', 'vitali', ref, '--',
                '01-canon/theorems/*.md', '05-knowledge/results/*.md',
                '07-reflections/*.md').splitlines()
    paths = sorted(n.split(':', 1)[1] for n in names
                   if n.split(':', 1)[1] != '05-knowledge/results/INDEX.md')
    records = []
    for path in paths:
        source = git('show', ref+':'+path)
        lines = source.splitlines()
        titles = [z[2:].strip() for z in lines if z.startswith('# ')]
        headings = [re.sub(r'^#+\s+', '', z) for z in lines if re.match(r'^#{2,6} ', z)]
        provenance = [z for z in lines if re.match(
            r'^\*{0,2}(?:Session|Source|Author|Arising from|Instance|Date)\b', z, re.I)]
        status = [z for z in lines if re.match(r'^\*{0,2}status\b', z, re.I)]
        contexts = [{'line': i+1, 'text': ' '.join(lines[max(0, i-1):i+2])}
                    for i, z in enumerate(lines) if 'vitali' in z.lower()]
        title = titles[0] if titles else path.rsplit('/', 1)[1]
        records.append({'path': path, 'title': title, 'lane': lane(path, title),
                        'provenance_as_written': provenance,
                        'status_as_written_not_revalidated': status,
                        'all_companion_headings': headings, 'vitali_contexts': contexts,
                        'theorem_ids': sorted(set(re.findall(r'THM-\d+', source))),
                        'markdown_links': re.findall(r'\]\(([^)]+)\)', source)})
    return {'base_commit': ref, 'scope': 'Literal Vitali matches in the three declared Markdown trees; INDEX excluded',
            'status': 'DISCOVERY INDEX; historical labels are not proof endorsements',
            'document_count': len(records),
            'directory_counts': dict(sorted(Counter(z['path'].split('/')[0] for z in records).items())),
            'lane_counts': dict(sorted(Counter(z['lane'] for z in records).items())),
            'records': records}


def markdown(data):
    out = ['# Vitali session and companion-topic discovery index', '',
           '**DISCOVERY INDEX, not a proof audit.** Every literal match in the',
           'declared trees at the pinned commit is listed. A document without an',
           'explicit session label is retained by filename; documents are not',
           'silently counted as distinct sessions. The JSON preserves every heading,',
           'Vitali context, source label, and linked theorem ID. Reading these',
           'metadata is not a claim that every historical proof was audited.', '',
           f"Pinned commit: `{data['base_commit']}`. Documents: {data['document_count']}.",
           f"Directory counts: `{json.dumps(data['directory_counts'], sort_keys=True)}`.", '',
           '[Current synthesis and corrections](vitali_atom_recovery_20261011.md);',
           '[full machine-readable index](vitali_recovery_index_20261011.json).', '']
    for group in sorted(data['lane_counts']):
        out += ['## '+group, '']
        for row in data['records']:
            if row['lane'] != group:
                continue
            target = ('../../'+row['path']) if not row['path'].startswith('05-knowledge/results/') else row['path'].split('/')[-1]
            title = row['title'].replace('[', '(').replace(']', ')')
            out += [f"- [{title}]({target})", '']
            prov = row['provenance_as_written']
            out += ['  Provenance as written: '+(' / '.join(prov[:3]) if prov else 'no explicit session/source line detected; use filename.'), '']
            heads = row['all_companion_headings']
            out += ['  Companion topics: '+('; '.join(heads[:8]) if heads else
                    'no section headings; the full nearby Vitali passages are in the JSON.')+
                    ('; further headings in JSON.' if len(heads) > 8 else ''), '']
    return '\n'.join(out).rstrip()+'\n'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--ref', default=BASE)
    ap.add_argument('--write', action='store_true')
    args = ap.parse_args()
    ref = git('rev-parse', args.ref).strip()
    data = make(ref)
    if args.write:
        dest = Path('05-knowledge/results')
        (dest/'vitali_recovery_index_20261011.json').write_text(json.dumps(data, indent=2, ensure_ascii=False)+'\n')
        (dest/'vitali_recovery_index_20261011.md').write_text(markdown(data))
    print(json.dumps({k: v for k, v in data.items() if k != 'records'}, indent=2))


if __name__ == '__main__':
    main()
