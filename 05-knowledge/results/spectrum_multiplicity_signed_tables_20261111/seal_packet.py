"""Seal this research packet and check its local Markdown destinations.

This checks artifact integrity and routes, not mathematical truth. Historical
dependencies may be outside the sparse checkout; their Git objects count as
present. No dependency is fetched or changed.
"""
from pathlib import Path
import hashlib
import json
import re
import subprocess
from urllib.parse import unquote

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUTPUTS = {'MANIFEST.json', 'LINK_CHECK.json'}


def main():
    paths = sorted(p for p in HERE.rglob('*') if p.is_file()
                   and p.name not in OUTPUTS and '__pycache__' not in p.parts)
    tracked = set(subprocess.check_output(
        ['git', 'ls-tree', '-r', '--name-only', 'HEAD'], cwd=ROOT,
        text=True).splitlines())
    links = []
    for p in paths:
        if p.suffix != '.md':
            continue
        prose = re.sub(r'```.*?```', '', p.read_text(), flags=re.S)
        prose = '\n'.join(line for line in prose.splitlines()
                          if not line.startswith(('    ', '\t')))
        prose = re.sub(r'`[^`]*`', '', prose)
        for destination in re.findall(r'\]\(([^\s)]+)(?:\s+[^)]*)?\)', prose):
            if destination.startswith(('http:', 'https:', '#', 'mailto:')):
                continue
            target = unquote(destination.split('#', 1)[0].strip('<>'))
            if not target:
                continue
            resolved = (p.parent / target).resolve()
            try:
                rel = str(resolved.relative_to(ROOT))
            except ValueError:
                rel = str(resolved)
            present = resolved.exists() or rel in tracked
            links.append({'source': str(p.relative_to(HERE)),
                          'destination': destination,
                          'repo_path': rel, 'present': present})
    failures = [x for x in links if not x['present']]
    report = {'scope': 'Local Markdown destinations; anchors not validated.',
              'checked': len(links), 'failures': failures, 'links': links}
    (HERE / 'LINK_CHECK.json').write_text(json.dumps(report, indent=2) + '\n')
    if failures:
        raise RuntimeError(json.dumps(failures, indent=2))
    manifest = {'scope': 'Owned packet bytes; proofs distinguish inherited, proved, and finite scopes.',
                'files': [{'path': str(p.relative_to(HERE)),
                           'bytes': len(p.read_bytes()),
                           'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}
                          for p in paths]}
    (HERE / 'MANIFEST.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(json.dumps({'status': 'PASS', 'files': len(paths), 'links': len(links)}))


if __name__ == '__main__':
    main()
