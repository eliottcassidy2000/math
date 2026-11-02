"""Run the declared finite checks normally and optimized, requiring identical receipts."""
from pathlib import Path
import hashlib
import json
import re
import subprocess
import sys

HERE = Path(__file__).resolve().parent
CASES = [('verify_stability.py', 'stability_checks.json'),
         ('verify_stability_320.py', 'stability_320_checks.json'),
         ('verify_stability_238.py', 'stability_238_checks.json'),
         ('verify_stability_56.py', 'stability_56_checks.json'),
         ('verify_clocks.py', 'clock_checks.json'),
         ('verify_carriers.py', 'carrier_checks.json')]


def main():
    rows, total = [], 0
    for script, receipt in CASES:
        normal = subprocess.run([sys.executable, '-B', str(HERE/script)],
                                cwd=HERE, check=True, capture_output=True, text=True)
        before = (HERE/receipt).read_bytes()
        optimized = subprocess.run([sys.executable, '-B', '-O', str(HERE/script)],
                                   cwd=HERE, check=True, capture_output=True, text=True)
        after = (HERE/receipt).read_bytes()
        if before != after or normal.stdout != optimized.stdout:
            raise RuntimeError('normal/optimized replay differs: '+script)
        record = json.loads(after)
        total += record['gates']
        rows.append({'script': script, 'receipt': receipt, 'gates': record['gates'],
                     'normal_optimized_identical': True,
                     'receipt_sha256': hashlib.sha256(after).hexdigest()})
    if total != 64626:
        raise RuntimeError('documented gate count changed')
    links = 0
    for path in HERE.glob('*.md'):
        prose = '\n'.join(line for line in path.read_text(encoding='utf-8').splitlines()
                          if not line.startswith('    '))
        for target in re.findall(r'\]\(([^)]+)\)', prose):
            if '://' in target or target.startswith('#'):
                continue
            local = (path.parent/target.split('#', 1)[0]).resolve()
            if not local.exists():
                # Existing very long Windows paths require the extended spelling.
                if sys.platform != 'win32' or not Path('\\\\?\\'+str(local)).exists():
                    raise RuntimeError('missing local link: '+str(path)+' -> '+target)
            links += 1
    report = {'status': 'PASS', 'local_gates': total, 'inherited_hostile_gates': 7,
              'local_links_checked': links, 'runs': rows}
    (HERE/'REPLAY.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    print(json.dumps({'status': 'PASS', 'local_gates': total, 'local_links_checked': links}))


if __name__ == '__main__':
    main()
