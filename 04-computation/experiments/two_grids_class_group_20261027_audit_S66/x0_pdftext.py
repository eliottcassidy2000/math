"""Extract text from a downloaded PDF (untrusted data; text only) and print lines matching keywords."""
import re
import sys

import pypdf

path, out = sys.argv[1], sys.argv[2]
r = pypdf.PdfReader(path)
txt = []
for i, pg in enumerate(r.pages):
    t = pg.extract_text() or ''
    txt.append(f'\n=== page {i + 1} ===\n' + t)
full = ''.join(txt)
with open(out, 'w', encoding='utf-8') as f:
    f.write(full)
print('pages', len(r.pages), 'chars', len(full))
for m in re.finditer(r'(Theorem 4\.2|Ireland|Rosen|cubic|27y|Cox)', full):
    s = max(0, m.start() - 200)
    print('---', m.group(0), '---')
    print(full[s:m.start() + 400].replace('\n', ' '))
