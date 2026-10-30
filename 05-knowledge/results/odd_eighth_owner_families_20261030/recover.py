"""Deterministic two-wave archaeology; occurrence counts are not truth labels."""
from pathlib import Path
import subprocess
import re
import json
import os

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]


def main():
    paths=subprocess.check_output(['git','ls-files','*.md'],cwd=ROOT).decode('utf-8').splitlines()
    rows=[]; eligible=0; read=0; failures=[]
    # Scan text, not generated JSON, frozen source replicas, or full chronology.
    for name in paths:
        if not name.startswith(('01-canon/theorems/','05-knowledge/results/','07-reflections/')): continue
        if any(x in name for x in ['/inputs/','/sources/','/frozen','/producer/']): continue
        eligible+=1
        path=ROOT/name
        if os.name=='nt': path=Path('\\\\?\\'+str(path.resolve()))
        try: text=path.read_text(encoding='utf-8')
        except (OSError,UnicodeError):
            failures.append(name)
            continue
        read+=1
        eighths=[k for k in [1,3,5,7] if re.search(fr'{k}\s*/\s*8\b|\\(?:t?frac)\{{{k}\}}\{{8\}}',text)]
        words=sorted(set(re.findall(r'odd.eighth|eighth.roots?|Boolean quartet|boundary.owner|primitive eighth|units modulo 8',text,re.I)))
        if len(eighths)>=3 or words:
            snippets=[{'line':i+1,'text':line[:300]} for i,line in enumerate(text.splitlines())
                      if re.search(r'(?:1|3|5|7)\s*/\s*8\b|eighth|Boolean quartet|boundary.owner',line,re.I)][:8]
            rows.append({'path':name,'fraction_numerators':eighths,'word_matches':words,'snippets':snippets})
    data={'revision':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT).decode().strip(),
          'tracked_markdown_inventory':len(paths),'eligible_documents':eligible,
          'read_documents':read,'read_failures':failures,'matched_documents':len(rows),'matches':rows}
    (HERE/'recovery.json').write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({k:v for k,v in data.items() if k not in ['matches','read_failures']}))


if __name__=='__main__':main()
