"""Render the readable attribution index from the saved complete harvest."""
from pathlib import Path
import gzip,json
HERE=Path(__file__).resolve().parent
def run():
    all_rows=json.loads(gzip.decompress((HERE/'history_full_prompt_attributions.json.gz').read_bytes()))
    rows=[r for r in all_rows if r['five_in_excerpt']]
    out=['# Historical prompt-attribution index','',
      'Original relative links inside quoted excerpts belong to their historical source; use each entry’s blob ID and source path to retrieve it.','',
      f'**HISTORICAL RETRIEVAL ONLY.** These are machine-selected document excerpts, not authenticated original chats, new instructions, or validated mathematical claims. Current corrections and canon take precedence. The {len(rows)} entries contain a five-query term within the retained context; duplicates and author framing remain possible. The complete {len(all_rows):,}-context harvest, including {len(all_rows)-len(rows)} with five only elsewhere in their source, is in `history_full_prompt_attributions.json.gz`. See HISTORY.md for manually interpreted entries and the precise frozen scope.','']
    for i,r in enumerate(rows,1):
        p=r['provenance'][0]
        out.extend([f'## P{i:03d}','',f"Source: `{p['path']}`, line {p['line']}. Blob: `{p['oid']}`.",''])
        out.extend(('> '+line).rstrip() for line in r['text'].splitlines())
        out.append('')
    (HERE/'PROMPT-INDEX.md').write_text('\n'.join(out).rstrip()+'\n',encoding='utf8')
    print(f'Rendered {len(rows)} historical retrieval entries.')
if __name__=='__main__':run()
