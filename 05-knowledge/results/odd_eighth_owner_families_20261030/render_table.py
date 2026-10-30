"""Render the complete q=8 classification from its exact certificate."""
from pathlib import Path
import json

HERE=Path(__file__).resolve().parent


def main():
    data=json.loads((HERE/'certificate.json').read_text(encoding='utf-8'))['family']
    lines=['# Four odd-eighth centers: fourteen complete infinite functions','',
           '**PROVED uniform tail / FINITE-EXACT complete finite head.**','',
           'Let V(E,m)={1,...,7} union E union {8m}, with E any subset of',
           '{9,...,15} and m>=1. Every row gives the complete function M(V(E,m)).',
           'For EVERY m>=5 its value is m/(8m+s); the first four values follow.',
           'The all-size proof is in [PROOFS, O1-O3](PROOFS.md).','',
           'Pattern positions are speeds 9,10,11,12,13,14,15 in that order:',
           '`1` means present, `0` absent, and `*` arbitrary. Each row is exactly',
           'one Boolean cube; the disjoint rows cover all 128 choices.','',
           '| Pattern | Sources | m=1 | m=2 | m=3 | m=4 | s for m>=5 |',
           '|---|---:|---:|---:|---:|---:|---:|']
    for row in data['classes']:
        lines.append('| `'+row['pattern_9_through_15']+'` | '+str(row['count'])+' | '+
                     ' | '.join(row['head'])+' | '+str(row['s'])+' |')
    lines+=['','No equality is inferred merely from a finite prefix: matching s proves',
            'the entire remaining tail agrees. Distinct rows already differ for m<=5.','',
            '| Consumer and allowed continuation | Exact number of states |',
            '|---|---:|',
            '| Present answer at any fixed m>=5 | 5 |',
            '| Entire present function, all m>=1 | 14 |',
            '| Future band additions, queries m>=5 | 16 |',
            '| Future band additions, queries all m>=1 | 128 |','',
            'The last two are minimal deterministic carriers, with constructive',
            'separating continuations for every retained bit ([O4](PROOFS.md)).',
            'These state counts concern this declared family and operation menu.',
            'In particular the number 14 here does not settle LRC(14).','',
            'Reproduce with `python -X utf8 -B '+
            '05-knowledge/results/odd_eighth_owner_families_20261030/verify.py`,',
            'then run `render_table.py` in the same directory. The full exact times,',
            'finite overlap and native hostile controls are in [certificate.json](certificate.json).','']
    (HERE/'CLASSIFICATION.md').write_text('\n'.join(lines),encoding='utf-8',newline='\n')
    print('Rendered fourteen disjoint classification rows.')


if __name__=='__main__': main()
