"""Draw all G25 edges; vertex order is one literal Hamiltonian path."""
from pathlib import Path
import json
HERE=Path(__file__).resolve().parent
def main():
    paths=json.loads((HERE/'paths.json').read_text())['orders'][24]['paths'];p=paths[0]
    row=json.loads((HERE/'verification.json').read_text())['orders'][24]
    pos={v:52+54*i for i,v in enumerate(p)};y=300
    role={'forced':('#182b49',4.5),'sometimes':('#176dab',3),'never':('#bd2b40',3)}
    svg=['<svg xmlns="http://www.w3.org/2000/svg" width="1400" height="640" viewBox="0 0 1400 640">',
         '<rect width="1400" height="640" fill="#f8fafc"/>',
         '<g font-family="Segoe UI,Arial,sans-serif" fill="#182b49">',
         '<text x="36" y="45" font-size="26" font-weight="600">Square-sum graph G₂₅: every edge, with its path-use count</text>',
         '<text x="36" y="75" font-size="17">Universe: all 10 undirected Hamiltonian paths, counted up to reversal. Labels are actual integers.</text>']
    path_edges={tuple(sorted(e)) for e in zip(p,p[1:])}
    chords=[r for r in row['edges'] if tuple(r['edge']) not in path_edges]
    for i,r in enumerate(chords):
        a,b=r['edge'];x1,x2=sorted((pos[a],pos[b]));color,width=role[r['status']]
        above=i%2==0;height=160 if (a,b)==(4,12) else min(195,70+0.15*(x2-x1));control=y+(-2*height if above else 2*height)
        dash=' stroke-dasharray="10 7"' if r['status']=='never' else ''
        svg.append(f'<path d="M{x1},{y} Q{(x1+x2)/2},{control} {x2},{y}" fill="none" stroke="{color}" stroke-width="{width}"{dash}/>')
        ty=y+(-height-8 if above else height+22)
        svg.append(f'<text x="{(x1+x2)/2}" y="{ty}" text-anchor="middle" fill="{color}" font-size="15">{a}–{b}: {r["uses"]}/10</text>')
    usage={tuple(r['edge']):r for r in row['edges']}
    for a,b in zip(p,p[1:]):
        r=usage[tuple(sorted((a,b)))];color,width=role[r['status']]
        svg.append(f'<line x1="{pos[a]}" y1="{y}" x2="{pos[b]}" y2="{y}" stroke="{color}" stroke-width="{width}"/>')
        if r['uses']!=10:svg.append(f'<text x="{(pos[a]+pos[b])/2}" y="{y+40}" text-anchor="middle" fill="{color}" font-size="13">{r["uses"]}/10</text>')
    for v,x in pos.items():
        svg.append(f'<circle cx="{x}" cy="{y}" r="19" fill="white" stroke="#657389" stroke-width="1.5"/><text x="{x}" y="{y+6}" text-anchor="middle" font-size="17" font-weight="600">{v}</text>')
    for x,status,label in [(40,'forced','Forced: 10/10 paths'),(390,'sometimes','Optional: 1–9/10 paths'),(780,'never','Never at size 25: 0/10 paths')]:
        color,width=role[status];dash=' stroke-dasharray="10 7"' if status=='never' else ''
        svg.append(f'<line x1="{x}" y1="565" x2="{x+60}" y2="565" stroke="{color}" stroke-width="{width}"{dash}/><text x="{x+72}" y="571" font-size="17">{label}</text>')
    svg.extend(['<text x="36" y="615" font-size="15">These colours encode fixed-size usage, not the proposed birth groups. Edge 1–15 is used at size 15, absent at size 25, and used again at size 27.</text>','</g></svg>'])
    (HERE/'G25_edge_roles.svg').write_text('\n'.join(svg)+'\n',encoding='utf-8')
if __name__=='__main__':main()
