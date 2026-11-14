from pathlib import Path
from branch.shell import compositions,flip_full,dp,triangles
from xml.etree import ElementTree
n=8
root=(0,1);gate=(6,57)
nodes={}
for w in compositions(n):
    if w!=w[::-1]:continue
    G=flip_full(w)
    point=(triangles(G),dp(G))
    if point!=root:nodes.setdefault(point,[]).append(w)
assert set(nodes)=={(6,65),(8,105),(10,117),(10,173),(10,225),(10,241)}
px=lambda c:94+54*c
py=lambda h:436-1.25*h
s=['<svg xmlns="http://www.w3.org/2000/svg" width="1070" height="500" viewBox="0 0 1070 500" role="img" aria-labelledby="t d">',
'<title id="t">The first two blue moves at order eight projected onto triangle count and Hamiltonian count</title>',
'<desc id="d">The transitive root goes to 57 paths and six triangles. Several blue second branches are shown. Two distinct tournaments with five-part or four-part layouts share the plotted coordinate with 173 paths and ten triangles.</desc>',
'<defs><style>text{font-family:system-ui,sans-serif;fill:#233657}.minor{fill:#62708c;font-size:12px}.label{font-size:14px}.b{font-weight:700}</style></defs>',
'<rect width="1070" height="500" fill="#f8faff"/><text class="b" x="35" y="55" font-size="26">A branching projection, not a planar embedding</text>',
'<text class="minor" x="36" y="83">Order 8 · blue moves only · ordinates are ordinary Hamiltonian path counts</text>']
for h in (0,50,100,150,200,250):
    y=py(h)
    s.append(f'<path d="M94 {y}H664" stroke="#dfe6f2"/><text class="minor" x="82" y="{y+4}" text-anchor="end">{h}</text>')
for c in range(11):
    x=px(c)
    s.append(f'<text class="minor" x="{x}" y="456" text-anchor="middle">{c}</text>')
s.append('<text class="minor" x="410" y="480" text-anchor="middle">cyclic triangles</text><path d="M418 112V436" stroke="#bac7de" stroke-dasharray="6 5"/>')
def line(a,b):
    return f'<path d="M{px(a[0])} {py(a[1])}L{px(b[0])} {py(b[1])}" stroke="#457bd1" stroke-width="2" opacity=".64" fill="none"/>'
s.append(line(root,gate))
for p in sorted(nodes):s.append(line(gate,p))
labels={root:(10,-12,'Transitive · 1'),gate:(11,20,'First blue mate · 57'),(6,65):(-11,-13,'65'),(8,105):(-12,-9,'105'),(10,117):(13,6,'117'),(10,173):(13,7,'173'),(10,225):(13,5,'225'),(10,241):(13,0,'241')}
for (c,h),(dx,dy,lab) in labels.items():
    x,y=px(c),py(h)
    s.append(f'<circle cx="{x}" cy="{y}" r="{8 if h==173 else 5.5}" fill="{"#e17f4b" if h==173 else "#215bbb"}" stroke="white" stroke-width="2"/>')
    s.append(f'<text class="label" x="{x+dx}" y="{y+dy}">{lab}</text>')
s+=['<rect x="725" y="128" width="309" height="168" rx="12" fill="#fff3eb" stroke="#f0cab2"/>',
'<text class="b" x="745" y="157" font-size="16">Two branches, one plotted point</text>',
'<text class="label" x="745" y="190">(1, 2, 2, 2, 1) and (2, 2, 2, 2)</text>',
'<text class="label" x="745" y="221">Both: 10 triangles, 173 paths,</text>',
'<text class="label" x="745" y="244">and the same score multiset.</text>',
'<text class="minor" x="745" y="272">Their triangle-incidence lists differ.</text>',
'<path d="M684 222L725 222" stroke="#e17f4b" stroke-width="2"/>',
'<text class="minor" x="725" y="350">Dashed line: first strong triangle level.</text>',
'<text class="minor" x="725" y="375">It contains a separate Boolean cube of</text>',
'<text class="minor" x="725" y="395">triangle moves, not all actual blue lines.</text></svg>']
p=Path(__file__).with_name('root_projection.svg')
p.write_text('\n'.join(s)+'\n');ElementTree.parse(p)
