"""Exact first-growth and boundary-parity controls; standard library."""
import hashlib,json
from pathlib import Path
from verify import T,bank,coords,need,section,signed,state,step
S=bank(T,3000)
for s in S:
    D,E,U,d=s;F=D-E+U
    need(F>D,'strict first growth')
    for r in range(7):
        v=state('O'*r,s)
        for c in 'ZJ':
            value=step(v,c)[0]
            need(value==(F if r==0 else D+2**r*U),'first changed D formula')
            need(value>=F,'growth gap')
for r in range(31):
    scale=3**r;s=state('Z'+'OJ'*r+'OOZ');v=state('OOOOJ'+'OJ'*r)
    need(s==tuple(scale*x for x in (15,1,23,1)),'scaled gap source')
    need(coords(s)[1]>0 and s[0]<signed(s)<s[0]-s[1]+s[2],'strict gap')
    need(v[0]==v[2]==signed(s),'balanced target independently constructed')
    q=state('OJJJJ',v);need(section(q) and signed(q)==signed(s),'replacement slice')
z=(-1,-7,-1,-7);formal={z:''};queue=[z];max_length=11
for s in queue:
    w=formal[s];D,E,U,d=s
    need(D%2==d%2==1 and (U-E)%2==0,'boundary parity')
    need(step(s,'O')!=step(s,'Z'),'no O/Z collision')
    if len(w)==max_length:continue
    for c in 'OZJ':
        q=step(s,c)
        if q not in formal:formal[q]=w+c;queue.append(q)
# Exact finite full-lattice controls of the three iff descriptions.
import itertools
lattice=0
for s in itertools.product(range(-3,4),repeat=4):
    D,E,U,d=s
    need((step(s,'O')==step(s,'J'))==(D==E==U==d),'O/J collision locus')
    need((step(s,'Z')==step(s,'J'))==(E==0 and d==D and U==2*D),'Z/J collision locus')
    need((step(s,'O')==step(s,'Z'))==(d==0 and E==U),'O/Z collision locus')
    lattice+=1
print(json.dumps(dict(status='PASS',physical_cap=3000,physical_profiles=len(S),growth_tail_lengths=[0,6],
    hostile_family_r=[0,30],formal_boundary_word_length=max_length,formal_boundary_states=len(formal),
    formal_lattice_controls=lattice,verifier_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    scope='Exact first-growth and formal parity statements; no coverage claim.'),indent=2,sort_keys=True))
