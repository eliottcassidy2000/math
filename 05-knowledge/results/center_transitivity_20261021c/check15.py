"""Small separate VT15 corollary: discovery is explicit, default is a consumer."""
import argparse
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
GATES = 0

def need(ok, why):
    global GATES
    GATES += 1
    if not ok:
        raise ValueError(why)

def source(mask):
    if type(mask) is not int or not 0 <= mask < 128:
        raise ValueError('literal C15 mask')
    s = {15-d if mask >> (d-1) & 1 else d for d in range(1,8)}
    return tuple(sum(1 << y for y in range(15) if (y-x)%15 in s) for x in range(15))

def read_word(text):
    if type(text) is not str or len(text) != 15 or set(text) != set('ABCDEFGHIJKLMNO'):
        raise ValueError('literal source permutation')
    return tuple(ord(c)-65 for c in text)

def key(rows,p):
    return tuple((rows[p[i]] >> p[j]) & 1 for i in range(15) for j in range(i+1,15))

def is_blue(rows,p):
    return all(((rows[p[i]] >> p[j]) & 1) != ((rows[p[14-i]] >> p[14-j]) & 1)
               for i in range(15) for j in range(i+1,15))

def audit(mask, supplied):
    rows=source(mask)
    if type(supplied) is not tuple or len(supplied) != 7:
        raise ValueError('exact seven-word source receipt')
    keys=set()
    for text in supplied:
        p=read_word(text)
        need(all(rows[a] >> b & 1 for a,b in zip(p,p[1:])), 'actual path')
        need(not is_blue(rows,p), 'literal blackness')
        keys.add(key(rows,p))
    need(len(keys)==7,'seven distinct complete adjacency keys')

def discover():
    """Bounded optional producer; a failed cap never proves exclusion."""
    frames=[]
    for mask in range(128):
        rows=source(mask)
        found={}; nodes=0
        def visit(p,used):
            nonlocal nodes
            nodes+=1
            if nodes>1000000:
                raise ValueError('discovery cap exhausted')
            if len(p)==15:
                if not is_blue(rows,p):
                    found.setdefault(key(rows,p),''.join(chr(65+x) for x in p))
                return
            for j in range(15):
                if not used>>j&1 and rows[p[-1]]>>j&1:
                    visit(p+(j,),used|1<<j)
                    if len(found)==7:
                        return
        visit((0,),1)
        if len(found)!=7:
            raise ValueError('missing certificate')
        frames.append(list(found.values()))
    (HERE/'inputs15.json').write_text(json.dumps({'schema':'C15-seven-black-v1','frames':frames},sort_keys=True,indent=2)+'\n',encoding='utf-8',newline='\n')

def run():
    global GATES
    GATES=0
    text=(HERE/'inputs15.json').read_text(encoding='utf-8')
    data=json.loads(text)
    need(type(data) is dict and set(data)=={'schema','frames'} and data['schema']=='C15-seven-black-v1','schema')
    need(type(data['frames']) is list and len(data['frames'])==128,'complete128-frame cover')
    for mask,words in enumerate(data['frames']):
        need(type(words) is list,'literal word-list input')
        audit(mask,tuple(words))
    first=tuple(data['frames'][0])
    for bad in (False,0.0,-1,128):
        try:audit(bad,first)
        except ValueError:need(True,'invalid source rejected')
        else:need(False,'invalid source accepted')
    for wrong in (first[:6],(first[0],)*7,(first[0][:-1]+first[0][0],)+first[1:]):
        try:audit(0,wrong)
        except ValueError:need(True,'invalid literal receipt rejected')
        else:need(False,'invalid literal receipt accepted')
    need(all((p-1)%q for p,q in ((3,5),(5,3))), 'pq cyclic-reduction guard')
    return {'status':'FINITE-EXACT all128 cyclicC15 frames; only vertex-transitive order15 excluded',
            'frames':128,'literal_words':896,'black_orbits_at_least':7,
            'input_lf_sha256':hashlib.sha256(text.encode()).hexdigest(),
            'checks':GATES,'remaining_order15_scope':'Any deficit-six center must be non-vertex-transitive; existence remains OPEN.'}

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--discover',action='store_true')
    p.add_argument('--verify',action='store_true')
    args=p.parse_args()
    if args.discover:
        discover()
    out=json.dumps(run(),sort_keys=True,indent=2)+'\n'
    target=HERE/'receipt15.json'
    if args.verify:
        need(target.read_text(encoding='utf-8')==out,'saved receipt differs')
    else:target.write_text(out,encoding='utf-8',newline='\n')
    print(f'PASS {GATES} checks; LF SHA256 {hashlib.sha256(out.encode()).hexdigest()}')

if __name__=='__main__':main()
