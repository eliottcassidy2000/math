#!/usr/bin/env python3
from pathlib import Path
from hashlib import sha256,sha1
import json,re,itertools,subprocess,difflib
from validate_package import read_package,need,RESULT,D,PUB

def main():
    old,oe,om,original=read_package('CP03-candidate');new,ne,nm,current=read_package('CP03-credit-final')
    extra=set(current)-set(original);need(not(set(original)-set(current)),'no removed original')
    expected={RESULT+'/spectrum/'+n for n in ['INCOMING_CREDIT.md','check_incoming_credit.py','INCOMING_ROUTE_MATCH.json','INCOMING_SPECTRUM_SOURCES.json']}
    need(extra==expected,'four precise additions')
    changed=[p for p in original if original[p]!=current[p]];need(changed==['00-navigation/S190.md'],'only navigation changed')
    delta=''.join(difflib.unified_diff(original[changed[0]].decode().splitlines(True),current[changed[0]].decode().splitlines(True),fromfile='CP03-original',tofile='CP03-shared-credit'))
    (D/'CREDIT_NAVIGATION_DELTA.txt').write_text(delta)
    source=json.loads((D/'INDEPENDENT_INCOMING_READS.json').read_text());retained=json.loads(current[RESULT+'/spectrum/INCOMING_SPECTRUM_SOURCES.json'])
    need(source['commit']==retained['commit'],'same incoming commit')
    for f,g in zip(source['files'],retained['files']):
        need(f['path']==g['path']and f['sha']==g['sha']and f['content']==g['content'],'independent immutable incoming read matches')
        raw=f['content'].encode();need(sha1(b'blob '+str(len(raw)).encode()+b'\0'+raw).hexdigest()==f['sha'],'incoming content git hash')
    out=source['files'][1]['content']
    initial=list(map(int,re.search(r'A out-rows \(root 0\): (.*?)   scores:',out).group(1).split()))
    lines=[line for line in out.splitlines()if 'below-start route'in line]
    target=lines[1];flips=[tuple(map(int,x))for x in re.findall(r'\{(\d+),(\d+)\}',target)]
    need(flips==[(0,1),(1,2),(0,3),(3,4)],'incoming second route parsed')
    route=[initial]
    for a,b in flips:
        r=route[-1][:];r[a]^=1<<b;r[b]^=1<<a;route.append(r)
    ours=json.loads(original[RESULT+'/spectrum/spectrum/receipt.json'])['patch'];need(len(ours)==len(route)==5,'five complete vertices')
    maps=[]
    for p in itertools.permutations(range(7)):
        if p[2]!=0:continue
        good=True
        for a,b in zip(ours,route):
            mapped=[0]*7
            for i in range(7):
                for j in range(7):
                    if a['rows'][i]>>j&1:mapped[p[i]]|=1<<p[j]
            if mapped!=b:good=False;break
        if good:maps.append(p)
    need(maps==[(2,1,0,3,6,4,5)],'exact unique simultaneous rooted transport')
    root=D/'credit_replay';need(not root.exists(),'fresh credit replay')
    for rel,raw in current.items():
        path=root/rel;path.parent.mkdir(parents=True,exist_ok=True);path.write_bytes(raw)
    script=root/RESULT/'spectrum/check_incoming_credit.py'
    normal=subprocess.check_output(['python3',str(script)],cwd=root)
    optimized=subprocess.check_output(['python3','-O',str(script)],cwd=root)
    need(normal==optimized==current[RESULT+'/spectrum/INCOMING_ROUTE_MATCH.json'],'portable exact incoming checker output')
    for f in new['source_files']:need(Path(f['source']).read_bytes()==current[f['path']],'delta source binding')
    result={'status':'PASS','candidate':'CP03-credit-final','package_sha256':sha256((PUB/'CP03-credit-final/package.json').read_bytes()).hexdigest(),'originals':len(current),'objects':len(ne),'added_paths':sorted(extra),'changed_paths':changed,'frozen_mathematical_originals_unchanged':62,'independently_fetched_incoming_commit':source['commit'],'unique_common_relabeling':maps[0],'normal_optimized_credit_replay':True,'inline_tree_test_sha256':sha256((D/'IDEMPOTENT_TREE_TEST.json').read_bytes()).hexdigest()}
    (D/'CREDIT_DELTA_RECEIPT.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2));print(delta)
if __name__=='__main__':main()
