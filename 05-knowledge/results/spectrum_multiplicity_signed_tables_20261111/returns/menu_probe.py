"""Complete two-newborn menus with an indegree-one marked newborn.

Every such extension has the marked newborn's sole predecessor either
in the old graph or at the other newborn. No general extension census.
"""
from pathlib import Path
from tempfile import TemporaryDirectory
import subprocess,json
from native_return import construct,birth,state,ends

HERE=Path(__file__).resolve().parent

def add_outcut(rows,cut):
    n=len(rows)
    return [x|(0 if cut>>i&1 else 1<<n) for i,x in enumerate(rows)]+[cut]

def main():
    cases=[]
    with TemporaryDirectory(prefix='native-return-menu-') as td:
        exe=Path(td)/'cut_menu'
        subprocess.run(['clang++','-O3','-std=c++17',str(HERE/'cut_menu.cpp'),'-o',str(exe)],check=True)
        for source,targets in ((construct(1)[1],[5993,5995]),(construct(1)[3],[17981,17983])):
            n=len(source)
            def run(rows,root,mode):
                data=' '.join(map(str,[mode,len(rows),root,*targets]+rows))+'\n'
                return json.loads(subprocess.check_output([str(exe)],input=data,text=True))
            menus=[];positive=[]
            r=run(source,n-1,1)
            menus.append({'predecessor':'other_newborn','result':r,'valid_hits':r['hits']})
            for hit in r['hits']:
                # mode1 cut is the first newborn's IN-neighbor set.
                I=hit['cut'];graph=add_outcut(source,((1<<n)-1)^I);graph=birth(graph,n)
                actual=state(graph,n+1)
                if actual[0]!=hit['H']:raise RuntimeError(('new predecessor recount',hit,actual))
                positive.append({'predecessor':'other_newborn','cut':I,'rows':graph,'root':n+1,'HESUVdK':actual})
            for pred in range(n):
                parent=birth(source,pred);r=run(parent,n,0)
                valid=[x for x in r['hits'] if x['retains_indegree_one']]
                menus.append({'predecessor':pred,'result':r,'valid_hits':valid})
                for hit in valid:
                    graph=add_outcut(parent,hit['cut']);actual=state(graph,n)
                    if actual[0]!=hit['H']:raise RuntimeError(('old predecessor recount',hit,actual))
                    if [i for i in range(len(graph)) if graph[i]>>n&1]!=[pred]:raise RuntimeError('indegree')
                    positive.append({'predecessor':pred,'cut':hit['cut'],'rows':graph,'root':n,'HESUVdK':actual})
            cases.append({'source_rows':source,'source_H':sum(ends(source)[-1]),'targets':targets,'legal_extensions':(n+1)*(1<<n),'menus':menus,'positive_literal_controls':positive})
    receipt={'scope':'For each named source, all two-newborn extensions whose specified marked newborn has indegree one. Its predecessor may be any old vertex or the other newborn. No constraint on the other newborn. Large intermediate H is pruned only by monotonicity.','cases':cases}
    (HERE/'menu_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps([{'source_H':c['source_H'],'targets':c['targets'],'legal_extensions':c['legal_extensions'],'hits':[{'H':x['HESUVdK'][0],'predecessor':x['predecessor'],'cut':x['cut']} for x in c['positive_literal_controls']]} for c in cases],indent=2))

if __name__=='__main__':main()
