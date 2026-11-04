from pathlib import Path
import itertools,json,subprocess,tempfile
from algebra import need,operator,response,strong,tables,birth,induced

HERE=Path(__file__).resolve().parent

def check(row,w,a):
    result=operator(row,w,a)
    need(tuple(result['predicted2'])==response(result['C2'],w),'two-birth +2 operator')
    need(tuple(result['predicted4'])==response(result['C4'],w),'two-birth +4 operator')
    need((result['H_X'],result['F_X'])==response(result['X'],w),'intermediate response')
    need(result['S_w_X']==tables(result['X'],w)[0][-1],'intermediate start response')
    need(result['S_a_2']>=result['h'] and result['L']+result['K']>=result['F'],'positive entry defects')
    return result

def main():
    checked=0
    for n in range(1,5):
        pairs=list(itertools.combinations(range(n),2))
        for code in range(1<<len(pairs)):
            rows=[0]*n
            for k,(a,b) in enumerate(pairs):
                if code>>k&1:rows[a]|=1<<b
                else:rows[b]|=1<<a
            # Identities apply to every tournament, not just strong ones.
            for w in range(n):
                for a in range(n):
                    result=check(rows,w,a);checked+=1
                    transported=operator(result['C2'],w,n)
                    h2,f2=result['predicted2']
                    need(transported['S_a_2']==h2,'transported auxiliary reverse readiness')
                    need(transported['L']+transported['K']==f2,'transported joint guard')
                    need(tuple(transported['predicted2'])==(3*h2,3*f2),'closed strong tripling coordinates')
                    same=operator(result['C2'],w,a)
                    last=operator(result['C2'],w,n+1)
                    deletion=induced(rows,[v for v in range(n) if v!=a])
                    dh=tables(deletion)[0][-1]
                    R=0 if w==a else response(deletion,w-int(w>a))[1]
                    G=result['L']+result['K']
                    need(same['S_a_2']==4*result['S_a_2']-dh,'old auxiliary start transport')
                    need(same['L']+same['K']==4*G-R,'old auxiliary joint transport')
                    need(last['S_a_2']==h2+result['h']-result['e_a'],'new y start transport')
                    need(last['L']+last['K']==f2+result['F']-result['J'],'new y joint transport')
    base=(8,1,3,6);first=check(base,0,1)
    prototypes=[('base',base),('child2',tuple(first['C2'])),('child4',tuple(first['C4']))]
    detail={}
    for name,rows in prototypes:
        need(strong(rows),'actual strong prototype')
        detail[name]=[check(rows,0,a) for a in range(len(rows))]
    menus={}
    with tempfile.TemporaryDirectory(prefix='two-birth-') as temp:
        exe=Path(temp)/'search.exe'
        audit=Path(temp)/'forward.exe'
        subprocess.run(['g++','-O3','-std=c++20',str(HERE/'birth_search.cpp'),'-o',str(exe)],check=True)
        subprocess.run(['g++','-O3','-std=c++20','-DFORWARD_AUDIT',str(HERE/'birth_search.cpp'),'-o',str(audit)],check=True)
        for name,rows in prototypes:
          for moving in (0,1):
            data=' '.join(map(str,[len(rows),0,moving,*rows]))+'\n'
            raw=subprocess.check_output([str(exe)],input=data,text=True)
            independent=subprocess.check_output([str(audit)],input=data,text=True)
            need(json.loads(raw)==json.loads(independent),'forward/deletion full-menu agreement')
            menu=json.loads(raw);menu.update(source_rows=rows,source_w=0)
            for hit in menu['hits']:
                x=birth(rows,hit['cut_x']);target=birth(x,hit['cut_y'])
                need(target==tuple(hit['rows']) and strong(target),'literal strong target')
                need(strong(x)==hit['intermediate_strong'],'intermediate strong receipt')
                need(response(target,hit['mark'])==(3*menu['source_H']+hit['offset'],3*menu['source_F']),'independent hit response')
            menus[name+('_moving' if moving else '')]=menu
    output=dict(status='FINITE_EXACT_WITH_GENERAL_IDENTITIES',formula_cases=checked,
                full_menu_independent_recurrences=['deletion','forward extension'],
                prototype_details=detail,two_birth_menus=menus)
    (HERE/'probe.json').write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps({'formula_cases':checked,'menus':{name:{'tested':m['tested'],'hits':len(m['hits']),
        'offsets':{str(d):sum(h['offset']==d for h in m['hits']) for d in (2,4)}} for name,m in menus.items()}},sort_keys=True))

if __name__=='__main__':main()
