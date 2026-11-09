"""audit D2 (THM-5056), item 3: the four seeds, every admissible w, the first PH(w, 2), and 513 = d(K_2 u T_171).
d by subset enumeration, parity by complete enumeration of Aut (own code, imported from d2build / d2base)."""
import importlib.util
def load(name):
    spec = importlib.util.spec_from_file_location(name, name + '.py'); m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m); return m
B = load('d2build'); Q = load('d2base')
def parity(G):
    auts = Q.all_auts(G); odd = sum(1 for g in auts if Q.rev_edges(G, g) % 2)
    return f'|Aut|={len(auts)} odd elements={odd} -> {"ODD" if odd else "EVEN"}'
names = {5: ['K_4 - e'], 1: ['K_2 u C_4'], 37: ['K_2 v (P_3 u E_2)', 'K_2 u T_55']}
for c in (5, 1, 37):
    for (G, tw, d0), nm in zip(B.seeds(c), names[c]):
        n = len(G); d = Q.dcount(G)
        ws = B.deg2(G, tw)
        print(f'{nm}: n={n} d={d} (claimed {d0}) offset={d - 2**(n-1)} twins{tw}={bool(B.adj_twins(G, *tw))} '
              f'admissible w={ws}; {parity(G)}')
        for w in ws:
            H = B.ph(G, w, 2); dH = Q.dcount(H)
            B.assert_odd_good(H, tw)
            print(f'   PH(w={w},2): n={len(H)} d={dH} = 2^{len(H)-1} + {dH - 2**(len(H)-1)}; odd-good; {parity(H)}')
T55 = B.threshold(55, 6); print('T_55 degrees:', [B.deg(T55, v) for v in range(6)], 'd =', Q.dcount(T55))
G513 = B.union(B.edges_of([(0, 1)], 2), B.threshold(171, 8))
print('K_2 u T_171: n =', len(G513), 'd =', Q.dcount(G513), 'L(513) =', (513).bit_length(), ';', parity(B.edges_of([(0,1)],2)), '(K_2 component)')
