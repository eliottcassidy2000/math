"""seeds.py - enumerate physical source states from the rooted edge t = (D,E,U,d) = (1,0,2,1) under the physical moves
O(D,E,U,d) = (D,0,2U,d), Z = (D-E+U,d,2U-d,d), J = (D-E+U,D-E,D+U,D) (complete source grammar, Session 198
bilinear-descent Theorem 1), with U <= UMAX; write the seed lists used by the confirmatory run. Seeds are written as
'h R e d' = 'D U E d'. A state is J-born if some word reaching it ends with J (its O/Z tail image is then not inside the
image of a shorter seed of the same word)."""
import sys, json
from pathlib import Path
HERE = Path(__file__).resolve().parent

def states(umax):
    seen = {(1, 0, 2, 1): ''}; q = [(1, 0, 2, 1)]
    for s in q:
        D, E, U, d = s
        for t, ch in (((D, 0, 2 * U, d), 'O'), ((D - E + U, d, 2 * U - d, d), 'Z'), ((D - E + U, D - E, D + U, D), 'J')):
            if t[2] <= umax and t not in seen:
                seen[t] = seen[s] + ch; q.append(t)
    return seen

def jborn(umax):
    """all states reachable by a word ending in J (any word, not only the BFS-first one)"""
    out = set(); seen = set([(1, 0, 2, 1)]); q = [(1, 0, 2, 1)]
    for s in q:
        D, E, U, d = s
        for t, ch in (((D, 0, 2 * U, d), 'O'), ((D - E + U, d, 2 * U - d, d), 'Z'), ((D - E + U, D - E, D + U, D), 'J')):
            if t[2] <= umax:
                if ch == 'J': out.add(t)
                if t not in seen: seen.add(t); q.append(t)
    return out

if __name__ == '__main__':
    jb = jborn(160)
    p1 = sorted(s for s in jb if 64 < s[2] <= 160 and s[3] > 1 and s[2] != 2 * s[3] and s[2] != s[0])
    bankA = [(1, 0, 2, 1)] + sorted(s for s in jb if s[2] <= 64 and s[3] > 1)
    # write LF bytes on every platform (the packet's receipts hash LF bytes)
    (HERE / 'p1_seeds.txt').write_bytes(''.join(f'{D} {U} {E} {d}\n' for D, E, U, d in p1).encode())
    (HERE / 'bankA.txt').write_bytes(''.join(f'{D} {U} {E} {d}\n' for D, E, U, d in bankA).encode())
    (HERE / 'root.txt').write_bytes(b'1 2 0 1\n')
    allst = states(4096)
    sym = sorted(s for s in allst if s[2] == 2 * s[3] and s[1] == 0)
    summary = {'p1_seeds': len(p1), 'bankA_seeds': len(bankA), 'states_U_le_4096': len(allst),
               'states_with_U_eq_2d_and_E_eq_0': [list(s) for s in sym]}
    (HERE / 'seeds_summary.json').write_bytes(json.dumps(summary, indent=1).encode())
    print(len(p1), 'P1 seeds;', len(bankA), 'bank A seeds;', len(allst), 'states with U <= 4096; symmetric:', sym)
