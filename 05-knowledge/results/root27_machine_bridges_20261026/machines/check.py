"""Native clocks and configuration-preserving tournament controls; no network."""
import itertools
import json


def need(ok, label):
    if not ok:
        raise RuntimeError(label)


def parse(code):
    rows = code.split('_')
    need(all(len(row) == 6 for row in rows), 'binary transition rows')
    table = {}
    for q, row in enumerate(rows):
        for bit in (0, 1):
            token = row[3 * bit:3 * bit + 3]
            if token == '---':
                table[q, bit] = None
            else:
                write, direction, target = int(token[0]), token[1], ord(token[2]) - 65
                need(write in (0, 1) and direction in 'LR' and 0 <= target < len(rows), 'valid rule')
                table[q, bit] = write, 1 if direction == 'R' else -1, target
    return table


BOUNCE = '1RB0LB_1LA0RA'
CUBIC = '1RB0LC_1LA1RB_0RA1LC'
EXP = '1RB0LC_1RC1RA_1LA0RB'
SANDWICH = '1RB1LC_1LA1RD_1RA0LC_1LB0RD'
IRREGULAR = '1RB0RC_1LC1RA_1RA1LD_0LC0LA'
# Transcribed from Leng's TuringMachine25.lean, SHA256 in README.
GOLD = ('1RW1RB_0LQ1RC_0LU1RD_0LU0LE_1LQ1LF_0RH0LF_0LI0LH_1RJ0LH_0RJ0LI_'
        '0RK0RV_1RK0RL_1RM1RL_0LO0LN_1LG1LN_0LP1LO_1RE1LS_1LR1LQ_0RJ0LY_'
        '1RT1LS_1RC1RT_0LX1LU_1LR1RV_1LX1RW_0RA1LX_---1RX')


def step(table, cfg):
    q, head, ones = cfg
    rule = table[q, int(head in ones)]
    if rule is None:
        return None
    write, direction, target = rule
    new = set(ones)
    if write:
        new.add(head)
    else:
        new.discard(head)
    return target, head + direction, frozenset(new)


def run_transcript(table, cfg, transcript):
    for token in transcript:
        need(cfg is not None, 'no internal halt')
        q, head, ones = cfg
        need(token == chr(65 + q) + str(int(head in ones)), 'literal read transcript')
        cfg = step(table, cfg)
    return cfg


def exp_word(n, published=False):
    word = ['A0', 'B0']
    for k in range(n):
        word += ['C1', 'B0'] * (2 * k) + ['C0', 'A1'] * (2 * k + 1)
    word += (['C0', 'B1'] if published else ['C1', 'B0']) * (2 * n)
    word += ['C0', 'A1'] * (2 * n) + ['C0'] + ['A0', 'B1'] * (2 * n + 1)
    return word


def colored_tournament(cfg):
    """Implicit transitive orientation; ordered cell colors (bit, head-state, origin)."""
    q, head, ones = cfg
    active = set(ones) | {head, 0}
    left, right = min(active) - 1, max(active) + 1
    return tuple((int(x in ones), q if x == head else None, x == 0)
                 for x in range(left, right + 1))


def decode(cells):
    origin = [i for i, cell in enumerate(cells) if cell[2]]
    heads = [(i, cell[1]) for i, cell in enumerate(cells) if cell[1] is not None]
    need(len(origin) == len(heads) == 1, 'marked origin and head')
    i, q = heads[0]
    ones = frozenset(j - origin[0] for j, cell in enumerate(cells) if cell[0])
    return q, i - origin[0], ones


def graph_step(table, cells):
    """Direct ordered-color rewrite, independent of native set-of-ones step."""
    cells = [list(cell) for cell in cells]
    i = next(i for i, cell in enumerate(cells) if cell[1] is not None)
    rule = table[cells[i][1], cells[i][0]]
    if rule is None:
        return None
    write, direction, target = rule
    cells[i][0], cells[i][1] = write, None
    j = i + direction
    need(0 <= j < len(cells), 'source blank guards suffice for one step')
    cells[j][1] = target
    active = [i for i, cell in enumerate(cells) if cell[0] or cell[1] is not None or cell[2]]
    lo, hi = min(active), max(active)
    result = [[0, None, False]] + cells[lo:hi + 1] + [[0, None, False]]
    return tuple(tuple(cell) for cell in result)


def main():
    clock_cases = 0
    for n in range(65):
        cfg = (0, 0, frozenset(range(-4 * n, 0, 2)))
        word = ['A0'] + ['B0', 'A1'] * (2 * n + 1) + ['B0'] + ['A0', 'B1'] * (2 * n + 2)
        end = run_transcript(parse(BOUNCE), cfg, word)
        need(len(word) == 8 * n + 8, 'bouncer native duration')
        need(end == (0, 2, frozenset(range(-4 * n - 2, 2, 2))), 'bouncer target')
        clock_cases += 1
        cfg = (0, 0, frozenset(range(1, n + 1)))
        word = []
        for k in range(n, 0, -1):
            word += ['A0'] + ['B1'] * k + ['B0', 'A1'] + ['C1'] * k + ['C0', 'A1', 'C0']
        word += ['A0', 'B0', 'A1', 'C0']
        end = run_transcript(parse(CUBIC), cfg, word)
        need(len(word) == n * n + 7 * n + 4, 'cubic native duration')
        need(end == (0, 0, frozenset(range(1, n + 2))), 'cubic target')
        clock_cases += 1
    for n in range(1, 65):
        word = exp_word(n)
        end = run_transcript(parse(EXP), (0, 0, frozenset(range(-2 * n + 2, 0))), word)
        need(len(word) == 4 * n * n + 10 * n + 5, 'exponential duration')
        need(end == (0, 2 * n + 3, frozenset(range(-2 * n + 1, 2 * n + 3))), 'exponential target')
        clock_cases += 1
    try:
        run_transcript(parse(EXP), (0, 0, frozenset()), exp_word(1, published=True))
    except RuntimeError as error:
        need(str(error) == 'literal read transcript', 'published transcript hostile')
    else:
        raise RuntimeError('wrong transcript must fail')
    # Two published claims fail separately: changing the pair does not repair +3.
    end = run_transcript(parse(EXP), (0, 0, frozenset()), exp_word(1))
    need(len(end[2]) == 6 and len(end[2]) != 2 * (2 * 1 + 3) - 2, 'wrong exponential recurrence')
    need(all(rule is not None for code in (SANDWICH, IRREGULAR, BOUNCE, CUBIC, EXP)
             for rule in parse(code).values()), 'syntactic nonhalting of all five examples')

    # Identical read-labelled target graphs and directions, opposite halting.
    halt = parse('1RB---_0LA0LA')
    loop = parse('0RB---_0LA0LA')
    need({k: None if v is None else v[1:] for k, v in halt.items()} ==
         {k: None if v is None else v[1:] for k, v in loop.items()}, 'same target/direction graph')
    blank = (0, 0, frozenset())
    a = blank
    for unused in range(3):
        a = step(halt, a)
    need(a is None and step(loop, step(loop, blank)) == blank, 'halting versus literal cycle')
    graph_cases = 0
    for table in (halt, loop):
        for mask in range(32):
            ones = frozenset(i - 2 for i in range(5) if mask >> i & 1)
            for q, head in itertools.product(range(2), range(-2, 3)):
                cfg = q, head, ones
                cells = colored_tournament(cfg)
                need(decode(cells) == cfg, 'literal encode/decode')
                raw = step(table, cfg)
                target = graph_step(table, cells)
                need((None if raw is None else colored_tournament(raw)) == target, 'one-step conjugacy')
                graph_cases += 1
    # Actual named-machine traces, checked at every native transition.
    for code in (BOUNCE, CUBIC, EXP, SANDWICH, IRREGULAR, GOLD):
        table, cfg = parse(code), blank
        for unused in range(1000):
            cells = colored_tournament(cfg)
            raw = step(table, cfg)
            need(graph_step(table, cells) == (None if raw is None else colored_tournament(raw)), 'named trace')
            graph_cases += 1
            if raw is None:
                break
            cfg = raw
    # Goldbach native checkpoint convention from the published Lean theorem cfg45.
    cfg, checkpoints, table = blank, [], parse(GOLD)
    for t in range(2000000):
        q, head, ones = cfg
        if q == 21 and head not in ones and len(ones) >= 4 and len(ones) % 2 == 0:
            if ones == frozenset(range(head - len(ones), head)):
                if not checkpoints or len(ones) > checkpoints[-1][1]:
                    checkpoints.append([t, len(ones)])
                    if len(checkpoints) == 10:
                        break
        cfg = step(table, cfg)
        need(cfg is not None, 'bounded Goldbach run')
    need(len(checkpoints) == 10 and checkpoints[0] == [45, 4], 'Goldbach stage checks')
    need([row[1] for row in checkpoints] == list(range(4, 24, 2)), 'consecutive even tests')
    out = {
        'status': 'PASS', 'clock_cases': clock_cases,
        'configuration_tournament_steps': graph_cases,
        'goldbach_native_checkpoints': checkpoints,
        'exponential_blank_times': [12 * 4 ** r - 18 * 2 ** r + r + 6 for r in range(7)],
        'exponential_blank_block_lengths': [6 * (2 ** r - 1) for r in range(7)],
        'hostiles': ['wrong published exponential pair', 'wrong exponential +3 recurrence',
                     'same control graph opposite halting'],
        'scope': 'finite native validation; all-order proofs and cited Lean boundary in README'
    }
    print(json.dumps(out, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
