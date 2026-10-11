"""Source-owned H/G/R credit receipts and their native-language boundary.

No production orbit discovery or universal ROOT claim.  Existing frozen
controllers and run receipts are authenticated rather than reclassified.
"""
from dataclasses import dataclass, replace
from fractions import Fraction as F
from itertools import product

import adaptive_credit_potential_20261004 as old
import collatz_projective_receipts_20261009 as runs
import paper_local_global_grounding_20261007 as native
import translation_phase_decoder_20261005 as five

KAPPA_R = F(57, 64)
CHECKS = 0


def need(ok, why):
    if not ok:
        raise ValueError(why)


def check(ok, why):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(why)


@dataclass(frozen=True)
class Action:
    letter: str
    run: int = 0


def audit_action(action):
    need(type(action) is Action and type(action.letter) is str
         and action.letter in ('H', 'G', 'R'), 'exact typed H/G/R action')
    old.integer(action.run)
    need(action.run >= 4 if action.letter == 'R' else action.run == 0,
         'R retains run >=4; H/G have canonical run zero')
    return action


def transition(current, action):
    """Return next state, credit change, sharp factor, and actual receipt."""
    old.positive_odd(current)
    audit_action(action)
    if action.letter == 'R':
        rec = runs.join(current, action.run)
        receipt = native.Receipt(current, rec.child, rec.source_word, rec.child_word)
        result = rec.child, 1, KAPPA_R, receipt
    else:
        following = old.action(current, action.letter)
        if action.letter == 'H':
            _, a = old.step(following)
            receipt = native.Receipt(current, following, old.HEAD+(a+2,), (a,))
            result = following, 2, old.KAPPA, receipt
        else:
            receipt = native.Receipt(current, following, (1, 2), ())
            result = following, -1, F(1), receipt
    native.verify_receipt(result[-1])
    return result


@dataclass(frozen=True)
class Episode:
    source: int
    actions: tuple
    states: tuple
    credits: tuple
    receipts: tuple


def execute(source, actions):
    """Zero-credit start. Every action and prefix is checked at this source."""
    old.positive_odd(source)
    need(type(actions) is tuple, 'ordered action tuple')
    current, credit, bound = source, 0, F(source+5)
    states, credits, receipts = [source], [0], []
    for action in actions:
        audit_action(action)
        need(action.letter != 'G' or credit >= 1, 'G requires available credit')
        following, award, factor, receipt = transition(current, action)
        prior = old.potential(current, credit)
        credit += award
        current = following
        energy = old.potential(current, credit)
        bound *= factor
        need(energy <= prior*factor and 116 <= energy <= bound < source+5,
             'native potential contract and immutable entry bound')
        need(current < source, 'strict immutable-source payment')
        states.append(current); credits.append(credit); receipts.append(receipt)
    return Episode(source, actions, tuple(states), tuple(credits), tuple(receipts))


def audit_episode(episode):
    need(type(episode) is Episode, 'typed Episode')
    old.positive_odd(episode.source)
    need(type(episode.states) is tuple and type(episode.credits) is tuple
         and type(episode.receipts) is tuple, 'exact cached tuples')
    for value in episode.states:
        old.positive_odd(value)
    for value in episode.credits:
        old.integer(value)
    for receipt in episode.receipts:
        native.verify_receipt(receipt)
    need(episode == execute(episode.source, episode.actions), 'receipt/cache provenance')
    return episode


def continue_episode(episode, actions):
    audit_episode(episode)
    need(type(actions) is tuple, 'ordered extension tuple')
    return execute(episode.source, episode.actions+actions)


def discharge(episode, supplied_root):
    """Ground only with the supplied final-state first-hit ROOT certificate."""
    audit_episode(episode)
    terminal = episode.states[-1]
    need(native.replay(terminal, supplied_root) == 1, 'supplied terminal ROOT word')
    final = native.Receipt(terminal, 1, supplied_root, ())
    path = tuple((rec, True) for rec in episode.receipts+(final,))
    _, word = native.authenticate_path(path, episode.source)
    return word


def funding_bound(source):
    """Upper bound on H+R count; total executed actions are at most threefold."""
    old.positive_odd(source)
    bound, count = F(source+5), 0
    while bound*old.KAPPA >= 116:
        bound *= old.KAPPA
        count += 1
    return count


def language_kind(episode):
    audit_episode(episode)
    letters = tuple(a.letter for a in episode.actions)
    if not letters:
        return 'EMPTY'
    if all(a == 'R' for a in letters):
        return 'RUN_ONLY'
    need(letters[0] == 'H' and all(a in ('H', 'G') for a in letters),
         'proved disjoint native episode languages')
    return 'HG_ONLY'


def bad(call):
    try:
        call()
    except ValueError:
        check(True, 'malformed/guard boundary rejected')
    else:
        raise RuntimeError('hostile accepted')


def main():
    h, g = Action('H'), Action('G')
    count = 0
    for length in range(1, 7):
        for letters in product('HG', repeat=length):
            word = ''.join(letters)
            spec = old.summary(word)
            if spec.minimum < 0:
                continue
            residue, modulus = old.cylinder(spec)
            for lift in (0, 1, 7):
                source = residue+lift*modulus
                episode = execute(source, tuple(Action(a) for a in letters))
                check(episode.states[-1] == old.verify_episode(source, word)[0],
                      'independent inherited H/G endpoint')
                check(language_kind(episode) == 'HG_ONLY', 'H/G language')
                check(sum(a.letter != 'G' for a in episode.actions) <= funding_bound(source),
                      'exact funding bound')
                for x in episode.states[1:]:
                    check(x % 9 in (3, 4) and runs.available(x) < 4,
                          'H/G cannot reach native R')
                count += 1

    run_count = 0
    for r in range(4, 17):
        for parameter in range(65):
            source, child, _ = runs.ray(r, parameter)
            episode = execute(source, (Action('R', r),))
            check(episode.states[-1] == child and language_kind(episode) == 'RUN_ONLY',
                  'source-owned run adapter')
            check(F(child+5, source+5) <= F(19, 24), 'uncredited sharp run ratio')
            check(old.potential(child, 1) <= KAPPA_R*(source+5), 'one-credit run ratio')
            check(child % 16 == 15, 'all old letters blocked after R')
            for letter, (_, _, _, residue, modulus) in five.LETTERS.items():
                check(child % modulus != residue, 'five-letter outgoing guard barrier')
            if r == 4:
                check(old.potential(child, 2)/(source+5) == 1+F(9, 16*(source+5)),
                      'two-credit potential failure, not legal RGG')
            run_count += 1
    check(old.potential(223, 1)/288 == KAPPA_R, 'sharp bound attained')
    check(old.potential(223, 2)/288 == F(513, 512), 'sharp two-credit hostile')
    for parameter in range(100):
        source = 13891+26244*parameter
        first = execute(source, (Action('R', 4),))
        both = continue_episode(first, (Action('R', 4),))
        check(both.states == (source, 10975+20736*parameter, 8671+16384*parameter),
              'same-original-source run return')
        check(both.credits == (0, 1, 2), 'stranded terminal credits retained')
        check(old.potential(both.states[-1], 2) <= (source+5)*KAPPA_R**2,
              'two-run potential composition')

    output_phases = {'H': 3, 'G': 4, 'A': 2, 'B': 5, 'L': 6}
    for letter, (r, a, b, residue, modulus) in five.LETTERS.items():
        p, q = 3**r, 1 << a
        check(p % 9 == 0 and (b*pow(q, -1, 9)) % 9 == output_phases[letter],
              'all-height old-letter output phase')
        check(residue % 16 == (7 if letter == 'B' else 11), 'outgoing binary guard class')
        for lift in range(64):
            x = residue+modulus*lift
            y = (p*x+b)//q
            check(y % 9 == output_phases[letter], 'independent native phase controls')
            if letter != 'G':
                check(runs.available(y) < 4, 'only old terminal G can enter R at all')
        # A preceding old-letter output is never8 mod9, so the following G
        # has v3(2G(x)+1)=2+v3(x+1)<=3, irrespective of binary admissibility.
        check(output_phases[letter] != 8, 'two-old-letter incoming obstruction')

    for t in range(128):
        x, n, m = 251+288*t, 283+324*t, 223+256*t
        check(old.action(x, 'G') == n and runs.join(n, 4).child == m,
              'isolated pre-funded incoming G/R ray')
        check(old.potential(m, 1) <= KAPPA_R*old.potential(x, 1),
              'credit1 entry G0 R1 stays within old source account')
        check(F(8*n-5, 9) == x and old.potential(x, 1) == n+5,
              'neutral inverse-G credit recycling')
        check(old.action(x, 'G') == n, 'literal C/G neutral cycle')

    # Supplied ROOT sidecars only: R uses the frozen223 word; H uses the
    # checked suffix of the inherited frozen27 word after its observed167.
    ep = execute(283, (Action('R', 4),))
    word283 = discharge(ep, runs.ROOT223)
    check(word283 == runs.discharge(runs.join(283, 4), runs.ROOT223), 'independent R splice')
    root111 = (1, 1, 1, 2, 1, 2, 1, 1, 1, 1, 1, 3, 1, 1, 1, 4, 2, 2, 4, 3, 1, 1, 5, 4)
    check(native.replay(111, root111) == 1, 'literal supplied H terminal ROOT')
    word155 = discharge(execute(155, (h,)), root111)
    check(native.replay(155, word155) == 1, 'H original-source ROOT export')
    check(discharge(execute(1, ()), ()) == (), 'ROOT identity boundary')
    for source in (True, 283.0, 0, 2):
        bad(lambda n=source: execute(n, ()))
    for action in (Action('R', True), Action('R', 4.0), Action('H', 1), Action('X'), 'R'):
        bad(lambda a=action: execute(283, (a,)))
    for source, acts in ((283, (g,)), (283, (Action('R', 4), g)),
                         (283, (Action('R', 4), h)), (155, (h, Action('R', 4)))):
        bad(lambda n=source, a=acts: execute(n, a))
    bad(lambda: execute(283, [Action('R', 4)]))
    bad(lambda: audit_episode(replace(ep, credits=(0, True))))
    bad(lambda: audit_episode(replace(ep, states=(283, 223.0))))
    bad(lambda: discharge(ep, runs.ROOT223+(2,)))
    bad(lambda: discharge(ep, (4,)))
    print('PROVED: H(+2), G(-1), R_r(+1) retain E=(x+5)(9/8)^credit and the original source.')
    print('Sharp R ratio57/64 at283->223; R4+2 gives1+9/[16(n+5)]>1, least513/512.')
    print('R outputs15mod16: G is not native, so this is not a legal RGG overspend example.')
    print('Native zero-credit episode language splits: H/G only or R only. No new mixed switching claim.')
    print('Old five-letter outputsmod9=(3,4,2,5,6); no two-old-letter prefix can enter R.')
    print('Isolated pre-funded ray:251+288t --G-->283+324t --R4-->223+256t.')
    print('Neutral inverse-G(+1)/G(-1) cycles preserve E forever; they are excluded from the finite-episode theorem.')
    print('H/G native funded controls:', count, '; R native controls:', run_count, '; repeated R rays:100.')
    print('Supplied ROOT exports:283 rank/cost', len(word283), sum(word283),
          ';155 rank/cost', len(word155), sum(word155))
    print('No new source-guard mass, universal coverage, or giant ROOT proof.')
    print('Checks:', CHECKS)


if __name__ == '__main__':
    main()
