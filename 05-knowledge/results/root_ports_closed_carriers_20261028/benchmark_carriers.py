#!/usr/bin/env python3
"""Paired production benchmark; no mathematical audit costs in the timing."""
from statistics import median
from time import perf_counter
import json
import gc
import check_carriers as new


def run(which):
    gc.collect()
    started = perf_counter(); parents = [()]
    for n in range(8):
        if which == 'baseline':
            parents, report = new.old.produce(parents, n, audit=False, selector='neighbor_degrees')
        else:
            parents, report = new.produce(parents, n)
    elapsed = perf_counter()-started
    new.need(len(parents) == 6880, 'benchmark scope')
    return elapsed


def main():
    samples = {'baseline': [], 'refined': []}
    for repetition in range(6):
        order = ('baseline', 'refined') if repetition % 2 == 0 else ('refined', 'baseline')
        for name in order:
            samples[name].append(run(name))
    medians = {name: median(xs) for name, xs in samples.items()}
    return dict(status='VERIFIED local timing, no asymptotic claim', max_order=8,
                representatives=6880, alternating_pairs=6, seconds=samples,
                median_seconds=medians, median_speed_ratio=medians['baseline']/medians['refined'],
                baseline='Session87 neighbor_degrees, production-only',
                new='degree prefilter plus stable selector and individualization/refinement')


if __name__ == '__main__':
    print(json.dumps(main(), indent=2, sort_keys=True))
