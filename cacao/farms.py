"""Pair both band rules within each Waslala tree and resample whole farms."""
import math
import random
from collections import Counter, defaultdict

from .chemistry import RULES, band, write_csv

SEED = 20261119
DRAWS = 10000


def quantile(sorted_values, p):
    at = p*(len(sorted_values)-1)
    a, b = math.floor(at), math.ceil(at)
    return sorted_values[a] + (sorted_values[b] - sorted_values[a])*(at-a)


def contrast(pairs, out):
    trees = [p for p in pairs if p.study == "Nicaragua"]
    by_basis = defaultdict(dict)
    for p in trees:
        by_basis[p.basis][p.unit] = {r: band(p.last, r) - band(p.first, r) for r in RULES}
    results, agreements = [], []
    for basis, values in by_basis.items():
        counts = Counter()
        farms = defaultdict(list)
        for p in trees:
            if p.basis != basis:
                continue
            four, three = (int(values[p.unit][r] != 0) for r in ("Osorio", "Pelaez"))
            farms[p.farm].append(four-three)
            counts[four, three] += 1
        blocks = [(sum(farms[f]), len(farms[f])) for f in sorted(farms, key=int)]
        total = sum(x for x, _ in blocks)
        n = sum(size for _, size in blocks)
        rng = random.Random(SEED)
        reps = []
        for _ in range(DRAWS):
            sample = rng.choices(blocks, k=len(blocks))
            reps.append(sum(x for x, _ in sample)/sum(size for _, size in sample))
        reps.sort()
        results.append({"ratio_basis": basis, "trees": n, "farms": len(blocks),
                        "both": counts[1, 1], "four_only": counts[1, 0],
                        "three_only": counts[0, 1], "neither": counts[0, 0],
                        "four_minus_three_pp": 100*total/n,
                        "farm_boot_low_pp": 100*quantile(reps, 0.025),
                        "farm_boot_high_pp": 100*quantile(reps, 0.975),
                        "bootstrap_draws": DRAWS, "seed": SEED})
    shared = by_basis["reported"].keys() & by_basis["component_means"].keys()
    for rule in RULES:
        counts = Counter()
        same_direction = 0
        for tree in shared:
            a = by_basis["reported"][tree][rule]
            b = by_basis["component_means"][tree][rule]
            counts[bool(a), bool(b)] += 1
            same_direction += int(a*b > 0)
        agreements.append({"rule": rule, "trees_with_both_bases": len(shared),
                           "both_changed": counts[True, True], "both_changed_same_direction": same_direction,
                           "reported_only": counts[True, False], "component_only": counts[False, True],
                           "neither": counts[False, False]})
    write_csv(out / "paired_rule_contrast.csv", results)
    write_csv(out / "ratio_definition_agreement.csv", agreements)
    return results, agreements
