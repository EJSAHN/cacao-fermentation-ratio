"""Tree-level change rates, whole-farm resampling and single-farm omissions."""
import math
import random
from collections import Counter, defaultdict
from decimal import Decimal

from .chemistry import RULES, band, read_csv, upper, write_csv
from .farms import DRAWS, SEED, quantile


def quantity(raw):
    if not raw:
        return None
    value = Decimal(raw)
    if not value.is_finite() or value <= 0:
        raise ValueError("Unexpected nonpositive or nonfinite concentration or ratio: " + raw)
    return value


def tree_data(path):
    pairs, excluded, ids = [], [], set()
    for raw in read_csv(path):
        tree = raw["tree"]
        if not tree or tree in ids or not raw["farm"]:
            raise ValueError("Missing tree/farm ID or repeated tree: " + tree)
        ids.add(tree)
        values = {key: quantity(raw[key]) for key in
                  ("n_tc", "f_tc", "n_theo", "n_caff", "f_theo", "f_caff")}
        comparisons = {"reported": (values["n_tc"], values["f_tc"])}
        for stage in ("n", "f"):
            t, c = values[stage + "_theo"], values[stage + "_caff"]
            values[stage + "_recalculated"] = t / c if t is not None and c is not None else None
        comparisons["component_means"] = (values["n_recalculated"], values["f_recalculated"])
        row = dict(raw)
        for method, (before, after) in comparisons.items():
            row[method + "_n"] = before
            row[method + "_f"] = after
            paired = before is not None and after is not None
            if not paired:
                excluded.append({"tree": tree, "farm": raw["farm"], "method": method,
                                 "missing": ",".join(s for s, v in [("N", before), ("F", after)]
                                                     if v is None)})
            row[method + "_log_change"] = math.log(after / before) if paired else None
            for rule in RULES:
                key = method + "_" + rule
                row[key + "_n_band"] = band(before, rule) if before is not None else None
                row[key + "_f_band"] = band(after, rule) if after is not None else None
                row[key + "_change"] = int(band(before, rule) != band(after, rule)) if paired else None
        pairs.append(row)
    return pairs, excluded


def farm_summary(rows, column, draws, seed, exponentiate=False):
    groups = defaultdict(list)
    for row in rows:
        if row[column] is not None:
            groups[row["farm"]].append(row[column])
    if len(groups) < 2:
        raise ValueError("At least two farms with paired measurements needed for " + column)
    farms = sorted(groups, key=int)
    blocks = [(len(groups[f]), sum(groups[f])) for f in farms]
    n = sum(size for size, _ in blocks)
    total = sum(value for _, value in blocks)
    generator = random.Random(seed)
    sampled = []
    for _ in range(draws):
        sample = generator.choices(blocks, k=len(farms))
        sampled.append(sum(v for _, v in sample) / sum(size for size, _ in sample))
    sampled.sort()
    transform = math.exp if exponentiate else float
    omissions = []
    for farm, (size, value) in zip(farms, blocks):
        omissions.append({"metric": column, "omitted_farm": farm, "trees_removed": size,
                          "trees_remaining": n-size, "estimate": transform((total-value)/(n-size))})
    omission_values = [r["estimate"] for r in omissions]
    summary = {"metric": column, "trees": n, "farms": len(farms),
               "changed": None if exponentiate else int(total),
               "farms_with_changes": None if exponentiate else sum(v > 0 for _, v in blocks),
               "estimate": transform(total/n),
               "cluster_ci_low": transform(quantile(sampled, 0.025)),
               "cluster_ci_high": transform(quantile(sampled, 0.975)),
               "leave_one_farm_min": min(omission_values), "leave_one_farm_max": max(omission_values)}
    return summary, omissions


def rule_details(trees):
    common = [r for r in trees if r["reported_log_change"] is not None]
    counts = Counter((r["reported_Osorio_change"], r["reported_Pelaez_change"]) for r in common)
    overlap = [{"osorio_changed": a, "pelaez_changed": b, "trees": counts[a, b]}
               for a in (0, 1) for b in (0, 1)]
    bounds, transitions = [], []
    for rule in RULES:
        for cut, side in RULES[rule]:
            before = [upper(r["reported_n"], cut, side) for r in common]
            after = [upper(r["reported_f"], cut, side) for r in common]
            up = sum((not a) and b for a, b in zip(before, after))
            down = sum(a and (not b) for a, b in zip(before, after))
            bounds.append({"rule": rule, "upper_band_begins": side+str(cut), "up": up,
                           "down": down, "any_crossing": up+down, "trees": len(common)})
        matrix = Counter((band(r["reported_n"], rule), band(r["reported_f"], rule)) for r in common)
        levels = range(1, len(RULES[rule]) + 2)
        transitions.extend({"rule": rule, "before_band": a, "after_band": b, "trees": matrix[a,b]}
                           for a in levels for b in levels)
    return overlap, bounds, transitions


def analyze_trees(data, out):
    trees, exclusions = tree_data(data / "waslala_s1_numeric.csv")
    summary, omissions = [], []
    for method in ("reported", "component_means"):
        metrics = [method + "_" + rule + "_change" for rule in RULES]
        metrics.append(method + "_log_change")
        for metric in metrics:
            result, leave = farm_summary(trees, metric, DRAWS, SEED, metric.endswith("_log_change"))
            summary.append(result)
            omissions.extend(leave)
    overlap, boundaries, transitions = rule_details(trees)
    write_csv(out / "tree_results.csv", trees)
    write_csv(out / "excluded_by_method.csv", exclusions, ["tree", "farm", "method", "missing"])
    write_csv(out / "farm_cluster_estimates.csv", summary)
    write_csv(out / "leave_one_farm_out.csv", omissions)
    write_csv(out / "rule_overlap.csv", overlap)
    write_csv(out / "boundary_crossings.csv", boundaries)
    write_csv(out / "transition_counts.csv", transitions)
