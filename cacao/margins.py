"""Threshold margins and component calibration stress tests."""
from collections import defaultdict
from fractions import Fraction
from statistics import median

from .chemistry import RULES, band, display, feasible, numeric, read_csv, upper, write_csv

PERCENT_GRID = sorted({Fraction(p) for p in range(21)} | {Fraction(5, 2), Fraction(25), Fraction(50)})
SECV_MULTIPLES = (Fraction(0), Fraction(1, 4), Fraction(1, 2), Fraction(1))


def limits(a, b):
    """Bounds on the change indicator when both endpoints can vary separately."""
    return int(a.isdisjoint(b)), int(len(a | b) > 1)


def cutoff_margins(x, y, rule):
    lo, hi = sorted((x, y))
    cuts = [cut for cut, side in RULES[rule] if upper(lo, cut, side) != upper(hi, cut, side)]
    init_distance, nearest = min((max(x/cut, cut/x), cut) for cut, _ in RULES[rule])
    # Every originally crossed cut must lose its ability to separate the two
    # perturbed ratios before both can share a band. Equality is tested by limits().
    erase = max((min(cut/lo, hi/cut) for cut in cuts), default=None)
    return cuts, nearest, 100*(init_distance-1), None if erase is None else 100*(erase-1)


def analyze_margins(pairs, out):
    details, groups = [], defaultdict(list)
    for p in pairs:
        for rule in RULES:
            first, last = band(p.first, rule), band(p.last, rule)
            cuts, nearest, distance, erase = cutoff_margins(p.first, p.last, rule)
            details.append({"study": p.study, "unit": p.unit, "unit_type": p.unit_type, "farm": p.farm,
                            "ratio_basis": p.basis, "role": p.role, "rule": rule,
                            "first_tc": display(p.first), "last_tc": display(p.last),
                            "first_band": first, "last_band": last,
                            "crossed_boundaries": "|".join(map(str, cuts)),
                            "first_nearest_boundary": nearest, "first_boundary_distance_pct": display(distance),
                            "erase_margin_pct_infimum": display(erase)})
            groups[(p.study, p.basis, p.role, rule)].append((p, first != last, distance, erase))
    summary, stress = [], []
    for (study, basis, role, rule), group in groups.items():
        fixed = {"study": study, "ratio_basis": basis, "role": role, "rule": rule, "units": len(group)}
        changers = [r for r in group if r[1]]
        summary.append({**fixed, "point_changed": len(changers),
                        "up": sum(p.last > p.first for p, _, _, _ in changers),
                        "down": sum(p.last < p.first for p, _, _, _ in changers),
                        "median_first_boundary_distance_pct_changed": display(
                            median(r[2] for r in changers) if changers else None),
                        "median_erase_margin_pct_changed": display(
                            median(r[3] for r in changers) if changers else None)})
        for percent in PERCENT_GRID:
            factor = 1 + percent/100
            minimum = maximum = 0
            for p, _, _, _ in group:
                a = feasible(p.first/factor, p.first*factor, rule)
                b = feasible(p.last/factor, p.last*factor, rule)
                small, large = limits(a, b)
                minimum += small
                maximum += large
            if not minimum <= len(changers) <= maximum <= len(group):
                raise ValueError("Inconsistent bounded change counts.")
            stress.append({**fixed, "factor_pct": display(percent), "factor": display(factor),
                           "point_changed": len(changers), "forced_changed": minimum,
                           "possible_changed": maximum})
    write_csv(out / "unit_margins.csv", details)
    write_csv(out / "endpoint_summary.csv", summary)
    write_csv(out / "ratio_stress.csv", stress)
    return summary, stress


def component_interval(t, c, te, ce):
    if t-te <= 0 or c-ce <= 0:
        raise ValueError("The chosen calibration envelope reaches zero for an analyte.")
    return (t-te)/(c+ce), (t+te)/(c-ce)


def calibration_stress(data, rows, out):
    scales = {r["analyte"]: numeric(r["secv_pct_air_dried_matter"])
              for r in read_csv(data / "nirs_calibration.csv")}
    included = [r for r in rows if all(r[k] for k in ("n_theo", "n_caff", "f_theo", "f_caff"))]
    summary, details = [], []
    for rule in RULES:
        for k in SECV_MULTIPLES:
            point = minimum = maximum = 0
            forced_trees = []
            for row in included:
                nt, nc, ft, fc = (numeric(row[key]) for key in ("n_theo", "n_caff", "f_theo", "f_caff"))
                te, ce = k*scales["theobromine"], k*scales["caffeine"]
                start = component_interval(nt, nc, te, ce)
                end = component_interval(ft, fc, te, ce)
                low, high = limits(feasible(*start, rule), feasible(*end, rule))
                changed = int(band(nt/nc, rule) != band(ft/fc, rule))
                if not low <= changed <= high:
                    raise ValueError("Inconsistent calibration envelope for " + row["tree"])
                if low:
                    forced_trees.append(row["tree"])
                point += changed
                minimum += low
                maximum += high
                details.append({"tree": row["tree"], "farm": row["farm"], "rule": rule,
                                "secv_multiple": display(k), "start_tc_low": display(start[0]),
                                "start_tc_high": display(start[1]), "end_tc_low": display(end[0]),
                                "end_tc_high": display(end[1]), "point_changed": changed,
                                "forced_changed": low, "possible_changed": high})
            summary.append({"study": "Nicaragua", "ratio_basis": "component_means", "rule": rule,
                            "secv_multiple": display(k), "units": len(included),
                            "point_changed": point, "forced_changed": minimum, "possible_changed": maximum,
                            "forced_tree_ids": "|".join(forced_trees)})
    write_csv(out / "nirs_stress.csv", summary)
    write_csv(out / "nirs_tree_detail.csv", details)
    return summary


def colombia_stress(data, out):
    """Propagate both Table 2 SDs and printed precision through the ratio."""
    sd_rows = read_csv(data / "colombia_table2_sd.csv")
    spreads = {(r["material"], int(r["day"])): (Fraction(r["theobromine_sd_mg_g"]),
                                                 Fraction(r["caffeine_sd_mg_g"])) for r in sd_rows}
    means = {(r["material"], int(r["day"])): (numeric(r["theobromine_mg_g"]), numeric(r["caffeine_mg_g"]))
             for r in read_csv(data / "component_means.csv") if r["study"] == "Colombia"}
    if len(spreads) != len(sd_rows) or means.keys() != spreads.keys():
        raise ValueError("Colombia SD and concentration records do not match.")
    units = sorted({name for name, day in means})
    rows, summary = [], []
    half_mean_unit, half_sd_unit = Fraction(5, 1000), Fraction(5, 10000)
    for rule in RULES:
        for mode in ("printed_rounding", "printed_rounding_plus_one_sd"):
            point = minimum = maximum = 0
            forced_units = []
            for unit in units:
                intervals, labels = [], []
                for day in (0, 6):
                    t, c = means[(unit, day)]
                    ts, cs = spreads[(unit, day)]
                    te = ce = half_mean_unit
                    if mode == "printed_rounding_plus_one_sd":
                        te += ts + half_sd_unit
                        ce += cs + half_sd_unit
                    intervals.append(component_interval(t, c, te, ce))
                    labels.append(band(t/c, rule))
                a, b = (feasible(*bounds, rule) for bounds in intervals)
                low, high = limits(a, b)
                changed = int(labels[0] != labels[1])
                if not low <= changed <= high:
                    raise ValueError("Inconsistent SD envelope in Colombia.")
                start_interval, end_interval = intervals
                rows.append({"rule": rule, "model": mode, "material": unit,
                             "start_tc_low": display(start_interval[0]), "start_tc_high": display(start_interval[1]),
                             "end_tc_low": display(end_interval[0]), "end_tc_high": display(end_interval[1]),
                             "start_feasible_bands": "|".join(map(str, sorted(a))),
                             "end_feasible_bands": "|".join(map(str, sorted(b))),
                             "point_changed": changed, "forced_changed": low, "possible_changed": high})
                if low:
                    forced_units.append(unit)
                point += changed
                minimum += low
                maximum += high
            summary.append({"study": "Colombia", "rule": rule, "model": mode, "units": len(units),
                            "point_changed": point, "forced_changed": minimum, "possible_changed": maximum,
                            "forced_units": "|".join(forced_units)})
    write_csv(out / "colombia_sd_detail.csv", rows)
    write_csv(out / "colombia_sd_stress.csv", summary)
    return summary
