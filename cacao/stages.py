"""Rounding intervals and stage-by-stage changes in published fermentation series."""
from collections import defaultdict
from dataclasses import dataclass
from fractions import Fraction

from .chemistry import RULES, band, display, read_csv, write_csv


def printed(raw):
    places = len(raw.partition(".")[2])
    step = Fraction(1, 2 * 10 ** places)
    center = Fraction(raw)
    if center <= step:
        raise ValueError("The rounding interval must be positive: " + raw)
    return center, center - step, center + step


@dataclass(frozen=True)
class Reading:
    day: int
    ratio: Fraction
    low: Fraction
    high: Fraction

    def bands(self, rule):
        return set(range(band(self.low, rule), band(self.high, rule) + 1))


def read_series(data):
    series = defaultdict(list)
    for row in read_csv(data / "ecuador_table1.csv"):
        key = ("Ecuador", row["material"], row["material_type"], "reported")
        for day in range(7):
            raw = row["day" + str(day)]
            if raw:
                series[key].append(Reading(day, *printed(raw)))
    for row in read_csv(data / "component_means.csv"):
        study = row["study"]
        key = (study, row["material"], row["material_type"])
        day = int(row["day"])
        t, t_low, t_high = printed(row["theobromine_mg_g"])
        c, c_low, c_high = printed(row["caffeine_mg_g"])
        series[key + ("component_means",)].append(Reading(day, t/c, t_low/c_high, t_high/c_low))
        if row["reported_tc"]:
            series[key + ("reported",)].append(Reading(day, *printed(row["reported_tc"])))
    for (study, material, _, _), readings in series.items():
        readings.sort(key=lambda r: r.day)
        days = tuple(r.day for r in readings)
        if study == "Ecuador":
            valid_days = (tuple(range(6)), tuple(range(7)))
        elif study == "Colombia":
            valid_days = ((0, 6),)
        elif study == "Venezuela":
            valid_days = ((0, 3),)
        else:
            raise ValueError("Unrecognized study: " + study)
        if days not in valid_days:
            raise ValueError("Missing or repeated day in " + study + " / " + material)
    return series


def change_bounds(sets):
    """One state per stage. A shared feasible band permits a constant path."""
    required = not set.intersection(*sets)
    possible = len(set.union(*sets)) > 1
    return int(required), int(possible)


def return_bounds(sets):
    """The path leaves its starting band and ends back in that band."""
    first, last = sets[0], sets[-1]
    common = first & last
    middle = sets[1:-1]
    possible = any(any(s - {b} for s in middle) for b in common)
    required = first == last and len(first) == 1 and any(first.isdisjoint(s) for s in middle)
    return int(required), int(possible)


def endpoint(readings, rule):
    first, last = readings[0], readings[-1]
    a, b = band(first.ratio, rule), band(last.ratio, rule)
    required, possible = change_bounds([first.bands(rule), last.bands(rule)])
    cuts = [str(cut) for cut, _ in RULES[rule][min(a, b)-1:max(a, b)-1]]
    return {"first_day": first.day, "last_day": last.day,
            "first_tc": display(first.ratio), "last_tc": display(last.ratio),
            "first_band": a, "last_band": b, "boundaries_crossed": "|".join(cuts),
            "direction": "up" if b > a else "down" if b < a else "same",
            "changed": int(a != b), "rounding_min_changed": required,
            "rounding_max_changed": possible}


def trajectory(readings, rule):
    baseline = band(readings[0].ratio, rule)
    states = [band(r.ratio, rule) for r in readings]
    sets = [r.bands(rule) for r in readings]
    away = [r.day for r, state in zip(readings, states) if state != baseline]
    different = bool(away)
    must_leave, may_leave = change_bounds(sets)
    must_return, may_return = return_bounds(sets)
    return {"days_observed": "|".join(str(r.day) for r in readings),
            "bands_by_day": "|".join(map(str, states)),
            "days_outside_initial_band": "|".join(map(str, away)),
            "adjacent_changes": sum(a != b for a, b in zip(states, states[1:])),
            "ever_changed": int(different), "rounding_min_ever": must_leave,
            "rounding_max_ever": may_leave,
            "back_in_initial_band_at_end": int(different and states[-1] == baseline),
            "rounding_min_back": must_return, "rounding_max_back": may_return}


def aggregate(rows, keys, count_fields):
    result = {}
    for row in rows:
        key = tuple(row[k] for k in keys)
        if key not in result:
            result[key] = {k: row[k] for k in keys}
            result[key].update({"units": 0, **{k: 0 for k in count_fields}})
        summary = result[key]
        summary["units"] += 1
        for field in count_fields:
            summary[field] += row[field]
    return list(result.values())


def analyze_external(data, out):
    stages, pairs, paths, same_day = [], [], [], []
    for (study, material, kind, basis), readings in read_series(data).items():
        role = "sensitivity" if (study, basis) == ("Venezuela", "component_means") else "primary"
        info = {"study": study, "material": material, "material_type": kind,
                "ratio_basis": basis, "role": role}
        for r in readings:
            stage = {**info, "day": r.day, "tc": display(r.ratio),
                     "rounding_low": display(r.low), "rounding_high": display(r.high)}
            for rule in RULES:
                stage[rule + "_band"] = band(r.ratio, rule)
                stage[rule + "_feasible_bands"] = "|".join(map(str, sorted(r.bands(rule))))
            stages.append(stage)
        for rule in RULES:
            row = {**info, "rule": rule}
            pairs.append({**row, **endpoint(readings, rule)})
            if len(readings) > 2:
                paths.append({**row, **trajectory(readings, rule)})
                for r in readings[1:]:
                    same_day.append({**row, "day": r.day, **endpoint([readings[0], r], rule)})
    keys = ("study", "ratio_basis", "role", "rule")
    fields = ("changed", "rounding_min_changed", "rounding_max_changed")
    endpoints = aggregate(pairs, keys, fields)
    time_summary = aggregate(paths, keys, ("ever_changed", "rounding_min_ever", "rounding_max_ever",
                                           "back_in_initial_band_at_end", "rounding_min_back",
                                           "rounding_max_back"))
    daily = sorted(aggregate(same_day, ("study", "rule", "day"), fields),
                   key=lambda r: (r["rule"], r["day"]))
    for filename, rows in (("external_stages.csv", stages), ("external_endpoints.csv", pairs),
                           ("external_summary.csv", endpoints), ("ecuador_paths.csv", paths),
                           ("ecuador_path_summary.csv", time_summary), ("ecuador_by_day.csv", daily)):
        write_csv(out / filename, rows)
    return endpoints, time_summary
