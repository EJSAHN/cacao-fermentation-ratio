"""Replay all Ecuador sampling schedules on the shared days 0 through 5."""
from collections import Counter
from fractions import Fraction
from itertools import combinations, product

from .chemistry import RULES, band, feasible, write_csv


def assignments(values, rule, model):
    if model == "printed_rounding":
        half_unit = Fraction(1, 20)
        intervals = [(v-half_unit, v+half_unit) for v in values]
    elif model == "ratio_factor_1.02":
        f = Fraction(102, 100)
        intervals = [(v/f, v*f) for v in values]
    else:
        raise ValueError("Unknown sensitivity model: " + model)
    return list(product(*(sorted(feasible(lo, hi, rule)) for lo, hi in intervals)))


def outcomes(path, days):
    daily = any(x != path[0] for x in path[1:])
    detected = any(path[d] != path[0] for d in days)
    return int(daily), int(detected), int(daily and not detected)


def analyze_schedules(series, out):
    group_rows, materials, extended = [], [], []
    for rule in RULES:
        for model in ("printed_rounding", "ratio_factor_1.02"):
            points = {name: tuple(band(v, rule) for v in values) for name, values in series.items()}
            options = {name: assignments(values[:6], rule, model) for name, values in series.items()}
            for name, values in series.items():
                if len(values) == 7:
                    full = assignments(values, rule, model)
                    def added(path):
                        return int(len(set(path[:6])) == 1 and path[-1] != path[0])
                    new_events = [added(path) for path in full]
                    extended.append({"material": name, "rule": rule, "model": model,
                                     "point_new_event_day6": added(points[name]),
                                     "min_new_event_day6": min(new_events), "max_new_event_day6": max(new_events)})
            for count in range(5):
                for extras in combinations(range(1, 5), count):
                    days = (0,) + extras + (5,)
                    keys = {"rule": rule, "model": model, "days": "|".join(map(str, days)),
                            "extra_visits": count}
                    totals = Counter()
                    for name, states in points.items():
                        point = outcomes(states[:6], days)
                        all_possible = [outcomes(path, days) for path in options[name]]
                        record = {}
                        for j, event in enumerate(("daily_event", "detected", "missed")):
                            choices = [opt[j] for opt in all_possible]
                            record["point_"+event] = point[j]
                            record["min_"+event] = min(choices)
                            record["max_"+event] = max(choices)
                        if point[0] != point[1]+point[2]:
                            raise ValueError("Daily events do not reconcile.")
                        totals.update(record)
                        materials.append({**keys, "material": name, **record})
                    group_rows.append({**keys, "materials": len(series), **dict(totals)})
    write_csv(out / "ecuador_schedules.csv", group_rows)
    write_csv(out / "ecuador_schedule_materials.csv", materials)
    write_csv(out / "ecuador_day6_extension.csv", extended)
    return group_rows
