"""Paired measurements and ordered intervals for the theobromine-to-caffeine ratio."""
import csv
from collections import defaultdict
from dataclasses import dataclass
from fractions import Fraction

RULES = {"Osorio": ((2, ">"), (6, ">="), (8, ">=")),
         "Pelaez": ((3, ">="), (9, ">"))}


def upper(value, cut, side):
    return value >= cut if side == ">=" else value > cut


def band(value, rule):
    return 1 + sum(upper(value, cut, side) for cut, side in RULES[rule])


def feasible(low, high, rule):
    if high < low or low <= 0:
        raise ValueError("Expected a positive, ordered T/C interval.")
    return set(range(band(low, rule), band(high, rule) + 1))


def read_csv(path):
    with path.open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def write_csv(path, rows, fields=None):
    with path.open("w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields or list(rows[0]), lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def numeric(raw):
    x = Fraction(raw)
    if x <= 0:
        raise ValueError("Expected a positive ratio or concentration.")
    return x


def display(q):
    return "" if q is None else format(float(q), ".10g")


@dataclass(frozen=True)
class Pair:
    study: str
    unit: str
    unit_type: str
    basis: str
    role: str
    first: Fraction
    last: Fraction
    farm: str = ""


def load_pairs(data):
    rows = read_csv(data / "waslala_s1_numeric.csv")
    if any(not r["tree"] or not r["farm"] for r in rows) or len({r["tree"] for r in rows}) != len(rows):
        raise ValueError("Missing identifier or repeated Waslala tree.")
    pairs = []
    for r in rows:
        for basis in ("reported", "component_means"):
            if basis == "reported":
                fields = ("n_tc", "f_tc")
                if not all(r[f] for f in fields):
                    continue
                x, y = (numeric(r[f]) for f in fields)
            else:
                fields = ("n_theo", "n_caff", "f_theo", "f_caff")
                if not all(r[f] for f in fields):
                    continue
                nt, nc, ft, fc = (numeric(r[f]) for f in fields)
                x, y = nt/nc, ft/fc
            pairs.append(Pair("Nicaragua", r["tree"], "tree", basis,
                              "primary" if basis == "reported" else "sensitivity", x, y, r["farm"]))
    series = {}
    for r in read_csv(data / "ecuador_table1.csv"):
        day_values = [(d, numeric(r["day"+str(d)])) for d in range(7) if r["day"+str(d)]]
        days = [d for d, _ in day_values]
        if days not in [list(range(6)), list(range(7))] or r["material"] in series:
            raise ValueError("Invalid Ecuador series for " + r["material"])
        values = tuple(v for d, v in day_values)
        series[r["material"]] = values
        pairs.append(Pair("Ecuador", r["material"], r["material_type"], "reported",
                          "primary", values[0], values[-1]))
    groups = defaultdict(dict)
    for r in read_csv(data / "component_means.csv"):
        key = r["study"], r["material"], r["material_type"]
        day = int(r["day"])
        if day in groups[key]:
            raise ValueError("Duplicate material and day: " + str(key))
        groups[key][day] = r
    for (study, material, kind), group in groups.items():
        days = sorted(group)
        if days != {"Colombia": [0, 6], "Venezuela": [0, 3]}[study]:
            raise ValueError("Unexpected fermentation stages for " + material)
        x, y = (numeric(group[d]["theobromine_mg_g"])/numeric(group[d]["caffeine_mg_g"]) for d in days)
        role = "sensitivity" if study == "Venezuela" else "primary"
        pairs.append(Pair(study, material, kind, "component_means", role, x, y))
        reported = [group[d]["reported_tc"] for d in days]
        if any(reported):
            if not all(reported):
                raise ValueError("Incomplete reported ratios for " + material)
            pairs.append(Pair(study, material, kind, "reported", "primary", *map(numeric, reported)))
    return rows, pairs, series
