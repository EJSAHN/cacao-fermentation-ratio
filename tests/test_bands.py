"""Regression tests for published cutoffs and sampling-day comparisons."""
import unittest
from fractions import Fraction
from pathlib import Path
from tempfile import TemporaryDirectory

from cacao.chemistry import band, feasible, load_pairs
from cacao.schedules import analyze_schedules
from cacao.margins import limits

DATA = Path(__file__).resolve().parents[1] / "data"


class IntervalTests(unittest.TestCase):
    def test_exact_cutoffs(self):
        for ratio, expected in [(2, 1), (6, 3), (8, 4)]:
            self.assertEqual(band(Fraction(ratio), "Osorio"), expected)
        for ratio in (3, 9):
            self.assertEqual(band(Fraction(ratio), "Pelaez"), 2)

    def test_straddling_the_boundary_allows_both_assignments(self):
        uncertain = feasible(Fraction("5.95"), Fraction("6.05"), "Osorio")
        below = feasible(Fraction("5.65"), Fraction("5.75"), "Osorio")
        self.assertEqual(limits(uncertain, below), (0, 1))

    def test_ecuador_cutoff_days_and_factor(self):
        _, _, series = load_pairs(DATA)
        upper_days = {"CCAT 4675": [1, 2, 3], "EET 103": [2, 5]}
        boundary = Fraction(6)
        for material, expected_days in upper_days.items():
            values = series[material][:6]
            self.assertEqual([d for d, r in enumerate(values) if band(r, "Osorio") == 3],
                             expected_days)
            peak = max(values)
            self.assertEqual(peak / boundary, Fraction(61, 60))
            self.assertEqual(band(peak / Fraction(61, 60), "Osorio"), 3)
            self.assertEqual(band(peak / Fraction("1.02"), "Osorio"), 2)

    def test_schedule_endpoints_and_intermediate_visits(self):
        _, _, series = load_pairs(DATA)
        with TemporaryDirectory() as folder:
            schedules = analyze_schedules(series, Path(folder))
        rounding = {row["days"]: row for row in schedules
                    if row["rule"] == "Osorio" and row["model"] == "printed_rounding"}
        self.assertEqual(len(rounding), 16)
        self.assertEqual(rounding["0|5"]["point_missed"], 1)
        self.assertEqual(rounding["0|2|3|5"]["min_detected"], 2)
        self.assertEqual(rounding["0|2|3|5"]["max_missed"], 0)


if __name__ == "__main__":
    unittest.main()
