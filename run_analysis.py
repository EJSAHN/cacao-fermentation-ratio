"""Classify paired cocoa bean samples across fermentation stages."""
import argparse
from pathlib import Path
from tempfile import TemporaryDirectory

from cacao.chemistry import load_pairs
from cacao.farms import contrast
from cacao.margins import analyze_margins, calibration_stress, colombia_stress
from cacao.schedules import analyze_schedules
from cacao.stages import analyze_external
from cacao.trees import analyze_trees


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=Path(__file__).resolve().parent / "data",
                        help="Directory containing the five input CSV files.")
    parser.add_argument("--output", type=Path, default=Path("results"),
                        help="New directory for the output tables (default: results).")
    args = parser.parse_args()
    output = args.output.expanduser().resolve()
    if output.exists():
        parser.error(f"Output already exists: {output}. Choose a new directory with --output.")
    output.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(dir=output.parent) as temporary:
        tables = Path(temporary) / "tables"
        tables.mkdir()
        rows, pairs, series = load_pairs(args.data)
        analyze_external(args.data, tables)
        analyze_trees(args.data, tables)
        analyze_margins(pairs, tables)
        calibration_stress(args.data, rows, tables)
        colombia_stress(args.data, tables)
        analyze_schedules(series, tables)
        contrast(pairs, tables)
        tables.rename(output)
    print(f"Wrote {len(list(output.glob('*.csv')))} tables to {output}")


if __name__ == "__main__":
    main()
