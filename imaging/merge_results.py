"""merge_results.py - Single-Cell-Results_*.csv eines Experiments zu Combined_Results.csv zusammenfuegen.

    python merge_results.py --exp-dir /prj/.../Data/WT/pH/6 [--results-subdir 03_results_v12]
    python merge_results.py --all /prj/.../Data [--results-subdir 03_results_v12]
"""
import argparse
from pathlib import Path
import pandas as pd


def merge_one(exp_dir: Path, results_subdir: str) -> int:
    d = exp_dir / results_subdir
    files = sorted(d.glob("Single-Cell-Results_*.csv"))
    if not files:
        return 0
    df = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
    df.to_csv(d / "Combined_Results.csv", index=False)
    print(f"{d}: {len(files)} Filme, {len(df)} Zeilen -> Combined_Results.csv")
    return len(files)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--exp-dir", type=Path); ap.add_argument("--all", type=Path)
    ap.add_argument("--results-subdir", default="03_results_v12"); a = ap.parse_args()
    if a.all:
        n = 0
        for d in sorted(a.all.glob("*/*/*")):
            if (d / a.results_subdir).is_dir():
                n += merge_one(d, a.results_subdir)
        print(f"{n} Filme zusammengefuehrt")
    elif a.exp_dir:
        merge_one(a.exp_dir, a.results_subdir)
    else:
        ap.error("--exp-dir oder --all")


if __name__ == "__main__":
    main()
