#!/usr/bin/env python3
"""Convert the agroramsse datasets to a plain text format the Rust code reads.

Mirrors agroramsse/utils.py load_dataset / load_dataset_3d exactly:
  * list/JSON datasets: duplicate points removed, every point gets value 1;
  * dict datasets (cali): keys are points, the dict value is the record value.
1D files (amazon-books, gowalla-1d) are one integer per line; agroramsse reads row[0]
  and removes duplicates.

Output format (<name>.pts):
  line 1: "<d> <N>"
  then N lines: "<x_1> ... <x_d> <value>"   (integers)
Usage: convert_datasets.py <agroramsse/datasets dir> <output dir>
"""
import json, os, pickle, sys, csv

SETS = {
    # name: (file, kind)
    "amazon-books": ("amazon-books.csv", "csv1d"),
    "gowalla-1d": ("gowalla-1d/5.0m-gowalla.csv", "csv1d"),
    "spitz-1024x1024": ("spitz-1024x1024.csv", "json"),
    "cali-1024x1024": ("cali-1024x1024.pickle", "pickle"),
    "gowalla_2d_50k": ("gowalla_2d_50k.pkl", "pickle"),
    "gowalla_2d_100k": ("gowalla_2d_100k.pkl", "pickle"),
    "synthetic_2d_1m": ("synthetic_2d_1m.pkl", "pickle"),
    "synthetic_2d_1m_sparse": ("synthetic_2d_1m_sparse.pkl", "pickle"),
    "synthetic_2d_1m-1024x1024": ("synthetic_2d_1m-1024x1024.pkl", "pickle"),
    "nh_64": ("nh_64.txt", "json"),
    "gowalla_3d_23k": ("gowalla_3d_23k.pkl", "pickle"),
    "synthetic_3d_1m_128": ("synthetic_3d_1m_128.pkl", "pickle"),
    "synthetic_3d_1m_256": ("synthetic_3d_1m_256.pkl", "pickle"),
}

def load(path, kind):
    if kind == "pickle":
        with open(path, "rb") as f:
            data = pickle.load(f)
    elif kind == "json":
        with open(path) as f:
            data = json.load(f)
    elif kind == "csv1d":
        data = []
        with open(path) as f:
            for row in csv.reader(f):
                if not row:
                    continue
                try:
                    data.append((int(float(row[0])),))
                except ValueError:
                    continue  # header
    if isinstance(data, dict):
        pts = {tuple(int(c) for c in k): int(v) for k, v in data.items()}
    else:
        pts = {}
        for p in data:
            p = tuple(int(c) for c in (p if isinstance(p, (list, tuple)) else [p]))
            pts[p] = 1
    return pts

def main(src, dst):
    os.makedirs(dst, exist_ok=True)
    for name, (fname, kind) in SETS.items():
        path = os.path.join(src, fname)
        if not os.path.exists(path):
            print(f"skip {name}: {path} missing")
            continue
        pts = load(path, kind)
        d = len(next(iter(pts)))
        with open(os.path.join(dst, name + ".pts"), "w") as f:
            f.write(f"{d} {len(pts)}\n")
            for p in sorted(pts):
                f.write(" ".join(map(str, p)) + f" {pts[p]}\n")
        print(f"{name}: d={d} N={len(pts)}")

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
