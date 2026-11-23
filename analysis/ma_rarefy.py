"""Rarefy every MicrobeAtlas sample to 1,000 reads, and samples with >= 5,000 reads to 5,000 (protocol v0.5).

Saves compressed sparse rows (sample -> OTU indices and counts) so F1-F10 can reuse them:
  data/microbeatlas/derived/rare_1000.npz, rare_5000.npz
Also adds the body-site field to meta.pkl as `body_site`.
Run from sampling-imbalance-biogeography-rules:  .venv/bin/python analysis/ma_rarefy.py
"""
import multiprocessing as mp
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
D = ROOT / "data" / "microbeatlas"
SEED = 20261121
G = {}


def init():
    with h5py.File(D / "samples-otus.97.metag.minfilter.biom", "r") as f:
        G["indptr"] = f["sample/matrix/indptr"][:]
        G["indices"] = f["sample/matrix/indices"][:]
        G["data"] = f["sample/matrix/data"][:]


def work(args):
    lo, hi, depth = args
    rng = np.random.default_rng([SEED, depth, lo])
    ip, ix, dt = G["indptr"], G["indices"], G["data"]
    keep, rows_idx, rows_cnt, lens = [], [], [], []
    for i in range(lo, hi):
        c = dt[ip[i]:ip[i + 1]].astype(np.int64)
        if c.sum() < depth:
            continue
        r = rng.multivariate_hypergeometric(c, depth)
        nz = r > 0
        keep.append(i)
        rows_idx.append(ix[ip[i]:ip[i + 1]][nz].astype(np.int32))
        rows_cnt.append(r[nz].astype(np.int16))
        lens.append(int(nz.sum()))
    return lo, keep, rows_idx, rows_cnt, lens


def run(depth, n):
    chunks = [(lo, min(lo + 20000, n), depth) for lo in range(0, n, 20000)]
    with mp.get_context("fork").Pool(16) as pool:
        parts = sorted(pool.map(work, chunks), key=lambda p: p[0])
    keep = np.concatenate([np.array(p[1], dtype=np.int64) for p in parts])
    lens = np.concatenate([np.array(p[4], dtype=np.int64) for p in parts])
    idx = np.concatenate([a for p in parts for a in p[2]])
    cnt = np.concatenate([a for p in parts for a in p[3]])
    indptr = np.concatenate([[0], np.cumsum(lens)])
    np.savez(D / "derived" / f"rare_{depth}.npz", sample_row=keep, indptr=indptr, otu=idx, count=cnt)
    print(f"depth {depth}: {len(keep)} samples, {len(idx)} nonzeros", flush=True)


def main():
    init()
    n = len(G["indptr"]) - 1
    for depth in (1000, 5000):
        run(depth, n)
    with h5py.File(D / "samples-otus.97.metag.minfilter.biom", "r") as f:
        ids = np.array([x.decode() for x in f["sample/ids"][:]])
        otus = np.array([x.decode() for x in f["observation/ids"][:]])
    np.save(D / "derived" / "sample_ids.npy", ids)
    np.save(D / "derived" / "otu_ids.npy", otus)
    env = pd.read_csv(D / "samples.env.info.tsv", sep="\t", usecols=[0, 2, 4], dtype=str)
    env.columns = ["map_sid", "body_site", "keywords"]
    meta = pd.read_pickle(D / "derived" / "meta.pkl")
    meta = meta.drop(columns=[c for c in ["body_site", "keywords"] if c in meta]).merge(env, on="map_sid", how="left")
    meta.to_pickle(D / "derived" / "meta.pkl")
    print("meta updated with body_site", flush=True)


if __name__ == "__main__":
    main()
