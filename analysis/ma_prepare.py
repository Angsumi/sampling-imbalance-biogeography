"""N1 MicrobeAtlas extension, step 1: build one metadata table per sample (protocol v0.4).

Joins pH/temperature, coordinates, environment labels, technology and project, and assigns a
country from coordinates (Natural Earth admin-0, 1:110m). Output: data/microbeatlas/derived/meta.pkl
Run from sampling-imbalance-biogeography-rules:  .venv/bin/python analysis/ma_prepare.py
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
import shapely
from shapely.geometry import shape

ROOT = Path(__file__).resolve().parents[1]
D = ROOT / "data" / "microbeatlas"
OUT = D / "derived"
OUT.mkdir(exist_ok=True)


def countries(lat, lon):
    gj = json.loads((ROOT / "data" / "naturalearth" / "ne_110m_admin_0_countries.geojson").read_text())
    geoms = [shape(f["geometry"]) for f in gj["features"]]
    names = np.array([f["properties"]["ADMIN"] for f in gj["features"]])
    tree = shapely.STRtree(geoms)
    pts = shapely.points(lon, lat)
    pi, gi = tree.query(pts, predicate="within")
    out = np.full(len(pts), "ocean", dtype=object)
    out[pi] = names[gi]
    return out


def main():
    env = pd.read_csv(D / "samples.env.info.tsv", sep="\t", usecols=[0, 1, 3, 6], dtype=str)
    env.columns = ["map_sid", "environments", "technology", "project"]
    env["sample"] = env.map_sid.str.split(".").str[-1]
    first = env.environments.fillna("").str.split("|").str[0]
    env["env_top"] = first.str.split(";").str[0].replace("", "unknown")
    env["env_2"] = first.replace("", "unknown")
    env["host_associated"] = env.env_top.isin(["animal", "plant"])

    m = pd.read_csv(D / "samples.measurements.tsv", sep="\t", usecols=[0, 1, 2])
    m = m[m.Meta_variable.isin(["ph", "temp_C"])]
    m["Value"] = pd.to_numeric(m.Value, errors="coerce")
    meas = m.pivot_table(index="SampleID", columns="Meta_variable", values="Value", aggfunc="first")
    meas = meas.rename(columns={"ph": "ph", "temp_C": "temp"})

    ll = pd.read_csv(D / "samples.info.latlon.parsed.tsv", sep="\t", usecols=[0, 1, 2],
                     na_values=["None"]).dropna()
    ll.columns = ["sample", "lat", "lon"]
    ll = ll[ll.lat.between(-90, 90) & ll.lon.between(-180, 180)].drop_duplicates("sample")
    ll["country"] = countries(ll.lat.to_numpy(), ll.lon.to_numpy())

    d = env.merge(meas, left_on="sample", right_index=True, how="left")
    d = d.merge(ll, on="sample", how="left")
    d["country"] = d["country"].fillna("no_coordinates")
    lat_c = np.floor(d.lat / 5).astype("Int64").astype(str)
    lon_c = np.floor(d.lon / 5).astype("Int64").astype(str)
    d["cell"] = (lat_c + "_" + lon_c).where(d.lat.notna(), "no_coordinates")
    d.to_pickle(OUT / "meta.pkl")
    print(len(d), "runs;", d["sample"].nunique(), "samples;", d.project.nunique(), "projects")
    print("with pH:", d.ph.notna().sum(), "| with temp:", d.temp.notna().sum(), "| with coords:", d.lat.notna().sum())
    vc = d[d.lat.notna()].country.value_counts()
    print("countries:", len(vc), "| top 10:", (vc.head(10) / vc.sum()).round(3).to_dict())


if __name__ == "__main__":
    main()
