"""Step 1: download the public input data (no account or API key needed).

Files go to the folder named by the RAW_DIR environment variable
(default: ./data_raw, which is git-ignored). Existing files are skipped.

Sources
- EPA TRI Basic Data Files, 2022-2024, Texas and Louisiana
- Census TIGER/Line 2022 census tracts, Texas (48) and Louisiana (22)
- CDC/ATSDR Social Vulnerability Index 2022, state-level CSVs
- EPA AQS pre-generated annual concentration by monitor files, 2022-2024
- EPA AirToxScreen 2020 block-level cancer risk by source group, Region 6
  (optional, about 240 MB; skip with SKIP_AIRTOX=1)
"""
import os
import sys
import time
from pathlib import Path

import requests

RAW = Path(os.environ.get("RAW_DIR", "data_raw")).resolve()
YEARS = [2022, 2023, 2024]

FILES = {}
for y in YEARS:
    for st in ("TX", "LA"):
        FILES[f"tri/tri_{y}_{st}.csv"] = (
            f"https://data.epa.gov/efservice/downloads/tri/mv_tri_basic_download/{y}_{st}/csv")
    FILES[f"aqs/annual_conc_by_monitor_{y}.zip"] = (
        f"https://aqs.epa.gov/aqsweb/airdata/annual_conc_by_monitor_{y}.zip")
for fips in ("48", "22"):
    FILES[f"tiger/tl_2022_{fips}_tract.zip"] = (
        f"https://www2.census.gov/geo/tiger/TIGER2022/TRACT/tl_2022_{fips}_tract.zip")
FILES["svi/SVI_2022_Texas.csv"] = "https://svi.cdc.gov/Documents/Data/2022/csv/states/Texas.csv"
FILES["svi/SVI_2022_Louisiana.csv"] = "https://svi.cdc.gov/Documents/Data/2022/csv/states/Louisiana.csv"
if os.environ.get("SKIP_AIRTOX") != "1":
    base = "https://gaftp.epa.gov/rtrmodeling_public/AirToxScreen/2020/Cancer/BySource/"
    for part in ("6a", "6b"):
        FILES[f"airtox/Region{part}_CancerRisk_by_block_srcgrp.xlsx"] = (
            f"{base}Region{part}_CancerRisk_by_block_srcgrp.xlsx")


def fetch(url, dest, tries=20):
    """Download with resume: the EPA gaftp server is slow and drops long transfers."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    for k in range(tries):
        try:
            have = tmp.stat().st_size if tmp.exists() else 0
            headers = {"User-Agent": "research-download/1.0"}
            if have:
                headers["Range"] = f"bytes={have}-"
            with requests.get(url, stream=True, timeout=300, headers=headers) as r:
                if r.status_code == 416:  # already complete
                    tmp.rename(dest)
                    return True
                r.raise_for_status()
                mode = "ab" if (have and r.status_code == 206) else "wb"
                with open(tmp, mode) as f:
                    for chunk in r.iter_content(1 << 20):
                        f.write(chunk)
            tmp.rename(dest)
            return True
        except Exception as e:  # network hiccup, resume on next attempt
            print(f"  attempt {k + 1} interrupted: {e}", flush=True)
            time.sleep(5)
    return False


def main():
    log = []
    for rel, url in FILES.items():
        dest = RAW / rel
        if dest.exists() and dest.stat().st_size > 0:
            print(f"skip {rel}")
        else:
            print(f"get  {rel}")
            if not fetch(url, dest):
                print(f"FAILED {url}", file=sys.stderr)
                continue
        acc = time.strftime("%Y-%m-%d", time.localtime(dest.stat().st_mtime))
        log.append(f"{rel},{url},{dest.stat().st_size},{acc}")
    text = "file,url,bytes,accessed\n" + "\n".join(log) + "\n"
    (RAW / "download_log.csv").write_text(text)
    data = Path(__file__).resolve().parents[1] / "data"
    data.mkdir(exist_ok=True)
    (data / "download_log.csv").write_text(text)
    print(f"done, {len(log)} files in {RAW}")


if __name__ == "__main__":
    main()
