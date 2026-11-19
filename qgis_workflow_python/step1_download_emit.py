"""
STEP 1 - Download the EMIT mineral map (EMIT L2B minerals) from NASA AppEEARS,
cropped to the box and orthorectified, in ENVI format (same as the manual Step 2).

Run:  python step1_download_emit.py
- Reads your Earthdata login from earthdata_login.txt next to this script
  (line 1 user name, line 2 password). Nothing is typed, nothing is printed.
- Looks up the exact EMIT mineral product name in AppEEARS's product list
  (prefers version 001, the version used in the paper).
- Submits one AppEEARS area request, waits until it is done (usually 10-60 min),
  then downloads the files into runs/<site>_<date>/emit/.
- Safe to stop and run again: it resumes the same request instead of submitting a new one.
- If the folder already holds EMIT "_min.img/_min.hdr" files (for example the ones you
  downloaded by hand), it does nothing. Delete them to force a new download.
"""
import os
import glob
import json
import time
from datetime import datetime

import requests

import config as C
from common import log, fail, ensure_dirs, earthdata_credentials

API = "https://appeears.earthdatacloud.nasa.gov/api"
TASK_FILE = lambda: os.path.join(C.DIR_EMIT, "appeears_task.json")


def aoi_geojson():
    w, s, e, n = C.BBOX
    return {"type": "FeatureCollection", "features": [{
        "type": "Feature", "properties": {"name": C.SITE},
        "geometry": {"type": "Polygon", "coordinates": [[[w, s], [e, s], [e, n], [w, n], [w, s]]]}}]}


def login():
    user, pw = earthdata_credentials()
    r = requests.post(f"{API}/login", auth=(user, pw), timeout=120)
    if r.status_code != 200:
        fail(f"AppEEARS login refused (HTTP {r.status_code}). Check earthdata_login.txt "
             f"and that you can sign in at https://appeears.earthdatacloud.nasa.gov")
    log("logged in to AppEEARS")
    return {"Authorization": "Bearer " + r.json()["token"]}


def find_emit_mineral_products():
    """Returns the AppEEARS product ids for EMIT L2B minerals, version 001 first."""
    r = requests.get(f"{API}/product", params={"limit": 2000}, timeout=120)
    r.raise_for_status()
    prods = r.json()
    hits = []
    for p in prods:
        pid = p.get("ProductAndVersion", "")
        text = (pid + " " + p.get("Product", "") + " " + p.get("Description", "")).upper()
        if "EMIT" in text and "MIN" in text and "L2B" in text.replace("_", "").replace(" ", ""):
            hits.append(pid)
    if not hits:
        emit = sorted(p.get("ProductAndVersion", "") for p in prods if "EMIT" in p.get("ProductAndVersion", "").upper())
        fail("AppEEARS lists no EMIT L2B mineral product. EMIT products it offers: " + (", ".join(emit) or "none"))
    hits.sort(key=lambda pid: (not pid.endswith(".001"), pid))
    log("EMIT mineral products in AppEEARS: " + ", ".join(hits))
    return hits


WANTED = ["group_1_band_depth", "group_1_mineral_id", "group_2_band_depth", "group_2_mineral_id"]


def product_layers(pid):
    """The four mineral layers the method needs. 'elev' and other layers AppEEARS lists but
    does not let you select are left out."""
    r = requests.get(f"{API}/product/{pid}", timeout=120)
    r.raise_for_status()
    info = r.json()
    log(f"{pid} layers listed: {', '.join(info)}")
    selectable = [k for k, v in info.items() if not isinstance(v, dict) or v.get("Available", True) is not False]
    layers = [k for k in WANTED if k in selectable]
    if len(layers) < 4:
        fail(f"{pid} does not offer all of {WANTED}; it offers {selectable}")
    return layers


def build_task():
    d = datetime.strptime(C.DATE, "%Y-%m-%d")
    dates = [{"startDate": d.strftime("%m-%d-%Y"), "endDate": d.strftime("%m-%d-%Y")}]
    name = f"{C.SITE}_emit_{C.DATE.replace('-', '')}"

    templates = [C.APPEEARS_TEMPLATE] + sorted(glob.glob(os.path.join(C.HERE, "*request*.json")))
    template = next((t for t in templates if os.path.exists(t)), None)
    if template:
        task = json.load(open(template, encoding="utf-8"))
        layers = [l for l in task.get("params", {}).get("layers", [])
                  if "EMIT" in l.get("product", "").upper() and l.get("layer") != "elev"]
        if layers:
            log(f"using the settings of {os.path.basename(template)} (product {layers[0]['product']})")
            task["task_name"] = name
            task["params"]["dates"] = dates
            task["params"]["geo"] = aoi_geojson()
            task["params"]["layers"] = layers
            task["params"].pop("additionalOptions", None)
            out = task["params"].setdefault("output", {"format": {"type": "envi"}, "projection": "geographic"})
            out.setdefault("additionalOptions", {})["orthorectify"] = True
            return task
        log(f"{os.path.basename(template)} has no EMIT layers - building the request from scratch")

    pids = find_emit_mineral_products()
    pid = pids[0]
    if not pid.endswith(".001"):
        log(f"WARNING: version 001 not offered; using {pid}. The paper used 001 (IDs 1-294 are identical; "
            f"the iron result can differ).")
    layers = product_layers(pid)
    log(f"requesting layers: {', '.join(layers)}")
    return {
        "task_type": "area",
        "task_name": name,
        "params": {
            "dates": dates,
            "layers": [{"product": pid, "layer": ln} for ln in layers],
            "geo": aoi_geojson(),
            "output": {"format": {"type": "envi"}, "projection": "geographic",
                       "additionalOptions": {"orthorectify": True}},
        },
    }


def task_status(task_id, hdr):
    r = requests.get(f"{API}/status/{task_id}", headers=hdr, timeout=120, allow_redirects=False)
    if r.status_code == 303:                         # AppEEARS answers 303 when the task is done
        return "done"
    if r.status_code == 404:
        return "unknown"
    r.raise_for_status()
    js = r.json()
    if isinstance(js, list):
        js = js[0] if js else {}
    return js.get("status", "processing"), js.get("progress")


def main():
    ensure_dirs()
    have = glob.glob(os.path.join(C.DIR_EMIT, "*_min.hdr"))
    if have:
        log(f"{len(have)} EMIT mineral file(s) already in {C.DIR_EMIT} - nothing to download.")
        log("STEP 1 DONE. Next: python step2_emit_abundance.py")
        return

    hdr = login()
    if os.path.exists(TASK_FILE()):
        task_id = json.load(open(TASK_FILE()))["task_id"]
        log(f"resuming AppEEARS request {task_id}")
    else:
        task = build_task()
        r = requests.post(f"{API}/task", json=task, headers=hdr, timeout=120)
        if r.status_code not in (200, 202):
            fail(f"AppEEARS refused the request (HTTP {r.status_code}): {r.text[:800]}")
        task_id = r.json()["task_id"]
        json.dump({"task_id": task_id, "submitted": time.strftime("%Y-%m-%d %H:%M"), "request": task},
                  open(TASK_FILE(), "w"), indent=1)
        log(f"request submitted: {task_id}. You can also follow it at https://appeears.earthdatacloud.nasa.gov/explore")

    t0 = time.time()
    while True:
        st = task_status(task_id, hdr)
        status, progress = (st, None) if isinstance(st, str) else st
        if status == "done":
            break
        if status == "error":
            fail(f"AppEEARS reports an error for request {task_id}. Open it in AppEEARS > Explore for details.")
        if status == "unknown":
            fail(f"AppEEARS does not know request {task_id}. Delete {TASK_FILE()} and run again.")
        if time.time() - t0 > C.MAX_WAIT_HOURS * 3600:
            fail(f"still '{status}' after {C.MAX_WAIT_HOURS} h. Run this step again later; it will resume.")
        log(f"AppEEARS status: {status}" + (f" ({progress})" if progress else "") + f", checking again in {C.POLL_SECONDS} s")
        time.sleep(C.POLL_SECONDS)

    r = requests.get(f"{API}/bundle/{task_id}", headers=hdr, timeout=120)
    r.raise_for_status()
    files = r.json()["files"]
    log(f"request done: {len(files)} files")
    for f in files:
        name = os.path.basename(f["file_name"])
        out = os.path.join(C.DIR_EMIT, name)
        if os.path.exists(out) and os.path.getsize(out) == f.get("file_size", -1):
            continue
        for attempt in range(3):
            try:
                with requests.get(f"{API}/bundle/{task_id}/{f['file_id']}", headers=hdr, stream=True,
                                  allow_redirects=True, timeout=600) as d:
                    d.raise_for_status()
                    with open(out + ".part", "wb") as fh:
                        for chunk in d.iter_content(1 << 20):
                            fh.write(chunk)
                os.replace(out + ".part", out)
                break
            except requests.RequestException as e:
                log(f"  retry {attempt + 1}/3 for {name}: {type(e).__name__}")
                time.sleep(15)
        else:
            fail(f"could not download {name}; run step 1 again (it resumes)")
        log(f"  downloaded {name}")
    if not glob.glob(os.path.join(C.DIR_EMIT, "*_min.hdr")):
        fail("no '_min' ENVI file in the download. Check the request in AppEEARS (format ENVI, Orthorectify).")
    log("STEP 1 DONE. Next: python step2_emit_abundance.py")


if __name__ == "__main__":
    main()
