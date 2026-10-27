"""Deterministic file writers so a rebuild is byte-identical to the committed files."""
import json
import os


def write_json(path, obj, indent=1):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    text = json.dumps(obj, indent=indent, sort_keys=True, ensure_ascii=True) + "\n"
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(text)


def write_text(path, text):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(text)


def write_bytes(path, raw):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "wb") as fh:
        fh.write(raw)


def write_cases(out_dir, group, cases):
    """cases: list of dicts with name, description, expected ('PASS' or 'FAIL'), and either obj or raw.

    Writes out_dir/<group>/<name>.json and out_dir/<group>/MANIFEST.json.
    """
    manifest = {}
    for case in cases:
        path = os.path.join(out_dir, group, case["name"] + ".json")
        if "raw" in case:
            write_bytes(path, case["raw"])
        else:
            write_json(path, case["obj"])
        manifest[case["name"]] = {"description": case["description"], "expected": case["expected"]}
    write_json(os.path.join(out_dir, group, "MANIFEST.json"), manifest)
