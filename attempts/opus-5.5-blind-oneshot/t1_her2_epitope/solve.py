import gzip, json, shlex
import numpy as np

path = "inputs/1N8Z.cif.gz"
with gzip.open(path, "rt") as f:
    lines = f.read().splitlines()

# locate _atom_site loop
cols = []
i = 0
n = len(lines)
start = None
while i < n:
    if lines[i].strip() == "loop_":
        j = i + 1
        if j < n and lines[j].strip().startswith("_atom_site."):
            while j < n and lines[j].strip().startswith("_atom_site."):
                cols.append(lines[j].strip().split()[0][len("_atom_site."):])
                j += 1
            start = j
            break
    i += 1
assert start is not None
ci = {c: k for k, c in enumerate(cols)}

def toks(line):
    if '"' in line or "'" in line:
        return shlex.split(line, posix=True)
    return line.split()

atoms = []
k = start
buf = []
while k < n:
    s = lines[k]
    st = s.strip()
    if st == "" :
        k += 1
        continue
    if st.startswith("#") or st.startswith("loop_") or st.startswith("_"):
        break
    buf.extend(toks(s))
    while len(buf) >= len(cols):
        atoms.append(buf[:len(cols)])
        buf = buf[len(cols):]
    k += 1

def g(a, name):
    return a[ci[name]]

fab = []
fabkey = []
her = []
herkey = []
for a in atoms:
    if "pdbx_PDB_model_num" in ci and g(a, "pdbx_PDB_model_num") != "1":
        continue
    if g(a, "label_entity_id") not in ("1", "2", "3"):
        continue
    if g(a, "type_symbol").upper() in ("H", "D"):
        continue
    if g(a, "label_alt_id") not in (".", "A"):
        continue
    ch = g(a, "auth_asym_id")
    ic = g(a, "pdbx_PDB_ins_code")
    if ic in ("?", "."):
        ic = ""
    resnum = int(g(a, "auth_seq_id"))
    resname = g(a, "auth_comp_id")
    xyz = (float(g(a, "Cartn_x")), float(g(a, "Cartn_y")), float(g(a, "Cartn_z")))
    if ch in ("A", "B"):
        fab.append(xyz); fabkey.append((ch, resnum, ic))
    elif ch == "C":
        her.append(xyz); herkey.append((ch, resnum, ic, resname))

fab = np.array(fab, dtype=np.float64)
her = np.array(her, dtype=np.float64)

# group HER2 atoms by residue
res_atoms = {}
res_order = []
for idx, key in enumerate(herkey):
    rk = key[:3]
    if rk not in res_atoms:
        res_atoms[rk] = {"name": key[3], "idx": []}
        res_order.append(rk)
    res_atoms[rk]["idx"].append(idx)

out = []
for rk in res_order:
    idx = res_atoms[rk]["idx"]
    d = np.sqrt(((her[idx][:, None, :] - fab[None, :, :]) ** 2).sum(-1))  # (na, nfab)
    dmin_per_fab = d.min(axis=0)
    mind = float(dmin_per_fab.min())
    if mind > 4.0:
        continue
    contact = sorted({fabkey[j][0] for j in np.nonzero(dmin_per_fab <= 4.0)[0]})
    cands = [fabkey[j] for j in np.nonzero(dmin_per_fab == dmin_per_fab.min())[0]]
    best = sorted(set(cands))[0]
    out.append({
        "chain": rk[0], "resnum": int(rk[1]), "icode": rk[2], "resname": res_atoms[rk]["name"],
        "contacted_by": contact,
        "closest": {"chain": best[0], "resnum": int(best[1]), "icode": best[2]},
        "min_distance": round(mind, 4),
    })

with open("answer.json", "w") as f:
    json.dump({"epitope": out}, f, indent=1)
print(len(out), "epitope residues")
