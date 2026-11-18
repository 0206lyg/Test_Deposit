"""Conservative parabolic inflow and matching exponential-PTT simple-shear data."""
import json
import math
from pathlib import Path
import re


def shear_factor(mode, gamma):
    """Positive simple-shear root f=exp(2*epsilon*(lambda*gamma/f)^2)."""
    target = 2 * mode["epsilon"] * (mode["lambda_s"] * gamma)**2
    lower, upper = 0., math.log1p(target)
    for _ in range(60):
        q = (lower + upper) / 2
        if q * math.exp(2 * q) < target:
            lower = q
        else:
            upper = q
    return math.exp((lower + upper) / 2)


def shear_state(mode, gamma):
    """C=[[1+2b²,b,0],[b,1,0],[0,0,1]], with b=lambda*gamma/f."""
    lam, eta = mode["lambda_s"], mode["eta_p_Pa_s"]
    b = lam * gamma / shear_factor(mode, gamma)
    # Orthonormal eigenbasis of the xy block of C. Obtain the small
    # eigenvalue from the determinant to avoid subtractive cancellation.
    angle = 0.5 * math.atan2(b, b*b)
    c, s = math.cos(angle), math.sin(angle)
    high = 1 + b*b + abs(b) * math.hypot(1, b)
    low = (1 + b*b) / high
    lh, ll = math.log(high), math.log(low)
    return {
        "tau": (2*b*b*eta/lam, b*eta/lam, 0, 0, 0, 0),
        "theta": (c*c*lh+s*s*ll, c*s*(lh-ll), 0, s*s*lh+c*c*ll, 0, 0),
        "eigVals": (high, 0, 0, 0, low, 0, 0, 0, 1),
        "eigVecs": (c, -s, 0, s, c, 0, 0, 0, 1),
    }


def foam_body(path):
    text = path.read_text()
    text = re.sub(r"/\*.*?\*/|//[^\n]*", "", text, flags=re.S)
    header = re.match(r"\s*FoamFile\s*\{([^}]*)\}", text)
    if not header or not re.search(r"\bformat\s+ascii\s*;", header[1]):
        raise ValueError(f"Expected an ASCII OpenFOAM mesh file at {path}")
    return text[header.end():].strip()


def mesh_list(path):
    body = foam_body(path)
    match = re.fullmatch(r"(\d+)\s*\((.*)\)\s*;?", body, flags=re.S)
    if not match:
        raise ValueError(f"Cannot read OpenFOAM list {path}")
    return int(match[1]), match[2]


def inlet_intervals(mesh):
    """Use actual patch face ordering after blockMesh, before decomposition."""
    boundary = foam_body(mesh / "boundary")
    patch = re.search(r"\binlet\s*\{([^}]*)\}", boundary)
    if not patch:
        raise ValueError("inlet patch missing")
    first = int(re.search(r"\bstartFace\s+(\d+);", patch[1])[1])
    count = int(re.search(r"\bnFaces\s+(\d+);", patch[1])[1])
    npoints, data = mesh_list(mesh / "points")
    pts = [tuple(map(float, p.split())) for p in re.findall(r"\(([^()]*)\)", data)]
    if len(pts) != npoints:
        raise ValueError("points count mismatch")
    nfaces, data = mesh_list(mesh / "faces")
    faces = re.findall(r"\d*\s*\(([^()]*)\)", data)
    if len(faces) != nfaces or first + count > nfaces or count == 0:
        raise ValueError("inlet face addressing or faces count invalid")
    intervals = []
    for face in faces[first:first+count]:
        ys = [pts[int(i)][1] for i in face.split()]
        a, b = min(ys), max(ys)
        if b <= a:
            raise ValueError("Inlet face must span a finite y interval")
        intervals.append((a, b))
    return intervals


def initialize(case):
    data = json.loads((case/"metadata.json").read_text())
    half_width = data["geometry_m"]["upstream_width"]/2
    mean = data["inlet_mean_velocity_m_s"]
    rows = inlet_intervals(case/"constant"/"polyMesh")
    fields = {"U": ("vector", [])}
    modes = data["fluid"].get("modes", [])
    for mode in modes:
        for key in ("tau", "theta", "eigVals", "eigVecs"):
            kind = "tensor" if key.startswith("eig") else "symmTensor"
            fields[key + mode["name"]] = (kind, [])
    flow = 0.
    for a, b in rows:
        # Exact face average of u=1.5*mean*(1-(y/half_width)^2).
        u = 1.5*mean*(1 - (a*a+a*b+b*b)/(3*half_width**2))
        fields["U"][1].append((u, 0, 0))
        flow += u*(b-a)
        gamma = -3*mean*(a+b)/(2*half_width**2)
        for mode in modes:
            for name, value in shear_state(mode, gamma).items():
                fields[name+mode["name"]][1].append(value)
    if not math.isclose(flow, data["q_half_per_depth_m2_s"], rel_tol=1e-8):
        raise ValueError(f"Discrete inlet flux {flow:g} differs from Q/(2h)")
    for name, (kind, values) in fields.items():
        contents = "\n".join("(" + " ".join(format(v, ".16g") for v in row) + ")" for row in values)
        p = case/"constant"/("inlet_"+name)
        p.write_text(f"value nonuniform List<{kind}>\n{len(values)}\n(\n{contents}\n);\n")
    print(f"Initialized {len(rows)} inlet faces at Q/(2h) = {flow:.12g} m2/s")
