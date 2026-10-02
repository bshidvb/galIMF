import re
from collections import OrderedDict

INPUT  = "yield_tables/agb_and_massive_stars_K10C20_lc18_r0_M.txt"
OUTPUT = "yield_tables/lc18_r0_M_recalculated.txt"
TARGET_Z = [0.02, 0.008, 0.004, 0.0001]
MASS_MIN = 10.0

def parse(filename):
    raw = open(filename).readlines()
    header, blocks = [], {}
    i = 0
    while i < len(raw) and not re.match(r'H Table:', raw[i].strip()):
        header.append(raw[i]); i += 1
    while i < len(raw):
        m = re.match(r'H Table:\s*\(M=([\d.]+),Z=([\deE.+-]+)\)', raw[i].strip())
        if not m: i += 1; continue
        M, Z = float(m.group(1)), float(m.group(2))
        lt = mf = ch = None; iso = OrderedDict(); i += 1
        while i < len(raw) and not re.match(r'H Table:', raw[i].strip()):
            l = raw[i]
            if 'Lifetime' in l: lt = float(l.split(':')[1])
            elif 'Mfinal'  in l: mf = float(l.split(':')[1])
            elif '&Isotopes' in l: ch = l
            elif l.strip().startswith('&') and ch:
                p = [x.strip() for x in l.split('&') if x.strip()]
                iso[p[0]] = {'y': float(p[1]), 'zn': p[-2], 'a': p[-1]}
            i += 1
        blocks[(M, Z)] = {'lt': lt, 'mf': mf, 'ch': ch, 'iso': iso}
    return header, blocks

def interp(blocks, M, Z_target):
    avail = sorted(Z for (Mm, Z) in blocks if abs(Mm-M) < 1e-9)
    if any(abs(Z_target-Z) < 1e-9 for Z in avail):
        return blocks[(M, next(Z for Z in avail if abs(Z-Z_target) < 1e-9))]
    Z_lo = max((Z for Z in avail if Z < Z_target), default=avail[0])
    Z_hi = min((Z for Z in avail if Z > Z_target), default=avail[-1])
    w = 0.0 if Z_lo == Z_hi else (Z_target - Z_lo) / (Z_hi - Z_lo)
    lo, hi = blocks[(M, Z_lo)], blocks[(M, Z_hi)]
    lp = lambda a, b: a + w * (b - a)
    new_iso = OrderedDict()
    for k in lo['iso']:
        if k in hi['iso']:
            new_iso[k] = {'y': lp(lo['iso'][k]['y'], hi['iso'][k]['y']),
                          'zn': lo['iso'][k]['zn'], 'a': lo['iso'][k]['a']}
    return {'lt': lp(lo['lt'], hi['lt']), 'mf': lp(lo['mf'], hi['mf']),
            'ch': lo['ch'], 'iso': new_iso}

def fmt_Z(Z):
    for val, s in [(0.0001,'0.0001'),(0.004,'0.004'),(0.008,'0.008'),(0.02,'0.02')]:
        if abs(Z-val) < 1e-9: return s
    return f"{Z:.4f}"

header, blocks = parse(INPUT)
masses = sorted({M for (M,Z) in blocks if M >= MASS_MIN})

with open(OUTPUT, 'w') as out:
    for line in header:
        out.write(f"H Number of metallicities: {len(TARGET_Z)}\n" if line.startswith('H Number') else line)
    for Z in TARGET_Z:
        for M in masses:
            b = interp(blocks, M, Z)
            M_str = str(M) if M != int(M) else f"{int(M)}.0"
            out.write(f"H Table: (M={M_str},Z={fmt_Z(Z)})\n")
            out.write(f"H Lifetime: {b['lt']:.3e}\n")
            out.write(f"H Mfinal: {b['mf']:.3e}\n")
            out.write(b['ch'])
            for iso, d in b['iso'].items():
                out.write(f"&{iso:<8}  &{d['y']:.3e}  &{d['zn']}  &{d['a']}\n")

print(f"Done → {OUTPUT}")