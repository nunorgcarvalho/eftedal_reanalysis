"""Logo variants for a dark background: white page -> transparent; the wordmark (right of `split`)
recoloured white with alpha from its ink density; the emblem keeps its colours, with only the white
connected to the image border made transparent (so white details inside the emblem stay white).
Usage: python ASHG/poster/make_dark_logos.py  (writes logos/*_dark.png)."""
from pathlib import Path
import numpy as np
from PIL import Image
from scipy import ndimage

HERE = Path(__file__).parent / "logos"   # logos are untracked (institutional marks)
LOGOS = {"hms.png": 214, "dfci.png": 160}   # first all-white column between emblem and wordmark

for name, split in LOGOS.items():
    rgb = np.asarray(Image.open(HERE / name).convert("RGB")).astype(float)
    h, w, _ = rgb.shape
    ink = 255 - rgb.min(axis=2)                     # 0 on white, high on any ink (incl. red)
    out = np.zeros((h, w, 4))
    # wordmark: white text, alpha = ink density (normalised to the darkest text)
    t = ink[:, split:]
    out[:, split:, :3] = 255
    out[:, split:, 3] = np.clip(t / np.percentile(t[t > 30], 90), 0, 1) * 255
    # emblem: original colours; transparent where near-white and connected to the border
    e = rgb[:, :split]
    white = e.min(axis=2) > 225
    lab, _ = ndimage.label(white)
    border = set(np.unique(np.concatenate([lab[0], lab[-1], lab[:, 0], lab[:, -1]]))) - {0}
    bg = np.isin(lab, list(border))
    alpha = np.where(bg, 0, 255).astype(float)
    # soften the edge: partially-white pixels next to the background get alpha from their ink
    edge = ndimage.binary_dilation(bg, iterations=2) & ~bg
    alpha[edge] = np.clip(ink[:, :split][edge] / 120, 0, 1) * 255
    out[:, :split, :3] = e
    out[:, :split, 3] = alpha
    Image.fromarray(out.astype(np.uint8), "RGBA").save(HERE / name.replace(".png", "_dark.png"))
    print("wrote", name.replace(".png", "_dark.png"))
