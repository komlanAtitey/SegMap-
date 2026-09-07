"""Generic array/column helpers."""
from .common import *  # noqa: F401,F403


def pick_first_present(
    columns: List[str],
     candidates: List[str]) -> Optional[str]:
    return next((c for c in candidates if c in columns),
                None)  # return first candidate present


def normalize01(img: np.ndarray) -> np.ndarray:  # normalize image robustly
    img = img.astype(np.float32, copy=False)  # ensure float32
    p1 = np.percentile(img, 1)  # 1st percentile
    p99 = np.percentile(img, 99)  # 99th percentile
    if p99 <= p1:  # degenerate case
        return np.zeros_like(img, dtype=np.float32)  # return zeros
    out = (img - p1) / (p99 - p1)  # scale
    return np.clip(out, 0.0, 1.0)  # clip to [0,1]


def downsample_mean(img: np.ndarray, factor: int) -> np.ndarray:  # block mean downsample
    if factor <= 1:  # no-op
        return img  # return original
    h, w = img.shape  # height,width
    nh, nw = h // factor, w // factor  # new dims
    if nh < 1 or nw < 1:  # too small
        return img  # return original
    img_c = img[: nh * factor, : nw * factor]  # crop to divisible
    img_rs = img_c.reshape(nh, factor, nw, factor)  # block view
    return img_rs.mean(axis=(1, 3))  # mean pool
