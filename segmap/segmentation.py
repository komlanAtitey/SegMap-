"""Tiled Cellpose segmentation: runner, multiprocessing workers, tiling, mask export."""
from .common import *  # noqa: F401,F403
from .utils import downsample_mean, normalize01


class CellposeRunner:

    def __init__(
        self,
        gpu: bool,
        model_type: str = "nuclei",
        diameter: float = 15,
        flow_threshold: float = 0.4,
        cellprob_threshold: float = 0.2,
        min_size: int = 15,
    ):

        if not CELLPOSE_AVAILABLE:
            raise RuntimeError(
                "Cellpose not installed. "
                "Install using: pip install cellpose"
            )

        self.gpu = bool(gpu)

        self.diameter = diameter
        self.flow_threshold = flow_threshold
        self.cellprob_threshold = cellprob_threshold
        self.min_size = int(min_size)

        # Build the model, actually honoring model_type (e.g. "nuclei" for DAPI).
        # Previously model_type was ignored, so the generic cyto model was used on
        # nuclei images, which misses dim/small nuclei. Try the requested built-in
        # model first, then fall back gracefully across Cellpose versions.
        mt = (model_type or "").strip()
        self.model = None
        ctors = []
        if mt:
            ctors += [
                lambda: models.CellposeModel(gpu=self.gpu, model_type=mt),
                lambda: models.CellposeModel(gpu=self.gpu, pretrained_model=mt),
                lambda: models.Cellpose(gpu=self.gpu, model_type=mt),
            ]
        ctors += [
            lambda: models.CellposeModel(gpu=self.gpu),
            lambda: models.Cellpose(gpu=self.gpu),
        ]
        last_err = None
        for ctor in ctors:
            try:
                self.model = ctor()
                break
            except Exception as e:  # try next constructor variant
                last_err = e
        if self.model is None:
            raise RuntimeError(f"Could not initialize Cellpose model: {last_err}")
        self.model_type = mt or "default"
        print(f"[CELLPOSE] model='{self.model_type}' diameter={self.diameter} "
              f"flow={self.flow_threshold} cellprob={self.cellprob_threshold} "
              f"min_size={self.min_size}", flush=True)

    def eval(self, img2d: np.ndarray) -> np.ndarray:

        img2d = img2d.astype(np.float32, copy=False)

        with warnings.catch_warnings():

            warnings.simplefilter("ignore")

            kw = dict(diameter=self.diameter,
                      flow_threshold=float(self.flow_threshold),
                      cellprob_threshold=float(self.cellprob_threshold),
                      min_size=int(self.min_size))
            try:
                masks, *_ = self.model.eval(img2d, **kw)
            except TypeError:
                try:
                    masks, *_ = self.model.eval(img2d, channels=[0, 0], **kw)
                except TypeError:
                    # very old/!min_size signature: drop min_size
                    kw.pop("min_size", None)
                    masks, *_ = self.model.eval(img2d, channels=[0, 0], **kw)

        return masks.astype(np.int32, copy=False)

    def eval_batch(self, imgs: List[np.ndarray]) -> List[np.ndarray]:
        """Evaluate a list of 2D images in one call.

        Cellpose batches a list of images internally, which on GPU is markedly
        faster than one model.eval() call per tile (it amortizes host<->device
        transfer and kernel-launch overhead). Returns a list of int32 masks in
        the same order as the input.
        """
        if len(imgs) == 0:
            return []
        imgs = [im.astype(np.float32, copy=False) for im in imgs]
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            kw = dict(diameter=self.diameter,
                      flow_threshold=float(self.flow_threshold),
                      cellprob_threshold=float(self.cellprob_threshold),
                      min_size=int(self.min_size))
            try:
                out = self.model.eval(imgs, **kw)
            except TypeError:
                try:
                    out = self.model.eval(imgs, channels=[0, 0], **kw)
                except TypeError:
                    kw.pop("min_size", None)
                    out = self.model.eval(imgs, channels=[0, 0], **kw)
        masks_list = out[0]  # (masks, flows, styles[, diams]); masks is a list for list input
        return [np.asarray(m).astype(np.int32, copy=False) for m in masks_list]


_WORKER_RUNNER: Optional[CellposeRunner] = None  # worker runner


_WORKER_DAPI_MEMMAP_PATH: Optional[str] = None  # memmap path


_WORKER_DAPI_SHAPE: Optional[Tuple[int, int]] = None  # memmap shape


def _worker_init(dapi_npy_path: str,
    shape: Tuple[int, int],
    gpu: bool,
    model_type: str,
    diameter: float = 15,
    flow_threshold: float = 0.4,
    cellprob_threshold: float = 0.2,
    min_size: int = 15):
    global _WORKER_RUNNER, _WORKER_DAPI_MEMMAP_PATH, _WORKER_DAPI_SHAPE  # globals
    _WORKER_RUNNER = CellposeRunner(  # build worker runner with the SAME params as the main run
        gpu=gpu,
        model_type=model_type,
        diameter=diameter,
        flow_threshold=flow_threshold,
        cellprob_threshold=cellprob_threshold,
        min_size=min_size,
    )
    _WORKER_DAPI_MEMMAP_PATH = dapi_npy_path  # store memmap path
    _WORKER_DAPI_SHAPE = shape  # store shape


def _worker_segment_tile(args_tuple):  # worker callable
    global _WORKER_RUNNER, _WORKER_DAPI_MEMMAP_PATH, _WORKER_DAPI_SHAPE  # globals
    (y0, x0, th, tw, tile_idx, tile_std, flow_threshold, cellprob_threshold,
     skip_threshold, overlap_half) = args_tuple  # unpack
    try:  # open memmap
        dapi_mm = np.load(
    _WORKER_DAPI_MEMMAP_PATH,
     mmap_mode='r')  # memmap load
    except Exception as e:  # error
        raise RuntimeError(
    # raise
    f"Worker cannot open memmap DAPI file: {_WORKER_DAPI_MEMMAP_PATH}: {e}")
    tile = np.asarray(dapi_mm[y0:y0 + th, x0:x0 + tw],
                      dtype=np.float32)  # read tile
    tile_norm = normalize01(tile)  # normalize
    if float(tile_norm.std()) < float(skip_threshold):  # skip if low std
        return {
    'tile_idx': tile_idx,
    'y0': y0,
    'x0': x0,
    'th': th,
    'tw': tw,
    'masks': None,
     'tile_std': float(tile_std)}  # return empty

    try:

        masks = _WORKER_RUNNER.eval(tile_norm)

    except Exception as e:

        print(
            f"[WORKER] Segmentation failed "
            f"(tile={tile_idx}, y0={y0}, x0={x0}): {e}",
            flush=True
        )

        return {
            'tile_idx': tile_idx,
            'y0': y0,
            'x0': x0,
            'th': th,
            'tw': tw,
            'masks': None,
            'tile_std': float(tile_std)
        }
    tmp = relabel_sequential(masks)
    masks0 = tmp[0].astype(np.int32, copy=False)

    return {
        'tile_idx': tile_idx,
        'y0': y0,
        'x0': x0,
        'th': th,
        'tw': tw,
        'masks': masks0,
        'tile_std': float(tile_std)
    }


def tile_cellpose_fullres_upgraded(
    dapi_full: np.ndarray,
    tile_size: int = 2048,
    overlap: int = 128,
    gpu: bool = False,
    flow_threshold: float = 0.4,
    cellprob_threshold: float = 0.2,
    cellpose_diameter: float = None,
    cpu_workers: int = 1,  # cpu workers
    skip_threshold: float = 0.01,  # skip threshold
    force_cpu_on_gpu_error: bool = False,  # fallback
    prioritize: bool = True,  # prioritize high-signal tiles
    downsample_preview_for_scan: int = 8,  # preview downsample
    tile_batch: int = 8,  # tile batch size (CPU sequential only)
    model_type: str = "nuclei",  # cellpose model (nuclei for DAPI)
    min_size: int = 15,  # min object size in px (lower -> keep smaller nuclei)
) -> np.ndarray:  # returns mask
    H, W = dapi_full.shape  # full dims
    tile = int(tile_size)  # tile int
    ov = int(overlap)  # overlap int
    step = int(max(1, tile - ov))  # step size
    xs = list(range(0, W, step))  # x starts
    ys = list(range(0, H, step))  # y starts
    jobs: List[Tuple[int, int]] = [(y0, x0) for y0 in ys for x0 in xs]  # job list
    print(f"[TILE-UPG] H={H} W={W} tile={tile} overlap={ov} step={step} n_tiles={len(jobs)}", flush=True)  # log

    preview = downsample_mean(dapi_full, downsample_preview_for_scan)  # preview image
    tile_info = []  # collect info
    t_idx = 0  # index
    for y0 in range(0, H, step):  # iterate grid
        for x0 in range(0, W, step):  # iterate grid
            th = int(min(tile, H - y0))  # tile height
            tw = int(min(tile, W - x0))  # tile width
            py0 = y0 // downsample_preview_for_scan  # preview y
            px0 = x0 // downsample_preview_for_scan  # preview x
            pth = max(1, th // downsample_preview_for_scan)  # preview th
            ptw = max(1, tw // downsample_preview_for_scan)  # preview tw
            win = preview[py0:py0 + pth, px0:px0 + ptw]  # preview window
            stdv = float(np.std(win)) if win.size else 0.0  # std
            tile_info.append((t_idx, y0, x0, th, tw, stdv))  # append info
            t_idx += 1  # increment

    if prioritize:  # sort by std desc
        tile_info = sorted(tile_info, key=lambda t: t[-1], reverse=True)  # sort

    tmpdir = tempfile.mkdtemp(prefix="dapi_memmap_")  # tmp dir
    dapi_npy = os.path.join(tmpdir, "dapi_memmap.npy")  # memmap path
    np.save(dapi_npy, dapi_full)  # save memmap

    global_mask = np.zeros((H, W), dtype=np.int32)  # global mask
    label_offset = 0  # label offset
    core_pad = ov // 2  # core padding

    def _core_slices(y0: int, x0: int, th: int, tw: int):  # compute slices
        y1 = min(H, y0 + th)  # y end
        x1 = min(W, x0 + tw)  # x end
        cy0 = y0 + (core_pad if y0 > 0 else 0)  # core y start
        cx0 = x0 + (core_pad if x0 > 0 else 0)  # core x start
        cy1 = y1 - (core_pad if y1 < H else 0)  # core y end
        cx1 = x1 - (core_pad if x1 < W else 0)  # core x end
        ly0 = (core_pad if y0 > 0 else 0)  # local y start
        lx0 = (core_pad if x0 > 0 else 0)  # local x start
        ly1 = ly0 + (cy1 - cy0)  # local y end
        lx1 = lx0 + (cx1 - cx0)  # local x end
        return slice(cy0, cy1), slice(cx0, cx1), slice(ly0, ly1), slice(lx0, lx1)  # return slices

    def _merge_core(masks, y0: int, x0: int, th: int, tw: int):  # shared merge of one tile into global mask
        nonlocal label_offset  # we update the running label offset
        if masks is None:  # nothing to merge
            return
        if int(masks.max()) == 0:  # empty tile
            return
        masks = relabel_sequential(masks)[0].astype(np.int32, copy=True)  # compact labels
        masks[masks > 0] += int(label_offset)  # offset to keep labels globally unique
        gy, gx, ly, lx = _core_slices(y0, x0, th, tw)  # core slices
        core = masks[ly, lx]  # tile core region
        gview = global_mask[gy, gx]  # corresponding global view
        core_labels = np.unique(core)  # tile labels
        core_labels = core_labels[core_labels > 0]  # drop background
        mapping = {}  # tile-label -> global-label (or None for new)
        for cl in core_labels:  # resolve overlaps with already-written cells
            mask_cl = (core == cl)  # this object's pixels
            overlap_vals, counts = np.unique(gview[mask_cl], return_counts=True)  # what's underneath
            nonzero_mask = overlap_vals != 0  # ignore background overlap
            if np.any(nonzero_mask):  # overlaps an existing cell
                candidates = overlap_vals[nonzero_mask]  # candidate global ids
                ccounts = counts[nonzero_mask]  # overlap sizes
                chosen = candidates[np.argmax(ccounts)]  # majority-overlap id
                mapping[int(cl)] = int(chosen)  # merge into it
            else:  # no overlap -> brand new cell
                mapping[int(cl)] = None  # mark new
        for cl, mapped in mapping.items():  # write results
            if mapped is None:  # new label only into empty pixels
                write_mask = (gview == 0) & (core == cl)  # avoid clobbering
                gview[write_mask] = cl  # write new id
            else:  # merge into existing id
                gview[core == cl] = mapped  # relabel
        global_mask[gy, gx] = gview  # commit back
        label_offset = int(global_mask.max())  # advance offset

    use_parallel = (not gpu) and (cpu_workers is not None) and int(cpu_workers) > 1  # decide parallel

    if use_parallel:  # multiprocessing path
        import multiprocessing as mp  # import mp
        print(f"[TILE-UPG] Using multiprocessing with {cpu_workers} workers", flush=True)  # log
        tasks = []  # tasks list
        for (t_idx, y0, x0, th, tw, stdv) in tile_info:  # build tasks
            tasks.append((y0, x0, th, tw, t_idx, stdv, flow_threshold, cellprob_threshold, skip_threshold, core_pad))  # append
        ctx = mp.get_context('spawn')  # spawn ctx
        pool = ctx.Pool(  # workers now receive the user's diameter/flow/cellprob
            processes=int(cpu_workers),
            initializer=_worker_init,
            initargs=(dapi_npy, (H, W), False, model_type,
                      cellpose_diameter, flow_threshold, cellprob_threshold, min_size),
        )  # pool
        try:  # iterate results
            for res in pool.imap_unordered(_worker_segment_tile, tasks, chunksize=1):  # unordered
                if res is None:  # skip
                    continue  # continue
                if res['masks'] is None:  # skip
                    continue  # continue
                y0 = res['y0']; x0 = res['x0']; th = res['th']; tw = res['tw']  # unpack
                _merge_core(res['masks'], y0, x0, th, tw)  # merge tile into global mask
        finally:  # cleanup
            pool.close(); pool.join()  # close pool
            try:  # remove memmap
                os.remove(dapi_npy)  # remove
            except Exception:  # ignore
                pass  # ignore
            shutil.rmtree(tmpdir, ignore_errors=True)  # remove tmpdir

    else:  # sequential path (GPU or single CPU)
        if gpu:  # prepare GPU runner if requested
            runner_gpu = CellposeRunner(
                gpu=True,
                model_type=model_type,
                diameter=cellpose_diameter,
                flow_threshold=flow_threshold,
                cellprob_threshold=cellprob_threshold,
                min_size=min_size,
            )
            runner_cpu = CellposeRunner(  # CPU fallback now uses the SAME params
                gpu=False,
                model_type=model_type,
                diameter=cellpose_diameter,
                flow_threshold=flow_threshold,
                cellprob_threshold=cellprob_threshold,
                min_size=min_size,
            )
            primary = runner_gpu  # primary runner
        else:  # CPU only
            runner_gpu = None  # no gpu runner
            runner_cpu = CellposeRunner(  # honor user diameter/flow/cellprob (was ignored before)
                gpu=False,
                model_type=model_type,
                diameter=cellpose_diameter,
                flow_threshold=flow_threshold,
                cellprob_threshold=cellprob_threshold,
                min_size=min_size,
            )
            primary = runner_cpu  # primary runner

        tile_batch = int(max(1, tile_batch))  # safe batch size

        # Process tiles in batches: batched model.eval amortizes overhead (big GPU win),
        # while merging is still done sequentially to keep label bookkeeping correct.
        for start in range(0, len(tile_info), tile_batch):  # iterate in chunks
            chunk = tile_info[start:start + tile_batch]  # this batch's tiles
            batch_imgs = []  # normalized tile images
            batch_meta = []  # (y0, x0, th, tw) per kept tile
            for (t_idx, y0, x0, th, tw, stdv) in chunk:  # prepare batch
                tile_img = dapi_full[y0:y0 + th, x0:x0 + tw].astype(np.float32, copy=False)  # crop
                tile_norm = normalize01(tile_img)  # normalize
                if float(tile_norm.std()) < float(skip_threshold):  # skip empty tiles
                    continue  # next
                batch_imgs.append(tile_norm)  # queue image
                batch_meta.append((y0, x0, th, tw))  # queue meta
            if not batch_imgs:  # whole batch was empty
                continue  # next chunk

            try:  # batched inference
                if len(batch_imgs) == 1:  # single image -> plain eval
                    masks_list = [primary.eval(batch_imgs[0])]  # one mask
                else:  # multiple -> batched eval
                    masks_list = primary.eval_batch(batch_imgs)  # list of masks
            except Exception as e:  # batch failed -> per-tile fallback
                print(
                    f"[TILE-UPG] WARNING: batch seg error near tile {start}: "
                    f"{type(e).__name__}: {e}",
                    flush=True
                )
                masks_list = []  # rebuild one-by-one
                for im in batch_imgs:  # per-tile retry
                    try:
                        masks_list.append(primary.eval(im))  # try primary
                    except Exception:  # primary failed
                        if (runner_gpu is not None) and bool(force_cpu_on_gpu_error):  # CPU fallback
                            try:
                                masks_list.append(runner_cpu.eval(im))  # cpu
                            except Exception:
                                masks_list.append(None)  # give up on this tile
                        else:
                            masks_list.append(None)  # give up on this tile

            for meta, masks in zip(batch_meta, masks_list):  # merge results sequentially
                if masks is None:  # skip failed tiles
                    continue  # next
                y0, x0, th, tw = meta  # unpack meta
                _merge_core(masks, y0, x0, th, tw)  # merge into global mask

    total_cells = int(global_mask.max())  # total labels
    print(f"[TILE-UPG] Tiling finished; total polygons (labels) = {total_cells}", flush=True)  # log
    return global_mask  # return mask


def save_segmentation_exports_fullres(mask_full: np.ndarray, fig_dir: str, boundary_width: int = 2) -> None:  # save exports
    b = find_boundaries(mask_full, mode="outer")  # find boundaries
    if int(boundary_width) > 1:  # thicken
        b = ndi.binary_dilation(b, iterations=int(boundary_width) - 1)  # dilate
    b_img = (b.astype(np.uint8) * 255)  # 0/255
    boundaries_png = os.path.join(fig_dir, "segmentation_boundaries_fullres.png")  # path
    Image.fromarray(b_img).save(boundaries_png)  # save
    print("Saved:", boundaries_png, flush=True)  # log
    rgb = (label2rgb(mask_full, bg_label=0, alpha=0.90) * 255).astype(np.uint8)  # colored labels
    rgb[b] = (0, 0, 0)  # boundaries black
    colored_png = os.path.join(fig_dir, "segmentation_colored_fullres.png")  # path
    Image.fromarray(rgb).save(colored_png)  # save
    print("Saved:", colored_png, flush=True)  # log
    props = measure.regionprops_table(mask_full, properties=("label", "centroid", "area", "bbox"))  # measure props
    cells_df = pd.DataFrame(props).rename(columns={"label": "cell_id", "centroid-1": "x_centroid_px", "centroid-0": "y_centroid_px", "area": "area_px", "bbox-1": "x_min_px", "bbox-0": "y_min_px", "bbox-3": "x_max_px", "bbox-2": "y_max_px"})  # df
    cells_csv = os.path.join(fig_dir, "cellpose_cells_table_fullres.csv")  # path
    cells_df.to_csv(cells_csv, index=False)  # save
    print("Saved:", cells_csv, flush=True)  # log
