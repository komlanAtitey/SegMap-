python -m segmap.run_segmap \
  --xenium_dir output-XETG00245__0058077__2_Day_Crush_DRGs__Naive_DRGs__9_Day_Transection_DRGs__20250806__180811 \
  --gpu \
  --tile_batch 16 \
  --no_tune_resolution \
  --resolution 0.8 \
  --cellpose_diameter 35 \
  --flow_threshold 0.4 \
  --cellprob_threshold 0.2 \
  --top_gene_cap 5050
