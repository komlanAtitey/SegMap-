python3 segmentation_sensitivity.py \
  --datasets AD=XeniumData_segmap/cellpose_pipeline_outputs_1/data \
             GBM=XeniumData_segmap/cellpose_pipeline_outputs_2/data \
             HEALTHY=XeniumData_segmap/cellpose_pipeline_outputs_3/data \
             LUNG=XeniumData_segmap/cellpose_pipeline_outputs_4/data \
  --outdir seg_sensitivity_out
