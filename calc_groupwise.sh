#!/usr/bin/zsh
mkdir -p groupwise_whole;
for trait (full_scores_whole/*.full_score*); do
  scdrs perform-downstream --h5ad-file ./converse_data/hs_whole_snn.h5ad \
  --score-file $trait \
  --out-folder groupwise_whole/ \
  --group-analysis cluster_lab \
  --flag-filter-data True \
  --flag-raw-count True;
done