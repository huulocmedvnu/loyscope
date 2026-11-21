# shared settings; override by exporting before running any step
export LOY_RESULTS=${LOY_RESULTS:-$PWD/results}        # all outputs
export LOY_DATA=${LOY_DATA:-$PWD/data}                 # temporary downloads (deleted after extraction)
export THREADS=${THREADS:-16}                           # rayon threads per job
export JOBS=${JOBS:-2}                                  # datasets processed in parallel
export BIN=${BIN:-$PWD/target/release}
export TSS=${TSS:-$PWD/resources/tss.tsv}
export RAYON_NUM_THREADS=$THREADS
mkdir -p "$LOY_RESULTS/extract" "$LOY_RESULTS/tables" "$LOY_RESULTS/logs" "$LOY_DATA"
