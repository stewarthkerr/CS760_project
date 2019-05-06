#!/bin/sh

expressions=("me","average")
timesteps=(2,3,4,5,6,7,8)
pids=()
threshold=0.02
IFS=','
for expression in ${expressions[@]}; do
    for timestep in ${timesteps[@]}; do
        #make EXPRESSION=$expression n_dbn_timesteps=$timestep analysis_at &
        cmd="R CMD BATCH '--args expression=\"$expression\" n_dbn_timepoints=$timestep n_reps=2 g_file=\"../data/Arab.Meristem/arabidopsis.meristem.modules.interactions.t${threshold}.tsv\" d_path=\"../output/Clustering.RDATA.RData\"' arab.meristem.analysis.R &"
        echo $cmd
        eval $cmd
        pids+=($!)
    done
done
wait ${pids[@]}
