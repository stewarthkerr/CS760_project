#!/bin/sh

expressions=("me","average")
timesteps=(2,3,4,5,6,7,8)
pids=()
threshold=0.02
IFS=','
for expression in ${expressions[@]}; do
    for timestep in ${timesteps[@]}; do
        cp "../output/arab.meristem.analysis.${expression}.${timestep}.RData" "../output/arab.meristem.analysis.atrm.${expression}.${timestep}.RData"
        cmd="R CMD BATCH '--args g_file=\"../data/ATRM/Regulations_in_ATRM.modules.interactions.t${threshold}.tsv\" save_file=\"../output/arab.meristem.analysis.atrm.${expression}.${timestep}.RData\"' arab.meristem.perf.R"
        echo $cmd
        eval $cmd
        pids+=($!)
    done
done
wait ${pids[@]}
