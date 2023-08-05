#!/bin/bash

# Set variables
dataset="ba"
query_set_sizes=(20 40 60 80 100)
program_names=("base")
k_s=(5 6 7 8 9 10 11 12 13 14 15)

# Loop over all combinations of query set sizes and program names
for k in ${k_s[@]}
do
    for query_set_size in ${query_set_sizes[@]}
    do
        for program_name in ${program_names[@]}
        do
            # Generate folder and log names
            log_folder=${dataset}/${program_name}
            output_log=${k}_${query_set_size}_output.log
            input_txt=${k}_${query_set_size}.txt

            # Empty log files
            echo "" > ${log_folder}/${output_log}
            echo "" > ${log_folder}/${input_txt}

            # Loop
            for i in {1..100}
            do
            ./${program_name} ${k} ${dataset}/${dataset} ${dataset}/Q/${dataset}_${query_set_size} "${i}" ${log_folder}/${input_txt} >> ${log_folder}/${output_log}
            done
        done
    done
done