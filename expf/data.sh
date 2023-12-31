#!/bin/bash

# Set variables
dataset="ep"
query_set_sizes=(20 40 60 80 100)
program_names=("ad" "ade")
k_s=(10)

# Generate folder and log names
output_log=${dataset}/query_output_sum.txt

# Empty log files
echo "" > ${output_log}

# Loop over all combinations of query set sizes and program names
for k in ${k_s[@]}
do
    for query_set_size in ${query_set_sizes[@]}
    do
        for program_name in ${program_names[@]}
        do
            # Construct the input file path
            input_file=${dataset}/${program_name}/${k}_${query_set_size}.txt
            
            # Construct the command
            cmd="python3 datacollect.py ${input_file}"

            # Echo the command to the console
            echo "Running: ${cmd}" >> ${output_log}



            # Execute the command and append output to the log file
            ${cmd} >> ${output_log}
        done
    done
done
