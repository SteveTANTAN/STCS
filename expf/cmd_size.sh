#!/bin/bash

# Set variables
dataset="ba"
# query_set_sizes=(20 40 60 80 100)
query_set_sizes=(60)
program_names=("ad" "ade" "base")
# k_s=(3 4 5 6 7 8 9 10 11 12 13 14 15)
k_s=(3 4 5 6 7)
# k_s=(4)

# Timeout in seconds (3 hours = 10800 seconds)
# timeout_duration=10800
timeout_duration=3600

# Memory threshold in KB (350GB)
memory_threshold=$((350 * 1024 * 1024))

# Loop over all combinations of query set sizes and program names

for query_set_size in ${query_set_sizes[@]}
do
    for k in ${k_s[@]}
    do
        for program_name in ${program_names[@]}
        do
            # Generate folder and log names
            log_folder=${dataset}/${program_name}
            output_log=${k}_${query_set_size}_output.log
            input_txt=${k}_${query_set_size}.txt

            # Ensure log folder exists
            mkdir -p ${log_folder}

            # Empty log files
            echo "" > ${log_folder}/${output_log}
            echo "" > ${log_folder}/${input_txt}

            # Loop
            for i in {1..100}
            do
            # Start the program in the background with a timeout
            timeout ${timeout_duration} ./${program_name} ${k} ${dataset}/${dataset} ${dataset}/Q/${dataset}_${query_set_size} "${i}" ${log_folder}/${input_txt} >> ${log_folder}/${output_log} &
            PID=$!
                            
            # Monitor the memory usage
            while [[ -e /proc/$PID ]]; do
                RSS=$(ps -p $PID -o rss=)
                                
                if [[ $RSS -gt $memory_threshold ]]; then
                    echo "Memory usage has exceeded 350GB for program ${program_name} with k=${k} and query_set_size=${query_set_size}. Killing the process." >> ${log_folder}/${output_log}
                    kill $PID
                    wait $PID
                    break
                fi
                                
                sleep 1
            done

            # Check for timeout
            wait $PID
            if [ $? -eq 124 ]; then
                echo "Execution of program ${program_name} with k=${k} and query_set_size=${query_set_size} timed out." >> ${log_folder}/${output_log}
            fi

            done
        done
    done
done
