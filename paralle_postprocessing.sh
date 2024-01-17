#!/bin/bash

# Define the parameters
cofs=(0.0 0.4 1.0 10.0)
#aspectRatio=(1.0 1.5 2.0 2.5 3.0)
phis=(0.5 0.6 0.7 0.8 0.9)
aspectRatio=(1.0 2.0 2.5 3.0)

cofs=(0.4)

phis=(0.5 0.6 0.7 0.8 0.9)

# Define the maximum number of parallel tasks
max_parallel_tasks=50

echo "Start of the loop"

# Loop through the values of -var flag for COF
for cof in "${cofs[@]}"
do
    # Loop through the values of -var flag for ap
    for ap in "${aspectRatio[@]}"
    do
        # Loop through the values of -var flag for phi
        for phi in "${phis[@]}"
        do
            # Create a job script for each combination of parameters
            job_script="post_process_job_c${cof}_a${ap}_phi${phi}.sh"
            echo "#!/bin/bash" > "$job_script"
            echo "python main.py -c $cof -a $ap -t phi -v $phi" >> "$job_script"
            chmod +x "$job_script"

            # Submit the job script to the background
            ./"$job_script" &
            
            # Limit the number of parallel tasks
            running_tasks=$(jobs -p | wc -l)
            while [ $running_tasks -ge $max_parallel_tasks ]; do
                sleep 1
                running_tasks=$(jobs -p | wc -l)
            done
        done
    done
done

# Wait for all background jobs to finish
wait

echo "All post-processing jobs completed."
