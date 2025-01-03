#!/bin/bash
executable=$1
inputFileTag=$2
outputFileTag=$3
datasetName=$4
process=$5
phoID=$6

# Path to the input text file containing the list of root file
input_file=$inputFileTag

# Path to the output file where the line will be written
execute_file="${input_file/.txt/}"_execute.txt
echo $execute_file
counter=1

# Read the input file line by line
while IFS= read -r line
do
    # Write the current line to execute.txt
    output_root="${outputFileTag/.root/}"_$counter.root
    echo
    echo $output_root
    echo /eos/home-c/cagrawal/SusySoftPhoAna/FR_estimation/SF_TnP/"${line##*/}" > "$execute_file"
    echo copying "$line"
    xrdcp -f "$line" /eos/home-c/cagrawal/SusySoftPhoAna/FR_estimation/SF_TnP/tmp_root_files/
    
    # Execute the other.sh script with execute.txt as an argument
    #./other.sh "$execute_file"
    echo EXECUTING "${line##*/}"
    ./$executable "$execute_file" $output_root $datasetName $process $phoID

    if [ $counter -eq 1 ]; then
	hadd -f "${outputFileTag/.root/}"_file$counter.root $output_root
    else
	hadd -f "${outputFileTag/.root/}"_file$counter.root $output_root "${outputFileTag/.root/}"_file$(expr $counter - 1).root
	rm "${outputFileTag/.root/}"_file$(expr $counter - 1).root
    fi
    rm $output_root

    ((counter++))
    remove="${line##*/}"
    echo removing /eos/home-c/cagrawal/SusySoftPhoAna/FR_estimation/SF_TnP/tmp_root_files/"$remove"
    rm /eos/home-c/cagrawal/SusySoftPhoAna/FR_estimation/SF_TnP/tmp_root_files/"$remove"

    # Clear the contents of execute.txt after execution
    > "$execute_file"
done < "$input_file"
echo removing "$execute_file"
rm $execute_file
