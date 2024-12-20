#!/bin/bash

# Check if the number of arguments is valid
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <number_of_iterations>"
    exit 1
fi

# If input is 1, no need to run further iterations
if [ "$1" -eq 1 ]; then
root -l  'GRINCH_LE_align_v2.C("zeros","try0")'
    exit 0
fi

# Always run the first command
root -l -q -b 'GRINCH_LE_align_v2.C("zeros","try0")'


# Loop through the given number of iterations
for ((i = 0; i < $1-1; i++)); do
    next=$((i + 1))
    if [ $i -eq $(($1 - 2)) ]; then
        root -l "GRINCH_LE_align_v2.C(\"try$i\",\"try$next\")"
    else
        root -l -q -b "GRINCH_LE_align_v2.C(\"try$i\",\"try$next\")"
    fi
done


