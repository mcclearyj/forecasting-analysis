#!/bin/bash

for i in {0..29}; do
    dir="r$i"
    file="$dir/pipe.log"
    if [ -f "$file" ]; then
        if grep -q "Finished pipeline run" "$file"; then
            echo "Phrase found in $file"
        else
            echo "Phrase not found in $file"
        fi
    else
        echo "$file not found"
    fi
done
