#!/bin/bash

count=0
for i in *.dat; do
    echo Standlone on file $i ...
    ./driver_conserv_to_prims $i &
done

wait
echo Finished!



