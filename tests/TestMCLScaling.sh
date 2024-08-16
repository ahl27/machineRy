#!/bin/sh

# This script requires either gtime (gnu-time, OSX) or time (Linux)
# Set the below alias accordingly
alias ftime="gtime -f '%E %M'"

for loopctr in 50 100 250 500 1000 2500 5000 10000 50000 100000 250000 500000 1000000 2500000 5000000 10000000
do
echo "$loopctr vertices and edges:"
# Create the graph
R -f tests/MemProfilingSetup.R --args ${1} $loopctr $loopctr >/dev/null

ftime ~/local/bin/mcl ${1} --abc -I 2.0 -o ./trashres.mcl -V all

done

# Elapsed time, max memory (KB)
# 50 vertices and edges:
# 0:00.01 1984
# 100 vertices and edges:
# 0:00.00 2208
# 250 vertices and edges:
# 0:00.00 2608
# 500 vertices and edges:
# 0:00.01 3200
# 1000 vertices and edges:
# 0:00.01 3184
# 2500 vertices and edges:
# 0:00.05 5008
# 5000 vertices and edges:
# 0:00.13 9840
# 10000 vertices and edges:
# 0:00.33 16384
# 50000 vertices and edges:
# 0:02.67 63072
# 100000 vertices and edges:
# 0:06.70 107472
# 250000 vertices and edges:
# 0:26.80 268304
# 500000 vertices and edges:
# 1:26.93 423056
# 1000000 vertices and edges:
# 5:03.19 797568
# 2500000 vertices and edges:
# 29:18.34 1332720
# 5000000 vertices and edges:
# 1:51:25 3174512
#10000000 vertices and edges:
# 7:15:45 4408864