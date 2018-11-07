CSE-6010 Assignment 2 part 1 K mean clustering

to compile:
cd CSE6010-A2
gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF src/KM.d -MT src/KM.o -o src/KM.o src/KM.c

the usage should be:
<directory_to_program> <input_file_path> <cluster_number> <output_file_path>