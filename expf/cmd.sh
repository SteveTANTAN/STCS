#!/bin/bash
echo "" > ba/ba_100_f/basline_output.log
echo "" > ba/ba_100_f/basline.txt
for i in {1..100}
do
   ./base 5 ba/ba ba/ba_100 "${i}" ba/ba_100_f/basline.txt >> ba/ba_100_f/basline_output.log
done