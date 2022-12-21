#!/bin/bash

a=$(cat $2 | cut -d ' ' -f1 | head -n 1)
a=$((a-2))
echo $a
parallel -j8 "./main {} < $2 > out/{1}_{2}_{3}" ::: `seq 0 $a` ::: `seq 0 $a` ::: $1

