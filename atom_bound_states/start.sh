#!/bin/bash
	read n
	i=0
	while [ $i -lt $n ]
	do
	gnuplot -e "out_file='test'; in_file='output.txt'; iterator=${i}" graphics.graph
	let "i=i+1"
	done