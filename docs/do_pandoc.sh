#!/bin/bash

for i in *.md
do
j=${i%%.md}
pandoc -f markdown -t html -s -o $j.html $i
done
