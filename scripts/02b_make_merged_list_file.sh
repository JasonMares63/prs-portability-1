#!/bin/bash

for i in {1..22}
do
  printf "data/merged/chr%s\n" $i >> data/merged_list.txt
done
