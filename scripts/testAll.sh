#!/bin/sh -e

for i in testMPI*.sh testXMP*.sh
do
  sh $i
done
