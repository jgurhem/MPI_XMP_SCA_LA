#!/bin/bash
format=clang-format
#use clang-format v8.0.0
# cmd to create .clang-format file :
#$format -dump-config -style=llvm > .clang-format

extensions=(c h cc)
for dir in sources utils
do
  cd $dir
  for elt in "${extensions[@]}"
  do
    find . -name "*.$elt" -print | xargs -i -t $format -style=file -i {}
  done
  cd -
done

