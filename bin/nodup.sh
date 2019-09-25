#!/bin/sh

#put lines in IFILE into OFILE checking for duplicates
IFILE=$1
OFILE=$2


# commented such that one can use it to remove duplicates lines from a file itself
#if [[ ! -f $OFILE ]]; then
# cat $IFILE > $OFILE
#else
 for L in `cat $IFILE`
 do
  if [[ ! `grep "$L" $OFILE` ]]; then
   echo $L >> $OFILE
  fi
 done
#fi
