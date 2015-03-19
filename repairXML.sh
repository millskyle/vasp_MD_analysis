#!/bin/bash


#requires xmllint program from the package  libxml2-utils:
#  do 
#      sudo apt-get install libxml2-utils
#
#  to get this package


#A hackish way to fix the incomplete XML from early-terminated VASP runs.  More elegant solution required
#This basically just closes all dangling XML tags



for tag in `xmllint $1 2>&1 > /dev/null | grep 'Premature end of data in tag' | sed 's|.*Premature end of data in tag ||g' | sed 's| line [0-9]*||g'`; do

echo "</$tag>" >> $1

done
