#!/usr/bin/env bash
if [ x"$1" == x ]; then

        echo "please specify a the split size"

        exit 1

fi

awk -v splitval=$1 'BEGIN{oldid=""; num=0; pieces = 0} {
	if ($1~">") {  
		split($0,a,"_"); id=a[2]; 
		if (oldid!=id) { 
			num++; oldid=id 
		} 
		if (num==splitval) {
			num = 0
			pieces++
		}
		outfile = "pieces_"pieces".fa"
	} 
	print $0 >outfile
}' Trinity.fasta 
