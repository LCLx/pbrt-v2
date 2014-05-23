#!/bin/bash
EXRTOTIFF="../bin/exrtotiff"
EXRDIFF="../bin/exrdiff"
PBRT="../bin/pbrt"
DATAFILE="hw4_data.txt"

#	part 2: occluder_shadow
for p in *.pbrt;
	do ../$PBRT $p 
done

#	compute the variance
if [ -f $DATAFILE ];
then
	rm $DATAFILE
fi

for p in *.exr;
	do ../$EXRDIFF $p ../occluder_shadow_ref.exr >> $DATAFILE
done

#	transform exr files into tiff files
for f in *.exr;
	do ../$EXRTOTIFF $f ${f/.exr/}.tiff;
done
