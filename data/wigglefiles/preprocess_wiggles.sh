#!/bin/bash
# here I make one DF from the wiggle file, making it easier to import it in pandas.

for i in $(ls ./*_no_PNA_3primeend_new.wig)
do
    NAME=${i%.wig}
    echo "${NAME}_mod.wig"
    grep -v "#" $i | grep -P "b\d{4}" > "${NAME}_mod.wig"
done

