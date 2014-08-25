for a in `ls *.LRT.result|sed -e 's/.LRT.result//'`;
do
        chi2 1 `cat $a.LRT.result`  | awk 'NF>1'|awk '{FS="=";print $6}'>$a.Chi2.result;
done;
