for a in `ls *.H0.result` ; do cat $a | grep lnL|awk '{FS=":";print $5}' >$a.Ln0; done;
for a in `ls *.H1.result` ; do cat $a | grep lnL|awk '{FS=":";print $5}' >$a.Ln1; done;
