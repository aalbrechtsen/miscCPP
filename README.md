# miscCPP
simple functions for c/cpp to test. 



 * readBinGZ
read large (>4Gb) binary gz files of unknown size 

test of dat with 755925680 doubles with bufSize 100000 doubles

> time ./readBinGz -f /home/albrecht/gitfuck/test/ceu8gl2.glf.gz 


readGZbin	reading failed nSites755925680 != nSites2 219054768
	[ALL done] walltime used =  30.00 sec
read 219054768 sites

readGZbin2		[ALL done] walltime used =  96.00 sec
read 755925680 sites

readGZbin3		[ALL done] walltime used =  64.00 sec
read 755925680 sites

readGZbin4		[ALL done] walltime used =  40.00 sec
read 755925680 sites

> time ./readBinGz -f /home/albrecht/gitfuck/test/ceu8gl2.glf.gz 
readGZbin	reading failed nSites755925680 != nSites2 219054768
	[ALL done] walltime used =  29.00 sec
read 219054768 sites

readGZbin2		[ALL done] walltime used =  46.00 sec
read 755925680 sites

readGZbin3		[ALL done] walltime used =  24.00 sec
read 755925680 sites

readGZbin4		[ALL done] walltime used =  24.00 sec




Smaller file (<4GB)
> time ./readBinGz -f /home/albrecht/gitfuck/test/tmpMulti.glf.gz

readGZbin		[ALL done] walltime used =  9.00 sec
read 168909600 sites

readGZbin2		[ALL done] walltime used =  10.00 sec
read 168909600 sites

readGZbin3		[ALL done] walltime used =  7.00 sec
read 168909600 sites

readGZbin4		[ALL done] walltime used =  11.00 sec
read 168909600 sites




** Conclusion


doubeling the size is most effective (single readin) (up to 100% ram use)
reading in the data twice is fairly effective (method 2) (+50% time)

hmmm. wierd..... behaivior 