# TCGAGenes

**Setup/install**

Make sure you have R as well as all of the packages loaded at the top of survivalgenes.R.

Download the repository and make sure survivalgenes.R is executable using the below code:
```
git clone https://github.com/lainemarrah/TCGAGenes
cd TCGAGenes
chmod u+x survivalgenes.R
```

**Usage**
```
Rscript survivalgenes.R [CANCER OF INTEREST]
```
Cancer of interest should be one of the TCGA studies listed in studies.txt. This will output a whitespace-separated file listing each of the genes found to be correlated with survival in that cancer from lowest to highest p-value below 0.01. This file also includes some statistical info used to arrive at that conclusion. Overall it's pretty fast, each cancer took about 5-15 on an 8 GB i5 Mac.

If you want to loop this to run every cancer simultaneously you can run the below. Fair warning that this did crash my computer when I ran it locally but it would probably be fine on HPC.
```
while read p; do Rscript survivalgenes.R "$p" &; done < studies.txt
```


