# syntenyPlotByR<br>
plot synteny genome alignment by mapper Mummer or LastZ.<br>
<br>
###########<br>
requriement<br>
###########<br>
R<br>
optparse # an R package <br>
<br>
###########<br>
plot mummer delta file<br>
###########<br>
There is a sample delta file and an output named "out.png" in project. user could try this by yourself.<br> 
<i>Rscript syntenyPlot.R -i &ltdelta.file&gt -o out.png</i><br>
<br>
###########<br>
Large delta file<br>
###########<br>
If the delta file is too large, it is very time-consuming for ploting it. We recommand format the delta file by our PERL script firstly. And then plot the format delta file.<br>
<i>perl nucmer2RPlotFormat.pl &ltdelta.file&gt &gtout.format</i><br>
<i>Rscript syntenyPlot.R -i &ltdelta.file&gt -c -o out.png</i><br>
<b>If the delta is formatted, the parameter -c should be setted.</b><br>
<br>
###########<br>
LastZ<br>
###########<br>
There are two steps to format LastZ into our plot format. Only support one chromosome in Reference right now.<br>
<i>perl maf2R_01.pl *.maf &gtout.maf_1</i><br>
<i>perl maf2R_02.pl *out.maf_1 &gtout.maf_2</i><br>
<i>Rscript syntenyPlot.R -i &ltout.maf_2&gt -c -o out.png</i><br>
<b>If the delta is formatted, the parameter -c should be setted.</b><br>

###########<br>
R Plot Format <br>
###########<br>
There are 11 columns in Our R plot Format. User could change any alignment format into our R plot format and then plot it by our R script. And there is a sample file named "out.format" in project.<br>
# refid <br>
# qryid <br>
# reflen <br>
# qrylen <br>
# refstart <br>
# refend <br>
# qrystart <br>
# qryend <br>
# refollen <br>
# belong ref start <br>
# belong ref <br>
<br><br>