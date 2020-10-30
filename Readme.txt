I've reworked some parts of the code, turned all the sweave files into vanilla .R files

I also added the files sent to me by Jiancang Zhuang. To get the 'inpoly' function to work you need to do run 'R CMD SHLIB poly.f' on the command line once. All the filepaths have already been modified by me to work








A semiparametric spatiotemporal Hawkes-type point process model with periodic background for crime data
J. Zhuang and J. Mateu
J. R. Statist. Soc. A, 182 (2019), 919 -- 942

Data 
The data are in the files

salida_puntos_llamada.csv: boundary of Castallon city center 
castallon_city_boundary.csv: boundary of Castallon city 
type1_crime_data.table: crime data file reported by Castallon Police during 2012 and 2013. 

 Code 
The R program in the .snw files is implemented to estimate the branching structure and background components in the space-time Hawkes model for describing crime behaviors in this study. All the components in the background and clustering parts are estimated nonparametrically by using the stochastic reconstruction method and the relaxation coefficients are estimated by the maximum likelihood method. It is not fully implemented into an R package at this current developing stage, but will be done if there is such a demand. 

R code in snw files. 
crime-excite.snw: main SWEAVE file to run in R 
crime-excite-functions.snw: snw file included by crime-excite.snw 
crime-excite-P1.snw: snw file included by crime-excite.snw 
crime-excite-P2.snw: snw file included by crime-excite.snw 
crime-excite-P3.snw: snw file included by crime-excite.snw 
crime-excite-P4-pn.snw: snw file included by crime-excite.snw 

Instructions for Use 
Requirements: R (3.4.2) with packages “foreach”, “doParallel”, “fields”, “polyCub”, and “spatstat” installed. 
To reproduce the results, first unzip CodeDataSubmittedToJASA-20171204.tar.gz in a temporary fold and then invoke R in command. 
$ tar zxvf CodeDataSubmittedTOJASA-20171204.tar.gz 
$ R 
> Sweave(‘crime-excite.snw’) 
> source(‘plot3d.r’) # reproduce Figure 2 
> source(‘work.r’) # to replace Figure 3 
> source(‘work3.r’) # to replace Figure 4. 
> source (‘residual.r’) # reproduce Figures 5 to 8 
> source(‘residual2.r’) # reproduce Figures 9 and 10 
> q() 
$ pdflatex crime-excite.tex # use latex to compile the output 
$ gv crime-excite.pdf # to see the output 

(1) Figure 1 and results in the first row in Table 1 are produced when running crime-excite.snw 
(2) plot3d.r: Figure 2 in the paper paper. 
(3) residual.r: Figures 5 to 8 in the paper 
(4) residual2.r: Figures 9 and 10 in the paper. 
(5) work.r: Figure 3 in the paper. 
(6) work3.r: Figure 4 in the paper. 

$ cd Nonperiodic 
$ R 
> Sweave(‘crime-excite.snw’) 
> q() 
$ pdflatex crime-excite.tex # use latex to compile the output 
$ gv crime-excite.pdf # to see the output 
(7) The second raw of the second row in Table 1 can produced by running Nonperiodic/crime-excite.snw 


Jiancang Zhuang
Institute of Statistical Mathematics
Research Organisation of Information 
 and Systems
Tokyo 190-8562
Japan

E-mail: zhuangjc@ism.ac.jp
