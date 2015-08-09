 
GOptions Reset=ALL;
GOptions LINESIZE=120 PAGESIZE=50 NoPrompt Vsize=6 Hsize=6 Horigin=1.2 Vorigin=2.5 FText=SwissX FTitle=SwissX HText=2 HTitle=2;

DATA datum;
     INFILE 'C:\Users\Administrator\Desktop\UHCL\2014_9\Multivariate\project\T7-7.dat' DLM='09'x;
     TITLE 'ANALYSIS OF 9.12 DATA';
     INPUT BL EM SF BS AFL LFF FFF ZST; 
run;

proc means Mean Std ;
  var BL EM SF BS AFL LFF FFF ZST;
  run;
proc corr pearson cov;
  var BL EM SF BS AFL LFF FFF ZST;
  run;
/*****simultaneous confidential interval 
if want a bonferroni confidential interval, Changing the C to
  ***********/
Proc IML;
  Use datum;
  Read ALL var{BL EM SF BS AFL LFF FFF ZST} into X; N=NROW(X); P=NCOL(X);
  close datum;
  One=SHAPE( 1,N,1);
  MeanVec=(One`*X)/N;
  *Calculate Mean Vector;
  M=REPEAT(MeanVec,N,1);
  Sigma=(X-M)`*(X-M)/(N-1);
  *Calculate Cov Matrix;
  *print MeanVec, Sigma;
  *C=(((N-1)*P)/(N-P))*FINV(0.95, P, N-P);
  *C=(tinv(1-.025/P,N-1))*(tinv(1-.025/P,N-1));
  *print N P C;
  Start SCI(Mean,S,a,level,P,N);*begin module named SCI;
  C=(((N-1)*P)/(N-P))*FINV(level, P, N-P); *simultaneous confidential interval;
  *C=(tinv(1-.025/P,N-1))*(tinv(1-.025/P,N-1)); *bonferroni confidential interval;
  Margin=sqrt(C*(a`*S*a)/N);*margine of error for T-sq simultaneous CI;
  LB=a`*Mean-Margin;
  UB=a`*Mean+Margin;
  Print LB UB;
  Finish;
  reset noprint;  * turns off auto printing;
  Mean=MeanVec`;
  S=Sigma;
  a={1 0 0 0 0 0 0 0}`;
  Run SCI(Mean,S,a,0.95,P,N);*execute module SCI for various coefficient vectors;
  a={0 1 0 0 0 0 0 0}`;
  Run SCI(Mean,S,a,0.95,P,N); 
  a={0 0 1 0 0 0 0 0}`;
  Run SCI(Mean,S,a,0.95,P,N); 
  a={0 0 0 1 0 0 0 0}`;
  Run SCI(Mean,S,a,0.95,P,N); 
  a={0 0 0 0 1 0 0 0}`;
  Run SCI(Mean,S,a,0.95,P,N); 
  a={0 0 0 0 0 1 0 0}`;
  Run SCI(Mean,S,a,0.95,P,N); 
  a={0 0 0 0 0 0 1 0}`;
  Run SCI(Mean,S,a,0.95,P,N); 
  a={0 0 0 0 0 0 0 1}`;
  Run SCI(Mean,S,a,0.95,P,N); 
quit;

ODS graphics on;
proc corr data=datum nomiss   plots=scatter(ellipse=CONFIDENCE nvar=4 alpha=0.05 0.01); *for new observations (ELLIPSE=PREDICTION) plots=matrix;
  var BL EM SF BS;
  run;
ODS graphics off; 
proc insight data=datum File;
Mult BL EM SF BS AFL LFF FFF ZST;
run; 

 /**end of confidential interval**********/
/********scatter plot*******/
ods html style=Journal; 
title'Scatter Plot Matrix'; 
proc sgscatter data=datum;
  matrix BL EM SF BS AFL LFF FFF ZST  /diagonal=(histogram normal);
run; 

proc sgscatter data=datum;
compare y=BL  x=(AFL LFF FFF ZST)/ellipse=(alpha=0.05 type=predicted);
title 'Scatter Plot Breaking Length (BL) ~ (z1-z4)';
title2 '-- with 95% prediction ellipse';
format BL dollar6.0;
label BL='BL';run;

proc sgscatter data=datum;
compare y=EM  x=(AFL LFF FFF ZST)/ellipse=(alpha=0.05 type=predicted);
title 'Scatter Plot Elastic Modulus (EM) ~ (z1-z4)';
title2 '-- with 95% prediction ellipse';
format BL dollar6.0;
label BL='EM';
run;
proc sgscatter data=datum;
compare y=SF  x=(AFL LFF FFF ZST)/ellipse=(alpha=0.05 type=predicted);
title 'Scatter Plot Stress at Failure (SF) ~ (z1-z4)';
title2 '-- with 95% prediction ellipse';
format BL dollar6.0;
label BL='SF';
run;
proc sgscatter data=datum;
compare y=BS  x=(AFL LFF FFF ZST)/ellipse=(alpha=0.05 type=predicted);
title 'Scatter Plot Burst Strength (BS) ~ (z1-z4)';
title2 '-- with 95% prediction ellipse';
run;

proc sgscatter data=datum;
compare y=BL  x=(EM SF BS)/ellipse=(alpha=0.05 type=predicted);
title 'Scatter Plot (BL) ~ (EM SF BS)';
title2 '-- with 95% prediction ellipse';
run;
proc sgscatter data=datum;
compare y=EM  x=(SF BS)/ellipse=(alpha=0.05 type=predicted);
title 'Scatter Plot (EM) ~ (SF BS)';
title2 '-- with 95% prediction ellipse';
run;
proc sgscatter data=datum;
plot  SF*BS /ellipse=(alpha=0.05 type=predicted);
title 'Scatter Plot SF~BS';
title2 '-- with 95% prediction ellipse';
run;
/********end of scatter plot*******/
/*Univariate Normality Tests per Variable*/
proc univariate data=datum OUTTABLE = NormaliltyTest NORMALTEST NOPRINT; 
   qqplot BL EM SF BS AFL LFF FFF ZST/square ctext=blue;
run;
 
PROC PRINT DATA = NormaliltyTest LABEL NOOBS; 
  var _VAR_ _MEAN_ _STD_ _NORMAL_ _PROBN_ ;
TITLE "Univariate Normality Tests per Variable";  
RUN;
/*standized the variable for outlier judgement*/
proc iml;
  USE datum;
  READ ALL var{BL EM SF BS AFL LFF FFF ZST} INTO x; 
  close datum;
  /* Standardize data: Assume no column has 0 variance */
  start stdMat(x);
    mean = mean(x);                        /* means for columns */
    cx = x - mean;                     /* center x to mean zero */
    std = std(x);                 /* standard deviation estimate*/
    y = cx / std(x);                    /* scaling to std dev 1 */
    return( y );
  finish stdMat;
  nm = {BLs EMs SFs BSs AFLs LFFs FFFs ZSTs};
  std = stdMat(x);
  print std[colname=nm label="Standardized Data"];
  create sdtied from std[colname=nm];
  append from std;
QUIT;
Data datum;
  set datum;
  set sdtied;
  if (abs(BLs)>64*4*0.006) then BLout=1;
  if (abs(EMs)>64*4*0.006) then EMout=1;
  if abs(SFs)>64*4*0.006 then SFout=1;
  if abs(BSs)>64*4*0.006 then BSout=1;
  if (abs(AFLs)>64*4*0.006) then AFLout=1;
  if (abs(LFFs)>64*4*0.006) then LFFout=1;
  if abs(FFFs)>64*4*0.006 then FFFout=1;
  if abs(ZSTs)>64*4*0.006 then ZSTout=1;
run;
proc print Data=datum; 
var AFL LFF FFF ZST AFLs LFFs FFFs ZSTs AFLout LFFout FFFout ZSTout;
run;
proc print Data=datum; 
var BL EM SF BS BLs EMs SFs BSs BLout EMout SFout BSout;
run;
proc sgscatter data=datum;
plot y=(BLs EMs SFs BSs)/ellipse=(alpha=0.05 type=predicted);

/* end of standized the variable */
/* chi-squaare plot for Paper Characteristics*/
PROC IML;
  USE datum;
  READ ALL var{BL EM SF BS} INTO X; 
  close datum;
  N=NROW(X);
  SUMX=X[+,];
  S=(T(X)*X-T(SUMX)*SUMX/N)/(N-1);
  XBAR=SUMX/N;
  Dist2=VECDIAG((X-(J(N,1)*XBAR))*INV(S)*T(X-(J(N,1)*XBAR)));
  ChiSquant=CINV((RANK(Dist2)/(N+1)),NCOL(X));
  CREATE CHIPLOT VAR {Dist2 ChiSquant};
  APPEND;
QUIT;
DATA CHIPLOT;
  set CHIPLOT;
  if Dist2>14.86 then name=_n_;
run;

PROC SGPLOT DATA=CHIPLOT ;
  SCATTER y=Dist2 x=ChiSquant / datalabel= name;
  refline 14.86; *chisquare(4,0.005);
  TITLE "ChiSquare plot of paper property";
RUN; 
/* chi-squaare plot for Pulp Fiber Characteristics*/
PROC IML;
  USE datum;
  READ ALL var{AFL LFF FFF ZST} INTO X; 
  close datum;
  N=NROW(X);
  SUMX=X[+,];
  S=(T(X)*X-T(SUMX)*SUMX/N)/(N-1);
  XBAR=SUMX/N;
  Dist2=VECDIAG((X-(J(N,1)*XBAR))*INV(S)*T(X-(J(N,1)*XBAR)));
  ChiSquant=CINV((RANK(Dist2)/(N+1)),NCOL(X));
  CREATE CHIPLOTFiber VAR {Dist2 ChiSquant};
  APPEND;
QUIT;
DATA CHIPLOTFiber;
  set CHIPLOTFiber;
  if Dist2>14.86 then name=_n_;
run;
proc print DATA=CHIPLOTFiber; run;
PROC SGPLOT DATA=CHIPLOTFiber;
  SCATTER y=Dist2 x=ChiSquant / datalabel= name;
  refline 14.86; *chisquare(4,0.005);
  TITLE "ChiSquare plot of Pulp Fiber Characteristics";
RUN; 
/***********end of chi-square plot*************/
/*reg*/
Proc reg data=datum;
  model BL EM SF BS = AFL LFF FFF ZST / cli xpx I;
  mtest /details print;
      run;

proc glm data = datum;
  class prog;
  model BL EM SF BS  
      = AFL LFF FFF ZST / solution ss3  ;
  manova h=_ALL_;
run;
quit;


* PRINCIPAL COMPONENT ANALYSIS  correlation;
Proc princomp data=datum  std out=pcresult 
plots(ncomp=2)=all n=4 prefix=prin;
var BL EM SF BS;
run;
PROC sgscatter DATA=pcresult;
compare y=prin1  x=(prin2 BL EM SF BS)/ellipse=(alpha=0.05 type=predicted);
title 'Scatter Plot of Prin1 ~ (prin2 BL EM SF BS)';
title2 '-- with 95% ellipse';
RUN;
* PRINCIPAL COMPONENT ANALYSIS  covariance;
Proc princomp data=datum cov std out=pcresult 
plots(ncomp=2)=all n=4 prefix=prin;
var BL EM SF BS;
run;
PROC sgscatter DATA=pcresult;
compare y=prin1  x=(prin2 BL EM SF BS)/ellipse=(alpha=0.05 type=predicted);
title 'Scatter Plot of Prin1 ~ (prin2 BL EM SF BS)';
title2 '-- with 95% ellipse';
RUN;

Proc princomp data=datum cov std out=pcresult 
plots(ncomp=2)=all n=4 prefix=prin;
var AFL LFF FFF ZST;
run;
PROC sgscatter DATA=pcresult;
compare y=prin1  x=(prin2 AFL LFF FFF ZST)/ellipse=(alpha=0.05 type=predicted);
title 'Scatter Plot of Prin1 ~ (prin2 AFL LFF FFF ZST)';
title2 '-- with 95% ellipse';
RUN;

*stdized principle analysis;
PROC standard  data=datum out=datum_t mean=0 std=1;
var BL EM SF BS AFL LFF FFF ZST;
run;
Proc princomp data=datum_t cov std out=pcresult
plots(ncomp=2)=all n=3;
var BL EM SF BS AFL LFF FFF ZST;
run;

/*or*/
* PRINCIPAL COMPONENT ANALYSIS;
PROC FACTOR data=datum OUT=FACTOR METHOD=PRINCIPAL N=7 ALL SCREE PREPLOT REORDER CORR;
     TITLE2 'PRINCIPAL COMPONENT ANALYSIS OF DATA';
*    PLOT THE SCORES;
PROC PLOT DATA=FACTOR;
     TITLE2 'PRINCIPAL COMPONENT SCORES';
     PLOT FACTOR1*FACTOR2;
     PLOT FACTOR1*FACTOR3;
     PLOT FACTOR2*FACTOR3;
     
* FACTOR ANALYSIS ON THE DATA for R (Heywood) ;
PROC FACTOR DATA=datum cov  OUT=FACTOR1 METHOD=prin PRIORS=SMC  N=2              
      ROTATE=VARIMAX  ALL SCREE PREPLOT PLOT REORDER;   
     VAR BL EM SF BS;
run;
Data FACTOR1; set FACTOR1;
if abs(FACTOR1)>2 then name=_n_;
run;

PROC SGPLOT DATA=FACTOR1;
  SCATTER y=FACTOR2 x=FACTOR1 / datalabel= name;
  refline 14.86; *chisquare(4,0.005);
  TITLE "FACTOR ANALYSIS with 2 factors for R";
RUN;  

* FACTOR ANALYSIS ON THE DATA for S  ;
PROC FACTOR DATA=datum  OUT=FACTOR1 METHOD=prin   N=2              
            ALL SCREE PREPLOT PLOT REORDER;   
     TITLE2 'FACTOR ANALYSIS (2) WITH PRIORS=SMC';
     VAR BL EM SF BS;
run;
Data FACTOR1; set FACTOR1;
if abs(FACTOR1)>2 then name=_n_;
run;
PROC SGPLOT DATA=FACTOR1;
  SCATTER y=FACTOR2 x=FACTOR1 / datalabel= name;
  refline 14.86; *chisquare(4,0.005);
  TITLE "FACTOR ANALYSIS with 2 factors for R";
RUN;  
* FACTOR ANALYSIS ON THE DATA for S with MLE ;
PROC FACTOR DATA=datum OUT=FACTOR1 METHOD=ML Heywood PRIORS=SMC N=2              
           ROTATE=VARIMAX ALL SCREE PREPLOT PLOT REORDER;   
     TITLE2 'FACTOR ANALYSIS (2) WITH MLE';
     VAR BL EM SF BS;
run; 
*repeat experiment for AFL LFF FFF ZST PRIORS=SMC ROTATE=VARIMAX;
PROC FACTOR DATA=datum  OUT=FACTOR1  METHOD=prin  N=2              
            ALL SCREE PREPLOT PLOT REORDER;   
     TITLE2 'FACTOR ANALYSIS (2) WITH PRIORS=SMC';
     VAR AFL LFF FFF ZST;
run; 
Data FACTOR1; set FACTOR1;
if abs(FACTOR1)>2 or abs(factor2)>1.5 then name=_n_;
run;
PROC SGPLOT DATA=FACTOR1;
  SCATTER y=FACTOR2 x=FACTOR1 / datalabel= name;
  TITLE "FACTOR ANALYSIS with 2 factors for Pulp Fiber Characteristic";
RUN; 

PROC FACTOR res DATA=datum OUT=FACTOR1_prin METHOD=prin PRIORS=MAX N=1              
           ROTATE=VARIMAX  PREPLOT PLOT ;   
     TITLE2 'TWO FACTORS ANALYSIS WITH Principle Components';
     VAR AFL LFF FFF ZST;
run;  
PROC FACTOR res  DATA=datum  OUT=FACTOR1_MLE METHOD=ML HEYWOOD PRIORS=MAX N=1              
            ROTATE=VARIMAX PREPLOT PLOT ;   
     TITLE2 'TWO FACTORS ANALYSIS WITH MLE';
     VAR AFL LFF FFF ZST;
run;  

Data FACTOR1_MLE;
set FACTOR1_MLE;
rename factor1=factor1mle ;
rename factor2=factor2mle ;
run;
Data factorcomp;
set FACTOR1_MLE ;
set FACTOR1_prin;
name=_n_;
run;
PROC sgscatter DATA=factorcomp;                                                         
     TITLE1 'PLOT OF 1st FACTOR PLOT of PRIN vs. MLE';                               
     PLOT factor1mle*factor1/ datalabel= name; run;
PROC sgscatter DATA=factorcomp;                                                         
     TITLE1 'PLOT OF 2nd FACTOR PLOT of PRIN vs. MLE';                               
     PLOT factor2mle*factor2/ datalabel= name; run;

*    PLOT THE FACTOR SCORED;   
PROC PLOT DATA=FACTOR1;                                                         
     TITLE3 'PLOT OF FACTOR SCORES';                               
     PLOT FACTOR2*FACTOR1;

*    FIND THE FINAL CORRELATION BETWEEN FACTORS AND VARIABLES;
PROC CORR DATA=FACTOR1;
     VAR  FACTOR1 FACTOR2 
          AFL LFF FFF ZST;
     TITLE3 'CORRELATION BETWEEN FACTOR AND VARIABLES';
     
* FACTOR ANALYSIS ON THE DATA METHOD=PRINCIPAL ROTATE=VARIMAX PRIORS=SMC AFL LFF FFF ZST;
PROC FACTOR DATA=datum OUT=FACTOR1 METHOD=PRINCIPAL PRIORS=MAX N=1              
            ALL SCREE PREPLOT PLOT REORDER;   
     TITLE2 'FACTOR ANALYSIS (3) WITH PRIORS=SMC';
     VAR BL EM SF BS ;
 RUN;    
*    PLOT THE FACTOR SCORED;   
PROC sgscatter DATA=FACTOR1;                                                         
     TITLE3 'PLOT OF FACTOR SCORES';                               
     PLOT FACTOR2*FACTOR1;

*    FIND THE FINAL CORRELATION BETWEEN FACTORS AND VARIABLES;
PROC CORR DATA=FACTOR1;
     VAR  FACTOR1 FACTOR2 
          BL EM SF BS AFL LFF FFF ZST;
     TITLE3 'CORRELATION BETWEEN FACTOR AND VARIABLES';


* FACTOR ANALYSIS ON THE DATA;
PROC FACTOR DATA=datum OUT=FACTOR1 METHOD=PRINCIPAL PRIORS=SMC N=2              
           ROTATE=VARIMAX ALL SCREE PREPLOT PLOT REORDER;   
     TITLE2 'FACTOR ANALYSIS (1) WITH PRIORS=SMC';
     VAR BL EM SF BS AFL LFF FFF ZST;
     
*    PLOT THE FACTOR SCORED;   
PROC PLOT DATA=FACTOR1;                                                         
     TITLE3 'PLOT OF FACTOR SCORES';                               
     PLOT FACTOR2*FACTOR1;

*    FIND THE FINAL CORRELATION BETWEEN FACTORS AND VARIABLES;
PROC CORR DATA=FACTOR1;
     VAR  FACTOR1 FACTOR2 
          BL EM SF BS AFL LFF FFF ZST;
     TITLE3 'CORRELATION BETWEEN FACTOR AND VARIABLES';
Run;


*Canonical analysis;
Proc cancorr corr data=datum_t vprefix=paper vname='BL-BS' wprefix=fiber wname='AFL-ZST';
var BL EM SF BS ;
with AFL LFF FFF ZST;
run;
proc import datafile= "C:\Users\Administrator\Desktop\UHCL\2014_9\Multivariate\ex7_3.xlsx" out=ex7_3 dbms=xlsx replace; 
proc print data =ex7_3;  run;
Proc cancorr simple corr data=ex7_3 (type=corr);
var X1-X9;
with Y1-Y5;
run;

