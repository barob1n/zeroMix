#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "defns.h"
#include "segy.h"
#include "segyIO_class.h"
#include <time.h>
#include <omp.h>

#define FLUSH while (getchar() !='\n');
//void velMix(float ***velTr, float ***velOutTr,int ns, int numTr,int iNdx, float***zeroTr, int *xNdxMn, int *xNdxMx, int iOrigin, int xOrigin, int rate);
void velMix(float ***velTr, float ***velOutTr,int ns, int numTr,int iNdx, float***zeroTr, int *xNdxMn, int *xNdxMx, int iOrigin, int xOrigin, int rate,int iline, int xline);
int slopeChk( float slopeOld, float slopeIXl, float dipMax);
int nearest (float input);
int boundCheck (int sampPlusWin, int slope, int ns);

/***********************************************************************
 * Function: main
 * Input/output:  input seismic file name, output dip file name, 
 * 		output dip file name, number of traces used in semblance calc,
 * 		dip as max ms/trace allowable, semblance window size
 * 	      
 * Descrip: Mixes velocity voulume along geologic dips.
 * ********************************************************************/
int main(int argc, char *argv[])
{
omp_set_dynamic(0);     // Explicitly disable dynamic teams
omp_set_num_threads(7);

FILE *finVel=NULL, *fzero=NULL, *foutVel=NULL;
int   fmt, i, j , itr, ierr, k;
int   ns;	/* Number of samples in a trace	*/
float dts;	/* Sample rate in seconds	*/
segy  thdr;	/* Segy trace header		*/
segy thdrI, thdrX, thdrS;
bhed bhI,bhX,bhS;
bhed  bh;	/* Segy file binary header	*/
char chdr[3202];	/* segy file character header */
int   endian=1; /* Which side of the egg	*/
int   status=0; /* Return status for segy calls, I should check*/
char  txthdr[3201];	/* The character file header*/

int choice;
char buff[80];

segy ***trGthHdr;
float *tempTrZ;
float *tempTrV;
float ***dipTrI;
float ***dipTrX;
float ***zeroTr;
float ***velTr;
float ***velOutTr;

int *xNdxMx; /*Keeps track of number of xline on each inline*/
int *xNdxMn; /*first xline on each inline*/
int statusR; /*Status of last segyRead*/
int statusW; /*Status of last segyWrite*/
int numTr; /*Radius of traces, about the trace being examined, used in
             the semblance-dip calculation*/
int maxDip; /*maximum allowed dip, IN SAMPLES PER TRACE*/
int window; /*Size of sembalnce window in number of samps*/
int shift; /*shifts starting sample point down for dip calc.  if you
				started at the 0th sample, the dip calc would run into
				* negative sample points*/
int iOrigin = 0; /*origin of inline grid*/
int xOrigin = 0; /*origin of xline grid*/
int iMax = 0;
int xMax = 0;
int line=0;
int iline = 0;
int xline = 0;
int ilineOld;			
int ilineCnt=0;
float thresh = 0;
int skip = 0;
int rate;
int ii;
int num10Percent;
int nsZ;

time_t now;

float temp;
endian = checkEndian();

thdr.iline = thdr.ep;
thdr.xline=thdr.cdpt;

	if(argc ==9){
/*Open fileS as readable/writable binary*/
		finVel = fopen(argv[1], "rb");
		if (finVel == NULL) {
			printf("Unable to open the Input file.  Please check the name.\n");
		return -1;
		}
		
		fzero = fopen(argv[2], "rb");
		if (fzero == NULL) {
			printf("Unable to open the semblance input file.  Please check the name.\n");
		return -1;
		}
			
		foutVel=fopen(argv[3],"r");
		if(foutVel!=NULL){
			printf("\noutput velocity file Seems to exist. ");
			printf("\nWould you like to overwrite it (-1 = yes)?");
			fgets(buff,80,stdin);
			sscanf(buff,"%d",&choice);
			if (choice == -1){
				fclose(foutVel);
				foutVel=fopen(argv[3],"wb");
			}else{
				fclose(foutVel);
				return 0;
			}
		choice=0;
		}else {
			foutVel=fopen(argv[3],"wb");
		}

		if(sscanf(argv[4],"%d", &numTr)!=1){
			printf("\nnumTr entered not valid. Exiting... ");
			return -1;
		}
		
		if(sscanf(argv[5], "%d", &iOrigin) != 1){
			printf("\niOrigin entered not valid. Exiting...");
			return -1;
		}
		if(sscanf(argv[6], "%d", &xOrigin) != 1){
			printf("\nxOrigin entered not valid. Exiting...");
			return -1;
		}
		
		if(sscanf(argv[7], "%d", &iMax) != 1){
			printf("\nxOrigin entered not valid. Exiting...");
			return -1;
		}
		if(sscanf(argv[8], "%d", &xMax) != 1){
			printf("\nxOrigin entered not valid. Exiting...");
			return -1;
		}		
		
	} else {
	fprintf(stderr,"\n***************************************************************************\n");
      fprintf(stderr,"Mixes velocities along geologic dip. Currently, the input dip files\n\n");
      fprintf(stderr,"are expected to be in ms/trace. The program rounds dip to nearest sample. \n\n");
      fprintf(stderr,"Program expects the following command line: \n ");
      fprintf(stderr,"sembDip <vel.sgy><zero.sgy><out.sgy><Num><iOrigin><xOrigin><iMax><xMax>\n\n");
      fprintf(stderr,"vel.sgy: Input velocity file. \n\n");
      fprintf(stderr,"zero.sgy: Output semblance  sgy filename.\n\n");
      fprintf(stderr,"Out.sgy: Output velocity  sgy filename.\n\n");
      fprintf(stderr,"Num: Number of traces for semblance given as  \n");
      fprintf(stderr,"     distance from trace being analysed. \n");
      fprintf(stderr,"     ex. numTr = 1 uses 9 traces. 3 by 3 block.\n");
      fprintf(stderr,"                 o--o--o\n");
      fprintf(stderr,"                 |  |  |\n");
      fprintf(stderr,"                 o--x--o\n");
      fprintf(stderr,"                 |  |  |\n");
      fprintf(stderr,"                 o--o--o\n\n");
      fprintf(stderr,"iOrigin: Starting inline number.\n\n");
      fprintf(stderr,"xOrigin: Starting xnline number.\n\n");
      fprintf(stderr,"iMax: Max iline. \n\n");
      fprintf(stderr,"xMax: Max xline.");
    fprintf(stderr, "\n**************************************************************************\n");
      return 0;
    } 

	/*print parameters. ask user to verify they are correct-option to abort*/
	printf("\nInput paramters are:");
	printf("\n\n  Input Velocity: %s", argv[1]);
	printf("\n  Input zero: %s", argv[2]);
	printf("\n  Output velocity:  %s",argv[3] );
	printf("\n  Trace radius: %d", numTr);
	printf("\n  inline origin: %d", iOrigin);
	printf("\n  xline origin: %d", xOrigin);
	printf("\n  iMax : %d", iMax);
	printf("\n  xMax: %d", xMax);
	
	printf("\n\nAre these correct? (-1: exit)");
	fgets(buff,80,stdin);
	sscanf(buff,"%d",&choice);
	if(choice == -1){ return 0;} 
	
	segyReadHeader(fzero, chdr, &bh, endian);
	nsZ=bh.hns;
	
	segyReadHeader(finVel, chdr, &bh, endian);
	ns=bh.hns;
	rate = bh.hdt;
	rate=rate/1000;
	
	if(nsZ !=ns){
		printf("\n\n Xline dip file's number of samps not matching input velocity file's number of samps!");
		return -1;
	}
	
	printf("\nNumber of samples: %d", ns);
	
	if(segyWriteHeader(foutVel, &chdr, &bh, endian)!=0){
		printf("\n Unable to write header for output data!! \n");
		return -1;
	}
	
	/*Converting numTr to the width of cube ex. 3x3*/
	numTr=(numTr)*2+1;
	
/* ******************************************************************* */
/* ******************************************************************* */
	/*allocating memory for the traces used for semblance calc/analysis*/
	velTr = (float***)calloc((numTr) , sizeof(float **));
    for (i=0; i<(numTr); ++i){
         velTr[i] = (float **)calloc((xMax-xOrigin+2),sizeof(float*));
         for ( j =0; j < (xMax-xOrigin+2); ++j){
			 velTr[i][j] = (float *) calloc((ns+1),sizeof(float));
		 }
	 }
	 
	 zeroTr = (float ***)calloc((numTr) , sizeof(float **));
     for (i=0; i<(numTr); ++i){
         zeroTr[i] = (float **)calloc((xMax-xOrigin+2),sizeof(float*));
         for ( j =0; j < (xMax-xOrigin+2); ++j){
			zeroTr[i][j] = (float *) calloc((ns+1),sizeof(float));
		 }
	 }
	 
	velOutTr = (float***)calloc((numTr) , sizeof(float **));
    for (i=0; i<(numTr); ++i){
         velOutTr[i] = (float **)calloc((xMax-xOrigin+2),sizeof(float*));
         for ( j =0; j < (xMax-xOrigin+2); ++j){
			 velOutTr[i][j] = (float *) calloc((ns+1),sizeof(float));
		 }
	 }
	 	
	 trGthHdr = (int ***)calloc((numTr) , sizeof(int **));
	 for (i=0; i<(numTr); ++i){
         trGthHdr[i] = (int **)calloc((xMax-xOrigin+2),sizeof(int*)); 
         for ( j =0; j < (xMax-xOrigin+2); ++j){
			 trGthHdr[i][j] = (int *) calloc((1),sizeof(segy));
		 }
	 }

	tempTrZ = (float *)calloc((ns+1),sizeof(float)); 
	tempTrV = (float *)calloc((ns+1),sizeof(float));
	 
	xNdxMx = (int *)calloc ((iMax - iOrigin +2 + numTr), sizeof(int)); //to track the starting and ending xlines
	xNdxMn = (int *)calloc ((iMax - iOrigin +2 + numTr), sizeof(int)); //to track the starting and ending xlines
/* ******************************************************************* */
/* ******************************************************************* */

	 /*Reading in frist Trace, I am currently reading in the vel trace last since it is the header I am most concerned with*/
	statusR = segyReadTrace(fzero, &bhS, &thdrS,tempTrZ, ns, endian);
	statusR = segyReadTrace(finVel, &bh, &thdr,tempTrV, ns, endian);
	
	
	
	/*im setting the iOrigin to the first inline read in...will be useful later*/
	if(thdr.iline != iOrigin){
		printf("\n\nFirst inline read does not match iOrigin: %d\n",iOrigin);
		printf("Setting iOrigin equal to first iline read in!: %d\n\n",thdr.iline);
		iOrigin=thdr.iline;
	}
	
	/*are the values sensibale?*/
	if (thdr.iline > iMax){
		printf("\n\nFirst iline read great than iMax.\n\n");
		return 0;
	}else if (thdr.xline - xOrigin < 0){
		printf("\n\nFirst xline read less than xOrigin.\n\n");
		return 0;
	}else if (thdr.xline > xMax){
		printf("\n\nFirst iline read great than iMax.\n\n");
		return 0;
	}

	/*save the xline values and set iline switch...really only need to set the Mn value*/
	xNdxMx[thdr.iline - iOrigin] = thdr.xline - xOrigin;//the first xline will be zero if it is at the xline axis
	xNdxMn[thdr.iline - iOrigin] = thdr.xline - xOrigin;//the first xline will be zero if it is at the xline axis
	ilineOld=thdr.iline;
	
	/*read in the first numtTr number of inlines.*/
	while ( thdr.iline - iOrigin < numTr  && statusR==0){
			
		memcpy(zeroTr[thdr.iline - iOrigin][thdr.xline -xOrigin], tempTrZ,(ns+1)*sizeof(float));
		memcpy(velTr[thdr.iline - iOrigin][thdr.xline -xOrigin], tempTrV,(ns+1)*sizeof(float));
		memcpy(trGthHdr[thdr.iline - iOrigin][thdr.xline -xOrigin], &thdr, sizeof(*trGthHdr[0][0]));
		
		statusR = segyReadTrace(fzero, &bhS, &thdrS,tempTrZ, ns, endian);
		statusR = segyReadTrace(finVel, &bh, &thdr,tempTrV, ns, endian);
		
		if (thdr.iline > iMax){
			printf("\n\nFound iline greater than iMax.\n\n");
			printf("\n\n iline: ", thdr.iline);
			return 0;
		}else if (thdr.iline < iOrigin){
			printf("\n\nFound iline less than iOrigin.\n\n");
			return 0;
		}else if (thdr.xline - xOrigin < 0){
			printf("\n\nFound xline  less than xOrigin.\n\n");
			return 0;
		}else if (thdr.xline > xMax){
			printf("\n\nFound xline greater than xMax.\n\n");
			return 0;
		}	
			
		/*if it gets to end of reading before having read in enough traces abort
		*else end o file and you have enough lines exit while loop */
		if(statusR!=0){
			if(thdr.iline - iOrigin < numTr -1 ){
				printf("\nNot enough ilines for a single execution calc!\n");
				return  0;
			}else{
					
				break; /*minimum number of lines met...exit loop*/
			}
		}
		    /*If the xline number read in for this inline is larger than the previously read in xline number, the
		     * this xline number becomes the new maximum saved as xNdxMx*/
			if ( (thdr.xline-xOrigin) > xNdxMx[thdr.iline-iOrigin]) xNdxMx[thdr.iline - iOrigin]=thdr.xline - xOrigin;
			if (ilineOld != thdr.iline){
				xNdxMn[thdr.iline-iOrigin]=thdr.xline-xOrigin;

				ilineOld=thdr.iline;
			}
	}

	/*copy the value read last, which caused while loop to exit, by the way 
	 * I could just compute ilineCnt from the difference between thdr.iline and iOrigin
	 * since this is a bit confusing*/
	ilineCnt=numTr/2;
	
	printf("\n\n ****************");
			time(&now);
			printf("\n Start Time: %s Percent Done: %d",  ctime(&now),0 );

	/*compute semblance and dips*/
	velMix(velTr, velOutTr, ns, numTr,ilineCnt,zeroTr,xNdxMn,xNdxMx,iOrigin,xOrigin,rate,thdr.iline,thdr.xline);
	
	/*outputing status*/
	num10Percent = (iMax-iOrigin + 1)/10;
	if((((*trGthHdr[numTr/2][xNdxMn[ilineCnt]]).iline) - iOrigin + 1)%num10Percent==0){
		printf("\n\n ****************");
		time(&now);
		printf("\n Finished Inline: %d Time: %s Percent Done: ", (*trGthHdr[numTr/2][xNdxMn[ilineCnt]]).iline, ctime(&now),((*trGthHdr[numTr/2][xNdxMn[ilineCnt]]).iline - iOrigin + 1)/num10Percent*10 );
	}

	/*spits out first numTr/2 +1  lines*/
	for (i=0;i<=numTr/2;++i){
		/*output first numTr traces for begining of inline*/
		for(j = xNdxMn[i]; j < xNdxMn[i] + numTr/2; ++j){
			if((*trGthHdr[i][j]).iline == 0){
				(*trGthHdr[i][j]).iline=(*trGthHdr[i][xNdxMn[ilineCnt]]).iline;
				(*trGthHdr[i][j]).xline=xOrigin + j;
			}
			segyWriteTrace(foutVel, trGthHdr[i][j], velOutTr[numTr/2][xNdxMn[i] + numTr/2], ns, endian);
		}
		
		//for (j = xNdxMn[i]; j <= xNdxMx[i]; ++j){
		for (j= xNdxMn[i]+numTr/2; j <= xNdxMx[i] - numTr/2; ++j){
			if((*trGthHdr[i][j]).iline == 0){
				(*trGthHdr[i][j]).iline=(*trGthHdr[i][xNdxMn[ilineCnt]]).iline;
				(*trGthHdr[i][j]).xline=xOrigin + j;
			}
			segyWriteTrace(foutVel, trGthHdr[i][j], velOutTr[numTr/2][j], ns, endian);
		}
	
		for(j = xNdxMx[i] - numTr/2 + 1; j <= xNdxMx[i]; ++j){
			if((*trGthHdr[i][j]).iline == 0){
				(*trGthHdr[i][j]).iline=(*trGthHdr[i][xNdxMn[ilineCnt]]).iline;
				(*trGthHdr[i][j]).xline=xOrigin + j;
			}
			segyWriteTrace(foutVel, trGthHdr[i][j], velOutTr[numTr/2][xNdxMx[i] - numTr/2], ns, endian);
		}
	}
	
   while (statusR == 0){ 
	   
		/*Shift all the traces back by one, meaning that
		* trGth[1] -> trGth[0], trGth[2] -> [1], and tempTr -> trGth[numTr-1]*/
		for (k =0 ;k<numTr-1;++k){
			for(j=0;j<xMax-xOrigin+1;++j){
				
				memcpy(zeroTr[k][j],zeroTr[k+1][j],(ns+1)*sizeof(float));
				memcpy(velTr[k][j],velTr[k+1][j],(ns+1)*sizeof(float));
				memcpy(trGthHdr[k][j],trGthHdr[k+1][j], sizeof(*trGthHdr[0][0]));
			}
		}

		/*setting last inline to zero before I read a new inline into it*/
		for (j =0; j < xMax-xOrigin +1;++j){
			
			memset(velTr[numTr-1][j],0,(ns+1)*sizeof(float));
			memset(trGthHdr[numTr-1][j],0,sizeof(*trGthHdr[0][0]));
			memset(zeroTr[numTr-1][j],0,(ns+1)*sizeof(float));
			memset(velOutTr[numTr/2][j],0,(ns+1)*sizeof(float));
		}
	
		//while(statusR ==0 && thdr.iline ==ilineOld){
		while (statusR==0 && thdr.iline ==ilineOld){	
			
			memcpy(zeroTr[numTr - 1][thdr.xline - xOrigin], tempTrZ, (ns + 1)*sizeof(float));
			memcpy(velTr[numTr - 1][thdr.xline - xOrigin], tempTrV, (ns + 1)*sizeof(float));
			memcpy(trGthHdr[numTr - 1][thdr.xline -xOrigin], &thdr, sizeof(*trGthHdr[0][0]));
			
			if ((thdr.xline - xOrigin) > xNdxMx[thdr.iline - iOrigin]) xNdxMx[thdr.iline - iOrigin] = thdr.xline - xOrigin;
			
			statusR = segyReadTrace(fzero, &bhS, &thdrS,tempTrZ, ns, endian);
			statusR = segyReadTrace(finVel, &bh, &thdr,tempTrV, ns, endian);
			
			if(statusR==0){
				if (thdr.iline > iMax){
					printf("\n\nFound iline greater than iMax.\n\n");
					return 0;
				}else if (thdr.iline < iOrigin){
					printf("\n\nFound iline less than iOrigin.\n\n");
					printf("\n\n iline: %d %d %d %d %d %d", thdr.iline, thdr.xline, status, xNdxMn[ilineCnt], xNdxMx[ilineCnt],ilineCnt);
					return 0;
				}else if (thdr.xline - xOrigin < 0){
					printf("\n\nFound xline  less than xOrigin.\n\n");
					return 0;
				}else if (thdr.xline > xMax){
					printf("\n\nFound xline greater than xMax.\n\n");
					return 0; 
				}	
			}		
		}

		if (statusR == 0){
			xNdxMn[thdr.iline - iOrigin] = thdr.xline - xOrigin;
			ilineOld = thdr.iline;
		}
			
		ilineCnt++;
		
		velMix(velTr, velOutTr, ns, numTr,ilineCnt,zeroTr,xNdxMn,xNdxMx,iOrigin,xOrigin,rate, thdr.iline,thdr.xline);
		
		if((((*trGthHdr[numTr/2][xNdxMn[ilineCnt]]).iline) - iOrigin + 1)%num10Percent==0){
			printf("\n\n ****************");
			time(&now);
			printf("\n Finished Inline: %d Time: %s Percent Done: %d", (*trGthHdr[numTr/2][xNdxMn[ilineCnt]]).iline, ctime(&now),((*trGthHdr[numTr/2][xNdxMn[ilineCnt]]).iline-iOrigin + 1)/num10Percent*10 );
			//printf("\\ Percent Done: 
			//printf("\n Time: %s", ctime(&now));
		}
		
		/*output current inline*/
		i=numTr/2;
		for(j = xNdxMn[ilineCnt]; j < xNdxMn[ilineCnt] + numTr/2; ++j){
			if((*trGthHdr[i][j]).iline == 0){
				(*trGthHdr[i][j]).iline=(*trGthHdr[i][xNdxMn[ilineCnt]]).iline;
				(*trGthHdr[i][j]).xline=xOrigin + j;
			}
			segyWriteTrace(foutVel, trGthHdr[i][j], velOutTr[numTr/2][xNdxMn[i] + numTr/2], ns, endian);
		}
		
		//for (j = xNdxMn[i]; j <= xNdxMx[i]; ++j){
		for (j= xNdxMn[ilineCnt]+numTr/2; j <= xNdxMx[ilineCnt] - numTr/2; ++j){
			if((*trGthHdr[i][j]).iline == 0){
				(*trGthHdr[i][j]).iline=(*trGthHdr[i][xNdxMn[ilineCnt]]).iline;
				(*trGthHdr[i][j]).xline=xOrigin + j;
			}
			segyWriteTrace(foutVel, trGthHdr[i][j], velOutTr[numTr/2][j], ns, endian);
		}
	
		for(j = xNdxMx[ilineCnt] - numTr/2 + 1; j <= xNdxMx[ilineCnt]; ++j){
			if((*trGthHdr[i][j]).iline == 0){
				(*trGthHdr[i][j]).iline=(*trGthHdr[i][xNdxMn[ilineCnt]]).iline;
				(*trGthHdr[i][j]).xline=xOrigin + j;
			}
			segyWriteTrace(foutVel, trGthHdr[i][j], velOutTr[numTr/2][xNdxMx[i] - numTr/2], ns, endian);
		}
	}
	
	for (i=1;i<=numTr/2;++i){
		/*output last numTr traces for begining of inline*/
		for(j = xNdxMn[ilineCnt + i]; j < xNdxMn[ilineCnt + i] + numTr/2; ++j){
			if((*trGthHdr[numTr/2 + i][j]).iline == 0){
				(*trGthHdr[numTr/2 + i][j]).iline=ilineCnt + i;
				(*trGthHdr[numTr/2 + i][j]).xline=xOrigin + j;
			}
			segyWriteTrace(foutVel, trGthHdr[numTr/2 + i][j], velOutTr[numTr/2][xNdxMn[ilineCnt + i ] + numTr/2], ns, endian);
		}
		
		//for (j = xNdxMn[i]; j <= xNdxMx[i]; ++j){
		for (j= xNdxMn[ilineCnt + i]+numTr/2; j <= xNdxMx[ilineCnt + i] - numTr/2; ++j){
			if((*trGthHdr[numTr/2 + i][j]).iline == 0){
				(*trGthHdr[ilineCnt + i][j]).iline=ilineCnt + i;
				(*trGthHdr[ilineCnt + i][j]).xline=xOrigin + j;
			}
			segyWriteTrace(foutVel, trGthHdr[numTr/2 + i][j], velOutTr[numTr/2][j], ns, endian);
		}
	
		for(j = xNdxMx[i] - numTr/2 + 1; j <= xNdxMx[i]; ++j){
			if((*trGthHdr[numTr/2 + i][j]).iline == 0){
				(*trGthHdr[ilineCnt + i][j]).iline=ilineCnt + i;
				(*trGthHdr[ilineCnt + i][j]).xline=xOrigin + j;
			}
			segyWriteTrace(foutVel, trGthHdr[numTr/2 + i][j], velOutTr[numTr/2][xNdxMx[ilineCnt + i] - numTr/2], ns, endian);
		}
	}


printf("\n\nDONE!\n\n");

free(zeroTr);
free(velTr);
free (velOutTr);
free(trGthHdr);
free(tempTrZ);
free(tempTrV);
free(xNdxMx);
free(xNdxMn);
return 0;
}
/*
************************************************************************
*/
void doMessage(char *str)
{
   fprintf(stderr, "%s\n", str);
}
/*
***********************************************************************
*/
/***********************************************************************
 * Function: Velocity Mixer
 * Input: Inline dip volume, Xline dip volume, semblance value volume,
 * 		  velocity volume, numTr, window, maxDip.s
 * 
 * Description:  Mixes velocities witing a radius of numTr along dips.
***********************************************************************/
void velMix(float ***velTr, float ***velOutTr,int ns, int numTr,int iNdx, float***zeroTr, int *xNdxMn, int *xNdxMx, int iOrigin, int xOrigin, int rate, int iline, int xline)
{
	int thread_id, nloops;	
//#pragma omp parallel private(thread_id, nloops)
//{	
	//velMix(velTr, velOutTr, ns, numTr,ilineCnt,window,maxDip,dipTrI,dipTrX,semTr,xNdxMn,xNdxMx,iOrigin,xOrigin,thresh,rate);
	int i,j,k,ii;
	int window =1;
	int samp = 0; //sample value
	int xStart=xNdxMn[iNdx];//-xOrigin;//index number of first xline.  if xNdxMn = xOrigin, then xStart =0 
	int xEnd = xNdxMx[iNdx];//-xOrigin; //xEnd indx number of last xline
	int xNdx; //current array value of xline being processed on this particular inline
	int shift = window/2 ; //shift is the starting point for the samples.
	int halfTr = numTr / 2; //should have called this radius
	int halfWin = window / 2; //hafwindow of samples - should hve called winradi
	float velSum=0; /*Current sum of velocities*/
	float sumVel; //sum of valid velocities
	int numVel; //number of valid velocities 

	//#pragma omp for
	/*walking down the inline.  Remeber, this will start index numTr/2*/
	for (xNdx=xStart+(halfTr);xNdx <= xEnd-(halfTr); ++xNdx){
		/*walking down the samples in the trace*/
			for (samp=0; samp<ns; ++samp){		
				numVel=0; 
				sumVel=0;		
				/*getting velocity at center of gather for each window*/
				if(zeroTr[halfTr][xNdx][samp] != 0){
				    sumVel+=velTr[halfTr][xNdx][samp];
				    ++numVel;
				  /* printf("\n\n value %f %d %d %d", zeroTr[halfTr][xNdx][samp],iline,xline,samp);
				    getchar();
				    getchar(); */
				/*going UP iline*/
				i=1;
				while (i<=halfTr && zeroTr[halfTr][xNdx+i][samp] !=0){
					sumVel+=velTr[halfTr][xNdx+i][samp];
					++numVel;
					++i;
				}
				
				/*going DOWN iline*/
				i=1;
				while (i<=halfTr && zeroTr[halfTr][xNdx-i][samp] !=0){
					sumVel+=velTr[halfTr][xNdx-i][samp];
					++numVel;
					++i;
				}
				
				/*going Right xline*/
				i=1;
				while (i<=halfTr && zeroTr[halfTr+i][xNdx][samp] !=0){
					sumVel+=velTr[halfTr+i][xNdx][samp];
					++numVel;
					++i;
				}
				/*going LEFT xline*/
				i=1;
				while (i<=halfTr && zeroTr[halfTr-i][xNdx][samp] !=0){
					sumVel+=velTr[halfTr-i][xNdx][samp];
					++numVel;
					++i;
				}
				/*going UP and RIGHT xline*/
				i=1;
				while (i<=halfTr && zeroTr[halfTr+i][xNdx+i][samp] !=0){
					sumVel+=velTr[halfTr+i][xNdx+i][samp];
					++numVel;
					++i;
				}
				/*going UP and LEFT xline*/
				i=1;
				while (i<=halfTr && zeroTr[halfTr-i][xNdx+i][samp] !=0){
					sumVel+=velTr[halfTr-i][xNdx+i][samp];
					++numVel;
					++i;
				}
				/*going DOWN and RIGHT xline*/
				i=1;
				while (i<=halfTr && zeroTr[halfTr+i][xNdx-i][samp] !=0){
					sumVel+=velTr[halfTr+i][xNdx-i][samp];
					++numVel;
					++i;
				}
				/*going DOWN and LEFT xline*/
				i=1;
				while (i<=halfTr && zeroTr[halfTr-i][xNdx-i][samp] !=0){
					sumVel+=velTr[halfTr-i][xNdx-i][samp];
					++numVel;
					++i;
				}
				
			} else{
				sumVel+=velTr[halfTr][xNdx][samp];
				++numVel;
			}
			if (numVel >0){
				velOutTr[halfTr][xNdx][samp]=sumVel/(numVel);
		    }
	}
}
//}
}
/***********************************************************************
 * Function: Checks slope for conflicting dips
 * Input: slopeOld - previously calculated slope, slopeIXl - most 
 * 	      recently computed slope, maxDip - largest acceptable value
 * 	      for the absolute value of the difference in the dips 
 * 		  
 * 
 * Description:  If the the dips are moving in different directions and
 * 	             are larger than specified value, return 0 else 1.
***********************************************************************/
int slopeChk( float slopeOld, float slopeIXl, float maxDip){
	if(slopeOld*slopeIXl < 0 &&  fabs(slopeOld -slopeIXl) > maxDip){
		return 0;
	}else {
		return 1;
	}
}

/***********************************************************************
 * Function: Nerest - finds nearest interger
 * Input: input- input floating point value
 * 		  
 * 
 * Description:  Adds +.5 if input is positve, subtracts .5 if input is
 * 	             negative.  Then performs integer truncation!
***********************************************************************/
int nearest (float input){
	if(input < 0){
		input = (int)(input - 0.5000);
	}else{
		input = (int)(input + 0.5000);
	}

return input;
}

int boundCheck (int sampPlusWin, int slope, int ns){
	if ((sampPlusWin - slope) >=0 && (sampPlusWin - slope) < ns){
		return 1;
	}else{
		return 0;
	}
}

	

