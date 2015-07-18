/*************************************************************************

   Program:    chisq
   File:       chisq.c
   
   Version:    V1.2
   Date:       09.02.94
   Function:   Do statistical analysis of seqan output
   
   Copyright:  (c) Dr. Andrew C. R. Martin, 1994
   Author:     Dr. Andrew C. R. Martin
   Address:    Biomolecular Structure and Modelling,
               University College London,
   Phone:      HOME: +44 (0)372 275775
   EMail:      JANET:    martin@uk.ac.ucl.bioc.bsm
               INTERNET: martin@bsm.bioc.ucl.ac.uk
                         amartin@scitec.adsp.sub.org
               UUCP:     ...cbmehq!cbmuk!scitec!amartin
               
**************************************************************************

   This program is not in the public domain, but it may be freely copied
   and distributed for no charge providing this header is included.
   The code may be modified as required, but any modifications must be
   documented so that the person responsible can be identified. If someone
   else breaks this code, I don't want to be blamed for code that does not
   work! The code may not be sold commercially without prior permission 
   from the author, although it may be given away free with commercial 
   products, providing it is made clear that this program is free and that 
   the source code is provided with the program.

**************************************************************************

   Description:
   ============

**************************************************************************

   Usage:
   ======

**************************************************************************

   Notes:
   ======
   
**************************************************************************

   Revision History:
   =================
   V1.0  03.02.94 Original
   V1.1  04.02.94 Fixed bug in ClearArray()
   V1.2  09.02.94 Added ChiSq calculation of individual data items


*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

/***********************************************************************/
/* Defines
*/
#ifndef TRUE
#define TRUE 1
#define FALSE 0
#endif

#define MAXAA      21
#define MINBIN     10
#define BIN_FIRST   1
#define BIN_SECOND  2
#define SMALL       ((double)1e-10)
#define TERMINATE(x) {                                            \
                         int i;                                   \
                         for(i=0; (x)[i]; i++)                    \
                         {                                        \
                            if((x)[i] == '\n') (x)[i] = '\0';     \
                            break;                                \
                         }                                        \
                      }


/***********************************************************************/
/* Type definitions
*/
typedef int BOOL;

/***********************************************************************/
/* Globals
*/
int  gData[MAXAA][MAXAA];
char gAAtab[]    = "ACDEFGHIKLMNPQRSTVWYB"; /* B is used for the bin   */
BOOL gWide       = FALSE,
     gIndividual = FALSE;
int  gMinBin     = MINBIN;

/***********************************************************************/
/* Prototypes
*/
int main(int argc, char **argv);
BOOL Initialise(void);
BOOL ParseCmdLine(int argc, char **argv, char *filename);
void Usage(void);
BOOL ProcessExample(FILE *fp);
void ClearArray(void);
void StoreData(char *buffer);
BOOL Lookup(char First, char Second, int *pos1, int *pos2);
void ProcessData(void);
char LookDown(int pos);
void PrintObsExpTable(double Expected[MAXAA][MAXAA]);
void CalcExpected(int *FirstTotal, int *SecondTotal, int NObs, 
                  double Expected[MAXAA][MAXAA]);
void BinResidues(int flag, int *totals);
void ShowTotals(int *FirstTotal, int *SecondTotal);
void PrintHeader(void);

/***********************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
   ChiSq main program

   03.02.94 Original   By: ACRM
*/
int main(int argc, char **argv)
{
   char filename[160];
   FILE *fp;
   
   if(Initialise())
   {
      PrintHeader();
      
      if(ParseCmdLine(argc, argv, filename))
      {
         if(!filename[0])
            fp = stdin;
         else
            fp=fopen(filename,"r");
         
         if(fp==NULL)
         {
            fprintf(stderr,"Unable to open input file %s\n",filename);
            exit(1);
         }
         else
         {
            while(ProcessExample(fp)) ;
         }
      }
      else
      {
         Usage();
      }
   }
   else
   {
      fprintf(stderr,"Unable to initialise\n");
   }
}

/***********************************************************************/
/*>BOOL Initialise(void)
   ---------------------
   Initialisation (dummy)

   03.02.94 Original   By: ACRM
*/
BOOL Initialise(void)
{
   return(TRUE);
}

/***********************************************************************/
/*>BOOL ParseCmdLine(int argc, char **argv, char *filename)
   --------------------------------------------------------
   Parse the command line getting switches and the filename.

   03.02.94 Original   By: ACRM
   09.02.94 Added -i
*/
BOOL ParseCmdLine(int argc, char **argv, char *filename)
{
   argc--;
   argv++;

   filename[0] = '\0';
      
   while(argc>0)
   {
      if(argv[0][0] == '-')
      {
         switch(argv[0][1])
         {
         case 'w':
            gWide = TRUE;
            break;
         case 'm':
            argv++; argc--;
            gMinBin = atoi(argv[0]);
            break;
         case 'i':
            gIndividual = TRUE;
            gWide       = TRUE;
            break;
         case 'h': case '?':
            Usage();
            exit(0);
         default:
            break;
         }
      }
      else
      {
         if(argc==1)
         {
            strcpy(filename, argv[0]);
            return(TRUE);
         }
         else
         {
            return(FALSE);
         }
      }

      argc--;
      argv++;
   }
}

/***********************************************************************/
/*>void Usage(void)
   ----------------
   Print usage message 

   03.02.94 Original   By: ACRM
   09.02.94 Added -i
*/
void Usage(void)
{
   printf("chisq V1.0 - A program to calculate Chi Squared from output of \
seqan\n");
   printf("Usage: chisq [-w] [-m <min>] [-h] [file.in]\n");
   printf("If an input file is not specified, input is read from stdin\n");
   printf("       -w Print results in wide format\n");
   printf("       -m Specify max frequency for binning (default: %d)\n",
          MINBIN);
   printf("       -i Show ChiSq on each item of data\n");
   printf("       -h/-? This help message\n");
}

/***********************************************************************/
/*>BOOL ProcessExample(FILE *fp)
   -----------------------------
   Read data from ffile for a single residue pair and call processing
   routines

   03.02.94 Original   By: ACRM
   09.02.94 Added separator line
*/
BOOL ProcessExample(FILE *fp)
{
   static BOOL FirstCall = TRUE;
   static char buffer[80];
      
   ClearArray();
   
   /* If not the first call, then display the buffer from the last go    */
   if(!FirstCall)
   {
      if(gWide)
         fprintf(stdout,"\n========================================\
=====================================================================\
======================\n");
      else
         fprintf(stdout,"\n========================================\
========================================\n");
      fprintf(stdout,"%s\n",buffer);
   }

   while(fgets(buffer,80,fp))
   {
      TERMINATE(buffer);

      if(!strncmp(buffer,"Pair:",5))
      {
         if(FirstCall)
         {
            /* If it's the first call, this is the first time a Pair line
               has been seen, so display it.
            */
            fprintf(stdout,"\n\n%s\n",buffer);
            FirstCall=FALSE;
         }
         else
         {
            /* Not the first call and a Pair line seen, so it's the end
               of this dataset.
            */
            ProcessData();
            return(TRUE);
         }
      }
      else if(buffer[2]==':')
      {
         /* It's a data line                                              */
         StoreData(buffer);
      }
   }
   
   ProcessData();
   return(FALSE);
}

/***********************************************************************/
/*>void ClearArray(void)
   ---------------------
   Clear the data array

   03.02.94 Original   By: ACRM
   04.02.94 Corrected count to MAXAA rather than 20
*/
void ClearArray(void)
{
   int i,j;
   
   for(i=0; i<MAXAA; i++)
      for(j=0; j<MAXAA; j++)
         gData[i][j] = 0;
}

/***********************************************************************/
/*>void StoreData(char *buffer)
   ----------------------------
   Store a line of data in the data array
   03.02.94 Original   By: ACRM
*/
void StoreData(char *buffer)
{
   char First,
        Second;
   int  i,
        ndata,
        pos1,
        pos2;

   First  = buffer[0];
   Second = buffer[1];

   for(i=0; i<strlen(buffer); i++)
   {
      if(buffer[i] == ',')
      {
         buffer[i] = '\0';
         break;
      }
   }
   
   ndata = atoi(buffer+3);

   if(Lookup(First,Second,&pos1,&pos2))
      gData[pos1][pos2] = ndata;
}

/***********************************************************************/
/*>BOOL Lookup(char First, char Second, int *pos1, int *pos2)
   ----------------------------------------------------------
   Look up the array positions for a residue pair

   03.02.94 Original   By: ACRM
*/
BOOL Lookup(char First, char Second, int *pos1, int *pos2)
{
   int i;
   
   *pos1 = *pos2 = (-1);

   for(i=0; gAAtab[i]; i++)
   {
      if(gAAtab[i]==First)  *pos1 = i;
      if(gAAtab[i]==Second) *pos2 = i;
   }
   
   if(*pos1 == (-1) || *pos2 == (-1))
      return(FALSE);
   
   return(TRUE);
}

/***********************************************************************/
/*>char LookDown(int pos)
   ----------------------
   Look up the residue type for a given array position

   03.02.94 Original   By: ACRM
*/
char LookDown(int pos)
{
   return(gAAtab[pos]);
}

/***********************************************************************/
/*>void ProcessData(void)
   ----------------------
   Process the data. Calculate expected values, perform binning and 
   recalcultlate. Calculate ChiSq

   03.02.94 Original   By: ACRM
*/
void ProcessData(void)
{
   int    FirstTotal[MAXAA],
          SecondTotal[MAXAA],
          NObs,
          NDoF,
          rows,
          cols, 
          i,
          j;
   
   double Expected[MAXAA][MAXAA],
          ChiSq;

   /* Clear the totals arrays                                          */
   for(i=0; i<MAXAA; i++)
   {
      FirstTotal[i]  = 0;
      SecondTotal[i] = 0;
      for(j=0; j<MAXAA; j++)
         Expected[i][j] = 0.0;
   }

   /* Sum the residue occurences for first position                   */
   for(i=0; i<MAXAA; i++)
      for(j=0; j<MAXAA; j++)
         FirstTotal[i] += gData[i][j];
      
   /* Sum the residue occurences for second position                  */
   for(j=0; j<MAXAA; j++)
      for(i=0; i<MAXAA; i++)
         SecondTotal[j] += gData[i][j];

   /* Calculate total number of observations                          */
   for(i=0, NObs=0; i<MAXAA; i++)
      NObs += FirstTotal[i];

   /* Calculate all the Expected values                              */
   CalcExpected(FirstTotal, SecondTotal, NObs, Expected);

   /* Print raw results                                              */
   printf("Raw results:\n============\n\n");
   printf("Number of observations: %d\n",NObs);
   ShowTotals(FirstTotal, SecondTotal);
   PrintObsExpTable(Expected);
   
   /* Now move all residues with <gMinBin occurences into the bins     */
   printf("\nThe following residues at the first position are now \
grouped:\n");
   BinResidues(BIN_FIRST,FirstTotal);
   printf("\nThe following residues at the second position are now \
grouped:\n");
   BinResidues(BIN_SECOND,SecondTotal);

   /* Recalculate all expected values                                  */
   CalcExpected(FirstTotal, SecondTotal, NObs, Expected);
   
   /* Show the binned results                                          */
   printf("\n\nBinned results:\n===============\n");
   ShowTotals(FirstTotal, SecondTotal);
   PrintObsExpTable(Expected);

   /* Calculate the ChiSq value.                                       */
   ChiSq = 0.0;
   for(i=0; i<MAXAA; i++)
   {
      for(j=0; j<MAXAA; j++)
      {
         if(FirstTotal[i] && SecondTotal[j])
            ChiSq += ((gData[i][j] - Expected[i][j]) * 
                      (gData[i][j] - Expected[i][j]) / 
                      Expected[i][j]);
      }
   }

   /* Calculate the number of degrees of freedom                      */
   for(i=0,rows=0; i<MAXAA; i++)
      if(FirstTotal[i]) rows++;
   for(i=0,cols=0; i<MAXAA; i++)
      if(SecondTotal[i]) cols++;
   NDoF = (rows-1) * (cols-1);

   /* Display these values                                            */
   printf("Chi Squared = %lf with %d degrees of freedom\n\n",ChiSq,NDoF);

}

/***********************************************************************/
/*>void ShowTotals(int *FirstTotal, int *SecondTotal)
   --------------------------------------------------
   Display total occurences of residue types

   03.02.94 Original   By: ACRM
*/
void ShowTotals(int *FirstTotal, int *SecondTotal)
{
   int i;
   
   printf("\nTotals at first position:\n=========================\n");
   for(i=0;i<MAXAA;i++)
      if(FirstTotal[i]) 
         printf("%c: %d\n",LookDown(i),FirstTotal[i]);

   printf("\nTotals at second position:\n==========================\n");
   for(i=0;i<MAXAA;i++)
      if(SecondTotal[i]) 
         printf("%c: %d\n",LookDown(i),SecondTotal[i]);
}

/***********************************************************************/
/*>void PrintObsExpTable(double Expected[MAXAA][MAXAA])
   ----------------------------------------------------
   Print table of observed and expected values

   03.02.94 Original   By: ACRM
   09.02.94 Added printing of individual ChiSq values
*/
void PrintObsExpTable(double Expected[MAXAA][MAXAA])
{
   int i,
       j;
   
   printf("\nObserved & expected values:\n===========================\n");
   if(gWide) printf("        A     C     D     E     F     G     H     \
I     K     L     M     N     P     Q     R     S     T     V     W     \
Y     B\n");
   else      printf("     A  C  D  E  F  G  H  I  K  L  M  N  P  Q  R  \
S  T  V  W  Y  B\n");
   for(i=0; i<MAXAA; i++)
   {
      printf("%c  ",LookDown(i));
      
      for(j=0; j<MAXAA; j++)
      {
         if(gWide) printf("%6d",gData[i][j]);
         else      printf("%3d",gData[i][j]);
      }
      printf("\n   ");

      for(j=0; j<MAXAA; j++)
      {
         if(gWide) printf("%6.1lf",Expected[i][j]);
         else      printf("%3d",(int)Expected[i][j]);
      }
      printf("\n   ");

      if(gIndividual)
      {
         for(j=0; j<MAXAA; j++)
         {
            double ChiSq;

            if(Expected[i][j] > SMALL)
            {
               ChiSq = (gData[i][j] - Expected[i][j]) *
                       (gData[i][j] - Expected[i][j]) /
                       Expected[i][j];
            
               printf("%6.1lf",ChiSq);
            }
            else
            {
               printf("   ---");
            }
         }
         printf("\n");
      }
      printf("\n");
   }
}


/***********************************************************************/
/*>void CalcExpected(int *FirstTotal, int *SecondTotal, int NObs, 
                     double Expected[MAXAA][MAXAA])
   ---------------------------------------------------------------
   Calculate the expected occurences

   03.02.94 Original   By: ACRM
*/
void CalcExpected(int *FirstTotal, int *SecondTotal, int NObs, 
                  double Expected[MAXAA][MAXAA])
{
   int i,
       j;
   
   for(i=0; i<MAXAA; i++)
   {
      for(j=0; j<MAXAA; j++)
      {
         Expected[i][j] = (double)FirstTotal[i]  * 
                          (double)SecondTotal[j] / 
                          (double)NObs;
      }
   }
}

/***********************************************************************/
/*>void BinResidues(int flag, int *totals)
   ---------------------------------------
   Place low frequency residues into a single bin

   03.02.94 Original   By: ACRM
*/
void BinResidues(int flag, int *totals)
{
   int i,
       j;
   
   for(i=0; i<MAXAA-1; i++)
   {
      if(totals[i] && totals[i] < gMinBin)
      {
         printf("%c ",LookDown(i));
         
         for(j=0; j<MAXAA; j++)
         {
            if(flag == BIN_FIRST)
            {
               gData[MAXAA-1][j] += gData[i][j];
               gData[i][j] = 0;
            }
            else
            {
               gData[j][MAXAA-1] += gData[j][i];
               gData[j][i] = 0;
            }
         }

         totals[MAXAA-1] += totals[i];
         totals[i] = 0;
      }
   }
   printf("\n");
}


/***********************************************************************/
/*>void PrintHeader(void)
   ----------------------
   Print a header.

   03.02.94 Original   By: ACRM
   09.02.94 Ammended for -i option
*/
void PrintHeader(void)
{
   printf("chisq V1.2 (c) 1994 Dr. Andrew C.R. Martin, UCL\n\n");
   printf("Takes results from seqan analysis and prints a contingency \
table containing\n");
   printf("observed and expected values for each residue pair together \
with the totals\n");
   printf("for each observed residue type. All residues occuring fewer \
than a minumum\n"); 
   printf("number of times (default: %d) are then merged into a single \
bin and the\n",MINBIN);
   printf("observed and expected values are reprinted using the bins. \
A value for\n");
   printf("Chi squared is then calculated and displayed with the number \
of degrees\n");
   printf("of freedom. Note that this value must be interpreted with \
care if any of\n");
   printf("the expected values are less than 5.0\n\n");
   printf("If the -i option has been specified, the third row of each \
entry shows\n");
   printf("the ChiSq value for this individual piece of data. This will \
have 1DoF\n");
}
