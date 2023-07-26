/*************************************************************************

   Program:    indirectrepeats
   File:       indirectrepeats.c
   
   Version:    V2.0
   Date:       07.03.13
   Function:   Identify indirect repeats in a FASTA file
   
   Copyright:  (c) UCL / Dr. Andrew C. R. Martin 2004-2013
   Author:     Dr. Andrew C. R. Martin
   Address:    Biomolecular Structure & Modelling Unit,
               Institute of Structural & Molecular Biology,
               University College,
               Gower Street,
               London.
               WC1E 6BT.
   EMail:      andrew.martin@ucl.ac.uk
               andrew@bioinf.org.uk
               
**************************************************************************

   This program is not in the public domain, but it may be copied
   according to the conditions laid out in the accompanying file
   COPYING.DOC

   The code may be modified as required, but any modifications must be
   documented so that the person responsible can be identified. If someone
   else breaks this code, I don't want to be blamed for code that does not
   work! 

   The code may not be sold commercially or included as part of a 
   commercial product except as described in the file COPYING.DOC.

**************************************************************************

   Description:
   ============

**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================
   V1.0   01.10.04 Original Perl implementation for direct repeats
   V1.1   19.02.13 Perl version modified for indirect repeats 
                   By: Rashmi Rajasabhai
   V2.0   07.03.13 C version   By: ACRM

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "bioplib/general.h"
#include "bioplib/macros.h"

/************************************************************************/
/* Defines and macros
*/
#define MAXBUFF    512
#define MAXAA      24
#define MAXPATLEN  360
#define MAXSEQ     100000

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/
int main(int argc, char **argv);
void SearchAllPatterns(FILE *in, BOOL exact, int minpat, int maxpat, 
                       BOOL verbose, BOOL quiet);
void SearchFileForPattern(char *pattern, FILE *in, BOOL exact, 
                          BOOL verbose, BOOL quiet);
BOOL CheckBounds(char *sequence, char *pattern, int offset);
int SearchSequenceForPattern(char *sequence, char *pattern, int offset);
BOOL GetFASTASequence(FILE *in, char *label, char *sequence);
void Usage(void);
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile,
                  int *minpat, int *maxpat, char *pattern, BOOL *verbose,
                  BOOL *quiet, BOOL *exact);

/************************************************************************/
int main(int argc, char **argv)
{
   char pattern[MAXPATLEN],
        InFile[MAXBUFF], OutFile[MAXBUFF];
   BOOL exact = TRUE,
        verbose = FALSE,
        quiet = FALSE;
   int  maxpat = 10, 
        minpat = 1;
   FILE *in  = stdin,
        *out = stdout;

   if(ParseCmdLine(argc, argv, InFile, OutFile, &minpat, &maxpat, 
                   pattern, &verbose, &quiet, &exact))
   {
      if(OpenStdFiles(InFile, OutFile, &in, &out))
      {
         if(in == stdin)
         {
            Usage();
         }
         else
         {
            if(pattern[0])
            {
               SearchFileForPattern(pattern, in, exact, verbose, quiet);
            }
            else
            {
               SearchAllPatterns(in, exact, minpat, maxpat, 
                                 verbose, quiet);
            }
         }
      }
   }
   else
   {
      Usage();
   }
   
   return(0);
}

/************************************************************************/
void SearchAllPatterns(FILE *in, BOOL exact, int minpat, int maxpat, 
                       BOOL verbose, BOOL quiet)
{
   char aa;
   int  i, j, k;
   char *letters = "ACDEFGHIKLMNPQRSTVWY",
      pat[MAXPATLEN];

   for(j=0; j<strlen(letters); j++)
   {
      aa = letters[j];
      for(i=minpat; i<=maxpat; i++)
      {
         for(k=0; k<i; k++)
         {
            pat[k*2] = aa;
            pat[k*2+1] = 'X';
         }
         pat[i*2-1] = '\0';

         fprintf(stdout, "Testing pattern '%s':\n", pat);
         fflush(stdout);
         
         SearchFileForPattern(pat, in, exact, verbose, quiet);
        }
    }
}

/************************************************************************/
void SearchFileForPattern(char *pattern, FILE *in, BOOL exact, 
                          BOOL verbose, BOOL quiet)
{
   char label[MAXBUFF];
   char sequence[MAXSEQ];
   BOOL ok, print;
   int count, offset, nseq=0;
   

   rewind(in);

   count = 0;
   while(1)
   {
      GetFASTASequence(in, label, sequence);
      if (label[0] == '\0') break;
      if(!quiet)
      {
         if(!(++nseq % 10000))
         {
            fprintf(stderr, "Processed %d sequences\n", nseq);
            fflush(stderr);
         }
      }
      
      ok     = FALSE;
      offset = 0;
      print  = FALSE;

      while((offset=SearchSequenceForPattern(sequence, pattern, offset))
            !=(-1))
      {
         ok = TRUE;
         if(exact)
         {
            ok = CheckBounds(sequence, pattern, offset);
         }
         if(ok)
         {
            count++;
            print = TRUE;
         }
         offset++;
      }
      if(verbose && print)
      {
         fprintf(stdout, "%s matches\n", label);
      }
   }
   fprintf(stdout, "Total matches: %d\n", count);
}

/************************************************************************/
/* takes a sequence, a pattern and the offset into the sequence
   where the pattern was found. Checks if this is a sub pattern
   i.e. the pattern continues before or after the identified 
   place
*/
BOOL CheckBounds(char *sequence, char *pattern, int offset)
{
   char ch;
   int  patlen, seqlen;

   ch     = pattern[0];
   patlen = strlen(pattern);
   seqlen = strlen(sequence);

   /* Check the N terminus                                              */
   if(offset >= 2)
   {
      if((sequence[offset-1] != ch) &&
         (sequence[offset-2] == ch))
      {
         return(FALSE);
      }
   }
   else if(offset == 1)
   {
      if(sequence[offset-1] == ch)
      {
         return(FALSE);
      }
   }
   
   /* Check the C terminus                                              */
   if((seqlen - (offset+patlen)) >= 2)
   {
      if((sequence[offset+patlen]   != ch) &&
         (sequence[offset+patlen+1] == ch))
      {
         return(FALSE);
      }
   }
   else if ((seqlen - (offset+patlen)) == 1)
   {
      if(sequence[offset+patlen] == ch)
      {
         return(FALSE);
      }
   }
   

   return(TRUE);
}

/************************************************************************/
/* Takes a sequence and pattern of the form AXAXAXA together with an
   offset into the sequence. It then searches for the pattern starting
   from the specified offset.
   If the pattern is found then it returns the new offset where the
   pattern was found.
   Returns (-1) when pattern not found
*/
int SearchSequenceForPattern(char *sequence, char *pattern, int offset)
{
   char *subseq, ch;
   int  patlen, nrepeat, seqlen, i, j;
   BOOL match;
   
   subseq = sequence+offset;
   seqlen  = strlen(subseq);

   patlen = strlen(pattern);
   nrepeat = (patlen-1)/2;
   ch      = pattern[0];
   
   for(i=0; i<(1+seqlen-patlen); i++)
   {
      match = TRUE;
      for(j=0; j<patlen; j+=2)
      {
         if(subseq[i+j] != ch)
         {
            match = FALSE;
            break;
         }
      }
      if(match)
      {
         for(j=1; j<patlen-1; j+=2)
         {
            if(subseq[i+j] == ch)
            {
               match = FALSE;
               break;
            }
         }
      }
      
      if(match)
      {
         return(offset+i);
      }
   }

   return(-1);
}

/************************************************************************/
BOOL GetFASTASequence(FILE *in, char *label, char *sequence)
{
   static char labelline[MAXBUFF] = "";
   char buffer[MAXBUFF];

   label[0] = '\0';
   sequence[0] = '\0';
   
   if(feof(in))
   {
      labelline[0] = '\0';
      return(FALSE);
   }
   
   while(fgets(buffer, MAXBUFF, in))
   {
      TERMINATE(buffer);
      if(buffer[0] == '>')
      {
         if(strlen(labelline))
         {
            strncpy(label, labelline, MAXBUFF);
            strncpy(labelline, buffer, MAXBUFF);
            UPPER(sequence);
            return(TRUE);
         }
         else
         {
            strncpy(labelline, buffer, MAXBUFF);
         }
      }
      else
      {
         strncat(sequence, buffer, MAXSEQ-strlen(sequence));
      }
   }
   strncpy(label, labelline, MAXBUFF);
   labelline[0] = '\0';
   return(TRUE);
}

/************************************************************************/
void Usage(void)
{
   fprintf(stderr,"indirectrepeats V2.0, (c) 2004-2013, \
Dr. Andrew C.R. Martin, UCL\n");
   fprintf(stderr,"\n");
   fprintf(stderr,"Usage: indirectrepeats [-x][-v][-q][-n minpat]\
[-m maxpat][-s pattern] file.faa [output]\n");
   fprintf(stderr,"       -x Do non-exact matching\n");
   fprintf(stderr,"       -v Verbose (report macthed sequences)\n");
   fprintf(stderr,"       -q Quiet (do not report progress)\n");
   fprintf(stderr,"       -n Minimum pattern length (default: 1)\n");
   fprintf(stderr,"       -m Maxmimum pattern length (default: 10)\n");
   fprintf(stderr,"       -s Specify a sequence pattern\n");
   fprintf(stderr,"\n");
   fprintf(stderr,"Searches a file for matches to a poly-amino acid \
sequence. By default,\n");
   fprintf(stderr,"works through each amino acid in turn and tests \
occurrences of n-mers\n");
   fprintf(stderr,"from 1..10 residues. You can override this with -s \
and just look\n");
   fprintf(stderr,"at the one specified pattern.\n");
   fprintf(stderr,"\n");
   fprintf(stderr,"By default, does exact matching. In other words, if \
you are looking for\n");
   fprintf(stderr,"the pattern 'AXA', then 'AXAXA' will NOT match. You \
can override this\n");
   fprintf(stderr,"with -x. Exact matching simply checks that the \
residues before and\n");
   fprintf(stderr,"after the pattern do not extend the pattern.\n");
   fprintf(stderr,"\n");

   exit(0);
}

/************************************************************************/
/*>BOOL ParseCmdLine()
   ----------------------------------------------------------------
   Input:   int    argc         Argument count
            char   **argv       Argument array
   Output:  
            
            char   *infile      Input file (or blank string)
            char   *outfile     Output file (or blank string)
            int    *minpat      
            int    *maxpat      
            char   *pattern
            BOOL   *verbose
            BOOL   *quiet
            BOOL   *exact
   Returns: BOOL                Success?

   Parse the command line

   01.06.09  Original   By: ACRM   
*/
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile,
                  int *minpat, int *maxpat, char *pattern, BOOL *verbose,
                  BOOL *quiet, BOOL *exact)
{
   argc--;
   argv++;

   infile[0] = outfile[0] = '\0';
   pattern[0] = '\0';
   *maxpat = 10;
   *minpat = 1;
   *verbose = FALSE;
   *exact = TRUE;
   *quiet = FALSE;
   
   while(argc)
   {
      if(argv[0][0] == '-')
      {
         if (argv [0][2]!='\0')
         {
           return(FALSE);
         }
         else
         {            
            switch(argv[0][1])
            {
            case 'h':
               return(FALSE);
               break;
            case 'v':
               *verbose = TRUE;
               break;
            case 'x':
               *exact = FALSE;
               break;
            case 'q':
               *quiet = TRUE;
               break;
            case 's':
               argv++;
               argc--;
               strncpy(pattern, argv[0], MAXPATLEN);
               break;
            case 'n':
               argv++;
               argc--;
               if(!sscanf(argv[0], "%d", minpat))
                  return(FALSE);
               break;
            case 'm':
               argv++;
               argc--;
               if(!sscanf(argv[0], "%d", maxpat))
                  return(FALSE);
               break;
            default:
               return(FALSE);
               break;
            }
         }         
      }
      else
      {
         /* Check that there are 1 or 2 arguments left                  */
         if(argc > 2)
            return(FALSE);
         
         /* Copy the first to infile                                    */
         if(argc)
         {
            strcpy(infile, argv[0]);
            argc--;
            argv++;
         }
         
         /* If there's another, copy it to outfile                      */
         if(argc)
         {
            strcpy(outfile, argv[0]);
            argc--;
            argv++;
         }
         
         return(TRUE);
      }
      
      argc--;
      argv++;
   }
   return(TRUE);
}
