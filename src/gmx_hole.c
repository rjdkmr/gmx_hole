/*
 * This file is part of gmx_hole
 *
 * Author: Rajendra Kumar
 * Copyright (C) 2014  Rajendra Kumar
 *
 * gmx_hole uses hole program for which documentation is available in the
 * following link:
 * http://www.csb.yale.edu/userguides/graphics/hole/doc/hole_d00.html
 * Please cite the original publication of hole:
 * O.S. Smart, J.M. Goodfellow and B.A. Wallace (1993)
 * The Pore Dimensions of Gramicidin A
 * Biophysical Journal 65:2455-2460.
 *
 * gmx_hole is a free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * gmx_hole is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with gmx_hole.  If not, see <http://www.gnu.org/licenses/>.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
 * TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *
 */

#include "typedefs.h"
#include "statutil.h"
#include "smalloc.h"
#include "do_fit.h"
#include "pbc.h"
#include "string2.h"
#include "copyrite.h"
#include "gmx_fatal.h"
#include "index.h"
#include "tpxio.h"
#include "rmpbc.h"
#include "vec.h"

#include "ExtractData.h"

void CopyRightMsg() {

    char *copyright[] = {
            "                                                                        ",
            "                  :-)  gmx_hole (-:                                     ",
            "                                                                        ",
            "               Author: Rajendra Kumar                                  ",
            "                                                                        ",
            "         Copyright (C) 2014  Rajendra Kumar                             ",
            "                                                                        ",
            "gmx_hole uses hole program for which documentation is available in the ",
            "following link:                                                         ",
            "http://www.csb.yale.edu/userguides/graphics/hole/doc/hole_d00.html      ",
            "Please cite the original publication of hole:                          ",
            "O.S. Smart, J.M. Goodfellow and B.A. Wallace (1993)                     ",
            "The Pore Dimensions of Gramicidin A                                     ",
            "Biophysical Journal 65:2455-2460.                                       ",
            "                                                                        ",
            "gmx_hole is a free software: you can redistribute it and/or modify      ",
            "it under the terms of the GNU General Public License as published by    ",
            "the Free Software Foundation, either version 3 of the License, or       ",
            "(at your option) any later version.                                     ",
            "                                                                        ",
            "gmx_hole is distributed in the hope that it will be useful,             ",
            "but WITHOUT ANY WARRANTY; without even the implied warranty of          ",
            "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           ",
            "GNU General Public License for more details.                            ",
            "                                                                        ",
            "You should have received a copy of the GNU General Public License       ",
            "along with gmx_hole.  If not, see <http://www.gnu.org/licenses/>.       ",
            "                                                                        ",
            "THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS     ",
            "\"AS IS\" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT     ",
            "LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR   ",
            "A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT    ",
            "OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,   ",
            "SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED",
            "TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR  ",
            "PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF  ",
            "LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING    ",
            "NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS      ",
            "SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.            ",
            "                                                                        ",
            "                           :-)  gmx_hole (-:                            ",
            "                                                                        ",
            "                                                                        "
    };
    int i = 0;
    char *str;
    for(i=0; i<43; i++) {
        str = strdup(copyright[i]);
        fprintf(stderr,"%s\n", str);
    }
}

int cat_pdb(int nframe, char *fn_input, FILE *fOut)	{
	int i = 0;
	int number=0;
	char **data=NULL;


	FILE *f_input;
	f_input = fopen(fn_input, "r");
	data = get_all_lines(f_input, &number);
	fclose(f_input);

	fprintf(fOut, "\nMODEL    %4d\n", nframe+1);
	for(i=0;i<number;i++)
		fprintf(fOut, "%s", data[i]);
	fprintf(fOut, "TER\n");
	fprintf(fOut, "ENDMDL\n");

	free(data);
	return 0;
}

int add_data_to_file(char *fn_input, FILE *fRad)	{
	int i = 0;
	FILE *f_input;
	int number=0;
	char **data=NULL;
	char **SplitData = NULL;
	gmx_bool bCenXYZ=FALSE;

	f_input = fopen(fn_input, "r");
	data = get_all_lines(f_input, &number);
	fclose(f_input);

	for(i=0;i<number;i++)	{
		if (data[i]==NULL)
			continue;

		if(strstr(data[i],"cenxyz.cvec      radius")!=NULL)	{
			bCenXYZ = TRUE;
			continue;
		}

		if(strstr(data[i],"Minimum radius found")!=NULL)	{
			bCenXYZ = FALSE;
			continue;
		}

		if(bCenXYZ)	{
			if(is_first_numerics(data[i])){
				SplitData = split_by_space(data[i]);
				fprintf(fRad,"%s   %s\n",SplitData[0],SplitData[1]);
				free(SplitData);
			}
		}
	}
	free(data);
	return 0;
}

int calculate_com(t_topology *top, int indsize, atom_id *index, rvec *x, rvec com)	{
    int d, i;
	real mass = 0;

	for(i=0;(i<indsize);i++) {
		mass += top->atoms.atom[index[i]].m;
	}

	for(d=0;(d<DIM);d++) {
  	  com[d]=0;
  	  for(i=0;(i<indsize);i++) {
  		  com[d] += x[index[i]][d] * top->atoms.atom[index[i]].m * 10;
  	  }
  	  com[d] /= mass;
	}

  return 0;
}

int main (int argc,char *argv[])	{
	  const char *desc[] = {
			  "To calculate channel radius using hole program. \"hole\" program should",
			  "be present in PATH variable. If -fit is enabled,",
			  "the molecule will be translated, rotated and centered at the origin. ",
			  "First, check with -pdb option for first frame that the hole program ",
			  "is able to identify channel correctly. Modify the values of -endrad, ",
			  "-sample,-cvect and -cpoint for correct identification of channel with ",
			  "-pdb option. The channel could be visualized by using this pdb file "
			  "in any visualization program."
	  };
	  real endrad = 5;
	  real sample = 0.5;
	  rvec cvect = {0, 0, 1}, cpoint = {999, 999, 999};
	  gmx_bool bFit=TRUE;
	  static char *radfile[] = {"simple.rad"};
	  output_env_t oenv;
	  t_pargs pa[] = {
			  { "-fit",    TRUE,  etBOOL, {&bFit},    "To fit structure" },
			  { "-endrad", FALSE, etREAL, {&endrad},  "radius above which the hole2 program regards a result as an indicating that the end of the pore has been reached" },
			  { "-sample", FALSE, etREAL, {&sample},  "The distance between the planes for the sphere centers." },
			  { "-cvect",  FALSE, etRVEC, {cvect},    "This specified a vector which lies in the direction of the channel/pore." },
			  { "-cpoint", FALSE, etRVEC, {cpoint},   "A point which lies within the channel. If not given, center of mass will be used." },
			  { "-rad",    FALSE, etSTR,  {&radfile}, "Name of the file specifying van der Waals radii for each atom." }
	    };

	  t_filenm   fnm[] = {
	     { efTRX, "-f",   NULL,      ffREAD },
	     { efTPS, NULL,   NULL,      ffREAD },
	     { efNDX, NULL,   NULL,      ffOPTRD },
	     { efDAT, "-o",  "radius", ffWRITE },
	     { efPDB, "-pdb",  "sphpdb", ffOPTWR }
	   };

#define NFILE asize(fnm)
   int npargs;
   CopyRightMsg();
   npargs = asize(pa);
   parse_common_args(&argc,argv, PCA_CAN_TIME | PCA_TIME_UNIT | PCA_BE_NICE ,
 		             NFILE,fnm,npargs,pa, asize(desc),desc,0,NULL,&oenv);


   FILE *fRad, *fOutPDB;

   t_trxstatus *status;
   t_topology top;
   char title[STRLEN];
   int        ePBC, natoms, nframe=0;
   real       t ;
   matrix     box;
   int 		  indsize, nfit;
   char       *grpnm=NULL,*fitname;
   atom_id    *index=NULL,*ifit=NULL;
   rvec       *xp, *x;
   char       prefix_name[32], pdbfile[32], hole_outfile[32], hole_outPDB[32];
   const char *fnOutPDB=NULL;
   char       hole_cmd[1024];
   FILE *tmpf;
   int i=0;
   gmx_bool bOutPDB=FALSE;

   //Reading topology
   read_tps_conf(ftp2fn(efTPS,NFILE,fnm),title,&top,&ePBC,&xp,NULL,box,FALSE);

   if(opt2fn_null("-pdb",NFILE,fnm)!=NULL)	{
	   fnOutPDB = opt2fn_null("-pdb",NFILE,fnm);
	   bOutPDB = TRUE;
   }

   if(bFit)	{
	   printf("\nChoose a group for the least squares fit\n");
	   get_index(&top.atoms,ftp2fn_null(efNDX,NFILE,fnm),1,&nfit,&ifit,&fitname);
	   if (nfit < 3)
		 gmx_fatal(FARGS,"Need >= 3 points to fit!\n");
   }

   real *w_rls=NULL;
   if(bFit)	{
 	   snew(w_rls,top.atoms.nr);
 	   for(i=0; (i<nfit); i++)
 		   w_rls[ifit[i]]=top.atoms.atom[ifit[i]].m;
    }


   //Getting index
   get_index(&top.atoms,ftp2fn_null(efNDX,NFILE,fnm),1,&indsize,&index,&grpnm);

   if (bFit)
	   reset_x(nfit,ifit,top.atoms.nr,NULL,xp,w_rls);
   natoms=read_first_x(oenv,&status,ftp2fn(efTRX,NFILE,fnm),&t,&x,box);
   if (bFit)	{
	   reset_x(nfit,ifit,top.atoms.nr,NULL,x,w_rls);
	   do_fit(natoms,w_rls,xp,x);
   }

   if( (cpoint[XX]==999) && (cpoint[YY]==999) && (cpoint[ZZ]==999) )
	   calculate_com(&top, indsize, index, x, cpoint);

   //Creating name for temporary PDB file
   strcpy(prefix_name,"ddXXXXXX");
   sprintf(prefix_name,"ddXXXXXX");
   gmx_tmpnam(prefix_name);
   remove(prefix_name);
   sprintf(pdbfile,"%s.pdb",prefix_name);

   if ((tmpf = fopen(pdbfile,"w")) == NULL)
 	  gmx_fatal(FARGS,"Can not open pdb file %s",pdbfile);
   fclose(tmpf);
   remove(pdbfile);

   //Creating name for output file of analyze program
   sprintf(hole_outfile,"%s.out",prefix_name);
   sprintf(hole_outPDB,"%s_sphere.pdb",prefix_name);


   fRad = ffopen(opt2fn("-o",NFILE,fnm),"w");
   fprintf(fRad,"#Axis \t Radius\n");

   if(bOutPDB)
	   fOutPDB = ffopen(fnOutPDB,"w");

   do	{
	   if (bFit)	{
		   reset_x(nfit,ifit,top.atoms.nr,NULL,x,w_rls);
		   do_fit(natoms,w_rls,xp,x);
	   }

	   //Creating variable for executing hole
	   if(bOutPDB)
		   sprintf(hole_cmd,"hole >%s <<EOF\ncoord %s\nradius %s\ncvect %2.2f %2.2f %2.2f\nsample %2.2f\nendrad %2.2f\ncpoint %2.2f %2.2f %2.2f\nsphpdb %s\nEOF", \
		   			   hole_outfile, pdbfile, radfile[0], cvect[XX], cvect[YY], cvect[ZZ], sample, endrad, cpoint[XX], cpoint[YY], cpoint[ZZ], hole_outPDB );
	   else
		   sprintf(hole_cmd,"hole >%s <<EOF\ncoord %s\nradius %s\ncvect %2.2f %2.2f %2.2f\nsample %2.2f\nendrad %2.2f\ncpoint %2.2f %2.2f %2.2f\nEOF", \
			   hole_outfile, pdbfile, radfile[0], cvect[XX], cvect[YY], cvect[ZZ], sample, endrad, cpoint[XX], cpoint[YY], cpoint[ZZ] );

	   tmpf = ffopen(pdbfile,"w");
	   write_pdbfile_indexed(tmpf,NULL,&top.atoms,x,ePBC,box,' ',-1,indsize,index,NULL,TRUE);
	   ffclose(tmpf);

	   if(0 != system(hole_cmd))
		   gmx_fatal(FARGS,"Failed to execute command: %s",hole_cmd);

	   fprintf(fRad,"\n# Time = %15.5f\n",t);
	   add_data_to_file(hole_outfile, fRad);

	   if(bOutPDB)
		   cat_pdb(nframe, hole_outPDB , fOutPDB);

	   remove(pdbfile);
	   remove(hole_outPDB);
	   nframe++;

   }while(read_next_x(oenv,status,&t,natoms,x,box));


   fprintf(stdout, "Thanks for using gmx_hole!!!\n");
   fprintf(stdout, "\n++++ PLEASE READ AND CITE THE FOLLOWING REFERENCE ++++\n");
   fprintf(stderr, "-------- -------- ------------------- -------- ----------\n");
   fprintf(stderr, "O.S. Smart, J.M. Goodfellow and B.A. Wallace (1993)\n");
   fprintf(stderr, "The Pore Dimensions of Gramicidin A\n");
   fprintf(stderr, "Biophysical Journal 65:2455-2460.\n");
   fprintf(stderr, "-------- -------- ------------------- -------- ----------\n");

   return 0;
}
