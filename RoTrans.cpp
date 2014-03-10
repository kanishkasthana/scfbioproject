//Made by Kanishk Asthana
#include<stdio.h>
#include<iostream.h>
#include<string.h>
#include<stdlib.h>
#include<math.h>

typedef  struct{

	char atom[5];
	char amino[4],atomtype[6];
	int aminoNo;
	int atomNo;
	char cI;
	double x;
	double y;
	double z;

}pdb;


int main(int argc,char *argv[])
{


	FILE *fpw,*fpd,*fout;
	fpw=fopen(argv[1],"r");
	char store[100];

	int cnt=0,totres_f1;
	pdb frag1[25000];
	//scanning information
	while(fgets(store,100,fpw)!=NULL)
	{

	 sscanf(store,"%s %d %s %s %d %lf %lf %lf",frag1[cnt].atom,&frag1[cnt].atomNo,frag1[cnt].atomtype,//
	 frag1[cnt].amino,&frag1[cnt].aminoNo,&frag1[cnt].x,&frag1[cnt].y,&frag1[cnt].z);

	 if(frag1[cnt].aminoNo==end_patch-1)
	  {
	   secAt=frag1[cnt].atomNo;
	  }
         
	 cnt++;
	}

       	fclose(fpw);


	totres_f1=cnt;
	cout<<"total number of atoms in fragment"<<totres_f1<<endl;



     //writing output of patching
      	fout=fopen(argv[2],"w");
	 for(int t=0;t<totres_f1;t++) //fragment 1 contains the updated patched coordinates
	{
		if(strlen(frag1[t].atom)==4)
			fprintf(fout,"ATOM  %5d %-4s %3s  %4d%12.3lf%8.3lf%8.3lf\n",frag1[t].atomNo,//
			frag1[t].atomtype,frag1[t].amino,frag1[t].aminoNo,frag1[t].x,frag1[t].y,frag1[t].z);
		else
			fprintf(fout,"ATOM  %5d  %-3s %3s  %4d%12.3lf%8.3lf%8.3lf\n",frag1[t].atomNo,//
			frag1[t].atomtype,frag1[t].amino,frag1[t].aminoNo,frag1[t].x,frag1[t].y,frag1[t].z);
	}
 
}



/*NOTE: "TER" AND "END" AT THE END OF FILES CAN BE INTERPETED AS LINES, LEADING
TWO EXTRA LINES BEING WRITTEN TO THE END OF FINAL PATCHED FILE*/
