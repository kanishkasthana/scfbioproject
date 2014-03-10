//Made by Kanishk Asthana with help from Ms.Priyanka Dhingra
//Program to patch two protein fragments
//Date: 30-April-2010
//Input: pdb1 pdb2 12 34 out(pdb1 is first fragment, pdb2 is second fragment,out is output file)
//to excute type this command: ./a.out <file1>.pdb <file2>.pdb 12 34 out
//Updated on 1 May 2010
//Updated on 2 May 2010
//Updated on 4 June 2010 (File I/O bugs resolved)
//Modified on 5 June 2010 1:45 AM(New equations(untested))
//Updated on 5 June 2010 7:33 PM

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

void print_atom(pdb);

int main(int argc,char *argv[])
{


	FILE *fpw,*fpd,*fout;
	fpw=fopen(argv[1],"r");
	fpd=fopen(argv[2],"r");
	char store[100],store2[100];
	int cnt=0,cnt2=0,totres_f1,totres_f2,secAt=0;
	int start_patch=atoi(argv[3]),end_patch=atoi(argv[4]);
	pdb frag1[50000],frag2[50000];

	cout<<"Start:"<<start_patch<<" End:"<<end_patch<<endl;
	//scanning information from the first fragment
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

       //scanning information from the second fragment
	while(fgets(store2,100,fpd)!=NULL)
	{

	 sscanf(store2,"%s %d %s %s %d %lf %lf %lf",frag2[cnt2].atom,&frag2[cnt2].atomNo,frag2[cnt2].atomtype,//
	 frag2[cnt2].amino,&frag2[cnt2].aminoNo,&frag2[cnt2].x,&frag2[cnt2].y,&frag2[cnt2].z);
	 cnt2++;
	}

	fclose(fpw);
	fclose(fpd);


	totres_f1=cnt; totres_f2=cnt2;
	cout<<"total number of atoms in two fragments "<<totres_f1<<"  "<<totres_f2<<endl;



 //First part: Calculating values of l,m and n and reference atoms


	double l,m,n,a,b,c;
	pdb p1,p2,p1d,p2d,ref,refd,temp;
	int len=(end_patch-start_patch)+1,i,f1=0,f2=0,f3=0;
	 for(i=0;i<cnt;i++)//Taking only Backbone atoms for calculations from first fragment
	 {
	  if(frag1[i].aminoNo==(start_patch +1) && f1==0)//Reference is first atom of second AA
	   {
	    ref=frag1[i];
	    f1++;
	   }
	  if(frag1[i].aminoNo==(start_patch +2) && f2==0)
	   {
	    p1=frag1[i];

	    f2++;
	   }
	  if(frag1[i].aminoNo==(start_patch +3) && f3==0)
	   {
	    p2=frag1[i];
	    f3++;
	   }

	  }

	  print_atom(ref);
	  print_atom(p1);
	  print_atom(p2);

	  f1=0;f2=0;f3=0;

	  for(i=0;i<cnt2;i++)//Taking only backbone atoms for calculations from second fragment
	 {
	  if(frag2[i].aminoNo==2 && f1==0)//Reference is first atom of second AA
	   {
	    refd=frag2[i];
	    f1++;
	   }
	  if(frag2[i].aminoNo==3 && f2==0)
	   {
	    p1d=frag2[i];
	    f2++;
	   }
	  if(frag2[i].aminoNo==4 && f3==0)
	   {
	    p2d=frag2[i];
	    f3++;
	   }

	 }

	   print_atom(refd);
	   print_atom(p1d);
	   print_atom(p2d);
	//Finding coordinates of p1,p2 wrt ref

	pdb r1,r2,r1d,r2d;

	r1=p1;r2=p2;r1d=p1d;r2d=p2d;

	r1.x=p1.x-ref.x;
	r1.y=p1.y-ref.y;
	r1.z=p1.z-ref.z;
	r2.x=p2.x-ref.x;
	r2.y=p2.y-ref.y;
	r2.z=p2.z-ref.z;
	r1d.x=p1d.x-refd.x;
	r1d.y=p1d.y-refd.y;
	r1d.z=p1d.z-refd.z;
	r2d.x=p2d.x-refd.x;
	r2d.y=p2d.y-refd.y;
	r2d.z=p2d.z-refd.z;

	double dx1,dx2,dy1,dy2,dz1,dz2;
	dx1=r1d.x-r1.x;
	dx2=r2d.x-r2.x;
	dy1=r1d.y-r1.y;
	dy2=r2d.y-r2.y;
	dz1=r1d.z-r1.z;
	dz2=r2d.z-r2.z;

	//Finding values of l,m and n
	if(dx1==0.0 && dx2==0.0 && dy1==0.0 && dy2==0.0 && dz1==0.0 && dz2==0.0)
	 {
	  l=0;
	  m=0;
	  n=1;
	 }
	else
	 {
	  a=((dy1*dz2)-(dy2*dz1));
	  b=((dz1*dx2)-(dz2*dx1));
	  c=((dx1*dy2)-(dx2*dy2));
	  l=a/(pow(((a*a)+(b*b)+(c*c)),0.5));
	  m=b/(pow(((a*a)+(b*b)+(c*c)),0.5));
	  n=c/(pow(((a*a)+(b*b)+(c*c)),0.5));
	 }
	//Finding values of (root 1-x^2) for l,m and n
	double ls,ms,ns;
	ls=pow((1-l*l),0.5);
	ms=pow((1-m*m),0.5);
	ns=pow((1-n*n),0.5);

	print_atom(r1);
	print_atom(r2);
	print_atom(r1d);
	print_atom(r2d);

   //Part two: Finding the rotation angle between the two sets of coordinates

    double xf,xfd,yf,yfd,zf,zfd,xs,xsd,ys,ysd,zs,zsd,cosphi,sinphi;
    //First rotation about X axis
    xf=r1.x;
    xfd=r1d.x;
    yf=((n*r1.y)-(m*r1.z))/(pow((n*n)+(m*m),0.5));
    yfd=((n*r1d.y)-(m*r1d.z))/(pow((n*n)+(m*m),0.5));
    zf=((m*r1.y)+(n*r1.z))/(pow((n*n)+(m*m),0.5));
    zfd=((m*r1d.y)+(n*r1d.z))/(pow((n*n)+(m*m),0.5));

    //Second rotation about Y axis
    xs=(ls*xf)-(l*zf);
    ys=yf;
    zs=(l*xf)+(ls*zf);
    xsd=(ls*xfd)-(l*zfd);
    ysd=yfd;
    zsd=(l*xfd)+(ls*zfd);
    //Angle between two fragments when axis of rotation is z axis

    cosphi=((xs*xsd)+(ys*ysd))/((xsd*xsd)+(ysd*ysd));
    sinphi=((xs*ysd)-(ys*xsd))/((xsd*xsd)+(ysd*ysd));

  //Part three: Finding the coordinates of the atoms of the second fragment wrt reference atom

    for(i=0;i<cnt2;i++)
    {
     temp=frag2[i];
     frag2[i].x=temp.x-refd.x;
     frag2[i].y=temp.y-refd.y;
     frag2[i].z=temp.z-refd.z;
    }

  //Part four: Rotation about reference atom for second fragment

    for(i=0;i<cnt2;i++)
    {
     double c1,c2,c3,ca,cb,x,y,z,xd,yd,zd,xdd,ydd,zdd;
     x=frag2[i].x;
     y=frag2[i].y;
     z=frag2[i].z;

     //Applying first rotation

     xd=x;
     yd=((n*y)-(m*z))/pow((m*m)+(n*n),0.5);
     zd=((m*y)+(n*z))/pow((m*m)+(n*n),0.5);

     //Applying second rotation

     xdd=(ls*xd)-(l*zd);
     ydd=yd;
     zdd=(l*xd)+(ls*zd);

     //Applying third rotation

     c1=(xdd*cosphi)+(ydd*sinphi);
     c2=(ydd*cosphi)-(xdd*sinphi);
     c3=zdd;

     //Moving frag2 elements to relative frag1 position
     frag2[i].x=((c1*pow((m*m)+(n*n),0.5))-c3)/((m*m)+(n*n)-l);
     ca=(c3-(l*frag2[i].x));
     cb=(c2*pow((m*m)+(n*n),0.5));
     frag2[i].y=((m*ca)+(n*cb))/((m*m)+(n*n));
     frag2[i].z=((n*ca)-(m*cb))/((m*m)+(n*n));

    }

  //Part five: Translation of atoms of second fragment to the first one
    for(i=0;i<cnt2;i++)
    {
     temp=frag2[i];
     frag2[i].x=(temp.x)+(ref.x);
     frag2[i].y=(temp.y)+(ref.y);
     frag2[i].z=(temp.z)+(ref.z);
     print_atom(frag2[i]);
    }


  //Part Six: Merging frag2 with frag1
   int op=0;

   for(i=0;i<cnt2;i++)
   {
    if((frag2[i].aminoNo<len))
     {
      op++;
     }
   }
//changing residue number of frag 2
  for(i=op;i<cnt2;i++)
  {
     frag2[i].aminoNo=frag2[i].aminoNo+(start_patch-1);
     frag2[i].atomNo=secAt+1;
     frag1[secAt]=frag2[i];
     secAt++;
   }
// Calculating new value of Total Residues

   totres_f1=secAt;

     //writing output of patching
	fout=fopen(argv[5],"w");
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

void print_atom(pdb p)
{
 cout<<p.atomNo<<" "<<p.aminoNo<<" "<<p.x<<" "<<p.y<<" "<<p.z<<endl;

}
/*NOTE: "TER" AND "END" AT THE END OF FILES CAN BE INTERPETED AS LINES, LEADING
TWO EXTRA LINES BEING WRITTEN TO THE END OF FINAL PATCHED FILE*/