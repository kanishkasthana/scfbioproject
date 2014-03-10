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
//Updated on 11 June 2010 9:44 PM(Translational Correction)
//Updated twice before 3 July 2010
//Updated on 6 & 7 July 2010(Clash removal and cosphi sinphi correction(pseudo rigid body))
//Updated on 8 July Rotation function defined
//Updated 11 July Bending defined
//Updated 15 July Correction of many mistakes
//Updated 16 July
//Updated 18,19,20,21,22,23,26,27 July
//Can't Seem to figure out why some proteins are still breaking up when I have done the sinphi and cosphi correction, Why? try to figure out...
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

void print_atom(pdb p);
int rotate(pdb frag1[],int totres_f1,int ang,int end_patch,int first_atm,int last_atm);
void bend(pdb frag1[],int totres_f1,int tweek,int end_patch);
int res_dif=0;
int tot_amino;
pdb p1_v,p2_v;
int swing(pdb frag1[],int totres_f1,int tweek,int end_patch,int first_atm,int last_atm);
void swing2(pdb frag1[],int totres_f1,int tweek,int end_patch);

int main(int argc,char *argv[])
{


	FILE *fpw,*fpd,*fout;
	fpw=fopen(argv[1],"r");
	fpd=fopen(argv[2],"r");
	char store[100],store2[100];
        double CA[3],C[3],O[3];

	int cnt=0,cnt2=0,totres_f1,totres_f2,secAt=0;
	int start_patch=atoi(argv[3]),end_patch=atoi(argv[4]);
	pdb frag1[25000],frag2[25000];

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


	double l,m,n,a,b,c,initialsum=0,finalsum=0;
	pdb p1,p2,p1d,p2d,ref,refd,temp;
	int len=(end_patch-start_patch)+1,i,f1=0,f2=0,f3=0;
	 for(i=0;i<cnt;i++)//Taking only Backbone atoms for calculations from first fragment
	 {
	  if(frag1[i].aminoNo==(end_patch) && f1==0)//Reference is first atom of second AA
	   {
	    if(strcmp(frag1[i].amino,"PRO")==0)
             { 
              ref=frag1[i];
	
	      p1=frag1[i-2];

	      p2=frag1[i+10];
	      f1++;
	     }
            else
             { 
              ref=frag1[i];
	
	      p1=frag1[i-2];

	      p2=frag1[i+2];
	      f1++;
	     }
           }

	  }

	  print_atom(ref);
	  print_atom(p1);
	  print_atom(p2);

	  f1=0;f2=0;f3=0;

	  for(i=0;i<cnt2;i++)//Taking only backbone atoms for calculations from second fragment
	 {
	  if(frag2[i].aminoNo==len && f1==0)//Reference is first atom of second AA
	   {
	    if(strcmp(frag2[i].amino,"PRO")==0)
             { 
              refd=frag2[i];
	
	      p1d=frag2[i-2];

	      p2d=frag2[i+10];
	      f1++;
	     }
            else
             { 
              refd=frag2[i];
	
	      p1d=frag2[i-2];

	      p2d=frag2[i+2];
	      f1++;
	     }
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
         cout<<"L2+M2+N2:"<<(l*l+m*m+n*n)<<endl;
	//Finding values of (root 1-x^2) for l,m and n
	double ls,ms,ns;
	ls=pow((1-l*l),0.5);
	ms=pow((1-m*m),0.5);
	ns=pow((1-n*n),0.5);
        cout<<l<<" "<<m<<" "<<n<<endl;
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
    //cout<<"These are the x,y,z:"<<endl;
    //cout<<ref.x+xf<<" "<<ref.y+yf<<" "<<ref.z+zf<<" "<<refd.x+xfd<<" "<<refd.y+yfd<<" "<<refd.z+zfd<<" "<<ref.x+xs<<" "<<ref.y+ys<<" "<<ref.z+zs<<" "<<refd.x+xsd<<" "<<refd.y+ysd<<" "<<refd.z+zsd<<endl;
    //Angle between two fragments when axis of rotation is z axis
     double checkcosphi,checksinphi;
    checkcosphi=((xs*xsd)+(ys*ysd))/((xsd*xsd)+(ysd*ysd));
    checksinphi=((xs*ysd)-(ys*xsd))/((xsd*xsd)+(ysd*ysd));
    
    if(checksinphi>=0 && checkcosphi>=0)
    { 
     cosphi=checkcosphi;
     sinphi=pow((1-(cosphi*cosphi)),0.5);
    }
    if(checksinphi<0 && checkcosphi<0)
    { 
     cosphi=checkcosphi;
     sinphi=-1*pow((1-(cosphi*cosphi)),0.5);
    }
    if(checksinphi<0 && checkcosphi>0)      //Doubt about where to put = sign here!
    { 
     sinphi=checksinphi;
     cosphi=pow((1-(sinphi*sinphi)),0.5);
    }
    if(checksinphi>0 && checkcosphi<0)
    { 
     sinphi=checksinphi;
     cosphi=-1*pow((1-(sinphi*sinphi)),0.5);
    }
    if(checksinphi==0 && checkcosphi>0)
    {
     sinphi=0;
     cosphi=1;
    }
    if(checksinphi==0 && checkcosphi<0)
    {
     sinphi=0;
     cosphi=-1;
    }
    if(checkcosphi==0 && checksinphi>0)
    {
     cosphi=0;
     sinphi=1;
    }
    if(checkcosphi==0 && checksinphi<0)
    {
     cosphi=0;
     sinphi=-1;
    }



    cout<<"COSPHI: "<<cosphi<<" SINPHI:"<<sinphi<<endl;
    cout<<"Cos2phi + sin2Phi:"<<(cosphi*cosphi+sinphi*sinphi)<<endl;

   //Finding initial distances wrt ref
   pdb sum_ref;
   double distance;
   sum_ref=frag2[0];
   for(i=0;i<cnt2;i++)
   {
    distance=pow(pow(frag2[i].x-sum_ref.x,2)+pow(frag2[i].y-sum_ref.y,2)+pow(frag2[i].z-sum_ref.z,2),0.5);
    initialsum=initialsum+distance;
   }


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
     frag2[i].x=(ls*c1+l*c3);
     ca=(c3-l*(ls*c1+l*c3));
     cb=(c2*ls);
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
 //Translational Correction due to non Ideal rigid body motion
 if(frag2[op-2].x!=frag1[secAt-2].x || frag2[op-2].y!=frag1[secAt-2].y || frag2[op-2].z!=frag1[secAt-2].z)
 {
   pdb corr;
   corr.x=frag2[op-2].x-frag1[secAt-2].x;
   corr.y=frag2[op-2].y-frag1[secAt-2].y;
   corr.z=frag2[op-2].z-frag1[secAt-2].z;
   for(i=0;i<cnt2;i++)
   {
    frag2[i].x=frag2[i].x-corr.x;
    frag2[i].y=frag2[i].y-corr.y;
    frag2[i].z=frag2[i].z-corr.z;
   }
 }
//Final Sum
   sum_ref=frag2[0];
   for(i=0;i<cnt2;i++)
   {
    distance=pow(pow(frag2[i].x-sum_ref.x,2)+pow(frag2[i].y-sum_ref.y,2)+pow(frag2[i].z-sum_ref.z,2),0.5);
    finalsum=finalsum+distance;
   }


if(finalsum>=0.95*initialsum && finalsum<=1.05*initialsum)
{ 
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
   tot_amino=frag1[totres_f1-1].aminoNo;
   cout<<"TOTAL AMINO ACIDS: "<<tot_amino<<endl;

//Part 7:Rotating Second Fragment and checking for clashes
 f1=0;
 f2=0;
 int last_atm,first_atm,clashes1=0,clashes2=0,ang=9,low1=0,pos1=0,low2=0,pos2=0,finalpos=0,finalclashes=0;
 for(i=0;i<totres_f1;i++)
   {
   if(frag1[i].aminoNo==end_patch && f1==0)
   {
    last_atm=frag1[i-1].atomNo;
    f1++;
   }
   if(frag1[i].aminoNo==end_patch+1 && f2==0)
   {
    first_atm=frag1[i].atomNo;
    f2++;
   }
   }
   cout<<"LAST ATOM:"<<last_atm<<endl;
   cout<<"FIRST ATOM:"<<first_atm<<endl;
  //Rotation
 pdb test1[25000],test2[25000],test3[25000],test4[25000];
  
 int j,bendclashes1=0,bendpos1=0,num1=20,num2=20;
 int bendclashes2=bendclashes1,bendpos2=bendpos1,finalbendpos=0,finalbendclashes=0,tweek=5,rotpos1=0,rotpos2=0,finalrotpos=0;
 //cout<<"INITIAL CLASHES:"<<bendclashes1<<endl;
   for(i=0;i<totres_f1;i++)
  {
   test1[i]=frag1[i];
   test2[i]=frag1[i];
  }
  low1=rotate(test1,totres_f1,0,end_patch,first_atm,last_atm);
  low2=low1;
  clashes1=low1;
  clashes2=low2;
  cout<<low1<<" 0"<<endl;
  for(i=0;i<=ang;i++)
  { 
   if(clashes1==0)
   {
    pos1=i;
    low1=clashes1;
    break;
   }
   if(clashes2==0)
   {
    pos2=i;
    low2=clashes2;
    break;
   }

   if(clashes1<low1)
   {
    low1=clashes1;
    pos1=i;
   }

   if(clashes2<low2)
   {
    low2=clashes2;
    pos2=i;
   }

   if(i<ang)
   {
    clashes1=rotate(test1,totres_f1,1*num2,end_patch,first_atm,last_atm);
    clashes2=rotate(test2,totres_f1,-1*num2,end_patch,first_atm,last_atm);
    cout<<clashes1<<" "<<(i+1)*num2<<endl;
    cout<<clashes2<<" "<<-1*(i+1)*num2<<endl;
   }
  }

  if(low1<low2)
  {
   finalclashes=low1;
   finalpos=pos1;
  }
  if(low2<low1)
  {
   finalclashes=low2;
   finalpos=(-1)*pos2;
  }
  if(low1==low2)
  {
   if(pos1<=pos2)
   {
    finalclashes=low1;
    finalpos=pos1;
   }
   else
   {
    finalclashes=low2;
    finalpos=(-1)*pos2;
   }
  
  }
 cout<<"FINAL: "<<finalclashes<<" "<<finalpos*num2<<endl;

 
 int right1=finalpos+1,left1=finalpos-1,center1=finalpos;
 int left2=finalpos-1,right2=finalpos+1,center2=finalpos; 
 int centerclashes1,rightclashes1,leftclashes1,centerclashes2,leftclashes2,rightclashes2;
 int check,k=0;

if(finalclashes!=0)
{
 rotpos1=finalpos;rotpos2=finalpos;bendpos1=0;bendpos2=0;
 bendclashes1=finalclashes;bendclashes2=finalclashes;
 
 for(k=1;k<=tweek;k++)
 {
 
  for(i=0;i<totres_f1;i++)
  {
   test1[i]=frag1[i];
   
   test3[i]=frag1[i];
   
  }
  left1=center1-1;right1=center1+1;left2=center2-1;right2=center2+1;
  bend(test1,totres_f1,k*num1,end_patch);
  bend(test3,totres_f1,(-1)*k*num1,end_patch);
  centerclashes1=rotate(test1,totres_f1,center1*num2,end_patch,first_atm,last_atm);

  for(i=0;i<totres_f1;i++)
  {
   test2[i]=test1[i];
  }

  leftclashes1=rotate(test2,totres_f1,-1*num2,end_patch,first_atm,last_atm); 
  rightclashes1=rotate(test1,totres_f1,num2,end_patch,first_atm,last_atm);
  centerclashes2=rotate(test3,totres_f1,center2*num2,end_patch,first_atm,last_atm);

  for(i=0;i<totres_f1;i++)
  {
   test4[i]=test3[i];
  }

  leftclashes2=rotate(test4,totres_f1,-1*num2,end_patch,first_atm,last_atm);
  rightclashes2=rotate(test3,totres_f1,num2,end_patch,first_atm,last_atm);
  cout<<"POSITIVE BEND:"<<k*num1<<endl;
  cout<<leftclashes1<<" "<<centerclashes1<<" "<<rightclashes1<<endl;
  
  if(leftclashes1<=centerclashes1)
  {
   check=rotate(test2,totres_f1,-1*num2,end_patch,first_atm,last_atm);
   cout<<"LEFT:"<<endl;
   cout<<check<<endl;
   while(check<leftclashes1)
   {
    left1=left1-1;
    leftclashes1=check;
    check=rotate(test2,totres_f1,-1*num2,end_patch,first_atm,last_atm);
    cout<<check<<endl;
   }
  } 
  cout<<"FINAL LEFT: "<<leftclashes1<<" "<<num2*left1<<endl;
  if(rightclashes1<=centerclashes1)
  {
   check=rotate(test1,totres_f1,num2,end_patch,first_atm,last_atm);
   cout<<"RIGHT:"<<endl;
   cout<<check<<endl;

   while(check<rightclashes1)
   {
    right1=right1+1;
    rightclashes1=check;
    check=rotate(test1,totres_f1,num2,end_patch,first_atm,last_atm);
    cout<<check<<endl;
   }
  }
  cout<<"FINAL RIGHT: "<<rightclashes1<<" "<<num2*right1<<endl;
  
  cout<<"NEGATIVE BEND:"<<-1*k*num1<<endl;
  cout<<leftclashes2<<" "<<centerclashes2<<" "<<rightclashes2<<endl;

  if(leftclashes2<=centerclashes2)
  {
   check=rotate(test4,totres_f1,-1*num2,end_patch,first_atm,last_atm);
   cout<<"LEFT:"<<endl;
   cout<<check<<endl;
   while(check<leftclashes2)
   {
    left2=left2-1;
    leftclashes2=check;
    check=rotate(test4,totres_f1,-1*num2,end_patch,first_atm,last_atm);
    cout<<check<<endl;
   }
  }
  cout<<"FINAL LEFT: "<<leftclashes2<<" "<<num2*left2<<endl; 
  if(rightclashes2<=centerclashes2)
  {
   check=rotate(test3,totres_f1,num2,end_patch,first_atm,last_atm);
   cout<<"RIGHT:"<<endl;
   cout<<check<<endl;
   while(check<rightclashes2)
   {
    right2=right2+1;
    rightclashes2=check;
    check=rotate(test3,totres_f1,num2,end_patch,first_atm,last_atm);
    cout<<check<<endl;
   }
  }
  cout<<"FINAL RIGHT: "<<rightclashes2<<" "<<num2*right2<<endl;

  if(rightclashes1<centerclashes1 || leftclashes1<centerclashes1)
 {
  if(rightclashes1<leftclashes1)
  {
   center1=right1;
   centerclashes1=rightclashes1;  
  } 

  if(leftclashes1<rightclashes1)
  {
   center1=left1;
   centerclashes1=leftclashes1;
  }
  if(leftclashes1==rightclashes1)
  {
   if((center1-left1)<=(right1-center1))
   {
    center1=left1;
    centerclashes1=leftclashes1;
   }
   
   if((right1-center1)<(center1-left1))
   {
    center1=right1;
    centerclashes1=rightclashes1;
   }
  }
 }
 cout<<"NEW POSITIVE CENTER:"<<centerclashes1<<" "<<center1*num2<<endl;
  if(rightclashes2<centerclashes2 || leftclashes2<centerclashes2)
 {
  if(rightclashes2<leftclashes2)
  {
   center2=right2;
   centerclashes2=rightclashes2;  
  } 

  if(leftclashes2<rightclashes2)
  {
   center2=left2;
   centerclashes2=leftclashes2;
  }
  if(leftclashes2==rightclashes2)
  {
   if((center2-left2)<=(right2-center2))
   {
    center2=left2;
    centerclashes2=leftclashes2;
   }
   
   if((right2-center2)<(center2-left2))
   {
    center2=right2;
    centerclashes2=rightclashes2;
   }
  }
 }
 cout<<"NEW NEGATIVE CENTER:"<<centerclashes2<<" "<<center2*num2<<endl;

 if(centerclashes1==0)
 {
   bendclashes1=centerclashes1;
   bendpos1=k;
   rotpos1=center1;
   break;
 }


 if(centerclashes2==0)
 {
  bendclashes2=centerclashes2;
  bendpos2=k;
  rotpos2=center2;
  break;
 }


 if(centerclashes1<bendclashes1)
 {
   bendclashes1=centerclashes1;
   bendpos1=k;
   rotpos1=center1;
 }


 if(centerclashes2<bendclashes2)
 {
  bendclashes2=centerclashes2;
  bendpos2=k;
  rotpos2=center2;
 }
 cout<<"MIN CLASHES: POS "<<bendclashes1<<" "<<bendpos1*num1<<" NEG "<<bendclashes2<<" "<<(-1)*bendpos2*num1<<endl;
}

 if(bendclashes1<bendclashes2)
 {
  finalbendclashes=bendclashes1;
  finalbendpos=bendpos1;
  finalrotpos=rotpos1;
 }
 if(bendclashes2<bendclashes1)
 {
  finalbendclashes=bendclashes2;
  finalbendpos=(-1)*bendpos2;
  finalrotpos=rotpos2;
 }

 if(bendclashes1==bendclashes2)
 {
  if(bendpos1<=bendpos2)
  {
   finalbendclashes=bendclashes1;
   finalbendpos=bendpos1;
   finalrotpos=rotpos1;
  }
  if(bendpos2<bendpos1)
  {
   finalbendclashes=bendclashes2;
   finalbendpos=(-1)*bendpos2;
   finalrotpos=rotpos2; 
  }
 }

} 

else
{
 finalbendclashes=finalclashes;
 finalbendpos=0;
 finalrotpos=finalpos;  
}
 cout<<"FINAL CLASHES:"<<finalbendclashes<<endl;
 cout<<"FINAL BEND:"<<finalbendpos*num1<<endl;
 cout<<"FINAL ROTATION:"<<finalrotpos*num2<<endl;
 bend(frag1,totres_f1,finalbendpos*num1,end_patch);
 rotate(frag1,totres_f1,finalrotpos*num2,end_patch,first_atm,last_atm);
 int an=end_patch+res_dif; 
   p1_v.x=0;p1_v.y=0;p1_v.z=0;p2_v.x=0;p2_v.y=0;p2_v.y=0;
   
   for(i=0;i<totres_f1;i++)
   {
    if(frag1[i].aminoNo==an)
    {
       ref=frag1[i];
       break; 
    }
   }
   for(i=0;i<totres_f1;i++)
   {
    if(frag1[i].atomNo<ref.atomNo)
     {
      p1_v.x=p1_v.x+(frag1[i].x-ref.x);
      p1_v.y=p1_v.y+(frag1[i].y-ref.y);
      p1_v.z=p1_v.z+(frag1[i].z-ref.z);
     }
    
    if(frag1[i].atomNo>ref.atomNo)
    {
      p2_v.x=p2_v.x+(frag1[i].x-ref.x);
      p2_v.y=p2_v.y+(frag1[i].y-ref.y);
      p2_v.z=p2_v.z+(frag1[i].z-ref.z);   
    }

   }
  for(i=0;i<totres_f1;i++)
  {
   test1[i]=frag1[i];
   test2[i]=frag1[i];
  }
 int num3=10,swingclashes1=0,swingclashes2=0,finalswingclashes=0;low1=0;low2=0;pos1=0;pos2=0;ang=5;
 low1=swing(test1,totres_f1,0,end_patch,first_atm,last_atm);
 low2=low1;
 swingclashes1=low1;
 swingclashes2=low2;
 cout<<low1<<" 0"<<endl;
  for(i=0;i<=ang;i++)
  { 
   if(swingclashes1==0)
   {
    pos1=i;
    low1=swingclashes1;
    break;
   }
   if(swingclashes1<low1)
   {
    low1=swingclashes1;
    pos1=i;
   }
   if(swingclashes2==0)
   {
    pos2=i;
    low2=swingclashes2;
    break;
   }
   if(swingclashes2<low2)
   {
    low2=swingclashes2;
    pos2=i;
   }

   if(i<ang)
   {
    swingclashes1=swing(test1,totres_f1,1*num3,end_patch,first_atm,last_atm);
    swingclashes2=swing(test2,totres_f1,(-1)*num3,end_patch,first_atm,last_atm);
    cout<<swingclashes1<<" "<<(i+1)*num3<<endl;
    cout<<swingclashes2<<" "<<(-1)*(i+1)*num3<<endl;
   }
  }
  if(low1<low2)
  {
   finalswingclashes=low1;
   finalpos=pos1;
  }
  if(low2<low1)
  {
   finalswingclashes=low2;
   finalpos=(-1)*pos2;
  }
  if(low1==low2)
  {
   if(pos1<=pos2)
   {
    finalswingclashes=low1;
    finalpos=pos1;
   }
   else
   {
    finalswingclashes=low2;
    finalpos=(-1)*pos2;
   }
  
  }

if(finalpos!=0)
{
  int modfinalpos=pow(finalpos*finalpos,0.5);
  int sign=finalpos/modfinalpos;
  for(i=1;i<=modfinalpos;i++)
  {
   swing2(frag1,totres_f1,sign*num3,end_patch);
  }
  //finalswingclashes=rotate(frag1,totres_f1,0,end_patch,first_atm,last_atm); 
}
  cout<<"FINAL SWING: "<<finalswingclashes<<" "<<finalpos*num3<<endl;

 


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
 else cout<<"FILE DEFORMED"<<endl;
}

void print_atom(pdb p)
{
 cout<<p.atomNo<<" "<<p.aminoNo<<" "<<p.x<<" "<<p.y<<" "<<p.z<<endl;

}

int rotate(pdb frag1[],int totres_f1,int ang,int end_patch,int first_atm,int last_atm)
  { 
   int an=end_patch+res_dif,i,f1=0,f2=0;
   double sinphi=sin(ang*3.1415926535898/180),cosphi=cos(ang*3.1415926535898/180),a,b,c,l,ls,m,ms,n,ns;
   pdb ref,temp;
   //Bond Selection for rotation   
   for(i=0;i<totres_f1;i++)
   {
    if(frag1[i].aminoNo==an && f1==0 && f2==0)
    {
     if(strcmp(frag1[i].amino,"PRO")==0)
      { 
       ref=frag1[i];
       a=frag1[i].x-frag1[i+10].x;
       b=frag1[i].y-frag1[i+10].y;
       c=frag1[i].z-frag1[i+10].z;
       f1++;
      }
      else
      { 
       ref=frag1[i];	
       a=frag1[i].x-frag1[i+2].x;
       b=frag1[i].y-frag1[i+2].y;
       c=frag1[i].z-frag1[i+2].z;
       f2++;
      }

    }
  }
  l=a/(pow(((a*a)+(b*b)+(c*c)),0.5));
  m=b/(pow(((a*a)+(b*b)+(c*c)),0.5));
  n=c/(pow(((a*a)+(b*b)+(c*c)),0.5));
  ls=pow((1-l*l),0.5);
  ms=pow((1-m*m),0.5);
  ns=pow((1-n*n),0.5);
  double c1,c2,c3,ca,cb,x,y,z,xd,yd,zd,xdd,ydd,zdd;
   for(i=0;i<totres_f1;i++)
    {
     if(frag1[i].atomNo>=ref.atomNo+2 && f2==1)
    {
     temp=frag1[i];
     frag1[i].x=temp.x-ref.x;
     frag1[i].y=temp.y-ref.y;
     frag1[i].z=temp.z-ref.z;
     x=frag1[i].x;
     y=frag1[i].y;
     z=frag1[i].z;

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
     frag1[i].x=(ls*c1+l*c3);
     ca=(c3-l*(ls*c1+l*c3));
     cb=(c2*ls);
     frag1[i].y=((m*ca)+(n*cb))/((m*m)+(n*n));
     frag1[i].z=((n*ca)-(m*cb))/((m*m)+(n*n));
     temp=frag1[i];
     frag1[i].x=(temp.x)+(ref.x);
     frag1[i].y=(temp.y)+(ref.y);
     frag1[i].z=(temp.z)+(ref.z);
     }
    if(frag1[i].atomNo>=ref.atomNo+12 && f1==1)
    {
     temp=frag1[i];
     frag1[i].x=temp.x-ref.x;
     frag1[i].y=temp.y-ref.y;
     frag1[i].z=temp.z-ref.z;
     x=frag1[i].x;
     y=frag1[i].y;
     z=frag1[i].z;

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
     frag1[i].x=(ls*c1+l*c3);
     ca=(c3-l*(ls*c1+l*c3));
     cb=(c2*ls);
     frag1[i].y=((m*ca)+(n*cb))/((m*m)+(n*n));
     frag1[i].z=((n*ca)-(m*cb))/((m*m)+(n*n));
     temp=frag1[i];
     frag1[i].x=(temp.x)+(ref.x);
     frag1[i].y=(temp.y)+(ref.y);
     frag1[i].z=(temp.z)+(ref.z);
     }
    }
   int clashes=0,j;f1=0;f2=0;
   double distance=0;
   
   for(i=0;i<totres_f1;i++)
  {
   if(frag1[i].atomNo<=last_atm)
   {
    for(j=first_atm-1;j<totres_f1;j++)
    {
     distance=pow(pow(frag1[i].x-frag1[j].x,2)+pow(frag1[i].y-frag1[j].y,2)+pow(frag1[i].z-frag1[j].z,2),0.5);
     if(distance<1.50)
     clashes++;
    }
   }
  }
   return clashes; 
}

void bend(pdb frag1[],int totres_f1,int tweek,int end_patch)
{
 int an=end_patch+res_dif,i,f1=0;
   tweek=tweek/3;
   double sinphi=sin(tweek*3.1415926535898/180),cosphi=cos(tweek*3.1415926535898/180),a,b,c,l,ls,m,ms,n,ns;
   pdb ref,p1,p2,temp;
   
   for(i=0;i<totres_f1;i++)
   {
    if(frag1[i].aminoNo==an && f1==0)
    {
      if(strcmp(frag1[i].amino,"PRO")==0)
      { 
       ref=frag1[i];
       p1=frag1[i-2];
       p2=frag1[i+10];
       f1++;
      }
      else
      { 
       ref=frag1[i];	
       p1=frag1[i-2];
       p2=frag1[i+2];
       f1++;
      }
    }
   }
  pdb r1=p1,r2=p2;
  r1.x=p1.x-ref.x;
  r1.y=p1.y-ref.y;
  r1.z=p1.z-ref.z;
  r2.x=p2.x-ref.x;
  r2.y=p2.y-ref.y;
  r2.z=p2.z-ref.z;
  //Cross product
  
  a=(r2.y*r1.z)-(r2.z*r1.y);

  b=(r2.z*r1.x)-(r2.x*r1.z);

  c=(r2.x*r1.y)-(r2.y*r1.x);
    
  l=a/(pow(((a*a)+(b*b)+(c*c)),0.5));
  m=b/(pow(((a*a)+(b*b)+(c*c)),0.5));
  n=c/(pow(((a*a)+(b*b)+(c*c)),0.5));
  ls=pow((1-l*l),0.5);
  ms=pow((1-m*m),0.5);
  ns=pow((1-n*n),0.5);
  
  int k=0,j=0;
  double c1,c2,c3,ca,cb,x,y,z,xd,yd,zd,xdd,ydd,zdd;
  for(j=0;j<totres_f1;j++)
   {
     if(frag1[j].aminoNo>=(an-0) && k<3)
    {
     if(strcmp(frag1[j].atomtype,"N")==0 || strcmp(frag1[j].atomtype,"CA")==0 || strcmp(frag1[j].atomtype,"C")==0)
     {
      ref=frag1[j];
      for(i=j;i<totres_f1;i++)
      {
       temp=frag1[i];
       frag1[i].x=temp.x-ref.x;
       frag1[i].y=temp.y-ref.y;
       frag1[i].z=temp.z-ref.z;
      
       x=frag1[i].x;
       y=frag1[i].y;
       z=frag1[i].z;
 
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
       frag1[i].x=(ls*c1+l*c3);
       ca=(c3-l*(ls*c1+l*c3));
       cb=(c2*ls);
       frag1[i].y=((m*ca)+(n*cb))/((m*m)+(n*n));
       frag1[i].z=((n*ca)-(m*cb))/((m*m)+(n*n));
       temp=frag1[i];
       frag1[i].x=(temp.x)+(ref.x);
       frag1[i].y=(temp.y)+(ref.y);
       frag1[i].z=(temp.z)+(ref.z); 
      }
      k++;
      //cout<<k<<" "<<frag1[j].atomtype<<" "<<frag1[j].aminoNo<<endl;
     }
    }
   } 
}

int swing(pdb frag1[],int totres_f1,int tweek,int end_patch,int first_atm,int last_atm)
{
 int an=end_patch+res_dif,i,f1=0,f2=0;
   tweek=tweek/3;
   double sinphi=sin(tweek*3.1415926535898/180),cosphi=cos(tweek*3.1415926535898/180),a,b,c,l,ls,m,ms,n,ns;
   pdb ref,temp;
   
  pdb r1=p1_v,r2=p2_v;

  //Cross product
  
  a=(r2.y*r1.z)-(r2.z*r1.y);

  b=(r2.z*r1.x)-(r2.x*r1.z);

  c=(r2.x*r1.y)-(r2.y*r1.x);
    
  l=a/(pow(((a*a)+(b*b)+(c*c)),0.5));
  m=b/(pow(((a*a)+(b*b)+(c*c)),0.5));
  n=c/(pow(((a*a)+(b*b)+(c*c)),0.5));
  ls=pow((1-l*l),0.5);
  ms=pow((1-m*m),0.5);
  ns=pow((1-n*n),0.5);
  
  int k=0,j=0;
  double c1,c2,c3,ca,cb,x,y,z,xd,yd,zd,xdd,ydd,zdd;
  for(j=0;j<totres_f1;j++)
   {
     if(frag1[j].aminoNo>=(an-0) && k<3)
    {
     if(strcmp(frag1[j].atomtype,"N")==0 || strcmp(frag1[j].atomtype,"CA")==0 || strcmp(frag1[j].atomtype,"C")==0)
     {
      ref=frag1[j];
      for(i=j;i<totres_f1;i++)
      {
       temp=frag1[i];
       frag1[i].x=temp.x-ref.x;
       frag1[i].y=temp.y-ref.y;
       frag1[i].z=temp.z-ref.z;
      
       x=frag1[i].x;
       y=frag1[i].y;
       z=frag1[i].z;
 
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
       frag1[i].x=(ls*c1+l*c3);
       ca=(c3-l*(ls*c1+l*c3));
       cb=(c2*ls);
       frag1[i].y=((m*ca)+(n*cb))/((m*m)+(n*n));
       frag1[i].z=((n*ca)-(m*cb))/((m*m)+(n*n));
       temp=frag1[i];
       frag1[i].x=(temp.x)+(ref.x);
       frag1[i].y=(temp.y)+(ref.y);
       frag1[i].z=(temp.z)+(ref.z); 
      }
      k++;
      //cout<<k<<" "<<frag1[j].atomtype<<" "<<frag1[j].aminoNo<<endl;
     }
    }
   } 
int clashes=0;j=0;f1=0;f2=0;
   double distance=0;
   
   for(i=0;i<totres_f1;i++)
  {
   if(frag1[i].atomNo<=last_atm)
   {
    for(j=first_atm-1;j<totres_f1;j++)
    {
     distance=pow(pow(frag1[i].x-frag1[j].x,2)+pow(frag1[i].y-frag1[j].y,2)+pow(frag1[i].z-frag1[j].z,2),0.5);
     if(distance<1.50)
     clashes++;
    }
   }
  }
   return clashes; 
}
void swing2(pdb frag1[],int totres_f1,int tweek,int end_patch)
{
 int an=end_patch+res_dif,i,f1=0,f2=0;
   tweek=tweek/3;
   double sinphi=sin(tweek*3.1415926535898/180),cosphi=cos(tweek*3.1415926535898/180),a,b,c,l,ls,m,ms,n,ns;
   pdb ref,temp;
   
  pdb r1=p1_v,r2=p2_v;

  //Cross product
  
  a=(r2.y*r1.z)-(r2.z*r1.y);

  b=(r2.z*r1.x)-(r2.x*r1.z);

  c=(r2.x*r1.y)-(r2.y*r1.x);
    
  l=a/(pow(((a*a)+(b*b)+(c*c)),0.5));
  m=b/(pow(((a*a)+(b*b)+(c*c)),0.5));
  n=c/(pow(((a*a)+(b*b)+(c*c)),0.5));
  ls=pow((1-l*l),0.5);
  ms=pow((1-m*m),0.5);
  ns=pow((1-n*n),0.5);
  
  int k=0,j=0;
  double c1,c2,c3,ca,cb,x,y,z,xd,yd,zd,xdd,ydd,zdd;
  for(j=0;j<totres_f1;j++)
   {
     if(frag1[j].aminoNo>=(an-0) && k<3)
    {
     if(strcmp(frag1[j].atomtype,"N")==0 || strcmp(frag1[j].atomtype,"CA")==0 || strcmp(frag1[j].atomtype,"C")==0)
     {
      ref=frag1[j];
      for(i=j;i<totres_f1;i++)
      {
       temp=frag1[i];
       frag1[i].x=temp.x-ref.x;
       frag1[i].y=temp.y-ref.y;
       frag1[i].z=temp.z-ref.z;
      
       x=frag1[i].x;
       y=frag1[i].y;
       z=frag1[i].z;
 
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
       frag1[i].x=(ls*c1+l*c3);
       ca=(c3-l*(ls*c1+l*c3));
       cb=(c2*ls);
       frag1[i].y=((m*ca)+(n*cb))/((m*m)+(n*n));
       frag1[i].z=((n*ca)-(m*cb))/((m*m)+(n*n));
       temp=frag1[i];
       frag1[i].x=(temp.x)+(ref.x);
       frag1[i].y=(temp.y)+(ref.y);
       frag1[i].z=(temp.z)+(ref.z); 
      }
      k++;
      //cout<<k<<" "<<frag1[j].atomtype<<" "<<frag1[j].aminoNo<<endl;
     }
    }
   } 
}


/*NOTE: "TER" AND "END" AT THE END OF FILES CAN BE INTERPETED AS LINES, LEADING
TWO EXTRA LINES BEING WRITTEN TO THE END OF FINAL PATCHED FILE*/
