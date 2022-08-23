{
  
  Double_t volt;//The desired value. 
  
  cout<<"Voltage increase: ";
  cin>>volt;
  
  Double_t start_volt[65]; //EPICS GUI uses 1 -> 64. So we need an extra array cell. 
  Double_t set_volt[65];


  start_volt[1] = 810; 
  start_volt[2] = 810; 
  start_volt[3] = 650;
  start_volt[4] = 650;
  start_volt[5] = 650;
  start_volt[6] = 650;
  start_volt[7] = 670;
  start_volt[8] = 670;
  start_volt[9] = 770;
  start_volt[10] = 700;
  start_volt[11] = 700;
  start_volt[12] = 700;
  
  start_volt[13] = 800;
  start_volt[14] = 800;
  start_volt[15] = 800;
  start_volt[16] = 820;
  start_volt[17] = 830;
  start_volt[18] = 830;
  start_volt[19] = 900;
  start_volt[20] = 900;
  start_volt[21] = 830;
  start_volt[22] = 830;
  start_volt[23] = 830;
  start_volt[24] = 830;

  start_volt[25] = 830;
  start_volt[26] = 850;
  start_volt[27] = 870;
  start_volt[28] = 870;
  start_volt[29] = 870;
  start_volt[30] = 880;
  start_volt[31] = 970;
  start_volt[32] = 970;
  start_volt[33] = 970;
  start_volt[34] = 970;
  start_volt[35] = 1000;
  start_volt[36] = 1000;

  start_volt[37] = 1050;
  start_volt[38] = 1050;
  start_volt[39] = 1200;
  start_volt[40] = 1200;
  start_volt[41] = 1250;
  start_volt[42] = 1250;
  start_volt[43] = 1150;
  start_volt[44] = 1250;
  start_volt[45] = 1150;
  start_volt[46] = 1150;
  start_volt[47] = 1200;
  start_volt[48] = 1150;

  start_volt[49] = 1150;
  start_volt[50] = 1150;
  start_volt[51] = 1200;
  start_volt[52] = 1100;
  start_volt[53] = 1200;
  start_volt[54] = 1200;
  start_volt[55] = 1300;
  start_volt[56] = 1150;
  start_volt[57] = 1300;
  start_volt[58] = 1200;
  start_volt[59] = 1350;
  start_volt[60] = 1300;

  start_volt[61] = 1200;//was 1200
  start_volt[62] = 1200;
  start_volt[63] = 1050;
  start_volt[64] = 1050;

  
  
  for (int i = 1; i<=64; i++){
    set_volt[i] = start_volt[i];//setting everything to default in case I don't want to change all of them in the next loop. 
  }

  for (int i = 3;i<=64;i++){
    set_volt[i] = start_volt[i] + volt;
  }


 

  char buffer[256];
  struct tm *sTm;
  time_t now = time (0);
  sTm = gmtime (&now);
  strftime(buffer, sizeof(buffer), "%a_%b_%d_%H_%M_%S", sTm);

    FILE *fp = fopen(Form("/adaqfs/home/a-onl/sbs/SBS_Replay/bbcal_macros/GRINCH_macros/HV_scripts/GRINCH_HV_settings_write_%1.0f_#_%s.set",volt,buffer),"w");

    cout << "written to "<< Form("/adaqfs/home/a-onl/sbs/SBS_Replay/bbcal_macros/GRINCH_macros/HV_scripts/GRINCH_HV_settings_write_%1.0f_#_%s.set",volt,buffer) <<endl;
 
 
 fprintf(fp,"# %s %s\n",buffer, "EDT 2021 HVFrame-0(rpi17:2001) demand voltage set(DV)" );
 fprintf(fp,"rpi17:2001 S10 DV %4.1f %4.1f %4.1f %4.1f 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0\n",
	 set_volt[61],
	 set_volt[62],
	 set_volt[63],
	 set_volt[64]);
 fprintf(fp,"rpi17:2001 S11 DV %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f\n",
	 set_volt[49],
	 set_volt[50],
	 set_volt[51],
	 set_volt[52],
	 set_volt[53],
	 set_volt[54],
	 set_volt[55],
	 set_volt[56],
	 set_volt[57],
	 set_volt[58],
	 set_volt[59],
	 set_volt[60]);
 fprintf(fp,"rpi17:2001 S12 DV %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f\n",
 	 set_volt[37],
	 set_volt[38],
	 set_volt[39],
	 set_volt[40],
	 set_volt[41],
	 set_volt[42],
	 set_volt[43],
	 set_volt[44],
	 set_volt[45],
	 set_volt[46],
	 set_volt[47],
	 set_volt[48]);
 fprintf(fp,"rpi17:2001 S13 DV %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f\n",
	 set_volt[25],
	 set_volt[26],
	 set_volt[27],
	 set_volt[28],
	 set_volt[29],
	 set_volt[30],
	 set_volt[31],
	 set_volt[32],
	 set_volt[33],
	 set_volt[34],
	 set_volt[35],
	 set_volt[36]);
 fprintf(fp,"rpi17:2001 S14 DV %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f\n",
	 set_volt[13],
	 set_volt[14],
	 set_volt[15],
	 set_volt[16],
	 set_volt[17],
	 set_volt[18],
	 set_volt[19],
	 set_volt[20],
	 set_volt[21],
	 set_volt[22],
	 set_volt[23],
	 set_volt[24]);
 fprintf(fp,"rpi17:2001 S15 DV %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f\n",
	 set_volt[1],
	 set_volt[2],
	 set_volt[3],
	 set_volt[4],
	 set_volt[5],
	 set_volt[6],
	 set_volt[7],
	 set_volt[8],
	 set_volt[9],
	 set_volt[10],
	 set_volt[11],
	 set_volt[12]);

 fclose(fp);


 /*
 for (Int_t i= 0; i<13; i++)
   {
     if(i==0) cout<<"\n\nMODULE S15"<<endl;
     if(i==8) cout<<"\nMODULE S14"<<endl;
     if(i>8)
       {
	 if(i==9)	 cout    <<"HV BOARD #"<<i<<"  0"<<endl;
	 cout    <<"HV BOARD #"<<i+1<<" "<<set_volt[i]<<endl;

       }
     else
       {
	 cout    <<"HV BOARD #"<<i<<" "<<set_volt[i]<<endl;
       }

   }
 */


}
