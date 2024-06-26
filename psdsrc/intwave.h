/* Macro to make a nice gradient pulse
 *	waveform is passed to macro
 *
 * rev 0	1/23/99	
 * rev 1	10/22/00	allows reverse load (negate wave)
 * rev 3	10/16/15	allows RF pulses and more flexible pw_ input
 */

short iw_buf[50000];
short iw_idx;

@pulsedef

INTWAVE(int_wgname, int_name, int_pos, int_amp, int_res, int_pw,  int_wave, int_dir, int_loggrd){

cv:{
	   float   a_$[int_name];
	   int    ia_$[int_name];
	   int    pw_$[int_name];
	   int   res_$[int_name];
}

insert: cvinit =>{
}

insert: predownload =>{
{
	float target;

	a_$[int_name] = $[int_amp];
	res_$[int_name] = $[int_res];
	gettarget(&target, $[int_wgname], &$[int_loggrd]);
	ia_$[int_name] = (a_$[int_name] / target) * MAX_PG_IAMP;
  	pw_$[int_name] = (int) res_$[int_name] * GRAD_UPDATE_TIME;
}
}

var:{

	    WF_PULSE $[int_name] = INITPULSE;
}

subst:{
  {
	  pulsename(&$[int_name],"$[int_name]");
	  createreserve(&$[int_name], $[int_wgname], res_$[int_name]);
	  if($[int_dir]>0) 		/* load normal */
		  for (iw_idx=0; iw_idx<res_$[int_name]; iw_idx++)
			  iw_buf[iw_idx] = $[int_wave][iw_idx];
	  else 			/* load reversed */
		  for (iw_idx=0; iw_idx<res_$[int_name]; iw_idx++)
			  iw_buf[iw_idx] = -$[int_wave][(res_$[int_name]-iw_idx-1)];
	  movewaveimm(iw_buf,&$[int_name],(int)0,res_$[int_name],TOHARDWARE);
	  setweos(EOS_DEAD,&$[int_name],res_$[int_name]-1);

	  createinstr(&$[int_name], (LONG)($[int_pos]) , pw_$[int_name], ia_$[int_name]);
	  if (($[int_wgname]==TYPRHO1)) 
	  {
		  addrfbits(&$[int_name], 0 , (LONG)($[int_pos]) , pw_$[int_name]);
	  }
  }
}
}




