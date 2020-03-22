#################################################################################
##
## Author:  Nat Goodman
## Created: 20-02-26
##          from misig/dat_confi.R created 19-09-23
##          from dat_ovrfx.R created 19-02-19
##          from dodata.R created 19-02-18
##          from ovrfx.R created 19-02-03 
##          from siglo.R 19-01-01
##          from repwr/R/repwr.R created 17-10-05 
##           and repwr/R/sim.R created 18-05-03
##
## Copyright (C) 2019 Nat Goodman.
## 
## Extract relevant data from
## "Health Insurance and Mortality: Experimental Evidence from Taxpayer Outreach"
## by Jacob Goldin, Ithai Z. Lurie, and Janet McCubbin, December 201
##
## This software is open source, distributed under the MIT License. See LICENSE
## file at https://github.com/natgoodman/NewPro/FDR/LICENSE 
##
#################################################################################
## ---- Extract relevant data from input tables ----
dat_nudge=function() {
  table1=read_table1();
  table2=read_table2();
  ## for reasons unknown, tables 1 and 2 report numbers of subjects differently
  ## hence this roundabout calculation
  n1.all=table1[39,'All.Included'];
  n1.treated=table1[39,'Treatment'];
  n1.control=table1[39,'Control'];
  n2.treated=round(unlist(table2[3,-1]*(n1.treated/n1.all)));
  n2.control=round(unlist(table2[3,-1]*(n1.control/n1.all)));
  n2.treated=round(unlist(n2.treated));
  n2.control=round(unlist(n2.control));
  n2.all=n2.treated+n2.control;
  rate.control=unlist(table2[2,-1]/100);
  rate.treated=rate.control+unlist((table2[1,-1]/100));
  ## rate.ct=rate.control/rate.treated;
  ## cvr.treated=n2.treated*(table2[1,-1]+table2[2,-1])/100;
  ## cvr.control=n2.control*table2[2,-1]/100;
  cvr.treated=n2.treated*rate.treated;
  cvr.control=n2.control*rate.control;
  cvr.treated=round(unlist(cvr.treated));
  cvr.control=round(unlist(cvr.control));
  cvr.all=cvr.treated+cvr.control;
  nudge=data.frame(n.all=n2.all,n.treated=n2.treated,n.control=n2.control,
                   cvr.all,cvr.treated,cvr.control,
                   rate.treated,rate.control);
  rownames(nudge)=cq(full,m0,m1_5,m6_10,m11,m0_10);
  save_data(nudge);
  invisible(nudge);
}
## simulate hypergeometric or binomial distribution. NOT USED in running code
## adapted from misig/R/dat.R:dosim_fixd
dosim_hyper=function(g1,g2,n,m=param(m.sim)) dosim('hyper',g1,g2,n,m);
dosim_binom=function(g1,g2,n,m=param(m.sim)) dosim('binom',g1,g2,n,m);
dosim=function(what,g1,g2,n,m=param(m.sim)) {
  if (n<=0) stop(paste('n must be positive, not',n));
  param(m1,verbose);
  vwhat=paste0('dosim_',what,':');    # used in verbose messages
  replace=if(what=='binom') TRUE else FALSE;
  cases=expand.grid(g1=g1,g2=g2,n=n,m=m);
  if (nrow(cases)>0)
    withrows(cases,case,{
      simdir=dirname_sim(what);
      if (!dir.exists(simdir)) dir.create(simdir);
      file=filename_sim(what,g1,g2,n,m,simdir=simdir);
      if (file.exists(file)&&(is.na(save)||!save)) {
        if (verbose) print(paste('>>>',vwhat,basename(file),'exists. skipping'));
        return();
      }
      n.all=g1+g2;
      more=m; i=1;
      while(more>0) {
        m1=min(m1,more);
        if (verbose) print(paste(sep=' ','>>>',vwhat,nvq(i,g1,g2,n)));
        sim=replicate(m,sample.int(n.all,n,replace=replace));
        ## convert n=1 case to 1-row matrix. from stackoverflow.com/questions/14614946. Thx!
        if (n==1) dim(sim)=c(1,m);
        save_sim_tmp(sim,i,g1,g2,n,m);
        more=more-m1; i=i+1;
      };
      ## consolidate subfiles into one
      sim=cat_sim(file=file);
    });
  return();
}
