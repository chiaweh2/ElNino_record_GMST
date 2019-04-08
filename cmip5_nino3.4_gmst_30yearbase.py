#!/usr/bin/env python


"""
Nino3.4 index V.S. GMST

The script extracts Nino3.4 and GMST time series from different CMIP5 models to
find El Nino event and record breaking GMST event. After the two events are
identify, the script find the record-breaking GMST following the El Nino.

For identify the two events seperately
- Nino3.4 monthly threshold +-0.5 for ENSO events
    (can be modified)
- Record high GMST is defined by the first year of each dataset
    (can be change to based on the hist. run if one has the GMST available)

Steps
---
1. Use monthly ONI to identify El Nino event (period)
    - ONI must pass threshold for more than 5 continuous months
2. Use yearly GMST to identify record breaking year
3. Pick record-breaking GMST happened within the El Nino period with 3-month-lag
4. Store the event in "record_elnino" and "record_gmst_yr"
5. Count the number of El Nino event relate to the record-breaking GMST
6. Count the total number of El Nino event
7. Count the number of record-breaking annual GMST related to El Nino
8. Count the total number of record-breaking annual GMST
---


"""
# to avoid x window forwarding
import matplotlib
matplotlib.use('Agg') # Must place before importing matplotlib.pyplot or pylab!

import sys
import os
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
import datetime
import scipy.stats
import matplotlib.transforms as mtransforms

#### Define the base directory
basedir=os.getcwd()
databasedir=basedir+'/data/'

#### read in the model names
model_info=pd.read_excel('%s/info/CMIP5_Models_Grid_Resolution_vERC.xlsx'\
                        %(basedir),header=0)
model_name=model_info['CMIP5 Model']
institude_name=model_info['Institute'].unique()
climate_sen=model_info['climate sensitivity']




#### function to find record-breaking GMST and El Nino
def record_gmst_induced_by_elnino(ds_stat,databasedir,\
        p,i,m,e,t,r,c,ens,v,filename_keyword2,var_name2,starttime,endtime,\
        basedir, exttext, OUTPUT_RECORD=True, OUTPUT_FIGURE=True):

    datadir=[]
    for comp in range(2):
        datadir.append(os.path.join(databasedir,"CMIP5",\
                                    p,i,m,e,t,r[comp],c[comp],ens,v[comp]))

    pickfile=[]
    if os.path.isdir(datadir[0]) and os.path.isdir(datadir[1]):
        for comp in range(2):
            for file in os.listdir(datadir[comp]):
                if file.endswith(".nc"):
                    if file.startswith(filename_keyword2[comp]) :
                        pickfile.append(file)
    else:
        return ds_stat

    if os.path.isdir(datadir[0]) and os.path.isdir(datadir[1]):
        print '---------------------------------------------'
        print datadir[0]
        print datadir[1]
        if len(pickfile)==2:
            try:
                ds_tos=xr.open_dataset(os.path.join(datadir[0],pickfile[0]))
                ds_tas=xr.open_dataset(os.path.join(datadir[1],pickfile[1]))
            except IOError:
                print "-----"
                print "file broken"
                print datafile
                print "-----"
                return ds_stat

            #### check time series length
            length=np.abs(starttime-endtime)*12.
            if np.float(len(ds_tas.time.values)/np.float(length))<0.5:
                print pickfile[1]
                print "tas time series too short (< 50%)"
                return ds_stat
            if np.float(len(ds_tos.time.values)/np.float(length))<0.5:
                print pickfile[0]
                print "tos time series too short (< 50%)"
                return ds_stat


            #### read the yearly gmst time series
            for file in os.listdir(datadir[1]):
                if file.startswith(filename_keyword2[2]):
                    ds_yr=xr.Dataset()
                    ds_tas_yr=xr.open_dataset(os.path.join(datadir[1],file))

            #### store yearly global mean time series in xr.Dataset
            yr_timestamp=[]
            for yy in range(len(ds_tas_yr.time)):
                yr_timestamp.append(\
                    datetime.datetime(ds_tas_yr['time.year'].values[yy],6,15))
            ds_yr['gmst']=xr.DataArray(ds_tas_yr[var_name2[2]].values,\
                                        coords={'time':yr_timestamp},\
                                        dims=['time'])



            #### Store yearly record hot event (based on historical run)
            record_yr=None
            if e not in ['historical']:
                histdatadir=os.path.join(databasedir,\
                                "CMIP5",p,i,m,"historical",t,r[1],c[1],ens,v[1])
                for file in os.listdir(histdatadir):
                    if file.endswith(".nc"):
                        if file.startswith("record_1850_1900"):
                            ds_tas_histrecord=\
                                xr.open_dataset(os.path.join(histdatadir,file))
                            record_yr=\
                                ds_tas_histrecord[var_name2[2]].max().values
            if not record_yr:
                print "No historical run available"
                print " => Using the first year as record"
                record_yr=ds_yr['gmst'].isel(time=range(1)).max().values
                if np.isnan(record_yr):
                    print " => first year NaN"
                    print " => Using the max in first 10 year as record"
                    record_yr=ds_yr['gmst'].isel(time=range(10)).max().values
            print "record gmst: %0.2f degree C "%record_yr

            ini_record=record_yr
            record_ind=[]
            for tas_ind in range(len(ds_yr['gmst'])):
                if ds_yr['gmst'].values[tas_ind] > record_yr:
                    record_yr=ds_yr['gmst'].values[tas_ind]
                    record_ind.append(tas_ind)
            da_gmst_record_yr=ds_yr['gmst'].copy()+np.nan
            da_gmst_record_yr[record_ind]=ds_yr['gmst'][record_ind]



            #### crop the dataset to the desired period
            # ds_tas, ds_tos, ds_tas_yr, ds_yr, da_gmst_record_yr
            crop_temp1=ds_tas[var_name2[1]][ds_tas['time.year']>=starttime].copy()
            crop_temp2=crop_temp1[crop_temp1['time.year']<=endtime].copy()
            ds_tas=None
            ds_tas=xr.Dataset()
            ds_tas[var_name2[1]]=crop_temp2.copy()
            crop_temp1=None
            crop_temp2=None

            crop_temp1=ds_tos[var_name2[0]][ds_tos['time.year']>=starttime].copy()
            crop_temp2=crop_temp1[crop_temp1['time.year']<=endtime].copy()
            ds_tos=None
            ds_tos=xr.Dataset()
            ds_tos[var_name2[0]]=crop_temp2.copy()
            ind_temp=ds_tos[var_name2[0]].notnull()
            mon_ind=np.arange(len(ds_tos[var_name2[0]]))[ind_temp]
            crop_temp1=None
            crop_temp2=None
            ind_temp=None

            crop_temp1=ds_tas_yr[var_name2[2]][ds_tas_yr['time.year']>=starttime].copy()
            crop_temp2=crop_temp1[crop_temp1['time.year']<=endtime].copy()
            ds_tas_yr=None
            ds_tas_yr=xr.Dataset()
            ds_tas_yr[var_name2[2]]=crop_temp2.copy()
            ind_temp=ds_tas_yr[var_name2[2]].notnull()
            yr_ind=np.arange(len(ds_tas_yr[var_name2[2]]))[ind_temp]
            crop_temp1=None
            crop_temp2=None
            ind_temp=None

            crop_temp1=ds_yr['gmst'][ds_yr['time.year']>=starttime].copy()
            crop_temp2=crop_temp1[crop_temp1['time.year']<=endtime].copy()
            ds_yr=None
            ds_yr=xr.Dataset()
            ds_yr['gmst']=crop_temp2.copy()
            crop_temp1=None
            crop_temp2=None

            if len(da_gmst_record_yr[da_gmst_record_yr['time.year']<starttime])>0 :
                ini_record=\
                da_gmst_record_yr[\
                        da_gmst_record_yr['time.year']<starttime\
                                  ].max().values
            crop_temp1=da_gmst_record_yr[da_gmst_record_yr['time.year']>=starttime].copy()
            # print da_gmst_record_yr[da_gmst_record_yr['time.year']<starttime].max()
            crop_temp2=crop_temp1[crop_temp1['time.year']<=endtime].copy()
            da_gmst_record_yr=None
            da_gmst_record_yr=crop_temp2.copy()
            crop_temp1=None
            crop_temp2=None

            # relocate the mean value of monthly array to yearly array
            mean_gmst=ds_tas_yr[var_name2[2]].mean().values
            mean_mon_gmst=ds_tas[var_name2[1]].mean().values
            corr_mean=-(mean_mon_gmst)+mean_gmst
            ds_tas[var_name2[1]].values+=corr_mean

            # # for after cropping the record change especially the case where the cropping is later than 2006
            # temp_record=ini_record
            # for file in os.listdir(datadir[1]):
            #     if file.startswith("gmtas_yr"):
            #         if starttime > xr.open_dataset(os.path.join(datadir[1],file))['time.year'].values[0]:
            #             ini_record=da_gmst_record_yr[da_gmst_record_yr['time.year']<starttime].max().values
            #             if not ini_record or np.isnan(ini_record):
            #                 ini_record=temp_record
            #
            # mask = np.ones(da_gmst_record_yr.values.shape,dtype=bool)
            # mask[yr_ind] = False
            # da_gmst_record_yr.values[mask]=np.nan


            #### remove trend based on CPC method of 30 year base
            da_tos_30year=ds_tos[var_name2[0]].copy()
            for ii in range(0,len(da_tos_30year.values),5*12):
                if ii < 15*12 :
                    da_tos_30year[ii:ii+5*12]=\
                         ds_tos[var_name2[0]][:ii+15*12].mean().values
                elif ii > len(da_tos_30year.values)-15*12:
                    da_tos_30year[ii:ii+5*12]=\
                         ds_tos[var_name2[0]][-15*12+ii:].mean().values
                else:
                    da_tos_30year[ii:ii+5*12]=\
                         ds_tos[var_name2[0]][ii-15*12:ii+15*12].mean().values
            ds_tos[var_name2[0]]-=da_tos_30year


            #### identify monthly El Nino event
            da_elnino=ds_tos[var_name2[0]].copy()
            da_elnino.values[da_elnino<elnino_crit_lower]+=np.nan


            #### counting total El Nino events
            da_elnino_crit=da_elnino.copy()+np.nan
            elnino_event_count=0
            event_length=0
            temp_ind=[]
            for kk in range(len(da_elnino)):
                if da_elnino.values[kk] > 0. :
                    event_length+=1
                    temp_ind.append(kk)
                else:
                    if event_length>=elnino_cont_mon:
                        if elnino_crit_max and elnino_crit_min:
                            # 3 month mean during El Nino period
                            da_temp=da_elnino[temp_ind].rolling(\
                                                dim={"time":3},\
                                                min_periods=3,\
                                                center=True).mean()
                            # 3 month mean determine catagory
                            if da_temp.max() < elnino_crit_max \
                                and da_temp.max() > elnino_crit_min :
                                elnino_event_count+=1
                                da_elnino_crit[temp_ind]=da_elnino[temp_ind]

                        elif elnino_crit_max is None \
                            and elnino_crit_min is None:
                            elnino_event_count+=1
                            da_elnino_crit[temp_ind]=da_elnino[temp_ind]

                        else:
                            sys.exit('please put both min and max El Nino \
                            Criterias or else put both as None')
                            return ds_stat
                    temp_ind=[]
                    event_length=0


            #### Calculate the likelhood of record breaking during El Nino
            # output number of El Nino event
            #  and the assoicated record breaking event
            record_elnino_num=da_elnino.copy()+np.nan
            record_elnino=da_elnino.copy()+np.nan
            record_gmst_yr=da_gmst_record_yr.copy()+np.nan
            record_gmst_yr_elnino_num=da_gmst_record_yr.copy()+np.nan

            #### storing the record breaking GMST and related El Nino event
            record_yr_elnino_event_count=0
            event_length=0
            for kk in range(len(da_elnino)-1-start_count_mon):
                record_flag=0
                if da_elnino_crit.values[kk] > 0. :
                    da_temp_e=da_elnino.isel(time=kk)
                    # year of 3 months after [kk]
                    da_temp_t=da_elnino.isel(time=kk+start_count_mon)
                    elnino_year=da_temp_e['time.year']
                    elnino_month=da_temp_e['time.month']
                    gmst_year=da_temp_t['time.year']
                    # index in da_gmst_record_yr
                    ind1=np.where(da_gmst_record_yr['time.year']==gmst_year)[0][0]

                    if elnino_year == gmst_year:
                        mon6mean=ds_tas[var_name2[1]][kk:kk+6].mean().values

                        if elnino_month < 6 \
                           and np.nanmax(np.abs(\
                                da_gmst_record_yr.isel(time=ind1).values))\
                                <= mon6mean:
                            # record the gmst record
                            record_flag=1
                            record_gmst_yr.values[ind1]=da_gmst_record_yr.values[ind1]
                            record_gmst_yr_elnino_num[ind1]=record_yr_elnino_event_count
                            # record the related elnino
                            record_elnino.values[kk]=da_elnino_crit.values[kk]
                            record_elnino_num[kk]=record_yr_elnino_event_count
                            event_length+=1
                        else:
                            record_flag=0

                    else :
                        if np.nanmax(np.abs(da_gmst_record_yr.isel(time=ind1))) > 0.:
                            # record the gmst record
                            record_flag=1
                            record_gmst_yr.values[ind1]=da_gmst_record_yr.values[ind1]
                            record_gmst_yr_elnino_num[ind1]=record_yr_elnino_event_count
                            # record the related elnino
                            record_elnino.values[kk]=da_elnino_crit.values[kk]
                            record_elnino_num[kk]=record_yr_elnino_event_count
                            event_length+=1
                        else:
                            record_flag=0
                else:
                    if event_length > 0 :
                        record_yr_elnino_event_count+=1
                    record_flag=0
                    event_length=0

            # likelihood define as 0 if no El Nino as all
            if elnino_event_count <1E-15:
                likelihood=0
            else:
                likelihood=float(record_yr_elnino_event_count)\
                            /float(elnino_event_count)*100
            print "likelihood: %0.2f %%"%(likelihood)
            print "number of total El Nino"
            print "- ",elnino_event_count
            print "number of El Nino related to record breaking GMST"
            print "- ",record_yr_elnino_event_count
            print "number of total record breaking GMST (Year)"
            print "- ",len(da_gmst_record_yr[da_gmst_record_yr.notnull()])
            print "number of record breaking GMST related to El Nino (Year)"
            print "- ",len(record_gmst_yr[record_gmst_yr.notnull()])
            # print "- ",len(np.unique(record_gmst_yr[record_gmst_yr.notnull()]['time.year'].values))


            #### storing the DataArray to Dataset
            ds_elnino_month=xr.Dataset()
            ds_elnino_month['record_elnino']=record_elnino
            ds_elnino_month['tot_elnino']=da_elnino_crit
            ds_elnino_month['record_elnino_index']=record_elnino_num
            ds_gmst_yr=xr.Dataset()
            ds_gmst_yr['record_gmst']=record_gmst_yr
            ds_gmst_yr['tot_gmst']=da_gmst_record_yr
            ds_gmst_yr['record_elnino_index']=record_gmst_yr_elnino_num


            #### Calculate the amplitude of record breaking event
            # amp of consecutive event of record breaking temp is cumulative
            # [previous year increase + this year increase]

            # all amp for all record-breaking event
            da_all_record_gmst_amp=da_gmst_record_yr.copy()+np.nan
            nonull_boolean_all_record=da_gmst_record_yr.notnull()

            all_record_gmst_amp=\
              da_gmst_record_yr[nonull_boolean_all_record].copy().values+np.nan
            all_record_gmst_amp[1:]=\
              da_gmst_record_yr[nonull_boolean_all_record].diff('time',n=1).values
            all_record_gmst_amp[0]=\
              da_gmst_record_yr[nonull_boolean_all_record][0]-ini_record
            da_all_record_gmst_amp[nonull_boolean_all_record]=\
              all_record_gmst_amp

            # amp for record-breaking event associated to elnino only
            da_record_gmst_amp=ds_gmst_yr.record_elnino_index.copy()+np.nan
            nonull_boolean=ds_gmst_yr['record_elnino_index'].notnull()

            da_record_gmst_amp[nonull_boolean]=\
              da_all_record_gmst_amp[nonull_boolean]
            ds_gmst_yr['record_gmst_amp']=da_record_gmst_amp

            # cumulate amp for same elnino
            namp=len(da_all_record_gmst_amp[nonull_boolean])
            for kk in range(1,len(nonull_boolean)):
                if nonull_boolean[kk] \
                 and np.abs(ds_gmst_yr['record_elnino_index'][kk-1].values\
                           -ds_gmst_yr['record_elnino_index'][kk].values)<1E-5:
                    ds_gmst_yr['record_gmst_amp'][kk].values+=\
                          ds_gmst_yr['record_gmst_amp'][kk-1].values

            #### calculate ensemble statistics of amplitude
            # calculate the error bar base on the number of standard deviation
            # the standard error is derived base on Students's T distribution
            # (due to small sample size)
            amp_mean=np.nanmean(ds_gmst_yr['record_gmst_amp'].values)
            amp_std=np.nanstd(ds_gmst_yr['record_gmst_amp'].values)

            stTconfint=0.95
            dof=namp-1
            alpha=1.0-stTconfint
            nstd=scipy.stats.t.ppf(1.0-(alpha/2.0),dof)
            conf=nstd*amp_std/np.sqrt(namp)
            if np.isnan(amp_mean):
                amp_mean=0
            if np.isnan(conf):
                conf=0


            if OUTPUT_RECORD is True:
                try :
                    ds_elnino_month.to_netcdf(path=\
                        os.path.join(datadir[0],\
                                    'record_simu%s_%s'%(exttext,pickfile[0])),\
                                     mode='w')
                    ds_gmst_yr.to_netcdf(path=\
                        os.path.join(datadir[1],\
                                    'record_simu%s_%s'%(exttext,pickfile[1])),\
                                     mode='w')
                except IOError:
                    os.remove(os.path.join(datadir[0],\
                                 'record_simu%s_%s'%(exttext,pickfile[0])))
                    os.remove(os.path.join(datadir[1],\
                                 'record_simu%s_%s'%(exttext,pickfile[1])))
                    ds_elnino_month.to_netcdf(path=\
                        os.path.join(datadir[0],\
                                    'record_simu%s_%s'%(exttext,pickfile[0])),\
                                     mode='w')
                    ds_gmst_yr.to_netcdf(path=\
                        os.path.join(datadir[1],\
                                    'record_simu%s_%s'%(exttext,pickfile[1])),\
                                     mode='w')

            if OUTPUT_FIGURE is True:

                #### plotting
                print "Plotting %s"%(m)
                fig=plt.figure(1,figsize=[1,1])
                ax1color='C0'
                ax2color='C1'
                ax1=fig.add_axes([0,0,10,5])

                # instantiate a second axes that shares the same x-axis
                ax2=ax1.twinx()
                ds_tos[var_name2[0]].plot(ax=ax1,label="Nino3.4",color=ax1color)
                ds_tas[var_name2[1]].plot(ax=ax2,label="GMST",color=ax2color)
                ds_yr['gmst'].plot(ax=ax2,label="GMST",color='tab:red')

                #### plot event
                record_yr_line=np.zeros(len(ds_yr['gmst']))+ini_record
                ax2.plot(ds_yr['gmst'].time.values,record_yr_line,\
                         color='tab:red',linestyle='dashed',alpha=0.5)
                record_gmst_yr.plot(ax=ax2,color='tab:red',marker='o',\
                         linestyle='none')
                if elnino_crit_max and elnino_crit_min:
                    elnino_low_line=np.zeros(len(ds_yr['gmst']))+elnino_crit_min
                    elnino_max_line=np.zeros(len(ds_yr['gmst']))+elnino_crit_max
                    ax1.plot(ds_yr['gmst'].time.values,elnino_low_line,\
                             color='tab:blue',linestyle='dashed',alpha=0.5)
                    ax1.plot(ds_yr['gmst'].time.values,elnino_max_line,\
                             color='tab:blue',linestyle='dashed',alpha=0.5)
                trans = mtransforms.blended_transform_factory\
                                         (ax1.transData, ax1.transAxes)
                ax1.fill_between(da_elnino_crit.time.values,0,1,\
                                 where=da_elnino_crit.notnull(),\
                                 facecolor='tab:blue', alpha=0.3,\
                                 transform=trans)

                #### setting the plotting format
                ax1.set_ylabel('Nino3.4 Index ($^o$C)',{'size':'16'},\
                               color=ax1color)
                ax1.tick_params(axis='y',labelsize=14,labelcolor=ax1color)
                ax2.set_ylabel('GMST ($^o$C)',{'size':'16'}, color=ax2color)
                ax2.tick_params(axis='y',labelsize=14,labelcolor=ax2color)
                ax1.grid(linestyle='dashed')
                ax1.set_xlabel('Year',{'size':'16'})
                ax1.set_title("%s %s"%(m,e),{'size':'24'},pad=24)
                xtick=[datetime.datetime(tt,1,15) \
                   for tt in np.arange(((endtime-starttime)/10)+1)*10+starttime]
                ax1.set_xlim(left=datetime.datetime(starttime,1,15),\
                            right=datetime.datetime(endtime,1,15))
                ax1.set_xticks(xtick)
                ax2.set_xticks(xtick)

                ax1.text(0.55,0.18 \
                        ,'Likelihood: %0.2f %%'%(likelihood)\
                        , fontdict={'size':'16'}, transform=ax1.transAxes)
                if conf > 0. :
                    ax1.text(0.55,0.1 \
                        ,'Amplitude: %0.2f $\pm$ %0.2f ($^o$C)'%(amp_mean,conf)\
                        , fontdict={'size':'16'}, transform=ax1.transAxes)
                else :
                    ax1.text(0.55,0.1 \
                        ,'Amplitude: %0.2f ($^o$C)'%(amp_mean)\
                        , fontdict={'size':'16'}, transform=ax1.transAxes)

                # export the figure
                fig.savefig('%s/figure/Nino3.4_GMST_%s_%s_event_yr%s.pdf'\
                            %(basedir,e,m,exttext), dpi=300, facecolor='w',\
                            edgecolor='w',orientation='portrait',\
                            papertype=None,format=None,transparent=False,\
                            bbox_inches="tight", pad_inches=None,
                            frameon=None)
                plt.close()


            iind=np.where(ds_stat['likelihood'].model_name.values==m)[0][0]
            jind=np.where(ds_stat['likelihood'].exp.values==e)[0][0]
            ds_stat['likelihood'][iind,jind]=likelihood
            ds_stat['amplitude'][iind,jind]=amp_mean
            ds_stat['amplitude_conf'][iind,jind]=conf
            ds_stat['total_elnino'][iind,jind]=elnino_event_count
            ds_stat['elnino_record'][iind,jind]=record_yr_elnino_event_count
            ds_stat['gmst_record'][iind,jind]=\
                            len(record_gmst_yr[record_gmst_yr.notnull()])
            ds_stat['total_gmst'][iind,jind]=\
                            len(da_gmst_record_yr[da_gmst_record_yr.notnull()])
            try :
                ds_stat['climate_sensitivity'][iind,jind]=\
                                       np.float(climate_sen[iind])
            except ValueError:
                ds_stat['climate_sensitivity'][iind,jind]=np.nan

            return ds_stat


def init(model,experiment):
    #### create Dataset for statistical result
    ds_stat=xr.Dataset()

    da_like=xr.DataArray(\
                        np.zeros([len(model),len(experiment)])+np.nan,\
                        coords={'model_name':model,'exp':experiment},\
                        dims={'model_name','exp'})
    da_ampl=xr.DataArray(\
                        np.zeros([len(model),len(experiment)])+np.nan,\
                        coords={'model_name':model,'exp':experiment},\
                        dims={'model_name','exp'})
    da_ampl_con=xr.DataArray(\
                        np.zeros([len(model),len(experiment)])+np.nan,\
                        coords={'model_name':model,'exp':experiment},\
                        dims={'model_name','exp'})
    da_tot_elnino_count=xr.DataArray(\
                        np.zeros([len(model),len(experiment)])+np.nan,\
                        coords={'model_name':model,'exp':experiment},\
                        dims={'model_name','exp'})
    da_record_elnino_count=xr.DataArray(\
                        np.zeros([len(model),len(experiment)])+np.nan,\
                        coords={'model_name':model,'exp':experiment},\
                        dims={'model_name','exp'})
    da_record_gmst_count=xr.DataArray(\
                        np.zeros([len(model),len(experiment)])+np.nan,\
                        coords={'model_name':model,'exp':experiment},\
                        dims={'model_name','exp'})
    da_tot_gmst_count=xr.DataArray(\
                        np.zeros([len(model),len(experiment)])+np.nan,\
                        coords={'model_name':model,'exp':experiment},\
                        dims={'model_name','exp'})
    da_climate_sen=xr.DataArray(\
                        np.zeros([len(model),len(experiment)])+np.nan,\
                        coords={'model_name':model,'exp':experiment},\
                        dims={'model_name','exp'})

    ds_stat['likelihood']=da_like
    ds_stat['amplitude']=da_ampl
    ds_stat['amplitude_conf']=da_ampl_con
    ds_stat['total_elnino']=da_tot_elnino_count
    ds_stat['elnino_record']=da_record_elnino_count
    ds_stat['gmst_record']=da_record_gmst_count
    ds_stat['total_gmst']=da_tot_gmst_count
    ds_stat['climate_sensitivity']=da_climate_sen

    return ds_stat


if __name__ == "__main__":

    #### setting
    starttime=2006
    endtime=2100
    elnino_crit_lower=0.5     # The threshold to determine all elnino event regardless of intensity (nino3.4 index in degree C)
    elnino_crit_max=None      # Upper threshold of nino3.4 index in degree C (if don't care about the intensity => set to None)
    elnino_crit_min=None      # Lower threshold of nino3.4 index in degree C (if don't care about the intensity => set to None)
    start_count_mon=3         # To be counted as record breaking GMST which is related to the El Nino event
                              #  the GMST must happen at least # month after the start of the El Nino event
    elnino_cont_mon=5         # To be counted as El Nino year, set the # continuous months which Nino3.4 over the threshold
    OUTPUT_RECORD=True        # Output the record ts data or not
    OUTPUT_STAT=True          # Output the stats data or not
    OUTPUT_FIGURE=True        # Output the figure or not

    #### to extend any text at the end of the file (filename)
    if elnino_crit_max and elnino_crit_min :
        exttext='_s%i_e%i_allnino%0.1f_smon%i_moncon%i_%0.1f_%0.1f.sen30yearbase'\
        %(starttime,endtime,elnino_crit_lower,start_count_mon,elnino_cont_mon,\
          elnino_crit_min,elnino_crit_max)
    else:
        exttext='_s%i_e%i_allnino%0.1f_smon%i_moncon%i.sen30yearbase'\
        %(starttime,endtime,elnino_crit_lower,start_count_mon,elnino_cont_mon)


    #### option setting (to find data in the corresponding directory)
    product=['output1']
    institute=institude_name
    model=model_name
    experiment=['rcp26','rcp45','rcp60','rcp85']

    #### fixed setting  (!!! order cannot be changed)
    time_frequency=['mon']
    realm2=['ocean','atmos']
    cmor_table2=['Omon','Amon']
    ensemble=['r1i1p1']
    variable2=['tos','tas']
    filename_keyword2=['nino3.4_tos','gmtas_tas','gmtas_yr']
    var_name2=['nino34_noclim','gmst_area_weighted_noclim','gmst_area_weighted']

    ds_stat=init(model,experiment)

    #### running seperate figure for exp.
    for e in experiment:
     for p in product:
      for i in institute:
       for m in model:
           t=time_frequency[0]
           r=realm2
           c=cmor_table2
           ens=ensemble[0]
           v=variable2
           ds_stat=record_gmst_induced_by_elnino(ds_stat,databasedir,\
             p,i,m,e,t,r,c,ens,v,filename_keyword2,var_name2,starttime,endtime,\
             basedir, exttext, OUTPUT_RECORD=OUTPUT_RECORD, \
             OUTPUT_FIGURE=OUTPUT_FIGURE)



if OUTPUT_STAT is True :
    try :
        ds_stat.to_netcdf('%s/Nino3.4_GMST_record_stat%s.nc'\
        %(databasedir,exttext),mode='w')
    except IOError:
        os.remove('%s/Nino3.4_GMST_record_stat%s.nc'
        %(databasedir,exttext))
        ds_stat.to_netcdf('%s/Nino3.4_GMST_record_stat%s.nc'\
        %(databasedir,exttext),mode='w')
