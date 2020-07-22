"""
    pyXSPEC.py  -  the pyXSPEC collection for X-IFU data analysis
    ---------------------------------------------------------------------------------
    Author: Mehdy Lefkir
    Modified by: V. Fioretti (INAF/OAS) valentina.fioretti@inaf.it
    ---------------------------------------------------------------------------------
    Dependencies:
    - python 2.7
    - numpy
    - scipy
    - matplotlib
    - PyPDF2
    - pyXSPEC running on Python 2.7
    ---------------------------------------------------------------------------------
    Example:
    import pyXIFU as px
    px.ml_get_errors('model.xcm','cstat',n_cores=8,level=2.706)
    """


import numpy as np
from xspec import *
import os
from scipy import interpolate  
import matplotlib.backends.backend_pdf
from PyPDF2 import PdfFileMerger, PdfFileReader
import sys
import shutil
import datetime
from plotting import ml_plots,ml_plotting_statistics_errors

Xset.chatter=0

def ml_interpolation_statistics_errors(init_value,x_graph,y_graph,hard_min_hit,level):
    """#--Spline interpolation--"""
    xnew=np.linspace(min(x_graph),max(x_graph),10000)
    f=interpolate.InterpolatedUnivariateSpline(x_graph,y_graph)
    yreduced=y_graph-level
    freduced=interpolate.InterpolatedUnivariateSpline(x_graph, yreduced)
    if hard_min_hit:
        return 0,np.abs(init_value-freduced.roots()[0]),xnew,f
    return init_value-freduced.roots()[0],freduced.roots()[1]-init_value,xnew,f
    
def ml_get_errors(filexcm,statistic,selection='all',blacklist=[''],n_cores=8,level=2.706,plot_statistic=True,interp_method="linear"):
    """Main function to evaluate errors of an XSPEC model.

    Parameters
    ----------
    filexcm : str
        Name of the .xcm XSPEC file to load both the model and the data.
    statistic : {'cstat', 'chi'}
        Statistic of the fit method.
    selection : 'all' or list of int
        Select the parameters used to get the errors with a list of parameter number.
    blacklist : list of strings
        Contains the list of the parameters's name to be frozen before the error computation.
        For instance if blacklist=['Sigma'] all parameters named 'Sigma' will be frozen and the error on those 
        parameters won't be computed.
        Default is ['']
    n_cores : float
        Number of cores to set for the XSPEC parallel variable.
    level : float
        Chi2 level to evaluate a confidence interval. Default is chi2 = 2.706 
        for a 90% confidence interval.
    plot_statistic : bool
        Used to plot or not the statistic graphs.
    interp_method : str
        Method for the interpolation of the statistic. "linear" uses a linear interpolation and the brentq method
        for a the root finding. "spline" uses a spline interpolation and find the roots with a scipy method in the spline class.
    
    +-----------+---------+---------+---------+---------+---------+---------+---------+---------+---------+
    |   sigma   |  1.00s  |  1.28s  |  1.64   |  1.96s  |  2.00s  |  2.58s  |  3.00s  |  3.29s  |  4.00s  |
    +-----------+---------+---------+---------+---------+---------+---------+---------+---------+---------+
    | conf_int  | 68.27%  | 80.00%  | 90.00%  | 95.00%  | 95.45%  | 99.00%  | 99.73%  | 99.90%  | 99.99%  |
    | level     | 1.000   | 1.642   | 2.706   | 3.841   | 4.000   | 6.635   | 9.000   | 10.828  | 16.000  |
    +-----------+---------+---------+---------+---------+---------+---------+---------+---------+---------+

    """
    filename=filexcm
    start=datetime.datetime.now()
    if not os.path.isdir(filename+"_plots"):
        os.mkdir(filename+"_plots")
    Xset.restore(filexcm+".xcm")
    try:
        if statistic=='cstat':
            Fit.statMethod='cstat'
        elif statistic=='chi':
            Fit.statMethod='chi'
        else : 
            raise ValueError('Wrong statistic method entered !')
    except ValueError:
         print "Please restart the script and change the statistic method"
    
    Fit.query='no'
    para_sigma=4.
    d=datetime.datetime.now()
    title="Title: "+filename+"Date: "+d.strftime("%c")
    head='para_nb name best_fit_value error_min error_max hard_min_hit hard_max_hit'
    dt = np.dtype([('para_nb',  np.int32, (1,)),('name', np.str_, 16) ,('best_fit_value', np.float64, (1,)),('error_min', np.float64, (1,)),('error_max', np.float64, (1,)),('hard_min_hit', np.bool, (1,)),('hard_max_hit', np.bool, (1,))])
    out_pdf = filename+'_plots.pdf'
    metadata={'Creator': 'Mehdy Lefkir', 'Author': 'Mehdy Lefkir', 'Title': 'Errors plots on model'}
    free_pars,to_be_frozen=[],[]
    for i in range(1,AllModels(1).nParameters+1):
        if AllModels(1)(i).link=='' and not AllModels(1)(i).frozen :
                free_pars.append(i)
        if not AllModels(1)(i).frozen and AllModels(1)(i).name in blacklist:
                to_be_frozen.append(i)
    if not selection=="all":
        free_pars=selection
        for i in free_pars:
            if AllModels(1)(i).frozen:
                print "<  WARNING  > : Parameter ",i,' ',AllModels(1)(i).name,' is frozen but in the selection ! It will be remove from the list.'
                selection.remove(i)
        free_pars=selection
    print "Free parameters =",free_pars
    print "Number of free parameters =",len(free_pars)
    print "Selected parameters :", selection
    print "Blacklisting all parameters with name :",blacklist
    print "Parameters blacklisted =",to_be_frozen
    for i in to_be_frozen:
        AllModels(1)(i).frozen=True
    Fit.nIterations=100 
    Fit.criticalDelta= 0.01
    Fit.perform()
    os.system("rm "+filexcm+".xcm"); Xset.save(filexcm+".xcm",info="a")
    Xset.parallel.leven = n_cores
    fitstatmin=Fit.statistic
    print "Fit statistic =",fitstatmin
    if os.path.isfile(filename+'_list.txt'):
        Array=np.loadtxt(filename+"_list.txt",dtype=dt)
        if len(Array["para_nb"].flatten())==len(free_pars): 
            if (Array["para_nb"].flatten()==free_pars).all():
                print '<  WARNING  > : Errors on this model were already computed !'
                print '<  INFO  > : Stopping script'
                return
            else: 
                j=free_pars.index(int(Array["para_nb"][-1]))+1
                print '<  INFO  > : Restarting steppar from parameter ',free_pars[j],' ',AllModels(1)(free_pars[j]).name
        else:
            j=free_pars.index(int(Array["para_nb"][-1]))+1
            print '<  INFO  > : Restarting steppar from parameter ',free_pars[j],' ',AllModels(1)(free_pars[j]).name
    else:
        Array=np.delete(np.zeros((1),dtype=dt),0)
        print '<  INFO  > : Initializing steppar'
        j=0      
    if plot_statistic :  mergedObject = PdfFileMerger()
        
    newbestfit=False
    while j<len(free_pars):
        Xset.restore(filexcm+".xcm")
        para_nb=free_pars[j] 
        if AllModels(1)(para_nb).link!='' or AllModels(1)(para_nb).frozen==True :
            print "<  INFO  > : Frozen parameter : I pass"
        else :  
            if plot_statistic : pdf = matplotlib.backends.backend_pdf.PdfPages(str(para_nb)+'.pdf')
            newbestfit=False
            cost_list,par_list=np.array([[],[]]),np.array([])       
            initial_value=AllModels(1)(para_nb).values[0]
            step_steppar=AllModels(1)(para_nb).sigma/para_sigma
            hardcap_hit=[False,False] # hardcap_hit[0]=hardcapmin, hardcap_hit[1]=hardcapmax
            if step_steppar <=0 and np.abs(initial_value - AllModels(1)(para_nb).values[2]) < 1e-8 :
                print "<  WARNING  > : Parameter pegged at the hard lower limit",par_value,AllModels(1)(para_nb).values[2]
                hardcap_hit[0]=True
            if step_steppar <=0 and np.abs(initial_value - AllModels(1)(para_nb).values[5]) < 1e-8 :
                print "<  WARNING  > : Parameter pegged at the hard upper limit"
                hardcap_hit[1]=True
            if step_steppar <= 0 : step_steppar=np.abs(AllModels(1)(para_nb).values[0]/10.)
            print "<  INFO  > : Starting steppar on parameter ",para_nb,' :',AllModels(1)(para_nb).name
            print "<  INFO  > : Initial value :", initial_value
            step_steppar_cur=step_steppar
            for par_dir in [-1,1]:               
                if par_dir==-1: index_cap=0
                elif par_dir==1: index_cap=1
                if index_cap:
                    cost_list=np.append(cost_list,0)
                    par_list=np.append(par_list,initial_value)
                    print "<  INFO  > : Initiating direction ==> right"
                else:
                    print "<  INFO  > : Initiating direction <== left"
                step=para_sigma
                dstat,n_fits=0,0
                step_steppar_cur=step_steppar 
                while not newbestfit and dstat<level+0.1 and not hardcap_hit[index_cap]:
                    Xset.restore(filexcm+".xcm")
                    par_value=initial_value+par_dir*step*step_steppar_cur
                    if index_cap==0:
                        if par_value < AllModels(1)(para_nb).values[2] :
                            hardcap_hit[0]=True
                            print "<  WARNING  > : ---Hard min hit---"
                    elif index_cap==1:
                        if par_value > AllModels(1)(para_nb).values[5] :
                            hardcap_hit[1]=True
                            print "<  WARNING  > : ---Hard max hit---"
                    if hardcap_hit[index_cap]==True:
                        print "<  WARNING  > : Hard cap hit, continue"
                    else :
                        AllModels(1).setPars({para_nb:par_value})
                        AllModels(1)(para_nb).frozen=True
                        Fit.perform()
                        n_fits+=1
                        dstat=Fit.statistic-fitstatmin
                        if dstat < -Fit.criticalDelta :
                            print "New miminum statistic found",dstat
                            AllModels(1)(para_nb).frozen=False
                            fitstatmin=Fit.statistic
                            os.system("rm "+filexcm+".xcm") ; Xset.save(filexcm+".xcm")
                            os.system("rm "+filename+'_list.txt')
                            plt.close('all')
                            Array=np.delete(np.zeros((1),dtype=dt),0)
                            newbestfit=True
                            j,n_fits=0,0 #with this line the whole procedure restarts completely from the first parameter when a new best fit is found!
                        else :                        
                            if dstat > level and n_fits<=2:
                                step_steppar_cur=step_steppar_cur/4.
                                n_fits,dstat=0,0
                                step=1
                                print "<  WARNING  > : Not enough points",step,step_steppar_cur,par_value
                                if index_cap==0:
                                    cost_list,par_list=np.array([[],[]]),np.array([])
                            else :
                                cost_list=np.append(cost_list,dstat) ; par_list=np.append(par_list,par_value)
                                print "<  STEP  > : ",int(step),par_value, "dstat=",dstat,initial_value-par_value
                                step=step+1
            if not newbestfit:
                #------ Finding the errors depends on the hard cap hit variable ------
                if hardcap_hit==[False,False]:
                    par_list, cost_list = zip(*sorted(zip(par_list, cost_list)))
                    par_list,cost_list=np.array(par_list),np.array(cost_list)
                    err_min,err_max,new_x,f=ml_interpolation_statistics_errors(initial_value,np.sort(par_list),cost_list,'None',level,interp_method)
                    if plot_statistic :
                        fig=ml_plotting_statistics_errors(new_x,f,err_min,err_max,par_list,cost_list,initial_value,para_nb,statistic,level,filexcm+".xcm",interp_method)
                        pdf.savefig(fig,bbox_inches='tight')
                        pdf.close()
                elif hardcap_hit==[True,True] :
                    err_min=initial_value-AllModels(1)(para_nb).values[2]
                    err_max=AllModels(1)(para_nb).values[5]-initial_value
                elif hardcap_hit[0]:
                    par_list, cost_list = zip(*sorted(zip(par_list, cost_list)))
                    par_list,cost_list=np.array(par_list),np.array(cost_list)
                    err_min=initial_value-AllModels(1)(para_nb).values[2]
                    err_max,new_x,f=ml_interpolation_statistics_errors(initial_value,np.sort(par_list),cost_list,hardcap_hit,level,'linear')
                    if plot_statistic :
                        fig=ml_plotting_statistics_errors(new_x,f,err_min,err_max,par_list,cost_list,initial_value,para_nb,statistic,level,filexcm+".xcm",interp_method)
                        pdf.savefig(fig,bbox_inches='tight')
                        pdf.close()
                elif hardcap_hit[1]:
                    par_list, cost_list = zip(*sorted(zip(par_list, cost_list)))
                    par_list,cost_list=np.array(par_list),np.array(cost_list)
                    err_min,new_x,f=ml_interpolation_statistics_errors(initial_value,np.sort(par_list),cost_list,hardcap_hit,level,'linear')
                    err_max=AllModels(1)(para_nb).values[5]-initial_value
                    if plot_statistic :
                        fig=ml_plotting_statistics_errors(new_x,f,err_min,err_max,par_list,cost_list,initial_value,para_nb,statistic,level,filexcm+".xcm",interp_method)
                        pdf.savefig(fig,bbox_inches='tight')
                        pdf.close()
                initial_value = np.array([initial_value])
                hardcap_hit_array = np.array([hardcap_hit[0]])
                val_cost = np.concatenate([initial_value, np.sort(par_list), cost_list, hardcap_hit_array])
                #val_cost=np.array([initial_value,np.sort(par_list),cost_list,hardcap_hit[0]])
                np.save(str(para_nb)+'_val_cost',val_cost,allow_pickle=True)
                Array=np.append(Array,np.array([(int(para_nb),AllModels(1)(para_nb).name,initial_value,err_min,err_max,hardcap_hit[0],hardcap_hit[1])],dtype=dt))                   
                AllModels(1)(para_nb).frozen=False
                if plot_statistic :
                    if os.path.exists(str(para_nb)+'.pdf'):
                        shutil.move(str(para_nb)+'.pdf',filename+"_plots/"+str(para_nb)+'.pdf')
                    elif os.path.exists(filename+"_plots/"+str(para_nb)+'.pdf'):
                        print "<  INFO  > :  The statistic plot was already moved"
                    else:
                        print "<  WARNING  > : The statistic plot was not found !"
                if os.path.exists(str(para_nb)+'_val_cost.npy'):
                    shutil.move(str(para_nb)+'_val_cost.npy',filename+"_plots/"+str(para_nb)+'_val_cost.npy')
                elif os.path.exists(filename+"_plots/"+str(para_nb)+'_val_cost.npy'):
                    print "<  INFO  > : The statistic array file was already moved"
                else:
                    print "<  WARNING  > : The statistic array file was not found !"
                #if plot_statistic : mergedObject.append(PdfFileReader(filename+"_plots/"+str(para_nb)+'.pdf', 'rb'))
        if not newbestfit:
            j=j+1
            info="Fit statMethod: "+str(Fit.statMethod)+" | "+"Fit statTest: "+str(Fit.statTest)+" | "+"Fit statistic: "+str(Fit.statistic)+" | "+"Fit DOF: "+str(Fit.dof)+" | "+"Confidence level: "+str(level)+" \n"
            np.savetxt(filename+'_list.txt',Array,header=info+head,fmt='%1.d %s %1.9f %1.9f %1.9f %1.d %1.d')
            print "Results :",Array
    #if plot_statistic : mergedObject.write(filename+"_error_plots.pdf")
    #convert_to_excel(filename,dt)
    end=datetime.datetime.now()
    print "<  INFO  > :  Finished in ",str(end-start)

"""
def convert_to_excel(filename,dtype):
    dt = np.dtype([('para_nb',  np.int32, (1,)),('name', np.string_, 16) ,('best_fit_value', np.float64, (1,)),('error_min', np.float64, (1,)),('error_max', np.float64, (1,))])
    df=pd.read_csv(filename+'_list.txt',header=None,delim_whitespace=True,comment="#",names=["para_nb","name","best_fit_value","error_min","error_max"])
    df.to_excel(filename+'.xlsx', 'Errors', index=False)
"""

if __name__ == '__main__':
    ml_get_errors()
