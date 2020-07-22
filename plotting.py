from matplotlib.ticker import ScalarFormatter
import matplotlib.ticker as ticker
from matplotlib.ticker import FuncFormatter
from matplotlib import rc
import numpy as np
import matplotlib.pyplot as plt
from xspec import *

rc('text', usetex=True)
rc('legend',fontsize=20)

def ml_plots(filename,expression,filexcm,x_min=0,x_max=0,sigma_res=3,redshift=0):
    Xset.chatter=0
    Xset.restore(filexcm)
    Plot.device='\null'
    Plot.redshift=redshift
    if redshift>0:
        rest_en='Rest ~'
    else:
        rest_en=''
    plt.rc('legend',fontsize=20)
    plt.rc('xtick',labelsize=20)
    plt.rc('ytick',labelsize=20)
    plt.rc('axes', titlesize=20)     # fontsize of the axes title
    plt.rc('axes', labelsize=20)
    plt.rc('text', usetex=True)
    plt.rcParams['xtick.major.pad']='8'
    plt.rcParams['ytick.major.pad']='8'
    #plt.rc( 'font', size=20, family="Times" )   # use a font with serifs
    if 'de' in expression:
        fig, axs = plt.subplots(2, sharex=True, sharey=False,clear=True, gridspec_kw={'hspace': 0,'height_ratios': [2*1.618, 1.618]},figsize=(12,9),tight_layout=True)
        upper_panel=axs[0]
        lower_panel=axs[1]
    else:
        fig, axs = plt.subplots(1, sharex=True, sharey=False,clear=True, gridspec_kw={'hspace': 0},figsize=(12,9),tight_layout=True)
        upper_panel=axs
    if 'eeufspec' in expression:
            Plot("eeufspec")
            flux_unit="\mathrm{keV}^2~ (\mathrm{Photons}~ \mathrm{cm}^{-2}\,\mathrm{s}^{-1}\,\mathrm{keV}^{-1})"
    elif 'eufspec' in expression:
            Plot("eufspec")
            flux_unit="\mathrm{keV}~ (\mathrm{Photons}~ \mathrm{cm}^{-2}\,\mathrm{s}^{-1}\,\mathrm{keV}^{-1})"
    elif 'ufspec' in expression:
            Plot("ufspec")
            flux_unit="\mathrm{Photons}~ \mathrm{cm}^{-2}\,\mathrm{s}^{-1}\,\mathrm{keV}^{-1}"
            
    Plot.xAxis=('keV')
    x = Plot.x()
    y = np.array(Plot.y())
    xErrs = Plot.xErr()
    yErrs = Plot.yErr()
    folded = np.array(Plot.model())
    linewidth=0.5
    if x_min!=0 or x_max!=0:
        upper_panel.set_xlim(x_min-1e-2,x_max+1e-2)
    else:
        upper_panel.set_xlim(np.min(AllData(1).energies)*(1+redshift),np.max(AllData(1).energies)*(1+redshift)+1e-2)

    upper_panel.set_xscale("log")
    upper_panel.set_yscale("log")
    upper_panel.xaxis.set_minor_formatter(ticker.ScalarFormatter(useMathText=True))
    upper_panel.xaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))    
    energies=np.array([0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,2,3,4,5,6,7,8,9,10,11,12])
    xticks=energies[np.intersect1d(np.where(energies>np.min(AllData(1).energies)*(1+redshift)),np.where(energies<np.max(AllData(1).energies)*(1+redshift)))]
    upper_panel.set_xticks([1,10])
    #upper_panel.set_xticks(, minor = True)
    upper_panel.xaxis.set_minor_formatter(ticker.ScalarFormatter(useMathText=True))

    upper_panel.set_xticklabels([r"$%.2f$" % i for i in [0.6,0.7,0.8,0.9,2,3,4,5,6,7,8,9]],minor=True)
    upper_panel.errorbar(x,y,xerr=xErrs,yerr=yErrs,fmt='none',elinewidth =linewidth, capsize=0,ecolor = 'k')
    upper_panel.step( x, folded,"blue",linewidth=2.)
    upper_panel.set_ylabel(r"$"+flux_unit+"$",labelpad=8)
    upper_panel.set_xlabel(r'$\mathrm{'+rest_en+' Energy~ (keV)}$',fontsize=20)

    upper_panel.tick_params(axis="both",which="both",direction="in",length=4,top=True,right=True)
    upper_panel.tick_params(axis="both",which="major",direction="in",length=10,top=True,right=True)

    if 'de' in expression: 
        Plot("delchi")
        x = Plot.x()
        y = np.array(Plot.y())
        xErrs = Plot.xErr()
        yErrs = Plot.yErr()
        lower_panel.errorbar(x,y,yerr=yErrs,xerr=xErrs,fmt='none',elinewidth =linewidth, capsize=0,markersize=2.,ecolor = 'k')
        lower_panel.axhline(0.,color='k')
        lower_panel.set_ylabel(r"$\sigma=\mathrm{(data-model)}$/$\mathrm{error}$",labelpad=35)
        lower_panel.set_ylim(-np.std(y)*sigma_res,np.std(y)*sigma_res)
        lower_panel.tick_params(axis="both",which="both",direction="in",length=4,top=True,right=True)
        lower_panel.set_xlabel(r'$\mathrm{'+rest_en+' Energy~ (keV)}$',fontsize=20)
    
    fig.savefig(filename+" "+expression+".pdf")    
    
def ml_plotting_statistics_errors(xnew,f,err_min,err_max,par_list,cost_list,initial_value,para_nb,statistic,level,filexcm,interp_method):
    Xset.chatter=0
    Xset.restore(filexcm)
    rc('legend',fontsize=14)
    rc('xtick',labelsize=14)
    rc('axes', titlesize=14)   # fontsize of the axes title  
    
    rc('axes', labelsize=14)
    if statistic=="chi":
        name_tag="\chi^2"
    elif statistic=="cstat":
        name_tag="cstat"
    if (err_max<0.01 and  err_max!=0 )or (err_min<0.01 and  err_min!=0) or ( initial_value<0.01 and not initial_value!=0):
        formatting='%.2e'
    else:
        formatting='%03.2f'

    if interp_method=="linear":
        ynew=f(xnew)+level
    elif interp_method=='spline':
        ynew=f(xnew)
    fig=plt.figure(figsize=(9,8))
    ax = fig.add_axes([0,0,1,1])
    plt.plot(xnew,ynew,"k--")
    plt.scatter(par_list,cost_list,marker='o',color='r')
    plt.axhline(y=level,label=r"$\Delta "+name_tag+"="+str(level)+"$",color='forestgreen')
    plt.text(0.5,0.9,r' Best fit value = $'+str(formatting %initial_value)+'^{+'+str(formatting %err_max)+'}_{-'+str(formatting %err_min)+'}$',fontsize=18, horizontalalignment='center',verticalalignment='center', transform=ax.transAxes)
    reform_name=AllModels(1)(para_nb).name
    if '_' in AllModels(1)(para_nb).name:
        reform_name=reform_name.replace('_',' ',10)
    plt.title("Parameter "+str(para_nb)+" : "+reform_name,fontsize=20)
    if AllModels(1)(para_nb).unit!='':
        if '^-' in AllModels(1)(para_nb).unit:
            reform_unit=AllModels(1)(para_nb).unit
            index_symbol=reform_unit.index('^')
            L=[]
            for s in reform_unit[index_symbol+2:].split():
                if  s.isdigit():
                    L.append(int(s))
                else:
                    break
            reform_unit=reform_unit.replace('^-','^{-',1)
            reform_unit=reform_unit.replace(str(L[0]),str(L[0])+'}')
            units_str=' $('+reform_unit+')$'
        elif '^' in AllModels(1)(para_nb).unit:
            reform_unit=AllModels(1)(para_nb).unit
            index_symbol=reform_unit.index('^')
            L=[]
            for s in reform_unit[index_symbol+1:].split():
                if  s.isdigit():
                    L.append(int(s))
                else:
                    break
            reform_unit=reform_unit.replace('^','^',1)
            reform_unit=reform_unit.replace(str(L[0]),'{'+str(L[0])+'}')
            units_str=' $('+reform_unit+')$'
        else:
            units_str=' ('+AllModels(1)(para_nb).unit+')'
    else:
        units_str=''
    plt.xlabel(r''+reform_name+units_str,fontsize=14)
    plt.ylabel(r"$\Delta "+name_tag+"$",fontsize=14)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.legend(fontsize=15)
    return fig
