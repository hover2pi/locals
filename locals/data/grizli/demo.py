def go():
    import matplotlib.pyplot as plt
    import numpy as np

    from grizli import multifit
    import seaborn as sns
    sns_colors = colors = sns.color_palette("cubehelix", 8)
    
    ################
    # Load the 2D spectrum object
    MW_EBV = 0.0148 # SFD dust map
    mb = multifit.MultiBeam('j100025+021651_04949.beams.fits', psf=True, MW_EBV=MW_EBV, group_name='j100025+021651', fcontam=0.5)
    
    ################
    # Fit Pickles templates to the observed spectrum
    pickles = np.load(os.getenv('GRIZLI')+'/templates/stars_pickles.npy')[0]
    
    N = len(pickles)
    chi2 = np.zeros(N)
    best_template = ''
    chimin = 1e30
    for i, t in enumerate(pickles):
        t_i = {}
        t_i[t] = pickles[t]
        out = mb.xfit_at_z(z=0, templates=t_i, fitter='nnls', fit_background=True, get_uncertainties=False, get_design_matrix=False, pscale=None, COEFF_SCALE=1e-19, get_components=False, huber_delta=4)
        chi2[i], coeffs = out[0], out[1]
        print('{0:30} {1:10.1f} ({2:30} {3:10.1f})'.format(t, chi2[i], best_template, chimin))
        if chi2[i] < chimin:
            chimin = chi2[i]
            best_template = t
    
    # Get normalized best-fit template
    t_i = {}
    t_i[best_template] = pickles[best_template]
    tfit = mb.template_at_z(z=0, templates=t_i, fit_background=True, fitter='nnls', fwhm=1400, get_uncertainties=2)
    
    ################
    # Drizzle 2D figure, with best-fit template as the model
    hdu, fig = mb.drizzle_grisms_and_PAs(fcontam=0.5, flambda=False, kernel='point', size=32, zfit=tfit)
    fig.savefig('{0}_{1:05d}.stack.png'.format(target, id))
    hdu.writeto('{0}_{1:05d}.stack.fits'.format(target, id), clobber=True)
    plt.close()
    
    ################
    # 1D figure, optimal extractions
    best_model = mb.get_flat_model([tfit['line1d'].wave, tfit['line1d'].flux])
    flat_model = mb.get_flat_model([tfit['line1d'].wave, tfit['line1d'].flux*0+1])
    bg = mb.get_flat_background(tfit['coeffs'])
    
    bin = 1 # relative to native resolution of the grism
    sp = mb.optimal_extract(mb.scif[mb.fit_mask] - bg, bin=bin)
    spm = mb.optimal_extract(best_model, bin=bin)
    spf = mb.optimal_extract(flat_model, bin=bin)
    
    # Can be multiple grisms, here just G141
    for g in sp:
        sn = sp[g]['flux']/sp[g]['err']
        clip = sn > 3
        clip = spf[g]['flux'] > 0.2*spf[g]['flux'].max()
        
        scale = 1
        
        plt.errorbar(sp[g]['wave'][clip]/1.e4, (sp[g]['flux']/spf[g]['flux']/scale)[clip]/1.e-19, (sp[g]['err']/spf[g]['flux']/scale)[clip]/1.e-19, marker='.', color='k', alpha=0.4, linestyle='None', elinewidth=0.5, label=g)
        plt.plot(sp[g]['wave']/1.e4, spm[g]['flux']/spf[g]['flux']/1.e-19, color=sns_colors[4], linewidth=2, alpha=0.8, zorder=10, label='Model, {0}'.format(best_template))
        #ymax = np.maximum(ymax, (spm[g]['flux']/spf[g]['flux']/1.e-19)[clip].max())
        #ymin = np.minimum(ymin, (spm[g]['flux']/spf[g]['flux']/1.e-19)[clip].min())
    
    plt.plot(tfit['line1d'].wave/1.e4, tfit['line1d'].flux/1.e-19, color='r', zorder=-10, alpha=0.4)
    
    plt.xlim(0.8, 1.9)
    plt.ylim(10,60)
    
    plt.xlabel(r'$\lambda$ ($\mu$m)')
    plt.ylabel(r'$f_\lambda$ ($10^{-19}$ erg/s/cm$^2$/$\mathrm{\AA}$)')
    plt.legend()
    plt.grid()
    plt.tight_layout(pad=0.1)
    plt.savefig('{0}_{1:05d}.sed.png'.format(mb.group_name, id))